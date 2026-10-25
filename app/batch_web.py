"""Batch screening from the browser.

The command-line runner exists for large panels and uses separate processes.
This is the smaller case: a panel pasted into the page, screened in the
background while the tab stays usable, with a row per molecule and a CSV.

It runs molecules one at a time in a single background thread rather than
forking workers. That is slower per molecule than ``run_batch.py``, and
deliberately so: a web request should not spawn a process pool inside the
server that outlives the request or competes with the interactive screen the
user is probably also running. For more than a few dozen molecules, the
command-line runner is the right tool, and the interface says so.
"""
from __future__ import annotations

import csv
import io
import threading
import time
import uuid
from datetime import datetime, timezone
from typing import Any

MAX_CONSTRUCTS = 40

_jobs: dict[str, dict[str, Any]] = {}
_lock = threading.Lock()


def _set(job_id: str, **changes) -> None:
    with _lock:
        job = _jobs.get(job_id)
        if job is not None:
            job.update(changes)


def get(job_id: str, include_records: bool = False) -> dict[str, Any] | None:
    """The job. Records are large and the page polls, so they are withheld
    unless a caller actually needs them."""
    with _lock:
        job = _jobs.get(job_id)
        if not job:
            return None
        out = dict(job)
    if not include_records:
        out.pop("records", None)
    return out


def record_for(job_id: str, name: str) -> dict | None:
    with _lock:
        job = _jobs.get(job_id)
        return (job or {}).get("records", {}).get(name)


def sequences_for(job_id: str, name: str) -> dict | None:
    with _lock:
        job = _jobs.get(job_id)
        return (job or {}).get("sequences", {}).get(name)


def _screen_one(name: str, chains_in: dict[str, str], options: dict) -> dict[str, Any]:
    """One molecule, through the same path the single screen uses."""
    from .analysis import clean_sequence, sequence_metrics
    from .numbering import number_chain
    from .humanness import framework_identity
    from .sequence_panel import annotate_chain, conclusion as panel_conclusion
    from .residue_table import build_residue_table
    from .molecule_flags import molecule_flags
    from .assessment import build_assessment
    from .reference_data import ReferenceDataMissing
    from .translation import resolve as resolve_sequence

    scheme, ph = options["scheme"], options["ph"]
    record: dict[str, Any] = {"scheme": scheme, "ph": ph, "chains": {}}
    for label, raw in chains_in.items():
        sequence = clean_sequence(resolve_sequence(raw)[0], allow_nucleotide=False)
        canonical = number_chain(sequence, "imgt")
        actual = canonical.get("chain_type")
        if (label == "H" and actual != "H") or (label == "L" and actual not in {"K", "L"}):
            raise ValueError(f"The {label} chain label does not match its numbered chain type ({actual}).")
        selected = canonical if scheme == "imgt" else number_chain(sequence, scheme)
        metrics = sequence_metrics(sequence, ph, canonical)
        record["chains"][label] = {
            "sequence": sequence, "numbering": selected, "numbering_imgt": canonical,
            "metrics": metrics, "humanness": framework_identity(sequence, canonical),
            "sequence_bands": annotate_chain(metrics, canonical.get("chain_type") or label)}
    record["sequence_conclusion"] = panel_conclusion(
        {k: v["sequence_bands"] for k, v in record["chains"].items()})

    if options.get("precedents", True):
        from .ptm_precedents import sequence_precedents
        from .clinical_frameworks import clinical_framework_precedents
        for key, fn in (("ptm_precedents", sequence_precedents),
                        ("clinical_framework_precedents", clinical_framework_precedents)):
            try:
                record[key] = fn(record["chains"], scheme)
            except ReferenceDataMissing as exc:
                record[key] = {"cases": [], "sites": [], "unavailable": True, "reason": str(exc)}
    else:
        record["ptm_precedents"] = {"cases": [], "unavailable": True}
        record["clinical_framework_precedents"] = {"sites": [], "unavailable": True}

    record["residue_table"] = build_residue_table(
        record["chains"], precedents=record.get("clinical_framework_precedents"))

    assembly = None
    if options.get("isotype"):
        from . import constant_regions
        try:
            assembly = constant_regions.assemble(
                record["chains"], heavy_isotype=options["isotype"])
        except constant_regions.ConstantRegionError as exc:
            # A full-length basis needs a heavy chain. Only the charge basis
            # depends on it, so a light-only construct falls back to the Fv
            # basis rather than dropping out of the panel entirely.
            record["charge_basis_note"] = (
                f"Full-length basis unavailable for this construct ({exc}). "
                "Charge descriptors are variable-domain only.")
    record["molecule_flags"] = molecule_flags(
        record["chains"], formulation_ph=options.get("formulation_ph", 6.0),
        assay_ph=ph, assembly=assembly)
    # DeepSP costs milliseconds per molecule where the structure route costs
    # seconds, so unlike the structure-derived descriptors it can run on every
    # construct in a panel. This is where it earns its place.
    from . import deepsp as deepsp_module
    if deepsp_module.available():
        try:
            record["deepsp"] = deepsp_module.analyze(
                (record["chains"].get("H") or {}).get("numbering_imgt"),
                (record["chains"].get("L") or {}).get("numbering_imgt"))
        except (ValueError, deepsp_module.DeepSPUnavailable) as exc:
            record["deepsp"] = {"status": "error", "reason": str(exc)}
    else:
        record["deepsp"] = {"status": "not_installed"}

    record["assessment"] = build_assessment(record)
    return record


def _row(name: str, record: dict) -> dict[str, Any]:
    table = record.get("residue_table") or {}
    verdicts = table.get("verdict_counts") or {}
    flags = record.get("molecule_flags") or {}
    # The browser batch is sequence-only. Focused assessment sites require
    # exposure and are empty here, even when the sequence has many motifs.
    sites = [r for r in table.get("rows", []) if r.get("liabilities")]
    cdr = [s for s in sites if str(s.get("imgt_region", "")).startswith("CDR")]
    framework = [s for s in sites if s.get("imgt_region") == "Framework"]
    row: dict[str, Any] = {
        "name": name,
        "cdr_liabilities": len(cdr),
        "framework_liabilities": len(framework),
        "unassigned_liabilities": len(sites) - len(cdr) - len(framework),
        "residues_notable": table.get("notable_rows"),
        "residues_fix": verdicts.get("evaluate", verdicts.get("fix")),
        "residues_watch": verdicts.get("watch"),
        "germline_fr_mismatches": (table.get("germline_counts") or {}).get("mismatch_framework"),
        "charge_basis": flags.get("basis"),
        "charge_basis_note": record.get("charge_basis_note"),
        "worst_flag": flags.get("worst_verdict"),
        "pI": (flags.get("fv") or {}).get("pi"),
        "net_charge": (flags.get("fv") or {}).get("net_charge_assay_ph"),
        "panel": (record.get("sequence_conclusion") or {}).get("state"),
    }
    for flag in flags.get("flags", []):
        # Keep the numeric descriptor sortable. The net_charge flag shares
        # its id with the descriptor and must not overwrite it with text.
        key = "flag_net_charge" if flag["id"] == "net_charge" else flag["id"]
        row[key] = flag.get("verdict")
    # DeepSP columns: the whole-Fv rollups, as percentiles rather than raw
    # values, because a percentile is comparable across a panel and the raw
    # SAP/SCM scale is not interpretable at a glance.
    deepsp = record.get("deepsp") or {}
    if deepsp.get("status") in {"ok", "unreliable"}:
        by_descriptor = {h["descriptor"]: h for h in deepsp.get("highlights", [])}
        for descriptor, column in (("SAP_pos_Fv", "deepsp_sap_pct"),
                                   ("SCM_neg_Fv", "deepsp_scm_neg_pct"),
                                   ("SCM_pos_Fv", "deepsp_scm_pos_pct")):
            row[column] = (by_descriptor.get(descriptor) or {}).get("percentile")
        row["deepsp_flags"] = len(deepsp.get("flagged", []))
        if deepsp.get("status") == "unreliable":
            row["deepsp_flags"] = None
    for label, chain in record.get("chains", {}).items():
        row[f"humanness_{label}"] = (chain.get("humanness") or {}).get("percent")
    return row


def _run(job_id: str, constructs: dict[str, dict[str, str]], options: dict) -> None:
    started = time.time()
    rows: list[dict] = []
    results: list[dict] = []
    # Records are kept so a row can be opened as a full report without
    # re-screening. Bounded by MAX_CONSTRUCTS, and the job is per server
    # session, so this does not grow without limit.
    records: dict[str, dict] = {}
    for index, (name, chains_in) in enumerate(constructs.items(), start=1):
        _set(job_id, current=name, done=index - 1)
        at = time.time()
        try:
            record = _screen_one(name, chains_in, options)
            rows.append(_row(name, record))
            records[name] = record
            results.append({"name": name, "status": "complete",
                            "seconds": round(time.time() - at, 1)})
        except Exception as exc:                     # one bad molecule must not stop the panel
            results.append({"name": name, "status": "failed",
                            "seconds": round(time.time() - at, 1),
                            "error": f"{type(exc).__name__}: {exc}"[:300]})
        _set(job_id, rows=list(rows), results=list(results), done=index,
             records=dict(records),
             sequences={n: dict(c) for n, c in constructs.items()})
    _set(job_id, status="complete", current=None,
         seconds=round(time.time() - started, 1),
         finished_utc=datetime.now(timezone.utc).isoformat())


def _looks_like_csv(text: str) -> bool:
    """True only for something that really is a delimited table."""
    first = next((l for l in text.splitlines() if l.strip()), "").lstrip('\ufeff')
    if first.lstrip().startswith(">"):
        return False
    if "," not in first and "\t" not in first:
        return False
    header = [c.strip().lower() for c in first.replace("\t", ",").split(",")]
    return any(c in {"name", "id", "construct", "antibody"} for c in header) and len(header) >= 2


def start(text: str, options: dict) -> dict[str, Any]:
    """Parse a pasted panel and begin screening it in the background."""
    from .batch_input import parse_fasta, parse_csv

    body = (text or "").strip()
    if not body:
        raise ValueError("Paste a panel first.")
    # Pick the reader by what the text looks like, not by whether it happens
    # to begin with '>'. A pasted panel of "label line, sequence line" is not
    # CSV, and sending it to the CSV reader produced "need at least a name
    # column and a heavy column", which describes a file the user never had.
    parsed = parse_csv(body) if _looks_like_csv(body) else parse_fasta(body)
    if parsed.problems:
        raise ValueError(parsed.report() if hasattr(parsed, "report") else
                         "; ".join(str(p) for p in parsed.problems))
    constructs = parsed.constructs
    if not constructs:
        raise ValueError("No constructs were found in that panel.")
    if len(constructs) > MAX_CONSTRUCTS:
        raise ValueError(
            f"{len(constructs)} constructs exceeds the {MAX_CONSTRUCTS} this page runs. "
            "Larger panels belong in the command-line runner, which uses separate worker "
            "processes: python run_batch.py panel.fasta --out results/")

    job_id = uuid.uuid4().hex
    with _lock:
        _jobs[job_id] = {
            "id": job_id, "status": "running", "total": len(constructs), "done": 0,
            "current": None, "rows": [], "results": [],
            "warnings": list(getattr(parsed, "warnings", []) or []),
            "records": {}, "sequences": {n: dict(c) for n, c in constructs.items()},
            "options": {k: v for k, v in options.items()},
            "started_utc": datetime.now(timezone.utc).isoformat(),
        }
    threading.Thread(target=_run, args=(job_id, constructs, options), daemon=True).start()
    return get(job_id)


def to_csv(job_id: str) -> str:
    job = get(job_id)
    rows = (job or {}).get("rows") or []
    if not rows:
        return ""
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()
