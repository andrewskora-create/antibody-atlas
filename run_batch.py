#!/usr/bin/env python
"""Screen a panel of antibodies from one multi-FASTA file.

    python run_batch.py panel.fasta --out results/
    python run_batch.py panel.fasta --out results/ --workers 4 --fold
    python run_batch.py panel.fasta --out results/ --resume        # default

Built for 20-100 molecules, where a run takes hours:

  * Input is validated in full before any work starts, so a typo on line 40
    fails in seconds rather than an hour in.
  * Each molecule is independent. One failure is recorded and the run
    continues; it never takes the batch down.
  * Progress is written to disk as it happens, so re-running skips what is
    already done. A crash at molecule 80 does not restart from 1.
  * Workers run in separate processes with their thread pools capped, because
    torch and OpenMM each try to use every core and together they thrash.

Outputs an HTML dashboard, a one-row-per-molecule CSV for ranking, and per
molecule a one-page brief and the four-page report.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import signal
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

STATE_NAME = "batch_state.json"


def input_options_fingerprint(chains: dict, options: dict) -> str:
    """Resume only the same sequences and scientific/output configuration."""
    payload = {"chains": chains, "options": {
        key: value for key, value in options.items()
        if key not in {"out", "workers", "threads", "resume"}}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def _outputs_exist(out_dir: Path, name: str, options: dict) -> bool:
    suffixes = ['.json']
    if options.get('fold'):
        suffixes.append('.pdb')
    if options.get('pdf'):
        suffixes.extend(['_report.pdf', '_brief.pdf'])
    return all((out_dir / f'{name}{suffix}').is_file() for suffix in suffixes)


# ----------------------------------------------------------------- worker
def _screen_one(args: tuple) -> dict:
    """Run one molecule in its own process. Returns a status record."""
    name, chains, options = args
    # Cap per-worker threading: several workers each grabbing every core is
    # slower than one, and on a laptop it stalls the machine.
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                     "NUMEXPR_NUM_THREADS", "TORCH_NUM_THREADS"):
        os.environ.setdefault(variable, str(options["threads"]))
    sys.path.insert(0, str(ROOT))
    started = time.time()
    out_dir = Path(options["out"])
    try:
        import torch
        torch.set_num_threads(options["threads"])
    except Exception:
        pass

    try:
        from app.analysis import sequence_metrics
        from app.numbering import number_chain
        from app.humanness import framework_identity
        from app.clinical_frameworks import clinical_framework_precedents
        from app.ptm_precedents import sequence_precedents
        from app.sequence_panel import annotate_chain, conclusion as panel_conclusion

        scheme, ph = options["scheme"], options["ph"]
        record: dict = {"scheme": scheme, "ph": ph, "chains": {}}
        for label, sequence in chains.items():
            canonical = number_chain(sequence, "imgt")
            selected = canonical if scheme == "imgt" else number_chain(sequence, scheme)
            metrics = sequence_metrics(sequence, ph, canonical)
            entry = {"sequence": sequence, "numbering": selected, "numbering_imgt": canonical,
                     "metrics": metrics, "humanness": framework_identity(sequence, canonical)}
            entry["sequence_bands"] = annotate_chain(metrics, canonical.get("chain_type") or label)
            record["chains"][label] = entry
        record["sequence_conclusion"] = panel_conclusion(
            {k: v["sequence_bands"] for k, v in record["chains"].items()})
        record["provenance"] = {
            "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "input_sha256": __import__("hashlib").sha256(
                "\n".join(f"{c}:{s}" for c, s in sorted(chains.items())).encode()).hexdigest(),
        }

        structure_image = None
        if options["fold"]:
            from ImmuneBuilder import ABodyBuilder2
            from app.structure_analysis import analyze_structure
            from app.structure_render import render_structure_png
            pdb_path = out_dir / f"{name}.pdb"
            ABodyBuilder2().predict(
                {("H" if k == "H" else "L"): v for k, v in chains.items()}).save(str(pdb_path))
            record["provenance"]["structure_pdb_sha256"] = __import__("hashlib").sha256(
                pdb_path.read_bytes()).hexdigest()
            record["structure_analysis"] = analyze_structure(
                pdb_path, {k: v["sequence"] for k, v in record["chains"].items()}, ph, scheme)
            structure_image = render_structure_png(pdb_path, record["structure_analysis"])
            if options["paratope"]:
                from app.paratope import predict, imgt_position_map
                record["paratope"] = predict(
                    pdb_path, heavy_chain_id="H",
                    light_chain_id="L" if "L" in chains else None,
                    imgt_to_position=imgt_position_map(record["chains"]))

        # Precedent corpora are not redistributed, so a fresh clone may not have
        # them. They are optional: the table, the flags and the liabilities do
        # not use them, so a batch run degrades rather than failing every row.
        from app.reference_data import ReferenceDataMissing
        for key, fn in (("ptm_precedents", sequence_precedents),
                        ("clinical_framework_precedents", clinical_framework_precedents)):
            if not options.get("precedents", True):
                record[key] = {"cases": [], "sites": [], "unavailable": True,
                               "reason": "Skipped: --no-precedents."}
                continue
            try:
                record[key] = fn(record["chains"], scheme)
            except ReferenceDataMissing as exc:
                record[key] = {"cases": [], "sites": [], "unavailable": True, "reason": str(exc)}

        # The consolidated per-residue table, with whatever layers this tier ran.
        from app.residue_table import build_residue_table
        record["residue_table"] = build_residue_table(
            record["chains"], structure=record.get("structure_analysis"),
            paratope=record.get("paratope"),
            precedents=record.get("clinical_framework_precedents"))

        # Molecule-level charge and hydrophobicity, on a full-length basis when
        # an isotype was requested.
        from app.molecule_flags import molecule_flags
        assembly = None
        if options.get("isotype"):
            from app import constant_regions
            assembly = constant_regions.assemble(
                record["chains"], heavy_isotype=options["isotype"],
                light_constant=options.get("light_constant"))
        record["molecule_flags"] = molecule_flags(
            record["chains"], formulation_ph=options.get("formulation_ph", 6.0),
            assay_ph=ph, assembly=assembly)

        if options["substitutions"]:
            from app.cdr_substitutions import propose
            record["substitutions"] = propose(
                record["chains"], structure=record.get("structure_analysis"),
                paratope=record.get("paratope"), ph=ph, max_per_site=options["max_per_site"])

        from app.assessment import build_assessment
        record["assessment"] = build_assessment(record)

        (out_dir / f"{name}.json").write_text(json.dumps(record, indent=2, default=str))
        if options["pdf"]:
            from app.report import make_report
            from app.report_brief import make_brief
            (out_dir / f"{name}_report.pdf").write_bytes(make_report(
                record, title=f"Developability report - {name}", structure_image=structure_image))
            (out_dir / f"{name}_brief.pdf").write_bytes(make_brief(
                record, title=name, detail_link=f"{name}_report.pdf", structure_image=structure_image))

        return {"name": name, "status": "complete", "seconds": round(time.time() - started, 1),
                "row": _row(name, record)}
    except Exception as exc:
        return {"name": name, "status": "failed", "seconds": round(time.time() - started, 1),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()[-1500:]}


def _row(name: str, record: dict) -> dict:
    """One row per molecule, for ranking."""
    assessment = record.get("assessment") or {}
    sites = assessment.get("sites", [])
    cdr = [s for s in sites if str(s.get("region", "")).startswith("CDR")]
    patches = {p["kind"]: p for p in (assessment.get("patch_risk") or {}).get("categories", [])}
    row = {"name": name,
           "cdr_sites_high": sum(s.get("priority") == "high" for s in cdr),
           "cdr_sites_review": sum(s.get("priority") == "review" for s in cdr),
           "framework_sites": len(sites) - len(cdr),
           "sequence_conclusion": (record.get("sequence_conclusion") or {}).get("state"),
           }
    for label, chain in record.get("chains", {}).items():
        metrics = chain.get("metrics", {})
        row[f"pI_{label}"] = metrics.get("isoelectric_point")
        row[f"charge_{label}"] = metrics.get("net_charge_at_ph")
        row[f"gravy_{label}"] = metrics.get("gravy")
        row[f"hydrophobic_windows_{label}"] = (metrics.get("aggregation_proxy") or {}).get("patch_count")
        row[f"humanness_{label}"] = (chain.get("humanness") or {}).get("percent")
    for kind in ("hydrophobic", "positive", "negative"):
        entry = patches.get(kind) or {}
        row[f"{kind}_patch_priority"] = entry.get("priority")
        row[f"{kind}_patch_largest_a2"] = entry.get("largest_area_a2")
    # The per-residue table, collapsed to the counts worth ranking on.
    table = record.get("residue_table") or {}
    verdicts = table.get("verdict_counts") or {}
    row["residues_notable"] = table.get("notable_rows")
    for name_ in ("preserve", "fix", "watch"):
        # Retain the existing CSV header for consumers of earlier releases.
        row[f"residues_{name_}"] = verdicts.get("evaluate" if name_ == "fix" else name_)
    row["germline_framework_mismatches"] = (table.get("germline_counts") or {}).get(
        "mismatch_framework")
    row["layers_missing"] = ";".join(table.get("missing_layers") or []) or None

    # Molecule-level flags: the worst verdict, plus the value behind each one.
    flags = record.get("molecule_flags") or {}
    row["molecule_basis"] = flags.get("basis")
    row["molecule_worst_flag"] = flags.get("worst_verdict")
    row["molecule_pI"] = (flags.get("fv") or {}).get("pi")
    row["molecule_net_charge"] = (flags.get("fv") or {}).get("net_charge_assay_ph")
    for flag in flags.get("flags", []):
        row[f"flag_{flag['id']}"] = flag.get("verdict")
        row[f"value_{flag['id']}"] = flag.get("value")
    row["tradeoffs"] = len(flags.get("tradeoffs") or []) or None

    paratope = (record.get("paratope") or {}).get("summary") or {}
    row["predicted_contacts"] = paratope.get("likely_contact_count")
    comparators = [s.get("comparators") or {} for s in sites
                   if not str(s.get("region", "")).startswith("CDR")]
    row["framework_approved_precedent"] = sum(c.get("approved", 0) for c in comparators)
    return row


# --------------------------------------------------------------- the queue
def load_state(out_dir: Path) -> dict:
    try:
        return json.loads((out_dir / STATE_NAME).read_text())
    except (OSError, ValueError):
        return {"molecules": {}}


def save_state(out_dir: Path, state: dict) -> None:
    # Written after every completion so a crash loses at most one molecule.
    temporary = out_dir / (STATE_NAME + ".tmp")
    temporary.write_text(json.dumps(state, indent=2))
    temporary.replace(out_dir / STATE_NAME)


def run_batch(panel, out_dir: Path, options: dict) -> dict:
    from app.batch_input import validate_name
    if not panel.ok:
        raise ValueError('Fix the panel input problems before screening.')
    for name in panel.constructs:
        validate_name(name)
    out_dir.mkdir(parents=True, exist_ok=True)
    state = load_state(out_dir) if options["resume"] else {"molecules": {}}
    state.setdefault("molecules", {})
    state["started_at_utc"] = state.get("started_at_utc") or datetime.now(timezone.utc).isoformat()
    state["options"] = {k: v for k, v in options.items() if k != "out"}
    fingerprints = {name: input_options_fingerprint(chains, options)
                    for name, chains in panel.constructs.items()}
    # Old records without a fingerprint are rerun once. Removed molecules and
    # mismatched configurations must not remain in the new dashboard/CSV.
    state['molecules'] = {name: result for name, result in state['molecules'].items()
                          if name in fingerprints and
                          result.get('input_options_sha256') == fingerprints[name]}
    done = {n for n, m in state["molecules"].items()
            if m.get("status") == "complete" and _outputs_exist(out_dir, n, options)}
    pending = [(n, c) for n, c in panel.constructs.items() if n not in done]
    if done:
        print(f"  resuming: {len(done)} already complete, {len(pending)} to run")

    if not pending:
        save_state(out_dir, state)
        return state

    jobs = [(name, chains, options) for name, chains in pending]
    workers = max(1, min(options["workers"], len(jobs)))
    print(f"  {len(jobs)} molecule(s) on {workers} worker(s)\n")

    completed = 0
    started_at = time.time()
    interrupted = False

    def on_interrupt(*_):
        nonlocal interrupted
        interrupted = True
        print("\n  interrupt received; finishing in-flight molecules then stopping. "
              "Re-run the same command to continue.")
    signal.signal(signal.SIGINT, on_interrupt)

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_screen_one, job): job[0] for job in jobs}
        for future in as_completed(futures):
            name = futures[future]
            try:
                result = future.result()
            except Exception as exc:   # a worker died outright
                result = {"name": name, "status": "failed",
                          "error": f"worker process failed: {type(exc).__name__}: {exc}"}
            completed += 1
            result['input_options_sha256'] = fingerprints[name]
            state["molecules"][name] = result
            save_state(out_dir, state)
            elapsed = time.time() - started_at
            rate = elapsed / completed
            remaining = (len(jobs) - completed) * rate
            mark = "ok  " if result["status"] == "complete" else "FAIL"
            detail = f"{result.get('seconds','?')}s" if result["status"] == "complete" \
                else str(result.get("error", ""))[:70]
            print(f"  [{completed}/{len(jobs)}] {mark} {name:<28} {detail}"
                  + (f"   ~{remaining/60:.0f} min left" if remaining > 60 and completed < len(jobs) else ""))
            if interrupted:
                for pending_future in futures:
                    pending_future.cancel()
                break

    state["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    save_state(out_dir, state)
    return state


# ------------------------------------------------------------- the outputs
def write_summary(out_dir: Path, state: dict) -> Path:
    rows = [m["row"] for m in state["molecules"].values()
            if m.get("status") == "complete" and m.get("row")]
    path = out_dir / "batch_summary.csv"
    if not rows:
        path.write_text("name\n")
        return path
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_dashboard(out_dir: Path, state: dict, panel) -> Path:
    from app.batch_dashboard import render
    path = out_dir / "index.html"
    path.write_text(render(state, panel))
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="multi-FASTA (or CSV) of the panel")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) // 3))
    parser.add_argument("--threads", type=int, default=2, help="threads per worker")
    parser.add_argument("--scheme", choices=("imgt", "kabat"), default="imgt")
    parser.add_argument("--ph", type=float, default=7.4,
                        help="pH for the net-charge descriptor")
    parser.add_argument("--formulation-ph", type=float, default=6.0,
                        help="pH the colloidal margin is judged at (default 6.0)")
    parser.add_argument("--isotype", choices=("IGHG1", "IGHG2", "IGHG3", "IGHG4"),
                        help="attach wild-type secreted constant regions, so pI and net charge "
                             "are computed on the full-length H2L2 rather than the Fv")
    parser.add_argument("--light-constant",
                        choices=("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7"),
                        help="light constant region (default follows the called chain type)")
    parser.add_argument("--fold", action="store_true", help="predict structure (slow; unlocks exposure, patches, paratope)")
    parser.add_argument("--no-paratope", action="store_true")
    parser.add_argument("--no-substitutions", action="store_true")
    parser.add_argument("--no-precedents", action="store_true",
                        help="skip the documented-PTM and clinical comparator lookups. These "
                             "dominate the sequence-only tier (about 5 s of the ~6 s warm cost "
                             "per molecule) and are the part a large triage run least needs.")
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--max-per-site", type=int, default=4)
    parser.add_argument("--no-resume", action="store_true", help="ignore previous progress and start over")
    parser.add_argument("--validate-only", action="store_true", help="check the input and stop")
    args = parser.parse_args()

    from app.batch_input import load
    print("Reading panel")
    panel = load(args.input)
    print(" ", panel.report().replace("\n", "\n  "))
    if not panel.ok:
        print("\nNothing was run. Fix the problems above and try again.")
        return 2
    if args.validate_only:
        print("\nInput is usable.")
        return 0

    out_dir = Path(args.out)
    options = {"out": str(out_dir), "scheme": args.scheme, "ph": args.ph, "fold": args.fold,
               "paratope": args.fold and not args.no_paratope,
               "substitutions": not args.no_substitutions, "pdf": not args.no_pdf,
               "precedents": not args.no_precedents,
               "max_per_site": args.max_per_site, "threads": max(1, args.threads),
               "workers": max(1, args.workers), "resume": not args.no_resume,
               "formulation_ph": args.formulation_ph, "isotype": args.isotype,
               "light_constant": args.light_constant}

    tier = "structure" if args.fold else "sequence-only"
    basis = (f"full-length {args.isotype}" if args.isotype else "variable domains only")
    print(f"\nScreening into {out_dir}")
    print(f"  tier: {tier}  |  charge basis: {basis}  |  workers: {max(1, args.workers)}")
    if not args.fold:
        print("  sequence-only: numbering, liabilities, germline, the per-residue table and the")
        print("  molecule flags. Exposure, surface patches and paratope need --fold.")
        print("  Each worker pays a one-off start-up of roughly a minute (imports, numbering")
        print("  model, reference corpora) before its first molecule, so throughput only")
        print("  settles after a few molecules per worker. Add --no-precedents for triage runs:")
        print("  the precedent lookups are most of the steady-state cost per molecule.")
    state = run_batch(panel, out_dir, options)

    summary = write_summary(out_dir, state)
    dashboard = write_dashboard(out_dir, state, panel)
    complete = sum(m.get("status") == "complete" for m in state["molecules"].values())
    failed = [m for m in state["molecules"].values() if m.get("status") != "complete"]
    print(f"\n{complete} complete, {len(failed)} failed.")
    for entry in failed[:10]:
        print(f"  FAILED {entry['name']}: {str(entry.get('error',''))[:110]}")
    if len(failed) > 10:
        print(f"  ... and {len(failed)-10} more; see {STATE_NAME}")
    print(f"\n  dashboard  {dashboard}")
    print(f"  ranking    {summary}")
    if failed:
        print("\nRe-run the same command to retry the failures; completed molecules are skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
