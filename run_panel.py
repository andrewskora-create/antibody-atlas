#!/usr/bin/env python
"""Headless panel runner for Antibody Atlas.

Runs the same engines as the localhost web application, without the browser, so a
defined panel of antibodies can be analysed in one unattended pass and the results
written to files. Intended for comparing a parent antibody against a set of
engineered variants, and for proposing conservative substitutions at CDR PTM sites.

Usage
-----
    python run_panel.py --panel panel.json --out results/ [options]

The panel file is JSON, one entry per construct::

    {
      "baseline": "MAB01_WT",
      "constructs": {
        "MAB01_WT":   {"H": "QVQL...", "L": "SYEL..."},
        "MAB01_N55A": {"H": "QVQL...", "L": "SYEL...",
                       "preserve": {"H": [52], "L": [50, 91]},
                       "notes": "free text carried into the output"}
      }
    }

``preserve`` lists one-based input positions that should be kept if possible --
normally known paratope contacts. Candidates that edit a preserved position are
ranked last rather than hidden.

Options
-------
    --fold                 predict a local Fv structure and run the surface and
                           structural PTM screens (slow: minutes per construct)
    --antiberty            run the bundled AntiBERTy model
    --propermab            run PROPERMAB descriptors where the local platform
                           supports them (see --check)
    --substitutions        propose conservative substitutions at PTM sites
    --all-regions          propose for framework sites too, not only CDR sites
    --pdf                  write the three-page PDF report per construct
    --scheme {imgt,kabat}  numbering convention (default imgt)
    --ph FLOAT             pH for charge descriptors (default 7.4)
    --max-per-site N       substitution candidates reported per site (default 5)
    --check                report engine availability and exit

Every value written is produced by the platform's own modules. Nothing is
estimated or filled in when an engine is unavailable; missing results are
recorded as unavailable with a reason.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import platform
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from app.analysis import clean_sequence, sequence_metrics            # noqa: E402
from app.humanness import framework_identity, REFERENCE_PATH         # noqa: E402
from app.numbering import number_chain                               # noqa: E402
from app.ptm_precedents import sequence_precedents                   # noqa: E402
from app.clinical_frameworks import clinical_framework_precedents    # noqa: E402
from app.source_provenance import source_records                     # noqa: E402

RUNNER_VERSION = "panel-runner-1.0"


# --- engine availability -----------------------------------------------------
def engine_report() -> dict:
    """What this machine can actually run. Checked rather than assumed."""
    def have(module: str) -> bool:
        return importlib.util.find_spec(module) is not None

    hmmscan = shutil.which("hmmscan")
    nanoshaper = shutil.which("NanoShaper") or shutil.which("nanoshaper")
    report = {
        "platform": f"{platform.system()} {platform.machine()}",
        "python": platform.python_version(),
        "numbering_anarcii": have("anarcii"),
        "structure_immunebuilder": have("ImmuneBuilder"),
        "legacy_anarci": have("anarci"),
        "hmmscan_on_path": bool(hmmscan),
        "hmmscan_path": hmmscan,
        "sasa_freesasa": have("freesasa"),
        "antiberty": have("antiberty"),
        "antiberty_weights_bundled": _antiberty_weights() is not None,
        "reportlab_pdf": have("reportlab"),
        "propermab_vendored": (ROOT / "vendor" / "propermab").exists(),
        "nanoshaper_on_path": bool(nanoshaper),
        "apbs_on_path": bool(shutil.which("apbs")),
    }
    notes = []
    if report["structure_immunebuilder"] and not report["hmmscan_on_path"]:
        notes.append(
            "ImmuneBuilder numbers its input through legacy ANARCI, which runs the external "
            "hmmscan binary. Without hmmscan on PATH, folding will fail. Install HMMER "
            "(macOS: brew install hmmer; Debian/Ubuntu: apt install hmmer) or build it from "
            "github.com/EddyRivasLab/hmmer.")
    if report["propermab_vendored"] and not report["nanoshaper_on_path"]:
        notes.append(
            "PROPERMAB is vendored but no NanoShaper binary is on PATH, so its surface-patch "
            "features cannot be calculated here. The bundled binary in the development "
            "release is Apple ARM64 only. Set PROPERMAB_NANOSHAPER to a compatible binary, or "
            "run PROPERMAB on the machine where it was built. Patch features will be reported "
            "as unavailable rather than zero-filled.")
    if not report["antiberty_weights_bundled"]:
        notes.append("AntiBERTy weights did not resolve. They download on first use, or set "
                     "ANTIBERTY_WEIGHTS_DIR to a local checkpoint directory.")
    report["notes"] = notes
    return report


# --- core analysis (mirrors the /api/analyze path, without HTTP) --------------
def _digest(sequences: dict[str, str]) -> str:
    payload = "\n".join(f"{k}:{sequences[k]}" for k in sorted(sequences))
    return hashlib.sha256(payload.encode()).hexdigest()


def _domain_boundary(seq: str, numbering: dict) -> dict:
    """Residues in the input that fall outside the numbered variable domain.

    These matter because the sequence descriptors (pI, net charge, GRAVY) are
    computed over the whole supplied chain, not over the numbered domain. A single
    stray constant-region residue -- commonly the first Cys-region Arg left on the
    end of a patent 'variable chain coding sequence' -- shifts light-chain net
    charge by a full unit and pI by several tenths. Reported, never silently
    removed; use --trim-to-domain to drop it deliberately.
    """
    start, end = numbering["domain_start"], numbering["domain_end"]
    lead, trail = seq[:start - 1], seq[end:]
    return {
        "domain_start": start, "domain_end": end,
        "input_length": len(seq), "domain_length": end - start + 1,
        "leading_residues": lead, "trailing_residues": trail,
        "residues_outside_domain": len(lead) + len(trail),
        "warning": (
            f"{len(lead) + len(trail)} residue(s) lie outside the numbered variable domain "
            f"(leading '{lead}', trailing '{trail}'). Sequence descriptors are computed over "
            f"the whole supplied chain, so these affect pI, net charge and GRAVY."
        ) if (lead or trail) else None,
    }


def analyse(name: str, heavy: str, light: str, scheme: str, ph: float,
            trim_to_domain: bool = False) -> dict:
    result = {
        "construct": name, "scheme": scheme, "ph": ph, "chains": {},
        "methods": {
            "numbering": "ANARCII with matching IMGT or Kabat CDR boundaries",
            "humanness": "Closest functional human V-gene framework identity (OGRDB)",
            "biophysical": "Biopython ProtParam and disclosed sequence-window heuristic",
            "ptm": "Sequence motifs and susceptible residues; no occupancy or reaction-rate prediction",
        },
    }
    for label, raw in (("H", heavy), ("L", light)):
        if not (raw or "").strip():
            continue
        seq = clean_sequence(raw)
        numbering = number_chain(seq, scheme)
        actual = numbering["chain_type"]
        if label == "H" and actual != "H":
            raise ValueError(f"{name}: the heavy sequence was classified as a light chain.")
        if label == "L" and actual not in {"K", "L"}:
            raise ValueError(f"{name}: the light sequence was classified as a heavy chain.")
        boundary = _domain_boundary(seq, numbering)
        if trim_to_domain and boundary["residues_outside_domain"]:
            seq = seq[boundary["domain_start"] - 1:boundary["domain_end"]]
            numbering = number_chain(seq, scheme)
            boundary["trimmed"] = True
            boundary["trimmed_to"] = len(seq)
        imgt = numbering if scheme == "imgt" else number_chain(seq, "imgt")
        result["chains"][label] = {
            "sequence": seq, "numbering": numbering, "numbering_imgt": imgt,
            "domain_boundary": boundary,
            "humanness": framework_identity(seq, imgt),
            "metrics": sequence_metrics(seq, ph, imgt),
        }
    if not result["chains"]:
        raise ValueError(f"{name}: no usable heavy or light sequence.")
    result["ptm_precedents"] = sequence_precedents(result["chains"], scheme)
    result["clinical_framework_precedents"] = clinical_framework_precedents(result["chains"], scheme)
    result["sources"] = source_records()
    sequences = {k: v["sequence"] for k, v in result["chains"].items()}
    result["provenance"] = {
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": _digest(sequences),
        "reference_sha256": hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest(),
        "reference_file": REFERENCE_PATH.name,
        "runner_version": RUNNER_VERSION,
        "variant_lineage": "Not supplied",
        "assay_provenance": "No measured assays supplied",
    }
    return result


def fold_and_screen(record: dict, out_dir: Path, name: str) -> dict:
    """Predict a local Fv model and run the surface/structural PTM screens."""
    from app.structure_analysis import analyze_structure, extract_trace
    sequences = {k: v["sequence"] for k, v in record["chains"].items()}
    status = {"status": "pending", "method": None, "pdb": None}
    try:
        if "L" not in sequences:
            from ImmuneBuilder import NanoBodyBuilder2
            predictor = NanoBodyBuilder2(model_ids=[1])
            status["method"] = "ImmuneBuilder NanoBodyBuilder2 (single model)"
        else:
            from ImmuneBuilder import ABodyBuilder2
            predictor = ABodyBuilder2(model_ids=[1])
            status["method"] = "ImmuneBuilder ABodyBuilder2 (single model)"
        model = predictor.predict(sequences)
        pdb_path = out_dir / f"{name}.pdb"
        model.save(str(pdb_path))
        if not pdb_path.exists() or pdb_path.stat().st_size == 0:
            raise RuntimeError("The model did not produce a PDB file.")
        status["pdb"] = str(pdb_path)
        record["provenance"]["structure_input_sha256"] = record["provenance"]["input_sha256"]
        record["provenance"]["structure_pdb_sha256"] = hashlib.sha256(pdb_path.read_bytes()).hexdigest()
        record["provenance"]["structure_model"] = status["method"]
        try:
            from importlib.metadata import version
            record["provenance"]["structure_model_version"] = version("ImmuneBuilder")
        except Exception:
            record["provenance"]["structure_model_version"] = "unknown"
        try:
            record["structure_analysis"] = analyze_structure(
                pdb_path, sequences, record["ph"], record["scheme"])
            status["status"] = "complete"
        except Exception as exc:
            record["structure_analysis"] = None
            status["status"] = "structure_screen_failed"
            status["error"] = f"Surface/PTM screen failed: {str(exc)[:300]}"
            try:
                record["structure_trace"] = extract_trace(pdb_path)
            except Exception:
                record["structure_trace"] = None
    except Exception as exc:
        status["status"] = "failed"
        status["error"] = str(exc)[:400]
    record["fold"] = status
    return record


def _antiberty_weights():
    """Where AntiBERTy's checkpoint resolves, or None.

    AntiBERTy resolves its own checkpoint: an explicit path, then
    ANTIBERTY_WEIGHTS_DIR, then a directory inside the installed package,
    then a Hugging Face download. download=False keeps this off the network,
    so a doctor report never triggers a 104 MB fetch.
    """
    try:
        from antiberty import get_weights
        return Path(get_weights(download=False))
    except Exception:
        return None


def run_antiberty(record: dict) -> dict:
    """AntiBERTy over the local checkpoint, without a model-hub call."""
    sequences = {k: v["sequence"] for k, v in record["chains"].items()}
    weights = _antiberty_weights()
    try:
        if weights is None:
            raise RuntimeError(
                "AntiBERTy weights did not resolve. They download on first use, or set "
                "ANTIBERTY_WEIGHTS_DIR to a local checkpoint directory.")
        if any(not 90 <= len(seq) <= 150 for seq in sequences.values()):
            raise ValueError("AntiBERTy expects 90-150 amino acids per variable domain.")
        from antiberty import AntiBERTyRunner
        from app.antiberty_analysis import analyze_antiberty
        from app.profiler_interpretation import annotate_antiberty
        runner = AntiBERTyRunner(device="cpu", checkpoint_path=str(weights))
        result = analyze_antiberty(runner, sequences, scheme=record["scheme"], batch_size=16)
        chains = {label: {"numbering": record["chains"][label]["numbering"]}
                  for label in sequences}
        result = annotate_antiberty(result, chains=chains)
        for label, chain in record["chains"].items():
            per_chain = (result or {}).get("chains", {}).get(label)
            if per_chain:
                chain["antiberty"] = per_chain
        record["antiberty"] = result
    except Exception as exc:
        record["antiberty"] = {"status": "unavailable", "reason": str(exc)[:300]}
    return record


def substitutions(record: dict, preserve: dict, cdr_only: bool, max_per_site: int) -> dict:
    from app.cdr_substitutions import propose, variant_specs
    preserve_sets = {k: set(v) for k, v in (preserve or {}).items()}
    proposals = propose(record["chains"], cdr_only=cdr_only, max_per_site=max_per_site,
                        structure=record.get("structure_analysis"), preserve=preserve_sets)
    proposals["variant_specs"] = variant_specs(proposals, top_n=1,
                                               name_prefix=record["construct"])
    record["cdr_substitutions"] = proposals
    return record


# --- outputs -----------------------------------------------------------------
def write_pdf(record: dict, out_dir: Path, name: str) -> str | None:
    try:
        from app.report import make_report
        from app.assessment import build_assessment
        payload = copy.deepcopy(record)
        payload["assessment"] = build_assessment(payload)
        path = out_dir / f"{name}.pdf"
        # The interactive app exports its NGL view; a headless run has no
        # browser, so render the same figure offline from the coordinates.
        image = None
        pdb_path = out_dir / f"{name}.pdb"
        if pdb_path.is_file():
            try:
                from app.structure_render import render_structure_png
                image = render_structure_png(pdb_path, record.get("structure_analysis"))
            except Exception as exc:
                print(f"  ! structure figure skipped for {name}: {str(exc)[:120]}", file=sys.stderr)
        data = make_report(payload, title=f"Developability report - {name}",
                           structure_image=image)
        path.write_bytes(data if isinstance(data, (bytes, bytearray)) else Path(data).read_bytes())
        return str(path)
    except Exception as exc:
        print(f"  ! PDF generation failed for {name}: {str(exc)[:200]}", file=sys.stderr)
        return None


def summary_rows(records: dict) -> list[dict]:
    rows = []
    for name, record in records.items():
        for label, chain in record["chains"].items():
            metrics = chain["metrics"]
            sites = metrics["ptm_sites"]
            numbering = chain["numbering"]
            regions = {r["sequence_position"]: r.get("region") for r in numbering["residues"]}
            cdr_sites = [s for s in sites if str(regions.get(s["start"], "")).startswith("CDR")]
            rows.append({
                "construct": name, "chain": label,
                "chain_type": numbering["chain_type"], "length": metrics["length"],
                "pI": metrics["isoelectric_point"],
                "charge_at_pH": metrics["net_charge_at_ph"],
                "gravy": metrics["gravy"],
                "residues_outside_v_domain": (chain.get("domain_boundary") or {}).get("residues_outside_domain"),
                "outside_domain_detail": ((chain.get("domain_boundary") or {}).get("trailing_residues") or "")
                                          + ("|lead:" + (chain.get("domain_boundary") or {}).get("leading_residues", "")
                                             if (chain.get("domain_boundary") or {}).get("leading_residues") else ""),
                "framework_identity_percent": (chain.get("humanness") or {}).get("percent"),
                "closest_v_allele": (chain.get("humanness") or {}).get("name"),
                "hydrophobic_patch_windows": metrics["aggregation_proxy"]["patch_count"],
                "ptm_sites_total": len(sites),
                "ptm_sites_in_cdr": len(cdr_sites),
                "cdr_deamidation": sum("deamidation" in s["kind"].lower() for s in cdr_sites),
                "cdr_isomerization": sum("isomerization" in s["kind"].lower() for s in cdr_sites),
                "cdr_oxidation": sum("oxidation" in s["kind"].lower() for s in cdr_sites),
                "cdr_glycosylation": sum("glycosylation" in s["kind"].lower() for s in cdr_sites),
                "cdr_site_list": "; ".join(
                    f"{label}{s['start']}{s['motif']}" for s in cdr_sites),
            })
    return rows


def substitution_rows(records: dict) -> list[dict]:
    rows = []
    for name, record in records.items():
        proposals = record.get("cdr_substitutions") or {}
        for chain in proposals.get("chains", {}).values():
            for site in chain["sites"]:
                for candidate in site["proposals"]:
                    germline = candidate["germline_precedent"]
                    clinical = candidate["clinical_precedent"]
                    rows.append({
                        "construct": name,
                        "chain": candidate["chain"],
                        "site_position": site["position"],
                        "site_imgt": site["imgt_number"],
                        "region": site["region"],
                        "chemistry": site["kind"],
                        "motif": site["motif"],
                        "exposure_relative_sasa": site["exposure"].get("relative_sidechain_sasa"),
                        "rank": candidate["rank"],
                        "mutation": candidate["mutation"],
                        "variant_spec": candidate["variant_spec"],
                        "motif_after": candidate["motif_after"],
                        "strategy": candidate["strategy"],
                        "conventional_choice": candidate["standard_choice"],
                        "conventional_note": candidate["standard_choice_note"] or "",
                        "edits_preserved_position": candidate["edits_preserved_position"],
                        "introduces_new_liability": not candidate["clean"],
                        "new_liabilities": "; ".join(
                            f"{i['kind']}@{i['position']}({i['motif']})" for i in candidate["introduces"]),
                        "removes_chemistries": "; ".join(candidate["removes_chemistries"]),
                        "blosum62": candidate["blosum62"],
                        "conservativeness": candidate["conservativeness"]["score"],
                        "physchem_distance": candidate["physicochemical"]["distance"],
                        "hydropathy_delta": candidate["physicochemical"]["hydropathy_delta"],
                        "charge_delta": candidate["physicochemical"]["charge_delta"],
                        "germline_percent": germline.get("percent"),
                        "germline_status": germline.get("status"),
                        "clinical_percent": clinical.get("percent"),
                        "clinical_approved_examples": "; ".join(clinical.get("approved_examples") or []),
                    })
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def compare_to_baseline(records: dict, baseline: str) -> dict:
    """Variant comparison against the named baseline, using the platform module."""
    from app.variant_comparison import compare_records
    if baseline not in records:
        return {"status": "unavailable", "reason": f"Baseline '{baseline}' is not in the panel."}
    out = {"baseline": baseline, "comparisons": {}}
    for name, record in records.items():
        if name == baseline:
            continue
        try:
            out["comparisons"][name] = compare_records(records[baseline], record)
        except Exception as exc:
            out["comparisons"][name] = {"status": "unavailable", "reason": str(exc)[:300]}
    return out


# --- entry point -------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel")
    parser.add_argument("--out", default="panel_results")
    parser.add_argument("--scheme", choices=["imgt", "kabat"], default="imgt")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--fold", action="store_true")
    parser.add_argument("--antiberty", action="store_true")
    parser.add_argument("--substitutions", action="store_true")
    parser.add_argument("--all-regions", action="store_true")
    parser.add_argument("--pdf", action="store_true")
    parser.add_argument("--max-per-site", type=int, default=5)
    parser.add_argument("--trim-to-domain", action="store_true",
                        help="drop residues outside the numbered variable domain before "
                             "computing descriptors (reported either way)")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        report = engine_report()
        print(json.dumps(report, indent=2))
        for note in report["notes"]:
            print(f"\nNOTE: {note}", file=sys.stderr)
        return 0

    if not args.panel:
        parser.error("--panel is required unless --check is given.")

    panel = json.loads(Path(args.panel).read_text())
    constructs = panel["constructs"]
    baseline = panel.get("baseline")
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    engines = engine_report()
    (out_dir / "engine_report.json").write_text(json.dumps(engines, indent=2))
    for note in engines["notes"]:
        print(f"NOTE: {note}\n", file=sys.stderr)

    records: dict[str, dict] = {}
    failures: dict[str, str] = {}
    for name, spec in constructs.items():
        print(f"[{name}] analysing...", flush=True)
        try:
            record = analyse(name, spec.get("H", ""), spec.get("L", ""), args.scheme, args.ph,
                             trim_to_domain=args.trim_to_domain)
            for label, chain in record["chains"].items():
                warning = (chain.get("domain_boundary") or {}).get("warning")
                if warning:
                    print(f"  ! {name} chain {label}: {warning}", file=sys.stderr)
            record["notes"] = spec.get("notes")
            if args.fold:
                print(f"[{name}] folding...", flush=True)
                record = fold_and_screen(record, out_dir, name)
                print(f"[{name}] fold: {record['fold']['status']}", flush=True)
            if args.antiberty:
                print(f"[{name}] antiberty...", flush=True)
                record = run_antiberty(record)
            if args.substitutions:
                record = substitutions(record, spec.get("preserve"),
                                       not args.all_regions, args.max_per_site)
            if args.pdf:
                record["pdf"] = write_pdf(record, out_dir, name)
            (out_dir / f"{name}.json").write_text(json.dumps(record, indent=2, default=str))
            records[name] = record
        except Exception as exc:
            failures[name] = f"{type(exc).__name__}: {exc}"
            print(f"[{name}] FAILED: {exc}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

    write_csv(out_dir / "panel_summary.csv", summary_rows(records))
    if args.substitutions:
        write_csv(out_dir / "substitution_candidates.csv", substitution_rows(records))
        specs = [line for r in records.values()
                 for line in (r.get("cdr_substitutions") or {}).get("variant_specs", [])]
        (out_dir / "variant_specs.txt").write_text("\n".join(specs) + ("\n" if specs else ""))
    if baseline:
        (out_dir / "comparison_to_baseline.json").write_text(
            json.dumps(compare_to_baseline(records, baseline), indent=2, default=str))

    manifest = {
        "runner_version": RUNNER_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "panel_file": str(Path(args.panel).resolve()),
        "scheme": args.scheme, "ph": args.ph,
        "options": {k: getattr(args, k) for k in
                    ("fold", "antiberty", "substitutions", "all_regions", "pdf",
                     "max_per_site", "trim_to_domain")},
        "constructs_analysed": sorted(records),
        "constructs_failed": failures,
        "engines": engines,
    }
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\n{len(records)} analysed, {len(failures)} failed -> {out_dir}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
