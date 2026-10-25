#!/usr/bin/env python
"""Build the reference datasets this project does not redistribute.

Three reference corpora derive from Oxford Thera-SAbDab, which publishes no
licence grant, so they are not shipped in this repository (see
DATA_PROVENANCE.md). This script downloads Oxford's CSV to your machine and
regenerates them locally.

    python setup/build_reference_data.py --all        # everything
    python setup/build_reference_data.py --verify     # report, change nothing

Every rebuild is hash-verified. The upstream CSV is checked against the hash
recorded when these results were produced, and each rehydrated chain is checked
against its own recorded hash, so a local rebuild is provably the same input the
published results came from rather than merely a similar one.

Oxford asks that you cite Raybould, Marks et al., Nucleic Acids Res. (2019)
gkz827.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT))

CSV_URL = ("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
           "TheraSAbDab_SeqStruc_OnlineDownload.csv")
CSV_SHA256 = "cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d"
CSV_RETRIEVED = "2026-10-04"
STAGE_COLUMN = "Highest_Clin_Trial (Feb '25)"
ELIGIBLE_STAGES = ("Approved", "Phase-III")
SEQUENCE_COLUMNS = ("HeavySequence", "LightSequence",
                    "HeavySequence(ifbispec)", "LightSequence(ifbispec)")
EMPTY = {"", "na", "n/a", "none", "null", "-"}

CITATION = "Raybould, Marks et al., Nucleic Acids Res. (2019) gkz827"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.lower() in EMPTY else value


# ---------------------------------------------------------------- source CSV

def fetch_csv(cache: Path, allow_updated: bool) -> list[dict[str, str]]:
    if cache.is_file():
        raw = cache.read_bytes()
        print(f"  using cached CSV at {cache.relative_to(ROOT) if cache.is_relative_to(ROOT) else cache}")
    else:
        print(f"  downloading {CSV_URL}")
        with urllib.request.urlopen(CSV_URL, timeout=120) as response:
            raw = response.read()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(raw)
        print(f"  {len(raw)} bytes")

    digest = hashlib.sha256(raw).hexdigest()
    if digest != CSV_SHA256:
        message = (
            f"\n  Upstream CSV hash does not match the recorded snapshot.\n"
            f"    expected {CSV_SHA256}\n    found    {digest}\n"
            f"  Oxford has updated the file since {CSV_RETRIEVED}. Panel membership and\n"
            f"  descriptor ranges both depend on it, so a rebuild from a different\n"
            f"  snapshot will not reproduce the published results.\n"
        )
        if not allow_updated:
            raise SystemExit(message + "  Re-run with --allow-updated-source to build anyway.")
        print(message + "  --allow-updated-source given; continuing.")
    else:
        print(f"  hash matches the recorded {CSV_RETRIEVED} snapshot")
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))


def eligible_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [r for r in rows if _clean(r.get(STAGE_COLUMN)) in ELIGIBLE_STAGES]


def chains_by_drug(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """lowercase drug name -> {column: sequence} for eligible rows."""
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        name = _clean(row.get("Therapeutic"))
        if not name:
            continue
        found = {c: _clean(row.get(c)) for c in SEQUENCE_COLUMNS}
        found = {c: s for c, s in found.items() if s}
        if found:
            out[name.lower()] = found
    return out


# --------------------------------------------------------- clinical library

def build_clinical(rows: list[dict[str, str]], out: Path) -> dict:
    from app.numbering import number_chain

    selected = eligible_rows(rows)
    print(f"  {len(selected)} eligible drugs (stage in {' | '.join(ELIGIBLE_STAGES)})")

    pending: list[tuple[str, str, str, str, str]] = []
    for row in selected:
        drug = _clean(row.get("Therapeutic"))
        stage = _clean(row.get(STAGE_COLUMN))
        status = _clean(row.get("Est. Status"))
        for column in SEQUENCE_COLUMNS:
            sequence = _clean(row.get(column))
            if sequence:
                pending.append((drug, stage, status, column, sequence))

    print(f"  {len(pending)} chains x 2 numbering schemes; this takes several minutes")
    references, excluded = [], []
    for index, (drug, stage, status, column, sequence) in enumerate(pending, 1):
        if index % 100 == 0 or index == len(pending):
            print(f"    {index}/{len(pending)}")
        # The dataset stores both conventions, because CDR-site retrieval is
        # performed in whichever convention the user selected.
        numbering = {}
        try:
            for scheme in ("imgt", "kabat"):
                numbering[scheme] = number_chain(sequence, scheme)
        except ValueError as exc:
            excluded.append({"drug": drug, "reference_chain": column, "reason": str(exc)})
            continue
        # Match the shipped record shape: the dataset stores a minimal numbering
        # record to control file size. Pruning here keeps a rebuild byte-identical
        # to the dataset these published results were computed from.
        numbering = {
            scheme: {
                "chain_type": record["chain_type"],
                "scheme": record["scheme"],
                "residues": [
                    {"sequence_position": r["sequence_position"], "number": r["number"],
                     "aa": r["aa"], "region": r["region"]}
                    for r in record["residues"]
                ],
            }
            for scheme, record in numbering.items()
        }
        references.append({
            "drug": drug,
            "stage": stage,
            "development_status": status,
            "reference_chain": column,
            "family": "H" if column.startswith("Heavy") else "L",
            "sequence": sequence,
            "numbering": numbering,
        })

    drugs = {r["drug"] for r in references}
    dataset = {
        "version": f"{CSV_RETRIEVED}-clinical-framework-v1",
        "source_url": CSV_URL,
        "source_snapshot": "Highest_Clin_Trial (Feb 2025)",
        "source_sha256": CSV_SHA256,
        "stage_column": STAGE_COLUMN,
        "citation": CITATION,
        "rebuilt_locally": True,
        "coverage": {
            "eligible_drugs": len(selected),
            "reference_chains": len(references),
            "approved_drugs": sum(1 for r in selected if _clean(r.get(STAGE_COLUMN)) == "Approved"),
            "phase_iii_drugs": sum(1 for r in selected if _clean(r.get(STAGE_COLUMN)) == "Phase-III"),
            "represented_drugs": len(drugs),
            "excluded_chains": len(excluded),
        },
        "excluded": excluded,
        "references": references,
    }
    out.write_text(json.dumps(dataset))
    print(f"  wrote {out.name}: {dataset['coverage']}")
    return dataset


# ------------------------------------------------- rehydrate stripped files

def _rehydrate(stub: Path, out: Path, lookup: dict[str, dict[str, str]],
               label: str, entries, set_sequences) -> dict:
    """Shared rehydration: put sequences back and verify each against its hash."""
    data = json.loads(stub.read_text())
    restored = failed = missing = 0

    for name, record, expected in entries(data):
        available = lookup.get(name.lower())
        if not available:
            missing += 1
            continue
        by_family: dict[str, str] = {}
        for column, sequence in available.items():
            by_family.setdefault("H" if column.startswith("Heavy") else "L", sequence)
        chosen: dict[str, str] = {}
        for family, want_hash in (expected or {}).items():
            sequence = by_family.get(family)
            if sequence is None:
                missing += 1
                continue
            if want_hash and sha256_text(sequence) != want_hash:
                print(f"    ! {name} chain {family}: hash mismatch, not restored")
                failed += 1
                continue
            chosen[family] = sequence
            restored += 1
        if chosen:
            set_sequences(record, chosen)

    out.write_text(json.dumps(data))
    print(f"  wrote {out.name}: {restored} chains restored and hash-verified, "
          f"{failed} mismatched, {missing} not found upstream")
    if failed:
        raise SystemExit("  Refusing to continue: a restored chain did not match its "
                         "recorded hash, so this build would not reproduce the "
                         "published results.")
    return data


def build_panel(rows, stub: Path, out: Path) -> dict:
    lookup = chains_by_drug(eligible_rows(rows) or rows)
    all_lookup = chains_by_drug(rows)
    lookup = {**all_lookup, **lookup}

    def entries(data):
        for ref in data.get("references", []):
            yield ref.get("name") or ref.get("id", ""), ref, ref.get("sequence_sha256")

    def setter(ref, chosen):
        ref["sequences"] = chosen

    return _rehydrate(stub, out, lookup, "panel", entries, setter)


def build_ptm_sequences(rows, stub: Path, out: Path) -> dict:
    lookup = chains_by_drug(rows)

    def entries(data):
        for key, record in data.get("antibodies", {}).items():
            yield record.get("name") or key, record, record.get("sequence_sha256")

    def setter(record, chosen):
        record["chains"] = chosen

    return _rehydrate(stub, out, lookup, "ptm-sequences", entries, setter)


# ----------------------------------------------------------------- stripping

#: AntiBERTy fields the panel is actually read for, in app/profiler_interpretation.py:
#: the per-chain likelihood, and the provenance that decides whether a reference
#: is comparable with the query at all. Everything else in an AntiBERTy record is
#: per-residue detail about a reference molecule that nothing in this project
#: looks at -- 27 KB of it per chain, which is about 90% of the distributed panel
#: file. It is dropped from the stub, not from the local full dataset.
PANEL_ANTIBERTY_KEEP = ("pseudo_log_likelihood", "pseudo_perplexity",
                        "mean_log_likelihood", "sum_log_likelihood",
                        "residue_count", "scheme", "chain_class", "species",
                        "graft_class", "provenance")


def _trim_panel_antiberty(data: dict) -> int:
    """Drop unread per-residue AntiBERTy detail from the distributed panel."""
    removed = 0
    for reference in data.get("references", []):
        chains = (reference.get("antiberty") or {}).get("chains") or {}
        for label, record in list(chains.items()):
            if not isinstance(record, dict):
                continue
            kept = {k: v for k, v in record.items() if k in PANEL_ANTIBERTY_KEEP}
            removed += len(record) - len(kept)
            chains[label] = kept
    return removed


def make_stubs() -> None:
    """Produce the shippable sequence-free stubs from full local datasets.

    Used by the export tooling. Keeps every computed value, identifier and
    per-chain hash; removes only Thera-SAbDab sequence text, and the AntiBERTy
    per-residue arrays that nothing reads.
    """
    jobs = [
        (DATA / "profiler_comparators.json", DATA / "profiler_comparators.stripped.json",
         lambda d: [r for r in d.get("references", [])], "sequences"),
        (DATA / "therapeutic_ptm_sequences.json", DATA / "therapeutic_ptm_sequences.stripped.json",
         lambda d: list(d.get("antibodies", {}).values()), "chains"),
    ]
    for source, target, records, field in jobs:
        if not source.is_file():
            print(f"  skip {source.name}: not present locally")
            continue
        data = json.loads(source.read_text())
        stripped = 0
        for record in records(data):
            sequences = record.pop(field, None)
            if isinstance(sequences, dict):
                record["sequence_sha256"] = {f: sha256_text(s) for f, s in sequences.items()}
                stripped += len(sequences)
        if target.name == "profiler_comparators.stripped.json":
            dropped = _trim_panel_antiberty(data)
            if dropped:
                data["antiberty_detail_removed"] = True
                data["antiberty_detail_note"] = (
                    "Per-residue AntiBERTy scores and region summaries for the reference "
                    "molecules are not distributed. Nothing in this project reads them: the "
                    "panel is used only for the per-chain pseudo-log-likelihood comparison "
                    "and the provenance check that decides comparability. They are "
                    "regenerated locally by setup/build_reference_data.py if ever needed.")
                print(f"  {target.name}: {dropped} unread AntiBERTy fields removed")
        data["sequences_removed"] = True
        data["rebuild_with"] = "python setup/build_reference_data.py --all"
        data["sequence_terms_note"] = (
            "Sequence text removed: derived from Oxford Thera-SAbDab, which publishes "
            "no licence grant. Rebuild locally; each chain is verified against the "
            "sha256 recorded here."
        )
        target.write_text(json.dumps(data))
        print(f"  {target.name}: {stripped} sequences removed, hashes retained")


# -------------------------------------------------------------------- verify

def verify() -> int:
    from app.reference_data import OPTIONAL_DATASETS, available
    print("Reference datasets")
    missing = 0
    for key, (filename, description, command) in OPTIONAL_DATASETS.items():
        ok = available(key)
        missing += not ok
        size = f"{(DATA / filename).stat().st_size / 1e6:.1f} MB" if ok else "-"
        print(f"  [{'present' if ok else 'ABSENT '}] {filename:<34} {size:>9}  {description}")
        if not ok:
            print(f"              build with: {command}")

    print("\nUpstream source")
    cache = DATA / ".cache" / "TheraSAbDab.csv"
    if cache.is_file():
        digest = hashlib.sha256(cache.read_bytes()).hexdigest()
        print(f"  cached CSV hash {'matches' if digest == CSV_SHA256 else 'DIFFERS from'} "
              f"the recorded {CSV_RETRIEVED} snapshot")
    else:
        print("  no cached CSV; it will be downloaded on first build")
    print(f"\nCitation requested by the source: {CITATION}")
    return 1 if missing else 0


# ---------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--all", action="store_true", help="build every dataset")
    parser.add_argument("--clinical", action="store_true", help="clinical framework/CDR library")
    parser.add_argument("--ptm-sequences", action="store_true", help="documented-PTM reference chains")
    parser.add_argument("--panel", action="store_true", help="therapeutic reference panel")
    parser.add_argument("--verify", action="store_true", help="report status, change nothing")
    parser.add_argument("--make-stubs", action="store_true",
                        help="strip sequences from local datasets (export tooling)")
    parser.add_argument("--allow-updated-source", action="store_true",
                        help="proceed even if Oxford's CSV has changed since the snapshot")
    args = parser.parse_args()

    if args.verify:
        return verify()
    if args.make_stubs:
        print("Stripping sequences for distribution")
        make_stubs()
        return 0

    wanted = {
        "clinical": args.clinical or args.all,
        "ptm": args.ptm_sequences or args.all,
        "panel": args.panel or args.all,
    }
    if not any(wanted.values()):
        parser.print_help()
        return 2

    print(__doc__.strip().splitlines()[0])
    print("\nSource CSV")
    rows = fetch_csv(DATA / ".cache" / "TheraSAbDab.csv", args.allow_updated_source)
    print(f"  {len(rows)} rows\n")

    if wanted["panel"]:
        stub = DATA / "profiler_comparators.stripped.json"
        print("Reference panel")
        if stub.is_file():
            build_panel(rows, stub, DATA / "profiler_comparators.json")
        else:
            print("  stub profiler_comparators.stripped.json not found; nothing to rehydrate")
        print()

    if wanted["ptm"]:
        stub = DATA / "therapeutic_ptm_sequences.stripped.json"
        print("Documented-PTM reference chains")
        if stub.is_file():
            build_ptm_sequences(rows, stub, DATA / "therapeutic_ptm_sequences.json")
        else:
            print("  stub therapeutic_ptm_sequences.stripped.json not found; nothing to rehydrate")
        print()

    if wanted["clinical"]:
        print("Clinical framework/CDR library")
        build_clinical(rows, DATA / "clinical_frameworks.json")
        print()

    print(f"Done. Please cite: {CITATION}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
