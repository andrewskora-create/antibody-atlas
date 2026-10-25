"""Auditable human V-gene framework identity by exact IMGT position."""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

REFERENCE_PATH = Path(__file__).resolve().parents[1] / "data" / "human_v_ogrdb_numbered.json"


@lru_cache(maxsize=1)
def _references() -> dict:
    return json.loads(REFERENCE_PATH.read_text())


def framework_identity(sequence: str, imgt_numbering: dict[str, Any]) -> dict[str, Any]:
    """Compare query FR1-FR3 amino acids at identical IMGT positions."""
    framework = [r for r in imgt_numbering["residues"]
                 if r["number_base"] <= 104 and r["region"] == "Framework"]
    if not framework:
        raise ValueError("No V-region framework was numbered.")
    query = {r["number"][1:]: r for r in framework}
    locus = {"H": "IGH", "K": "IGK", "L": "IGL"}[imgt_numbering["chain_type"]]
    candidates = []
    for ref in _references()["references"]:
        if ref["locus"] != locus:
            continue
        positions = ref["positions"]
        covered = sum(position in positions for position in query)
        if covered / len(query) < 0.7:
            continue
        mismatches = [{"number": query[position]["number"],
                       "query_aa": residue["aa"], "reference_aa": positions.get(position)}
                      for position, residue in query.items()
                      if positions.get(position) != residue["aa"]]
        matches = len(query) - len(mismatches)
        candidates.append({"name": ref["name"], "percent": round(matches / len(query) * 100, 1),
                           "matches": matches, "framework_positions": len(query),
                           "aligned_positions": covered,
                           "coverage_percent": round(covered / len(query) * 100, 1),
                           "mismatches": mismatches})
    if not candidates:
        raise ValueError("No human V reference covered enough of the numbered framework.")
    candidates.sort(key=lambda x: (-x["matches"], -x["aligned_positions"], x["name"]))
    best = candidates[0]
    audit = candidates[:5]
    for family in ("IGHV3-48", "IGHV3-69-1") if locus == "IGH" else ():
        family_candidate = next((c for c in candidates if c["name"].startswith(family + "*")), None)
        if family_candidate and family_candidate["name"] not in {c["name"] for c in audit}:
            audit.append(family_candidate)
    return {
        "name": "Human V-gene framework identity",
        "percent": best["percent"], "matches": best["matches"],
        "framework_positions": best["framework_positions"],
        "aligned_positions": best["aligned_positions"],
        "coverage_percent": best["coverage_percent"],
        "nearest_human_v_gene": best["name"], "locus": locus,
        "candidate_alleles": audit,
        "input_length": len(sequence),
        "input_sha256_prefix": hashlib.sha256(sequence.encode()).hexdigest()[:16],
        "definition": "Exact amino-acid identity at shared IMGT FR1-FR3 positions against numbered functional human V-gene alleles. Query framework positions absent from a reference count as mismatches. CDRs and FR4 excluded.",
        "reference_source": _references()["source"],
        "reference_url": _references()["source_url"],
        "reference_release": _references()["release_versions"][locus],
        "caveat": "Framework identity is a similarity measure, not a validated immunogenicity or clinical ADA prediction. Other tools may use different reference sets, regions, and denominators.",
    }
