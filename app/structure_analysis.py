"""Local, residue-level surface screening for an ImmuneBuilder Fv model.

This is an independent, disclosed FreeSASA/contact-graph method. It does not
implement PROPERMAB's licensed mesh, APBS, or trained developability models.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import freesasa
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

from .analysis import ptm_sites
from .numbering import number_chain


HYDROPHOBIC = set("AFILMVWCY")
# Relative visual weights only. They are not transfer free energies or
# calibrated aggregation probabilities.
HYDROPHOBIC_DISPLAY_WEIGHT = {
    "A": 0.35, "C": 0.45, "F": 0.95, "I": 0.90, "L": 0.85,
    "M": 0.65, "V": 0.80, "W": 1.00, "Y": 0.75,
}
PATCH_CONFIG = {
    "hydrophobic": {"contact_angstrom": 8.0, "min_area_a2": 60.0, "min_members": 2},
    "positive": {"contact_angstrom": 8.0, "min_area_a2": 40.0, "min_members": 2},
    "negative": {"contact_angstrom": 8.0, "min_area_a2": 40.0, "min_members": 2},
}


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def _sidechain_charge(aa: str, ph: float) -> float:
    positive_pka = {"K": 10.5, "R": 12.5, "H": 6.0}
    negative_pka = {"D": 3.9, "E": 4.1}
    if aa in positive_pka:
        return 1.0 / (1.0 + 10 ** (ph - positive_pka[aa]))
    if aa in negative_pka:
        return -1.0 / (1.0 + 10 ** (negative_pka[aa] - ph))
    return 0.0


def _centroid(residue) -> tuple[float, float, float]:
    sidechain = [atom.coord for atom in residue if atom.element != "H"
                 and atom.name not in {"N", "CA", "C", "O", "OXT"}]
    atoms = sidechain or [residue["CA"].coord]
    return tuple(float(sum(a[i] for a in atoms) / len(atoms)) for i in range(3))


def _map_model_positions(expected: str, observed: str) -> dict[int, int]:
    """Map only identical modeled residues to one-based input positions."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(expected, observed)[0]
    return {int(model_i) + 1: int(input_i) + 1
            for input_i, model_i in zip(*alignment.indices)
            if input_i >= 0 and model_i >= 0 and expected[input_i] == observed[model_i]}


def extract_trace(pdb_path: str | Path) -> dict[str, Any]:
    """Return model coordinates without relying on sequence-to-model mapping."""
    model = PDBParser(QUIET=True).get_structure("fv", str(pdb_path))[0]
    residues = []
    for chain in model:
        for ordinal, residue in enumerate((r for r in chain if r.id[0] == " " and "CA" in r), 1):
            residues.append({"chain": chain.id, "sequence_position": ordinal,
                             "aa": protein_letters_3to1.get(residue.resname, "X"),
                             "ca_xyz": [round(float(x), 3) for x in residue["CA"].coord]})
    return {"residues": residues, "patches": [], "ptm_sites": []}


def _components(candidates: list[dict[str, Any]], cutoff: float) -> list[list[dict[str, Any]]]:
    remaining = set(range(len(candidates)))
    groups = []
    while remaining:
        stack = [remaining.pop()]
        component = []
        while stack:
            i = stack.pop()
            component.append(candidates[i])
            neighbors = {j for j in remaining if _distance(candidates[i]["sidechain_xyz"],
                                                            candidates[j]["sidechain_xyz"]) <= cutoff}
            stack.extend(neighbors)
            remaining.difference_update(neighbors)
        groups.append(component)
    return groups


def _site_priority(site: dict[str, Any], residue: dict[str, Any], scheme: str) -> dict[str, Any]:
    """Evidence-informed triage, deliberately not a calibrated PTM predictor.

    Exposure bins are display heuristics. Published thresholds use different
    atom selections, normalizations, dynamics and stress conditions.
    """
    kind = site["kind"]
    motif = site["motif"]
    raw_exposure = residue.get("relative_sidechain_sasa")
    exposure = raw_exposure if isinstance(raw_exposure, (int, float)) and math.isfinite(raw_exposure) and raw_exposure >= 0 else None
    exposure_band = "unknown" if exposure is None else "low" if exposure < 0.10 else "intermediate" if exposure < 0.30 else "high"
    strong_motif = kind in {"Asparagine deamidation motif", "Aspartate isomerization motif"} and motif in {"NG", "NS", "DG", "DS"}
    if kind == "N-linked glycosylation sequon":
        priority = "review"
        evidence_note = "Sequon requires occupancy assessment; glycosylation occurs during biosynthesis, so native-state burial does not remove the flag."
    elif kind in {'Asp-Pro cleavage motif', 'Noncanonical cysteine'}:
        priority = 'review'
        evidence_note = site['interpretation'] + ' This review flag is retained regardless of native-state exposure.'
    elif kind == 'N-terminal pyroglutamate formation':
        priority = 'lower'
        evidence_note = site['interpretation'] + ' Track as characterization context; no generic adverse functional effect is assigned.'
    elif strong_motif:
        priority = "high" if exposure is not None and exposure >= 0.30 else "review"
        evidence_note = "Sequence liability retained even when buried; exposure raises screening priority but does not predict modification rate."
    else:
        priority = "review" if exposure is None or exposure >= 0.10 else "lower"
        evidence_note = "Exposure contributes to susceptibility, but chemistry, local environment, dynamics and stress conditions are also required; exposure alone does not establish high risk."
    # A framework PTM is not a high-priority finding. The chemistry is the same
    # as in a CDR, but the functional consequence is not: framework motifs are
    # carried by approved products, and the ones that matter are a CMC question
    # rather than a binding one. High is reserved for CDR sites; framework sites
    # are capped at review and routed to the approved/Phase-III comparator search.
    framework_capped = False
    if priority == "high" and not residue["region"].startswith("CDR"):
        priority = "review"
        framework_capped = True
        evidence_note += (" Priority is capped below high because the site is outside a CDR; "
                          "framework liabilities are assessed against approved and Phase-III "
                          "comparators and handled as a CMC question.")
    risk_label = {"high": "High screening priority", "review": "Review", "lower": "Lower structural concern"}[priority]
    if exposure is None:
        risk_label = "Unknown structural risk — review"
    points = {"high": 4, "review": 3, "lower": 1}[priority]
    impact = "CDR: binding impact requires testing" if residue["region"].startswith("CDR") else "Functional impact unknown"
    reasons = [f"relative side-chain SASA {exposure:.1%}" if exposure is not None else "Surface exposure unavailable; structural risk is unknown"]
    if residue["region"].startswith("CDR"):
        reasons.append(f"in {residue['region']}")
    if strong_motif:
        reasons.append("sequence liability retained independently of burial")
    if framework_capped:
        reasons.append("outside a CDR, so not rated high; see approved/Phase-III comparators")
    if kind == "N-linked glycosylation sequon":
        reasons.append("sequon presence does not establish glycan occupancy")
    return {**site, "chain": residue["chain"], "residue": residue["aa"],
            "number": residue["number"], "region": residue["region"],
            "scheme": scheme, "sequence_position": site["start"],
            "relative_sidechain_sasa": round(exposure, 3) if exposure is not None else None,
            "sidechain_sasa_a2": round(residue["sidechain_sasa_a2"], 1) if residue.get("sidechain_sasa_a2") is not None and math.isfinite(residue["sidechain_sasa_a2"]) else None,
            "exposure_percent": round(exposure * 100, 1) if exposure is not None else None,
            "exposure_band": exposure_band, "risk_label": risk_label,
            "evidence_note": evidence_note,
            "priority_points": points, "priority": priority,
            "framework_priority_capped": framework_capped,
            "modification_likelihood": risk_label + "; modification probability not estimated",
            "functional_impact": impact, "evidence_status": "Predicted only",
            "reasons": reasons}


def _group_physical_sites(sites: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Count a residue once while retaining its distinct chemistry flags."""
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for site in sites:
        grouped.setdefault((site["chain"], site["start"]), []).append(site)
    result = []
    for members in grouped.values():
        lead = max(members, key=lambda s: s.get("priority_points", -1))
        item = dict(lead)
        item["precedent_ids"] = list(dict.fromkeys(pid for s in members for pid in s.get("precedent_ids", [])))
        item["flags"] = [{"kind": s["kind"], "motif": s["motif"],
                          "interpretation": s.get("interpretation", ""),
                          "priority": s.get("priority"),
                          "priority_points": s.get("priority_points"),
                          "exposure_assessed": isinstance(s.get("relative_sidechain_sasa"), (int, float))
                          and math.isfinite(s["relative_sidechain_sasa"])
                          and s["relative_sidechain_sasa"] >= 0,
                          "chemistry": "glycosylation" if s["kind"] == "N-linked glycosylation sequon"
                          else "deamidation" if "deamidation" in s["kind"].lower()
                          else "isomerization" if "isomerization" in s["kind"].lower()
                          else "oxidation" if "oxidation" in s["kind"].lower()
                          else "cleavage" if s["kind"] == "Asp-Pro cleavage motif"
                          else "pyroglutamate" if s["kind"] == "N-terminal pyroglutamate formation"
                          else "cysteine"} for s in members]
        item["kind"] = " + ".join(dict.fromkeys(s["kind"] for s in members))
        result.append(item)
    return result


def summarize_site_priorities(sites: list[dict[str, Any]],
                             unassessed_sites: list[dict[str, Any]]) -> dict[str, Any]:
    """Refresh summary fields after either single-model or ensemble scoring."""
    counts = {level: sum(s.get("priority") == level for s in sites)
              for level in ("high", "review", "lower")}
    call = "Follow-up recommended" if counts["high"] or counts["review"] >= 3 else \
           "Review sites" if counts["review"] or counts["lower"] >= 5 else "No priority site identified"
    evidence_status = "Incomplete mapping" if unassessed_sites else "Predicted only"
    missing_exposure = any(not isinstance(s.get("exposure_percent"), (int, float))
                           or not math.isfinite(s["exposure_percent"]) for s in sites)
    if missing_exposure:
        evidence_status = "Incomplete exposure data"
    if unassessed_sites or missing_exposure:
        call = "Incomplete evidence — follow-up recommended" if counts["high"] or counts["review"] >= 3 else "Incomplete evidence — review required"
    return {"call": call, "counts": counts, "unassessed_sites": len(unassessed_sites),
            "evidence_status": evidence_status}


def analyze_structure(pdb_path: str | Path, sequences: dict[str, str], ph: float = 7.4,
                      scheme: str = "imgt") -> dict[str, Any]:
    """Analyze a paired Fv or heavy-only predicted PDB without network calls."""
    if not 0 <= ph <= 14:
        raise ValueError("pH must be between 0 and 14.")
    if scheme not in {"imgt", "kabat"}:
        raise ValueError("Numbering scheme must be IMGT or Kabat.")
    pdb_path = Path(pdb_path)
    sasa_structure = freesasa.Structure(str(pdb_path))
    sasa = freesasa.calc(sasa_structure).residueAreas()
    model = PDBParser(QUIET=True).get_structure("fv", str(pdb_path))[0]
    residues: list[dict[str, Any]] = []
    by_position: dict[tuple[str, int], dict[str, Any]] = {}
    chain_mapping = {}
    mapping_warnings = []
    for chain_id, expected_seq in sequences.items():
        if chain_id not in model:
            chain_mapping[chain_id] = {"input_residues": len(expected_seq), "model_residues": 0,
                                       "mapped_residues": 0, "coverage_percent": 0.0}
            mapping_warnings.append(f"{chain_id} chain is absent from the model; its PTM sites were not structurally assessed.")
            continue
        chain_residues = [r for r in model[chain_id] if r.id[0] == " " and "CA" in r]
        observed_seq = "".join(protein_letters_3to1.get(r.resname, "X") for r in chain_residues)
        position_map = _map_model_positions(expected_seq, observed_seq)
        coverage = len(position_map) / len(expected_seq)
        chain_mapping[chain_id] = {"input_residues": len(expected_seq),
                                   "model_residues": len(observed_seq),
                                   "mapped_residues": len(position_map),
                                   "coverage_percent": round(coverage * 100, 1)}
        usable = coverage >= 0.90 and len(position_map) / max(1, len(observed_seq)) >= 0.90
        if not usable:
            mapping_warnings.append(f"{chain_id} model maps fewer than 90% of residues; its surface and PTM calls were withheld.")
            position_map = {}
        elif len(position_map) != len(expected_seq) or len(position_map) != len(observed_seq):
            mapping_warnings.append(f"{chain_id} model maps {len(position_map)}/{len(expected_seq)} input residues; unmatched sites were not structurally scored.")
        numbered = number_chain(expected_seq, scheme) if usable else None
        numbering = {r["sequence_position"]: r for r in numbered["residues"]} if numbered else {}
        for model_position, residue in enumerate(chain_residues, 1):
            seq_position = position_map.get(model_position)
            aa = observed_seq[model_position - 1]
            index, insertion = residue.id[1], residue.id[2].strip()
            area_key = f"{index}{insertion}"
            area = sasa.get(chain_id, {}).get(area_key)
            if area is None:
                raise ValueError(f"FreeSASA omitted residue {chain_id}{area_key}.")
            annotated = numbering.get(seq_position, {})
            charge = _sidechain_charge(aa, ph)
            item = {
                "chain": chain_id, "sequence_position": seq_position,
                "scheme": scheme,
                "number": annotated.get("number", f"{chain_id}{area_key}"),
                "aa": aa, "region": annotated.get("region", "Unassigned"),
                "pdb_number": area_key,
                "ca_xyz": [round(float(x), 3) for x in residue["CA"].coord],
                "sidechain_xyz": _centroid(residue),
                "total_sasa_a2": round(area.total, 1),
                "sidechain_sasa_a2": round(area.sideChain, 1),
                # Glycine has no side chain. Any other missing measurement
                # remains unknown, rather than falsely implying burial.
                "relative_sidechain_sasa": round(area.relativeSideChain, 3)
                if math.isfinite(area.relativeSideChain) else 0.0 if aa == "G" else None,
                "sidechain_charge_at_ph": round(charge, 3),
            }
            residues.append(item)
            if seq_position is not None:
                by_position[(chain_id, seq_position)] = item

    if not by_position:
        raise ValueError("No model chain maps reliably to the input; the PDB is still available for inspection.")

    cdr_residues = [r for r in residues if r["region"].startswith("CDR")]
    patches: list[dict[str, Any]] = []
    for kind, config in PATCH_CONFIG.items():
        if kind == "hydrophobic":
            candidates = [r for r in residues if r["sequence_position"] is not None and r["aa"] in HYDROPHOBIC
                          and r["relative_sidechain_sasa"] is not None
                          and r["relative_sidechain_sasa"] >= 0.20
                          and r["sidechain_sasa_a2"] >= 12]
        else:
            sign = 1 if kind == "positive" else -1
            candidates = [r for r in residues if r["sequence_position"] is not None and r["sidechain_charge_at_ph"] * sign >= 0.5
                          and r["relative_sidechain_sasa"] is not None
                          and r["relative_sidechain_sasa"] >= 0.15
                          and r["sidechain_sasa_a2"] >= 10]
        for component in _components(candidates, config["contact_angstrom"]):
            area_sum = sum(r["sidechain_sasa_a2"] for r in component)
            if len(component) < config["min_members"] or area_sum < config["min_area_a2"]:
                continue
            near_cdr = any(_distance(member["sidechain_xyz"], cdr["sidechain_xyz"]) <= 8.0
                           for member in component for cdr in cdr_residues)
            patches.append({
                "id": f"{kind}-{len(patches)+1}", "kind": kind,
                "member_count": len(component), "exposed_sidechain_sasa_sum_a2": round(area_sum, 1),
                "near_cdr": near_cdr,
                "residues": [{"chain": r["chain"], "sequence_position": r["sequence_position"],
                              "number": r["number"], "aa": r["aa"], "region": r["region"]}
                             for r in sorted(component, key=lambda x: (x["chain"], x["sequence_position"]))],
            })
    patches.sort(key=lambda p: (p["kind"], -p["exposed_sidechain_sasa_sum_a2"]))

    # A continuous display signal combines each exposed side chain with the
    # area of its qualifying cluster. A residue outside a cluster can still
    # look cool-to-mildly warm; large exposed clusters approach the hot end.
    for residue in residues:
        exposure = min(1.0, residue["sidechain_sasa_a2"] / 80.0)
        charge = residue["sidechain_charge_at_ph"]
        residue["surface_heat"] = {
            "hydrophobic": round(0.32 * HYDROPHOBIC_DISPLAY_WEIGHT.get(residue["aa"], 0) * exposure, 3),
            "positive": round(0.32 * max(charge, 0) * exposure, 3),
            "negative": round(0.32 * max(-charge, 0) * exposure, 3),
        }
    by_model_position = {(r["chain"], r["sequence_position"]): r for r in residues
                         if r["sequence_position"] is not None}
    for patch in patches:
        cluster_signal = min(1.0, 0.35 + 0.65 * patch["exposed_sidechain_sasa_sum_a2"] / 230.0)
        for member in patch["residues"]:
            residue = by_model_position.get((member["chain"], member["sequence_position"]))
            if residue is not None:
                kind = patch["kind"]
                residue["surface_heat"][kind] = round(max(residue["surface_heat"][kind], cluster_signal), 3)

    sites = []
    unassessed_sites = []
    for chain_id, sequence in sequences.items():
        try:
            input_regions = {r["sequence_position"]: r for r in number_chain(sequence, scheme)["residues"]}
        except ValueError:
            input_regions = {}
        try:
            imgt_numbering = number_chain(sequence, 'imgt')
        except ValueError:
            imgt_numbering = None
        for site in ptm_sites(sequence, imgt_numbering):
            residue = by_position.get((chain_id, site["start"]))
            motif_mapped = all((chain_id, pos) in by_position
                               for pos in range(site["start"], site["end"] + 1))
            if residue and motif_mapped:
                sites.append(_site_priority(site, residue, scheme))
            else:
                input_residue = input_regions.get(site["start"], {})
                unassessed_sites.append({"chain": chain_id, **site,
                                         "region": input_residue.get("region", "Unassigned"),
                                         "number": input_residue.get("number", f"{chain_id}:{site['start']}"),
                                         "scheme": scheme, "sequence_position": site["start"],
                                         "relative_sidechain_sasa": None, "sidechain_sasa_a2": None,
                                         "exposure_percent": None, "exposure_band": "unknown",
                                         "risk_label": "Unknown structural risk — unassessed",
                                         "evidence_note": "Sequence flag only; missing coordinates cannot establish low PTM risk.",
                                         "reason": "No reliable coordinate mapping for this site or full motif."})
    # Documented-PTM precedent comes from a corpus this project cannot
    # redistribute, so a fresh clone does not have it. Exposure, priority and
    # grouping do not depend on it, so its absence leaves precedent_ids empty
    # rather than failing the structural analysis.
    from .ptm_precedents import sequence_precedents
    from .reference_data import ReferenceDataMissing
    try:
        precedents = sequence_precedents(
            {chain: {"sequence": seq} for chain, seq in sequences.items()}, scheme)
    except ReferenceDataMissing as exc:
        precedents = {"cases": [], "unavailable": True, "reason": str(exc)}
    for site in sites + unassessed_sites:
        site["precedent_ids"] = [case["id"] for case in precedents["cases"]
            if case["chemistry"] == site["kind"] and any(
                m["query_chain"] == site["chain"] and m["query_sequence_position"] == site["start"]
                for m in case["matches"])]
    sites = _group_physical_sites(sites)
    unassessed_sites = _group_physical_sites(unassessed_sites)
    sites.sort(key=lambda s: (-s["priority_points"], s["chain"], s["start"]))
    mapping_coverage = {chain: metrics["coverage_percent"] for chain, metrics in chain_mapping.items()}
    return {
        "method": "FreeSASA per-residue exposure + exposed side-chain contact graph",
        "scheme": scheme,
        "ph": ph,
        "residues": residues,
        "patches": patches,
        "patch_counts": {kind: sum(p["kind"] == kind for p in patches)
                         for kind in PATCH_CONFIG},
        "ptm_sites": sites,
        "unassessed_ptm_sites": unassessed_sites,
        "chain_mapping": chain_mapping,
        "mapping_coverage_percent": mapping_coverage,
        "mapping_warnings": mapping_warnings,
        "overall_ptm_risk": {
            **summarize_site_priorities(sites, unassessed_sites),
            "calibration_note": "Uncalibrated means these screening categories have not been fitted and independently validated against measured PTM frequencies or rates under specified conditions. They do not estimate percent modified, probability, shelf life or clinical risk.",
            "definition": "Rule-based follow-up screen; one physical residue is counted once. Relative side-chain SASA bins are heuristic: low <10%, intermediate 10–<30%, high >=30%; these are not universal PTM thresholds. NG/NS deamidation and DG/DS isomerization are high priority at >=30% exposure and review otherwise, including buried sites. Other Asn/Asp motifs and Met/Trp oxidation are review at >=10% and lower structural concern below 10%; exposure alone never makes oxidation high priority. Glycosylation sequons, Asp-Pro cleavage motifs and noncanonical cysteines require review regardless of native-state exposure. Additional cysteines are not assumed unpaired. N-terminal Q/E cyclization is a conditional characterization flag with lower default priority, not a generic adverse effect; the input must represent a free mature N terminus. Missing exposure is unknown and requires review; unmapped sites are unassessed separately. CDR location indicates potential functional impact separately. Overall call: follow-up recommended for any high site or >=3 review sites; review sites for any review site or >=5 lower sites; otherwise no priority site identified. This is triage, not total PTM probability.",
            "exposure_definition": "Relative side-chain SASA from FreeSASA, expressed as a percentage of its residue-specific reference side-chain area, using the default 1.4 Å water probe. This is a single predicted conformation, not percent of molecules modified; values can exceed 100% for conformations exceeding the reference area.",
        },
        "limitations": [
            "Single predicted structure; surface exposure and side-chain placement may change with conformation and formulation.",
            "Patch area is the sum of member side-chain SASA, not a triangulated surface-mesh area or PROPERMAB value.",
            "Fixed residue pKa approximations are used for charge at the selected pH; no APBS electrostatics calculation.",
            "Folded-state exposure does not establish PTM occupancy, chemical rate, glycosylation, or clinical developability.",
        ],
        "patch_definition": PATCH_CONFIG,
        "surface_heat_definition": "Relative blue-to-red display signal only: exposed residue-class contribution (0-0.32) plus qualifying cluster-area signal (up to 1). Each mode is independently scaled; colors are not measured aggregation, electrostatic potential, or PK risk.",
    }
