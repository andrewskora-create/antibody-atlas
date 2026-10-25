"""Evidence-bounded comparison of two already analyzed antibody variants.

This module performs no inference and never uses one variant's structure for another.
All positions refer to the one-based input sequence, not an aligned numbering scheme.
"""
from __future__ import annotations

import hashlib
import math
from typing import Any


def _digest(chains: dict[str, dict]) -> str:
    payload = "\n".join(f"{label}:{chains[label]['sequence']}" for label in sorted(chains))
    return hashlib.sha256(payload.encode()).hexdigest()


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _numbering(chain: dict) -> dict[int, dict]:
    return {r["sequence_position"]: r for r in chain.get("numbering", {}).get("residues", [])
            if isinstance(r.get("sequence_position"), int)}


def _metric(left: Any, right: Any, reason: str = "Both values are required.") -> dict:
    if not _finite(left) or not _finite(right):
        return {"baseline": None, "candidate": None, "delta": None,
                "status": "unavailable", "reason": reason}
    return {"baseline": left, "candidate": right,
            "delta": round(right - left, 6), "status": "available"}


def _humanness(chain_a: dict, chain_b: dict, record_a: dict, record_b: dict) -> dict:
    a, b = chain_a.get("humanness"), chain_b.get("humanness")
    if not isinstance(a, dict) or not isinstance(b, dict):
        return _metric(None, None, "Humanness result is missing for one or both variants.")
    provenance_keys = ("name", "locus", "reference_source", "reference_release", "definition")
    ref_a = record_a.get("provenance", {}).get("reference_sha256")
    ref_b = record_b.get("provenance", {}).get("reference_sha256")
    if any(not a.get(key) or a.get(key) != b.get(key) for key in provenance_keys) or (
        not ref_a or ref_a != ref_b
    ):
        return _metric(None, None, "Humanness references or definitions differ.")
    return _metric(a.get("percent"), b.get("percent"), "Humanness percent is missing.")


def _antiberty(chain_a: dict, chain_b: dict) -> dict:
    a, b = chain_a.get("antiberty"), chain_b.get("antiberty")
    if not isinstance(a, dict) or not isinstance(b, dict):
        return _metric(None, None, "AntiBERTy was not run for both variants.")
    pa, pb = a.get("provenance"), b.get("provenance")
    if not isinstance(pa, dict) or not isinstance(pb, dict) or not pa or pa != pb or a.get("scheme") != b.get("scheme"):
        return _metric(None, None, "AntiBERTy checkpoint, method, version, or numbering differs.")
    return _metric(a.get("mean_log_likelihood"), b.get("mean_log_likelihood"),
                   "AntiBERTy mean log likelihood is missing.")


def _sequence_sites(record: dict) -> dict[tuple[str, int, str], dict]:
    result = {}
    for label, chain in record["chains"].items():
        for site in chain.get("metrics", {}).get("ptm_sites", []):
            position, kind = site.get("start"), site.get("kind")
            if not isinstance(position, int) or not kind:
                continue
            key = (label, position, kind)
            if key in result:
                raise ValueError(f"Duplicate PTM event: {label} {position} {kind}.")
            result[key] = site
    return result


def _structure_ready(record: dict) -> tuple[bool, str]:
    structure = record.get("structure_analysis")
    provenance = record.get("provenance", {})
    if not isinstance(structure, dict):
        return False, "Structure analysis is missing for one or both variants."
    digest = _digest(record["chains"])
    if not provenance.get("input_sha256") or provenance.get("input_sha256") != digest:
        return False, "Analysis sequence digest is absent or does not match the input."
    if provenance.get("structure_input_sha256") != digest:
        return False, "Structure input digest is absent or does not match this variant."
    if not all(provenance.get(key) for key in
               ("structure_pdb_sha256", "structure_model", "structure_model_version")):
        return False, "Structure model identity or PDB digest is missing."
    if structure.get("scheme") != record.get("scheme") or structure.get("ph") != record.get("ph"):
        return False, "Structure analysis numbering or pH differs from its sequence analysis."
    mapping = structure.get("chain_mapping")
    if not isinstance(mapping, dict) or set(mapping) != set(record["chains"]):
        return False, "Structure chain mapping does not match the analyzed chain set."
    return True, ""


def _structural_sites(structure: dict) -> dict[tuple[str, int, str], dict]:
    result = {}
    for state, collection in (("mapped", structure.get("ptm_sites", [])),
                              ("unassessed", structure.get("unassessed_ptm_sites", []))):
        for site in collection:
            chain, position = site.get("chain"), site.get("sequence_position")
            if not chain or not isinstance(position, int):
                continue
            flags = site.get("flags") or [{"kind": site.get("kind")}]
            for flag in flags:
                kind = flag.get("kind")
                if kind:
                    result[(chain, position, kind)] = {
                        "state": state,
                        "priority": flag.get("priority", site.get("priority") if len(flags) == 1 else None),
                        "exposure_band": site.get("exposure_band"),
                        "exposure_percent": site.get("exposure_percent"),
                        "grouped_flags": len(flags) > 1,
                    }
    return result


_PM_SEQUENCE_KEYS = {"n_charged_res_fv", "vh_charge", "vl_charge", "fv_charge", "fv_csp"}
_PM_PATCH_KEYS = {"hyd_patch_area", "pos_patch_area", "neg_patch_area",
                  "hyd_patch_area_cdr", "pos_patch_area_cdr", "neg_patch_area_cdr", "heiden_score"}
_PM_STRUCTURE_KEYS = {"hyd_asa", "hph_asa", "aromatic_asa", "cdr_length", "aromatic_cdr",
                      "exposed_aromatic", "net_charge", "exposed_net_charge", "net_charge_cdr",
                      "fv_chml", "exposed_fv_chml", "dipole_moment", "scm",
                      "ann_index_pos", "ann_index_neg", "ann_index_aro",
                      "ripley_k_pos", "ripley_k_neg", "ripley_k_aro"}


def _propermab(a: dict, b: dict) -> dict:
    result = {"status": "unavailable", "reason": None, "features": [], "unavailable": []}
    pm_a, pm_b = a.get("propermab"), b.get("propermab")
    if not isinstance(pm_a, dict) or not isinstance(pm_b, dict):
        result["reason"] = "PROPERMAB results are missing for one or both variants."
        return result
    digest_a, digest_b = _digest(a["chains"]), _digest(b["chains"])
    prov_a, prov_b = a.get("provenance", {}), b.get("provenance", {})
    if any(prov.get("input_sha256") != digest or prov.get("propermab_input_sha256") != digest
           for prov, digest in ((prov_a, digest_a), (prov_b, digest_b))):
        result["reason"] = "PROPERMAB input identity is absent or differs from the analyzed sequence."
        return result
    compatibility = ("source_commit", "comparison_protocol", "adapter_version", "numbering_scheme", "ph", "scope")
    if (any(not pm_a.get(key) or pm_a.get(key) != pm_b.get(key) for key in compatibility) or
            pm_a.get("model_checkpoint") != pm_b.get("model_checkpoint")):
        result["reason"] = "PROPERMAB source, protocol, adapter, numbering, pH, or scope differs."
        return result
    def feature_map(pm: dict) -> dict[str, dict]:
        mapped = {}
        for f in pm.get("features", []):
            key = f.get("key")
            if key and key not in mapped:
                mapped[key] = f
        return mapped
    fa, fb = feature_map(pm_a), feature_map(pm_b)
    result["status"] = "available"
    for key in sorted(set(fa) | set(fb)):
        stage = ("sequence" if key in _PM_SEQUENCE_KEYS else
                 "surface_patch" if key in _PM_PATCH_KEYS else
                 "structure" if key in _PM_STRUCTURE_KEYS else "unknown")
        f1, f2 = fa.get(key), fb.get(key)
        item = {"key": key, "name": (f1 or f2).get("label", key),
                "unit": (f1 or f2).get("unit", ""), "stage": stage}
        reason = None
        if stage == "unknown":
            reason = "Descriptor stage is unrecognized."
        elif not f1 or not f2:
            reason = "Descriptor is missing for one variant."
        elif f1.get("unit", "") != f2.get("unit", ""):
            reason = "Descriptor units differ."
        elif stage != "sequence":
            owned = all(
                prov.get("propermab_structure_pdb_sha256") and
                prov.get("propermab_structure_pdb_sha256") == prov.get("structure_pdb_sha256") == pm.get("input_pdb_sha256")
                for prov, pm in ((prov_a, pm_a), (prov_b, pm_b)))
            if pm_a.get('structure_partial') or pm_b.get('structure_partial'):
                reason = 'Structural comparisons are withheld for terminal-trimmed models.'
            elif not owned:
                reason = "Both structure descriptors require their own matched input PDBs."
            elif stage == "surface_patch" and (
                not pm_a.get("patch_pipeline") or pm_a.get("patch_pipeline") != pm_b.get("patch_pipeline") or
                not pm_a.get("patch_conditions") or pm_a.get("patch_conditions") != pm_b.get("patch_conditions")
            ):
                reason = "Surface patch pipelines or conditions differ."
        metric = _metric(f1.get("value") if f1 else None, f2.get("value") if f2 else None)
        if not reason and metric["status"] == "unavailable":
            reason = "Descriptor is nonnumeric or missing."
        if reason:
            result["unavailable"].append({**item, "reason": reason})
        else:
            result["features"].append({**item, "baseline": metric["baseline"],
                                       "candidate": metric["candidate"], "delta": metric["delta"]})
    for side, pm in (("Original", pm_a), ("Candidate", pm_b)):
        for missing in pm.get("unavailable", []):
            result["unavailable"].append({"key": "", "name": f"{side}: {missing.get('feature', 'Descriptor group')}",
                "unit": "", "stage": "unavailable", "reason": missing.get("reason", "Not calculated.")})
    if not result["features"]:
        result["status"] = "unavailable"
        result["reason"] = "No compatible numeric descriptors are available for both inputs."
    return result


def _transition(a: dict | None, b: dict | None, key: str) -> dict:
    left = a.get(key) if a else None
    right = b.get(key) if b else None
    if not a or not b or a["state"] != "mapped" or b["state"] != "mapped" or left is None or right is None:
        return {"status": "unavailable", "baseline": None, "candidate": None,
                "reason": "PTM site was not mapped with this descriptor in both structures."}
    return {"status": "changed" if left != right else "unchanged",
            "baseline": left, "candidate": right}


def _molecule_comparison(baseline: dict, candidate: dict, isotype: str = "fv",
                         light_constant: str | None = None,
                         clip_c_terminal_lysine: bool = True) -> dict:
    """Same physical-chain format and constants for both sequence descriptors."""
    from .constant_regions import assemble, ConstantRegionError
    from .molecule_flags import _molecular_charge, _molecular_pi
    if isotype not in {"fv", "IGHG1", "IGHG2", "IGHG3", "IGHG4"}:
        raise ValueError("Choose fv or a supported human IgG subclass.")
    assemblies = []
    components = []
    ph = baseline["ph"]
    try:
        for record in (baseline, candidate):
            chains = record["chains"]
            if isotype == "fv":
                components.append({label: (chain["sequence"], 1) for label, chain in chains.items()})
            else:
                assembly = assemble(chains, heavy_isotype=isotype, light_constant=light_constant,
                                    clip_c_terminal_lysine=clip_c_terminal_lysine)
                assemblies.append(assembly)
                # A mutation must not silently switch the reference light constant.
                if assembly.get("light"):
                    light_constant = assembly["light"]["constant"]
                components.append({label: (assembly[key]["sequence"], assembly[key]["copies"])
                                   for label, key in (("H", "heavy"), ("L", "light")) if assembly.get(key)})
    except ConstantRegionError as exc:
        return {"status": "unavailable", "reason": str(exc), "isotype": isotype,
                "ph": ph, "metrics": {}, "per_chain": {}}
    def values(parts):
        return {"net_charge_at_ph": round(_molecular_charge(list(parts.values()), ph), 2),
                "isoelectric_point": _molecular_pi(list(parts.values()))}
    before, after = map(values, components)
    per_chain = {}
    for label in sorted(components[0]):
        a, b = components[0][label], components[1][label]
        va, vb = values({label: (a[0], 1)}), values({label: (b[0], 1)})
        per_chain[label] = {"copies": a[1], "metrics": {key: _metric(va[key], vb[key]) for key in va}}
    basis = ("Variable domains: one copy of each supplied chain." if isotype == "fv" else
             f"Full-length {assemblies[0]['format']}, {isotype}, " +
             (f"{light_constant} light constant, " if light_constant else "") +
             ("heavy C-terminal lysine clipped." if assemblies[0]['c_terminal_lysine_clipped'] else
              "heavy C-terminal lysine retained."))
    return {"status": "available", "isotype": isotype, "ph": ph,
            "format": "Fv" if isotype == "fv" and len(components[0]) == 2 else
                      "Single variable domain" if isotype == "fv" else assemblies[0]["format"],
            "basis": basis + " Charge sums physical chains, with unmodified termini; pI is the zero of that summed charge.",
            "metrics": {key: _metric(before[key], after[key]) for key in before},
            "per_chain": per_chain,
            "evidence": "Sequence calculation; not an experimental measurement or electrostatic-potential calculation.",
            "structure_scope": "Variable domains only; full-length sequence assembly does not produce a full-length structure.",
            "assumptions": assemblies[0]["assumptions"] if assemblies else
                           ["Constant regions excluded; sequence-only ionization model with unmodified chain termini."]}


def _structure_selection(record: dict, label: str, position: int) -> dict | None:
    ready, _ = _structure_ready(record)
    if not ready:
        return None
    matches = [r for r in record["structure_analysis"].get("residues", [])
               if r.get("chain") == label and r.get("sequence_position") == position and r.get("pdb_number")]
    if len(matches) != 1:
        return None
    return {"chain": label, "residue": str(matches[0]["pdb_number"])}


def compare_records(baseline: dict, candidate: dict, *, isotype: str = "fv",
                    light_constant: str | None = None,
                    clip_c_terminal_lysine: bool = True) -> dict:
    """Compare analyzed records. Reject chain, pH, scheme, or length mismatch.

    Raises ValueError for unsupported input/indels. Optional methods become
    explicitly unavailable instead of yielding a surrogate or combined score.
    """
    if not isinstance(baseline, dict) or not isinstance(candidate, dict):
        raise ValueError("Provide two analyzed records.")
    chains_a, chains_b = baseline.get("chains"), candidate.get("chains")
    if not isinstance(chains_a, dict) or not chains_a or not isinstance(chains_b, dict):
        raise ValueError("Both records require analyzed chains.")
    if set(chains_a) != set(chains_b):
        raise ValueError("Variants must contain the same chain set.")
    if baseline.get("scheme") != candidate.get("scheme") or baseline.get("scheme") not in {"imgt", "kabat"}:
        raise ValueError("Variants require the same IMGT or Kabat numbering scheme.")
    if not _finite(baseline.get("ph")) or baseline.get("ph") != candidate.get("ph"):
        raise ValueError("Variants require the same analysis pH.")
    for label in chains_a:
        seq_a, seq_b = chains_a[label].get("sequence"), chains_b[label].get("sequence")
        if not isinstance(seq_a, str) or not isinstance(seq_b, str) or not seq_a or not seq_b:
            raise ValueError(f"{label} chain sequence is missing.")
        if len(seq_a) != len(seq_b):
            raise ValueError(f"{label} chain lengths differ; insertions and deletions are unsupported in this release.")

    ready_a, why_a = _structure_ready(baseline)
    ready_b, why_b = _structure_ready(candidate)
    structure_a, structure_b = baseline.get("structure_analysis"), candidate.get("structure_analysis")
    pa, pb = baseline.get("provenance", {}), candidate.get("provenance", {})
    structure_reason = why_a or why_b
    if ready_a and ready_b:
        if any(pa.get(key) != pb.get(key) for key in ("structure_model", "structure_model_version")) or (
            structure_a.get("method") != structure_b.get("method") or
            structure_a.get("patch_definition") != structure_b.get("patch_definition")
        ):
            structure_reason = "Structure models or analysis definitions differ."
    structure_available = ready_a and ready_b and not structure_reason
    structural_sites_a = _structural_sites(structure_a) if structure_available else {}
    structural_sites_b = _structural_sites(structure_b) if structure_available else {}

    mutations = []
    chains = {}
    for label in sorted(chains_a):
        chain_a, chain_b = chains_a[label], chains_b[label]
        numbered_a, numbered_b = _numbering(chain_a), _numbering(chain_b)
        changes = []
        for position, (aa_a, aa_b) in enumerate(zip(chain_a["sequence"], chain_b["sequence"]), 1):
            if aa_a != aa_b:
                na, nb = numbered_a.get(position, {}), numbered_b.get(position, {})
                changes.append({"chain": label, "position": position, "from": aa_a, "to": aa_b,
                                "baseline_number": na.get("number"), "candidate_number": nb.get("number"),
                                "baseline_region": na.get("region"), "candidate_region": nb.get("region"),
                                "baseline_structure": _structure_selection(baseline, label, position),
                                "candidate_structure": _structure_selection(candidate, label, position)})
        mutations.extend(changes)
        ma, mb = chain_a.get("metrics", {}), chain_b.get("metrics", {})
        chains[label] = {"mutations": changes, "metrics": {
            key: _metric(ma.get(key), mb.get(key), f"{key} is missing for one or both variants.")
            for key in ("isoelectric_point", "net_charge_at_ph", "gravy")}}
        chains[label]["metrics"]["humanness_percent"] = _humanness(chain_a, chain_b, baseline, candidate)
        chains[label]["metrics"]["antiberty_mean_log_likelihood"] = _antiberty(chain_a, chain_b)

    engineering_notes = []
    if any(m['from'] == 'C' or m['to'] == 'C' for m in mutations):
        engineering_notes.append('Cysteine changed: verify canonical and additional disulfide connectivity and free thiols; sequence scoring cannot establish whether this edit preserves folding.')
    if any(str(m.get('baseline_region', '')).startswith('CDR') for m in mutations):
        engineering_notes.append('CDR substitution: target binding and functional activity may change even if a PTM motif is removed.')
    if len(mutations) > 1:
        engineering_notes.append('Combination: compare component single substitutions first; effects may interact and are not additive.')
    if any(not m.get('baseline_region') or m.get('baseline_region') == 'Unassigned' for m in mutations):
        engineering_notes.append('Some substituted residues lack a reliable region assignment; functional context is unknown.')
    sites_a, sites_b = _sequence_sites(baseline), _sequence_sites(candidate)
    def describe(key: tuple[str, int, str], site: dict) -> dict:
        return {"chain": key[0], "position": key[1], "kind": key[2], "motif": site.get("motif")}
    removed = [describe(k, sites_a[k]) for k in sorted(sites_a.keys() - sites_b.keys())]
    added = [describe(k, sites_b[k]) for k in sorted(sites_b.keys() - sites_a.keys())]
    retained = []
    for key in sorted(sites_a.keys() & sites_b.keys()):
        sa, sb = structural_sites_a.get(key), structural_sites_b.get(key)
        retained.append({"chain": key[0], "position": key[1], "kind": key[2],
                         "baseline_motif": sites_a[key].get("motif"),
                         "candidate_motif": sites_b[key].get("motif"),
                         "priority_transition": _transition(sa, sb, "priority"),
                         "exposure_transition": _transition(sa, sb, "exposure_band"),
                         "exposure_percent_transition": _transition(sa, sb, "exposure_percent")})

    structure = {"status": "available" if structure_available else "unavailable",
                 "reason": None if structure_available else structure_reason,
                 "patch_counts": {}, "patch_areas": {}, "mapping_coverage_percent": {}}
    if structure_available:
        full_mapping = all(
            all(s.get("input_residues") == s.get("mapped_residues") == s.get("model_residues") ==
                len(record["chains"][label]["sequence"]) for label, s in structure["chain_mapping"].items())
            for record, structure in ((baseline, structure_a), (candidate, structure_b)))
        if full_mapping:
            for kind in sorted(set(structure_a.get("patch_counts", {})) | set(structure_b.get("patch_counts", {}))):
                structure["patch_counts"][kind] = _metric(structure_a.get("patch_counts", {}).get(kind),
                                                           structure_b.get("patch_counts", {}).get(kind))
            for kind in sorted(set(p.get("kind") for p in structure_a.get("patches", [])) |
                               set(p.get("kind") for p in structure_b.get("patches", []))):
                if kind is None:
                    continue
                sums = []
                for source in (structure_a, structure_b):
                    values = [p.get("exposed_sidechain_sasa_sum_a2") for p in source.get("patches", [])
                              if p.get("kind") == kind]
                    sums.append(round(sum(values), 3) if all(_finite(v) for v in values) else None)
                structure["patch_areas"][kind] = _metric(*sums,
                    reason="Patch area is missing or nonnumeric for one or both variants.")
        else:
            structure["aggregate_reason"] = "Patch count changes withheld because one or more modeled chains are not fully mapped."
        for label in sorted(chains_a):
            structure["mapping_coverage_percent"][label] = _metric(
                structure_a.get("mapping_coverage_percent", {}).get(label),
                structure_b.get("mapping_coverage_percent", {}).get(label))
    return {
        "baseline_analysis_id": baseline.get("analysis_id"),
        "candidate_analysis_id": candidate.get("analysis_id"),
        "scheme": baseline["scheme"], "ph": baseline["ph"],
        "mutations": mutations, "mutation_count": len(mutations), "chains": chains,
        "engineering_notes": engineering_notes,
        "molecule": _molecule_comparison(baseline, candidate, isotype, light_constant, clip_c_terminal_lysine),
        "free_energy": {"status": "unavailable", "reason": "No folding-stability or binding free-energy method was run. Geometry and raw force-field energies are not ΔΔG."},
        "ptm": {"added": added, "removed": removed, "retained": retained},
        "structure": structure, "propermab": _propermab(baseline, candidate),
        "conclusion": (f"{len(mutations)} substitution(s); {len(added)} added and {len(removed)} removed "
                       "sequence PTM flag(s). Predicted/calculated descriptor changes require binding and "
                       "stability validation; no combined rank, pass/fail call, or clinical risk estimate."),
        "limitations": [
            "PTM motifs are sequence flags, not measured modification occupancy or rates.",
            "Predicted structure reflects one modeled conformation; unmapped PTM sites remain unknown.",
            "Humanness and AntiBERTy outputs are similarity or model-compatibility descriptors, not immunogenicity predictions.",
            "Substitutions are compared at input positions; insertions and deletions are unsupported.",
        ],
        "next_validation": [
            "Confirm target binding and functional activity for each variant.",
            "Measure expression, purity, stability, and aggregation under relevant formulation conditions.",
            "Experimentally assess PTM sites and liabilities prioritized for this use case.",
        ],
    }
