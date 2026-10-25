"""Enumerate and rank conservative substitutions at PTM-liable positions.

Scope and intent
----------------
This module answers one question: *given a PTM motif that the existing screen has
flagged, which single-residue substitutions remove that chemistry while changing
the molecule as little as possible?* It is aimed at CDR sites, where an arbitrary
"textbook" fix is most likely to cost binding, but it will run on framework sites
too.

It does not predict binding, potency, expression, or aggregation. It ranks a blend of residue similarity and clinical reference frequency,
with new-liability and preserved-position safeguards. Choosing among the ranked options remains an experimental decision.

Method
------
For each flagged motif, the editable positions are derived from the chemistry
(for example, a deamidation NG motif can be addressed at the Asn or at the +1
Gly). At each editable position all 19 alternative residues are generated, then:

  1. **Filtered for efficacy.** The candidate chain is rescanned with
     ``analysis.ptm_sites``. A substitution is kept only if it actually removes
     the target chemistry at that site.
  2. **Checked for collateral damage.** The same rescan reports any *new* PTM
     chemistry the substitution introduces anywhere in the chain. This catches
     the self-defeating fixes that a static substitution table misses -- for
     example NG -> NA, which still matches the slow-deamidation NT/NQ/NA class
     in several screening conventions, or D -> N next to Gly, which trades an
     isomerization motif for a deamidation motif.
  3. **Shortlisted by residue similarity and clinical reference frequency**
     with equal weight after preserved-position and new-liability safeguards.
     Four components are reported separately (germline is only a tiebreaker):
       - BLOSUM62 exchange score;
       - a physicochemical distance over hydropathy, side-chain volume, charge,
         polarity and aromaticity;
       - the residue's frequency at the same IMGT position across functional
         human V alleles (OGRDB), which indicates whether the substitution is
         something human germline actually does there;
       - the residue's frequency at the same IMGT position across the local
         clinical-stage chain library, which indicates whether approved or
         Phase-III molecules carry it there.

Germline and clinical frequencies are *observational*. A residue being common at
a position does not establish that it is tolerated in this antibody, and absence
is not evidence against it. CDR3 is largely outside V-gene coverage, so germline
frequency is usually unavailable there and is reported as such rather than as
zero.

Nothing here is applied automatically. The output includes a variant
specification string that can be pasted into the existing variant comparison, so
that every proposal is re-analysed by the same engines as any hand-entered
candidate.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any

from .analysis import ptm_sites
from .developability_delta import mutation_effect
from .ptm_precedents import _canonical

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
GERMLINE_PATH = DATA_DIR / "human_v_ogrdb_numbered.json"
CLINICAL_PATH = DATA_DIR / "clinical_frameworks.json"

AA = "ACDEFGHIKLMNPQRSTVWY"

LIMITATION = (
    "Substitutions are shortlisted by residue similarity and clinical reference frequency, "
    "after preserved-position and new-liability safeguards. "
    "They are not predictions of retained antigen binding, potency, expression, "
    "thermostability or aggregation behaviour. Germline and clinical frequencies are "
    "observational: a common residue is not established as tolerated in this antibody, "
    "and an absent residue is not established as harmful. Every candidate requires "
    "binding and activity testing, and the underlying liability should be confirmed by "
    "site-resolved peptide mapping before any residue is changed."
)

# --- physicochemical reference values ---------------------------------------
# Kyte-Doolittle hydropathy; side-chain volume (A^3, Zamyatnin); formal side-chain
# charge at neutral pH; polarity and aromaticity as booleans.
_HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
    "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
    "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
_VOLUME = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1,
    "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1,
    "P": 112.7, "Q": 143.8, "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0,
    "W": 227.8, "Y": 193.6,
}
_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
_POLAR = set("DEKRHNQSTY")
_AROMATIC = set("FWYH")


@lru_cache(maxsize=1)
def _blosum62():
    """BLOSUM62 as a plain dict; falls back to None if Biopython is unavailable."""
    try:
        from Bio.Align import substitution_matrices
    except Exception:
        return None
    matrix = substitution_matrices.load("BLOSUM62")
    table: dict[tuple[str, str], int] = {}
    for a in AA:
        for b in AA:
            try:
                table[(a, b)] = int(matrix[a, b])
            except Exception:
                pass
    return table or None


def _physchem_distance(old: str, new: str) -> dict[str, Any]:
    """Normalised 0-1 distance plus the per-property deltas that produced it.

    Each property is scaled by its range across the 20 residues so no single
    axis dominates. Reported alongside the components, not as a bare number.
    """
    hyd = abs(_HYDROPATHY[old] - _HYDROPATHY[new]) / 9.0          # range -4.5..4.5
    vol = abs(_VOLUME[old] - _VOLUME[new]) / 167.7                # range 60.1..227.8
    chg = abs(_CHARGE.get(old, 0.0) - _CHARGE.get(new, 0.0)) / 2.0
    pol = 0.0 if (old in _POLAR) == (new in _POLAR) else 1.0
    aro = 0.0 if (old in _AROMATIC) == (new in _AROMATIC) else 1.0
    distance = (hyd + vol + chg + pol + aro) / 5.0
    return {
        "distance": round(distance, 4),
        "hydropathy_delta": round(_HYDROPATHY[new] - _HYDROPATHY[old], 2),
        "volume_delta_a3": round(_VOLUME[new] - _VOLUME[old], 1),
        "charge_delta": round(_CHARGE.get(new, 0.0) - _CHARGE.get(old, 0.0), 2),
        "polarity_changed": bool(pol),
        "aromaticity_changed": bool(aro),
    }


# --- positional precedent ----------------------------------------------------
@lru_cache(maxsize=1)
def _germline_index() -> dict[str, Any]:
    """Residue counts per canonical IMGT position per locus, from OGRDB."""
    data = json.loads(GERMLINE_PATH.read_text())
    index: dict[str, dict[str, Counter]] = {}
    totals: dict[str, dict[str, int]] = {}
    for ref in data["references"]:
        if ref.get("functional") is False:
            continue
        locus = ref["locus"]
        per_locus = index.setdefault(locus, {})
        per_total = totals.setdefault(locus, {})
        for key, aa in ref["positions"].items():
            per_locus.setdefault(key, Counter())[aa] += 1
            per_total[key] = per_total.get(key, 0) + 1
    return {
        "counts": index,
        "totals": totals,
        "source": data["source"],
        "source_url": data["source_url"],
        "release_versions": data.get("release_versions", {}),
    }


@lru_cache(maxsize=1)
def _clinical_index() -> dict[str, Any]:
    """Residue counts per canonical IMGT position per chain type, from the local
    clinical-stage library. Uses the numbering already stored in the dataset."""
    from .reference_data import load
    data = load('clinical_frameworks')
    index: dict[str, dict[str, Counter]] = {}
    totals: dict[str, dict[str, int]] = {}
    drugs: dict[str, dict[str, dict[str, list]]] = {}
    for ref in data["references"]:
        numbering = (ref.get("numbering") or {}).get("imgt")
        if not numbering:
            continue
        chain_type = numbering.get("chain_type")
        if not chain_type:
            continue
        per_type = index.setdefault(chain_type, {})
        per_total = totals.setdefault(chain_type, {})
        per_drug = drugs.setdefault(chain_type, {})
        for residue in numbering["residues"]:
            key = _canonical(residue)
            aa = residue["aa"]
            per_type.setdefault(key, Counter())[aa] += 1
            per_total[key] = per_total.get(key, 0) + 1
            if ref.get("stage") == "Approved":
                per_drug.setdefault(key, {}).setdefault(aa, []).append(ref["drug"])
    return {
        "counts": index, "totals": totals, "approved_examples": drugs,
        "source_url": data.get("source_url"), "snapshot": data.get("source_snapshot"),
        "coverage": data.get("coverage", {}),
    }


_LOCUS = {"H": "IGH", "K": "IGK", "L": "IGL"}


def _germline_frequency(chain_type: str, imgt_key: str, aa: str) -> dict[str, Any]:
    index = _germline_index()
    locus = _LOCUS.get(chain_type)
    counts = index["counts"].get(locus, {}).get(imgt_key)
    total = index["totals"].get(locus, {}).get(imgt_key, 0)
    if not locus or not counts or not total:
        return {"status": "unavailable", "percent": None, "observed": None, "total": 0,
                "reason": "Position is outside local functional human V-allele coverage "
                          "(expected for most of CDR3 and FR4).",
                "source": index["source"], "release": index["release_versions"].get(locus)}
    return {"status": "available", "percent": round(100 * counts.get(aa, 0) / total, 1),
            "observed": counts.get(aa, 0), "total": total,
            "most_common": counts.most_common(3),
            "source": index["source"], "release": index["release_versions"].get(locus)}


def _clinical_reference_note() -> dict[str, Any]:
    """Provenance for the clinical column, or why it is absent."""
    from .reference_data import ReferenceDataMissing
    try:
        index = _clinical_index()
    except ReferenceDataMissing as exc:
        return {"status": "unavailable", "source_url": None, "snapshot": None,
                "coverage": {}, "reason": str(exc).split("\n")[0],
                "build_command": next((line.strip() for line in str(exc).splitlines()
                                       if line.strip().startswith("python ")), None)}
    return {"status": "available",
            **{k: index[k] for k in ("source_url", "snapshot", "coverage")}}


def _clinical_frequency(chain_type: str, imgt_key: str, aa: str) -> dict[str, Any]:
    # The clinical-stage library is not redistributed, so a fresh clone does
    # not have it until it is built. It informs one column of the ranking;
    # conservativeness, germline precedent and the developability delta do not
    # use it. Losing the whole feature over one missing column would be the
    # wrong trade, so its absence degrades to "unavailable" exactly as an
    # uncovered position does.
    from .reference_data import ReferenceDataMissing
    try:
        index = _clinical_index()
    except ReferenceDataMissing as exc:
        return {"status": "unavailable", "percent": None, "observed": None, "total": 0,
                "reason": str(exc).split("\n")[0], "snapshot": None,
                "build_command": next((line.strip() for line in str(exc).splitlines()
                                       if line.strip().startswith("python ")), None)}
    counts = index["counts"].get(chain_type, {}).get(imgt_key)
    total = index["totals"].get(chain_type, {}).get(imgt_key, 0)
    if not counts or not total:
        return {"status": "unavailable", "percent": None, "observed": None, "total": 0,
                "reason": "Position is not represented in the local clinical-stage library.",
                "snapshot": index["snapshot"]}
    examples = index["approved_examples"].get(chain_type, {}).get(imgt_key, {}).get(aa, [])
    return {"status": "available", "percent": round(100 * counts.get(aa, 0) / total, 1),
            "observed": counts.get(aa, 0), "total": total,
            "most_common": counts.most_common(3),
            "approved_examples": sorted(set(examples))[:5],
            "approved_example_count": len(set(examples)),
            "snapshot": index["snapshot"]}


# --- chemistry definitions ---------------------------------------------------
# Which offsets within a flagged motif can be edited to remove that chemistry,
# and a short description of what editing each offset accomplishes. The actual
# replacement residues are enumerated and filtered, never taken from a table.
_EDITABLE: dict[str, Any] = {
    "Asparagine deamidation motif": {
        0: "Replace the labile Asn.",
        1: "Keep the Asn and break the motif at the +1 partner.",
    },
    "Aspartate isomerization motif": {
        0: "Replace the isomerization-prone Asp.",
        1: "Keep the Asp and break the motif at the +1 partner.",
    },
    "N-linked glycosylation sequon": {
        0: "Replace the sequon Asn; also removes deamidation at the same residue.",
        2: "Break the sequon at the +2 Ser/Thr hydroxyl.",
    },
    "Oxidation-susceptible residue": {
        0: "Replace the oxidation-susceptible residue.",
    },
    "Asp-Pro cleavage motif": {
        0: "Replace the Asp of the acid-labile Asp-Pro bond.",
        1: "Replace the Pro and retain the Asp.",
    },
    "Noncanonical cysteine": {
        0: "Replace the noncanonical cysteine.",
    },
}

# Chemistries deliberately not offered a substitution.
_MONITOR_ONLY = {
    "N-terminal pyroglutamate formation":
        "N-terminal pyroglutamate is normally controlled by process and analytics rather "
        "than by engineering, and the flag is conditional on this being the true free "
        "mature N-terminus. No substitution is proposed.",
}

# Severity used only to order which sites are offered first.
_PRIORITY_CHEMISTRY = {
    "Asparagine deamidation motif": 3,
    "Aspartate isomerization motif": 3,
    "N-linked glycosylation sequon": 3,
    "Noncanonical cysteine": 2,
    "Asp-Pro cleavage motif": 2,
    "Oxidation-susceptible residue": 1,
}


def _site_key(site: dict) -> tuple:
    return (site["kind"], site["start"], site["motif"])


def _region_at(numbering: dict, position: int) -> str | None:
    for residue in numbering.get("residues", []):
        if residue["sequence_position"] == position:
            return residue.get("region")
    return None


def _imgt_key_at(canonical: dict, position: int) -> str | None:
    for residue in canonical.get("residues", []):
        if residue["sequence_position"] == position:
            return _canonical(residue)
    return None


def _exposure_at(structure: dict | None, chain: str, position: int) -> dict[str, Any]:
    """Join exposure by input sequence position, including non-PTM residues.

    PTM rows may carry conservative ensemble exposure, so valid values there
    take precedence over the primary-model residue table. A motif need not
    exist for a mapped residue's surface exposure to be useful for engineering.
    """
    if not isinstance(structure, dict):
        return {"status": "unavailable", "reason": "No structure analysis in this record."}
    for collection, source in (("ptm_sites", "structural_ptm"), ("residues", "mapped_residue")):
        for site in structure.get(collection, []):
            if site.get("chain") != chain or site.get("sequence_position") != position:
                continue
            value = site.get("relative_sidechain_sasa")
            if (isinstance(value, (int, float)) and not isinstance(value, bool)
                    and math.isfinite(value) and value >= 0):
                return {"status": "available", "relative_sidechain_sasa": round(value, 3),
                        "exposure_band": "low" if value < .10 else "intermediate" if value < .30 else "high",
                        "source": source, "pdb_number": site.get("pdb_number")}
    return {"status": "unavailable",
            "reason": "No finite exposure value is available for this mapped chain and sequence position."}


def propose_for_chain(
    sequence: str,
    canonical: dict,
    selected: dict,
    *,
    chain_label: str = "H",
    cdr_only: bool = True,
    max_per_site: int = 5,
    structure: dict | None = None,
    preserve: set[int] | None = None,
    all_chains: dict[str, str] | None = None,
    paratope: dict | None = None,
    ph: float = 7.4,
) -> dict[str, Any]:
    """Propose ranked substitutions for the PTM sites of one chain.

    ``canonical`` must be IMGT numbering; ``selected`` is the user's chosen
    convention and decides which sites count as CDR. With ``cdr_only`` the
    framework sites are still counted and reported, but no substitution is
    enumerated for them.
    """
    baseline_sites = ptm_sites(sequence, canonical)
    baseline_keys = {_site_key(s) for s in baseline_sites}
    blosum = _blosum62()
    chain_type = canonical.get("chain_type")

    out_sites: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for site in baseline_sites:
        kind, start, motif = site["kind"], site["start"], site["motif"]
        region = _region_at(selected, start)
        in_cdr = bool(region and region.startswith("CDR"))

        if kind in _MONITOR_ONLY:
            skipped.append({"chain": chain_label, "position": start, "kind": kind,
                            "motif": motif, "region": region,
                            "reason": _MONITOR_ONLY[kind]})
            continue
        if cdr_only and not in_cdr:
            skipped.append({"chain": chain_label, "position": start, "kind": kind,
                            "motif": motif, "region": region,
                            "reason": "Framework site; substitution enumeration limited to "
                                      "CDR sites by the cdr_only setting."})
            continue
        editable = _EDITABLE.get(kind)
        if not editable:
            skipped.append({"chain": chain_label, "position": start, "kind": kind,
                            "motif": motif, "region": region,
                            "reason": "No substitution strategy is defined for this chemistry."})
            continue

        candidates: list[dict[str, Any]] = []
        for offset, strategy in editable.items():
            position = start + offset
            if position < 1 or position > len(sequence):
                continue
            original = sequence[position - 1]
            edit_region = _region_at(selected, position)
            imgt_key = _imgt_key_at(canonical, position)

            for replacement in AA:
                if replacement == original:
                    continue
                mutated = sequence[:position - 1] + replacement + sequence[position:]
                # Renumbering every candidate would be prohibitive and the length is
                # unchanged, so the baseline numbering is reused for the rescan. This
                # keeps the noncanonical-cysteine test anchored to the same positions.
                rescanned = ptm_sites(mutated, canonical)
                rescanned_keys = {_site_key(s) for s in rescanned}

                if _site_key(site) in rescanned_keys:
                    continue  # target chemistry survives; not a fix

                introduced = [s for s in rescanned if _site_key(s) not in baseline_keys]
                removed = [s for s in baseline_sites if _site_key(s) not in rescanned_keys]

                physchem = _physchem_distance(original, replacement)
                germline = _germline_frequency(chain_type, imgt_key, replacement) if imgt_key else \
                    {"status": "unavailable", "percent": None,
                     "reason": "Residue could not be mapped to an IMGT position."}
                clinical = _clinical_frequency(chain_type, imgt_key, replacement) if imgt_key else \
                    {"status": "unavailable", "percent": None,
                     "reason": "Residue could not be mapped to an IMGT position."}

                exposure = _exposure_at(structure, chain_label, position)
                sasa = exposure.get("relative_sidechain_sasa")
                try:
                    delta = mutation_effect(
                        all_chains or {chain_label: sequence}, chain_label, position,
                        replacement, ph=ph, relative_sidechain_sasa=sasa,
                        in_cdr=bool(edit_region and str(edit_region).startswith("CDR")))
                except (ValueError, KeyError) as exc:
                    delta = {"status": "unavailable", "reason": str(exc)}

                candidates.append({
                    "chain": chain_label,
                    "edit_position": position,
                    "edit_imgt_number": imgt_key and f"{chain_type}{imgt_key}",
                    "edit_region": edit_region,
                    "original_aa": original,
                    "replacement_aa": replacement,
                    "mutation": f"{original}{position}{replacement}",
                    "variant_spec": f"{chain_label}:{original}{position}{replacement}",
                    "strategy": strategy,
                    "motif_before": motif,
                    "motif_after": mutated[start - 1:start - 1 + len(motif)],
                    "blosum62": blosum.get((original, replacement)) if blosum else None,
                    "physicochemical": physchem,
                    "germline_precedent": germline,
                    "clinical_precedent": clinical,
                    "removes_chemistries": sorted({s["kind"] for s in removed}),
                    "removes_site_count": len(removed),
                    "introduces": [
                        {"kind": s["kind"], "position": s["start"], "motif": s["motif"],
                         "region": _region_at(selected, s["start"])}
                        for s in introduced
                    ],
                    "clean": not introduced,
                    "developability_delta": delta,
                    "paratope": paratope_context(paratope, chain_label, position),
                })

        ranked = _rank(candidates, kind, preserve)
        selected_proposals = _select(ranked, max_per_site)
        out_sites.append({
            "chain": chain_label,
            "position": start,
            "imgt_number": (lambda k: k and f"{chain_type}{k}")(_imgt_key_at(canonical, start)),
            "kind": kind,
            "motif": motif,
            "region": region,
            "in_cdr": in_cdr,
            "exposure": _exposure_at(structure, chain_label, start),
            "chemistry_priority": _PRIORITY_CHEMISTRY.get(kind, 0),
            "candidate_count": len(ranked),
            "clean_candidate_count": sum(c["clean"] for c in ranked),
            "standard_choice_count": sum(c["standard_choice"] for c in ranked),
            "proposals": selected_proposals,
            "withheld_candidate_count": max(0, len(ranked) - len(selected_proposals)),
        })

    out_sites.sort(key=lambda s: (not s["in_cdr"], -s["chemistry_priority"], s["position"]))
    return {
        "chain": chain_label,
        "chain_type": chain_type,
        "scheme": selected.get("scheme"),
        "cdr_only": cdr_only,
        "preserved_positions": sorted(preserve) if preserve else [],
        "sites": out_sites,
        "skipped_sites": skipped,
        "ranking_method": _RANKING_NOTE,
        "limitation": LIMITATION,
        "references": {
            "germline": {k: _germline_index()[k] for k in ("source", "source_url", "release_versions")},
            "clinical": _clinical_reference_note(),
            "substitution_matrix": "BLOSUM62 (Biopython)" if blosum else "unavailable",
        },
    }


_RANKING_NOTE = (
    "Candidates are ordered by: (1) avoids caller-preserved positions; (2) introduces "
    "no new PTM chemistry; (3) higher selection score, an equal-weight blend of residue "
    "similarity and clinical reference frequency divided by 100. Residue similarity is "
    "the mean of normalised BLOSUM62 and physicochemical similarity; if BLOSUM62 is "
    "unavailable it uses physicochemical similarity alone. If clinical coverage is "
    "unavailable, selection uses residue similarity alone and reports that fallback; "
    "an observed clinical frequency of zero remains zero, not missing. Germline "
    "frequency breaks ties, followed by the developability penalty. The top five are "
    "a modeling shortlist, not a prediction of folding stability or retained binding. "
    "Clinical frequency counts reference chain records covering the same IMGT position "
    "and chain type, not successful tests of this mutation in this antibody. "
    "Germline and clinical frequencies are observational; component values and weights "
    "are reported for review."
)

# The conventional first-choice substitution for each chemistry, as used in
# biopharmaceutical sequence optimisation. These are annotated and guaranteed a
# place in the returned set so that a well-established option is never displaced
# by a higher-frequency but less conservative alternative. Being conventional is
# not evidence that it is correct for a given antibody.
_STANDARD_FIXES: dict[str, dict[str, str]] = {
    "Asparagine deamidation motif": {
        "N->Q": "Isosteric replacement of the labile Asn; the conventional deamidation fix.",
        "N->A": "Removes the Asn side-chain amide; small and frequently tolerated in CDRs.",
        "N->S": "Removes the amide while retaining a small polar side chain.",
    },
    "Aspartate isomerization motif": {
        "D->E": "Charge-conservative replacement of the isomerization-prone Asp.",
        "S->A": "Breaks the motif at the +1 hydroxyl and retains the Asp.",
        "T->A": "Breaks the motif at the +1 hydroxyl and retains the Asp.",
        "G->A": "Breaks the motif at the +1 Gly and retains the Asp.",
    },
    "N-linked glycosylation sequon": {
        "N->Q": "Removes the sequon Asn and any deamidation at the same residue.",
        "S->A": "Breaks the sequon at the +2 hydroxyl.",
        "T->A": "Breaks the sequon at the +2 hydroxyl.",
    },
    "Oxidation-susceptible residue": {
        "M->L": "Isosteric hydrophobic replacement removing the oxidisable thioether.",
        "W->F": "Retains aromaticity while removing the readily oxidised indole.",
    },
    "Asp-Pro cleavage motif": {
        "D->E": "Charge-conservative; removes the acid-labile Asp side of the bond.",
        "P->A": "Removes the proline and retains the Asp.",
    },
    "Noncanonical cysteine": {
        "C->S": "Isosteric replacement removing the thiol.",
        "C->A": "Removes the thiol with a minimal side chain.",
    },
}

# BLOSUM62 values span roughly -4..11; normalise onto 0..1 for blending.
_BLOSUM_MIN, _BLOSUM_MAX = -4.0, 11.0


def _precedent_score(candidate: dict) -> float | None:
    values = [candidate[key]["percent"] for key in ("germline_precedent", "clinical_precedent")
              if candidate[key].get("percent") is not None]
    return sum(values) / len(values) if values else None


def _conservativeness(candidate: dict) -> dict[str, Any]:
    """Blend normalised BLOSUM62 with the complement of physicochemical distance."""
    blosum = candidate.get("blosum62")
    distance = candidate["physicochemical"]["distance"]
    similarity = 1.0 - min(1.0, max(0.0, distance))
    if blosum is None:
        score = similarity
        components = {"blosum_normalised": None, "physicochemical_similarity": round(similarity, 4)}
    else:
        blosum_norm = (max(_BLOSUM_MIN, min(_BLOSUM_MAX, float(blosum))) - _BLOSUM_MIN) / \
                      (_BLOSUM_MAX - _BLOSUM_MIN)
        score = (blosum_norm + similarity) / 2.0
        components = {"blosum_normalised": round(blosum_norm, 4),
                      "physicochemical_similarity": round(similarity, 4)}
    components["score"] = round(score, 4)
    return components



def _selection_score(candidate: dict) -> dict[str, Any]:
    """An explicit shortlist heuristic, never a structural stability estimate.

    Missing clinical coverage is not evidence of absence. Available zero counts
    do contribute zero; absent, malformed, or nonfinite frequencies do not.
    """
    similarity = candidate["conservativeness"]["score"]
    clinical = candidate.get("clinical_precedent") or {}
    percent = clinical.get("percent")
    available = (clinical.get("status") == "available"
                 and isinstance(percent, (int, float)) and not isinstance(percent, bool)
                 and math.isfinite(percent) and 0 <= percent <= 100)
    frequency = float(percent) / 100 if available else None
    return {
        "score": round((similarity + frequency) / 2 if available else similarity, 4),
        "residue_similarity": similarity,
        "clinical_frequency": frequency,
        "weights": {"residue_similarity": .5 if available else 1.0,
                    "clinical_frequency": .5 if available else 0.0},
        "basis": "similarity_clinical_blend" if available else "similarity_only",
        "reason": ("Equal-weight residue similarity and clinical reference frequency."
                   if available else "Clinical frequency unavailable; residue similarity only."),
    }


def _standard_annotation(kind: str, original: str, replacement: str) -> str | None:
    return _STANDARD_FIXES.get(kind, {}).get(f"{original}->{replacement}")


def _rank(candidates: list[dict], kind: str, preserve: set[int] | None = None) -> list[dict]:
    """Order candidates. ``preserve`` holds input positions the caller wants kept.

    A position is normally preserved because it is a known paratope contact, which
    this module has no way to determine itself. Candidates that edit a preserved
    position are ranked last rather than removed, so the trade-off stays visible.
    """
    preserve = preserve or set()
    for candidate in candidates:
        candidate["conservativeness"] = _conservativeness(candidate)
        candidate["selection_score"] = _selection_score(candidate)
        note = _standard_annotation(kind, candidate["original_aa"], candidate["replacement_aa"])
        candidate["standard_choice"] = note is not None
        candidate["standard_choice_note"] = note
        candidate["edits_preserved_position"] = candidate["edit_position"] in preserve
        candidate["developability_penalty"] = _delta_penalty(candidate)

    def key(candidate: dict):
        germline = candidate["germline_precedent"].get("percent")
        return (
            candidate["edits_preserved_position"],
            not candidate["clean"],
            -candidate["selection_score"]["score"],
            -(germline if germline is not None else -1.0),
            candidate["developability_penalty"],
            candidate["replacement_aa"],
        )

    ranked = sorted(candidates, key=key)
    for index, candidate in enumerate(ranked, 1):
        precedent = _precedent_score(candidate)
        candidate["rank"] = index
        candidate["combined_precedent_percent"] = round(precedent, 1) if precedent is not None else None
    return ranked


def _select(ranked: list[dict], max_per_site: int) -> list[dict]:
    """Top-ranked candidates, with any conventional choice guaranteed inclusion.

    A standard choice that falls outside the top slice is appended rather than
    replacing a better-ranked candidate, so the list can exceed ``max_per_site``
    by the number of conventional options recovered.
    """
    selected = ranked[:max_per_site]
    chosen = {c["mutation"] for c in selected}
    recovered = [c for c in ranked[max_per_site:]
                 if c["standard_choice"] and c["mutation"] not in chosen]
    for candidate in recovered:
        candidate["included_as"] = "conventional choice recovered below the top ranks"
    return selected + recovered


def propose(
    chains: dict[str, dict],
    *,
    cdr_only: bool = True,
    max_per_site: int = 5,
    structure: dict | None = None,
    preserve: dict[str, set[int]] | None = None,
    paratope: dict | None = None,
    auto_preserve_paratope: bool = True,
    ph: float = 7.4,
) -> dict[str, Any]:
    """Run the proposer across an analysed chain set.

    ``chains`` follows the analysis record shape: ``{"H": {"sequence": ...,
    "numbering": ..., "numbering_imgt": ...}}``. A chain whose IMGT numbering is
    missing is reported in ``skipped_chains`` rather than silently dropped.
    """
    sequences = {label: (chain["sequence"] if isinstance(chain, dict) else chain)
                 for label, chain in chains.items()}
    resolved_preserve, preserve_source = _resolve_preserve(preserve, paratope,
                                                           auto_preserve_paratope)
    results: dict[str, Any] = {"chains": {}, "skipped_chains": [],
                               "limitation": LIMITATION, "ranking_method": _RANKING_NOTE,
                               "preserved_positions": {k: sorted(v) for k, v in
                                                       resolved_preserve.items()},
                               "preserved_positions_source": preserve_source,
                               "paratope": (paratope or {}).get("summary")
                               if isinstance(paratope, dict) else None}
    for label, chain in chains.items():
        sequence = chain["sequence"] if isinstance(chain, dict) else chain
        canonical = chain.get("numbering_imgt") or chain.get("numbering") if isinstance(chain, dict) else None
        selected = chain.get("numbering") if isinstance(chain, dict) else None
        if not canonical or canonical.get("scheme") != "imgt":
            results["skipped_chains"].append(
                {"chain": label, "reason": "Canonical IMGT numbering is required and was not available."})
            continue
        results["chains"][label] = propose_for_chain(
            sequence, canonical, selected or canonical, chain_label=label,
            cdr_only=cdr_only, max_per_site=max_per_site, structure=structure,
            preserve=resolved_preserve.get(label), all_chains=sequences,
            paratope=paratope, ph=ph)
    totals = [site for chain in results["chains"].values() for site in chain["sites"]]
    results["summary"] = {
        "sites_with_proposals": sum(bool(s["proposals"]) for s in totals),
        "cdr_sites": sum(s["in_cdr"] for s in totals),
        "sites_without_clean_candidate": sum(not s["clean_candidate_count"] for s in totals),
        "total_candidates": sum(s["candidate_count"] for s in totals),
    }
    return results


def variant_specs(proposals: dict[str, Any], *, top_n: int = 1,
                  name_prefix: str = "auto") -> list[str]:
    """Render top-ranked proposals as variant-comparison specification lines.

    Each line is a single substitution so the comparison attributes any change to
    one edit. Combinations are intentionally not auto-generated: the existing
    comparison runs supplied singles before combinations, and stacking edits
    without testing each one first confounds the readout.
    """
    lines: list[str] = []
    for chain in proposals.get("chains", {}).values():
        for site in chain["sites"]:
            for candidate in site["proposals"][:top_n]:
                label = f"{name_prefix}_{candidate['chain']}{candidate['mutation']}"
                lines.append(f"{label} | {candidate['variant_spec']}")
    return lines


def _delta_penalty(candidate: dict) -> float:
    """The bounded surface-consequence tie-break, or 0 when not computable."""
    delta = candidate.get("developability_delta")
    if isinstance(delta, dict):
        penalty = delta.get("penalty")
        if isinstance(penalty, dict) and isinstance(penalty.get("score"), (int, float)):
            return float(penalty["score"])
    return 0.0


def paratope_context(paratope: dict | None, chain: str, position: int) -> dict[str, Any]:
    """Predicted antigen-contact context for the residue being edited."""
    if not isinstance(paratope, dict):
        return {"status": "unavailable",
                "reason": "No paratope prediction was supplied for this run."}
    residues = (paratope.get("by_chain") or {}).get(chain) or {}
    record = residues.get(position) or residues.get(str(position))
    if not record:
        return {"status": "unavailable",
                "reason": "Residue was not covered by the paratope prediction."}
    return {
        "status": "available",
        "method": paratope.get("method"),
        "probability": record.get("probability"),
        "band": record.get("band"),
        "likely_contact": record.get("likely_contact"),
        "note": "Predicted antigen contact, not a measured epitope. A likely "
                "contact position is reported so an equivalent fix elsewhere can "
                "be preferred; it is not a prohibition.",
    }


def _resolve_preserve(
    preserve: dict[str, set[int]] | None,
    paratope: dict | None,
    auto: bool,
) -> tuple[dict[str, set[int]], dict[str, Any]]:
    """Combine user-specified preserved positions with predicted paratope contacts.

    A user-supplied set always wins for that chain: a measured or known contact
    is better evidence than a prediction. Where the user said nothing and a
    prediction exists, likely contacts are preserved automatically, which pushes
    candidates that edit them to the end of the ranking rather than hiding them.
    """
    resolved = {label: set(positions) for label, positions in (preserve or {}).items()}
    source: dict[str, Any] = {
        label: "supplied by caller" for label in resolved
    }
    if auto and isinstance(paratope, dict):
        for label, residues in (paratope.get("by_chain") or {}).items():
            if label in resolved and resolved[label]:
                source[label] = "supplied by caller; paratope prediction not merged"
                continue
            contacts = {int(position) for position, record in residues.items()
                        if record.get("likely_contact")}
            if contacts:
                resolved[label] = contacts
                source[label] = (f"predicted paratope ({paratope.get('method')}), "
                                 f"{len(contacts)} likely-contact positions")
    return resolved, source


def propose_at_position(
    chains: dict[str, dict],
    chain_label: str,
    position: int,
    *,
    max_candidates: int = 8,
    structure: dict | None = None,
    paratope: dict | None = None,
    preserve: set[int] | None = None,
    ph: float = 7.4,
) -> dict[str, Any]:
    """Rank substitutions at a position the user chose, liability or not.

    The automatic path only enumerates positions where a chemical liability was
    detected. A user often wants the same analysis somewhere else -- a residue
    they suspect, one an assay implicated, or a position they are considering
    for an unrelated reason. The scoring is identical: residue similarity and
    clinical frequency are equally weighted after safety ordering, with germline
    precedent and developability consequence as tiebreakers. Predicted contact
    is reported separately.

    The one difference is the test each candidate has to pass. A liability fix
    must remove its chemistry; here there is no chemistry to remove, so a
    candidate is judged on what it *introduces* instead. ``clean`` therefore
    means "creates no new liability", not "removed the old one".
    """
    chain = chains.get(chain_label)
    if chain is None:
        raise ValueError(f"No chain {chain_label!r} was supplied. "
                         f"Available: {', '.join(sorted(chains)) or 'none'}.")
    sequence = chain["sequence"]
    if not 1 <= position <= len(sequence):
        raise ValueError(f"Position {position} is outside chain {chain_label}, "
                         f"which has {len(sequence)} residues.")

    canonical = chain.get("numbering_imgt") or chain.get("numbering")
    selected = chain.get("numbering") or canonical
    if not canonical or canonical.get("scheme") != "imgt":
        raise ValueError("Canonical IMGT numbering is required to rank substitutions.")

    original = sequence[position - 1]
    chain_type = canonical.get("chain_type")
    imgt_key = _imgt_key_at(canonical, position)
    region = _region_at(selected, position)
    baseline_sites = ptm_sites(sequence, canonical)
    baseline_keys = {_site_key(s) for s in baseline_sites}
    blosum = _blosum62()
    exposure = _exposure_at(structure, chain_label, position)
    sasa = exposure.get("relative_sidechain_sasa")

    candidates: list[dict[str, Any]] = []
    for replacement in AA:
        if replacement == original:
            continue
        mutated = sequence[:position - 1] + replacement + sequence[position:]
        rescanned = ptm_sites(mutated, canonical)
        rescanned_keys = {_site_key(s) for s in rescanned}
        introduced = [s for s in rescanned if _site_key(s) not in baseline_keys]
        removed = [s for s in baseline_sites if _site_key(s) not in rescanned_keys]

        try:
            delta = mutation_effect(
                {k: v["sequence"] for k, v in chains.items()}, chain_label, position,
                replacement, ph=ph, relative_sidechain_sasa=sasa,
                in_cdr=bool(region and str(region).startswith("CDR")))
        except (ValueError, KeyError) as exc:
            delta = {"status": "unavailable", "reason": str(exc)}

        candidates.append({
            "chain": chain_label,
            "edit_position": position,
            "edit_imgt_number": imgt_key and f"{chain_type}{imgt_key}",
            "edit_region": region,
            "original_aa": original,
            "replacement_aa": replacement,
            "mutation": f"{original}{position}{replacement}",
            "variant_spec": f"{chain_label}:{original}{position}{replacement}",
            "strategy": "User-requested position; no target chemistry to remove.",
            "blosum62": blosum.get((original, replacement)) if blosum else None,
            "physicochemical": _physchem_distance(original, replacement),
            "germline_precedent": (_germline_frequency(chain_type, imgt_key, replacement)
                                   if imgt_key else
                                   {"status": "unavailable", "percent": None,
                                    "reason": "Residue could not be mapped to an IMGT position."}),
            "clinical_precedent": (_clinical_frequency(chain_type, imgt_key, replacement)
                                   if imgt_key else
                                   {"status": "unavailable", "percent": None,
                                    "reason": "Residue could not be mapped to an IMGT position."}),
            "removes_chemistries": sorted({s["kind"] for s in removed}),
            "removes_site_count": len(removed),
            "introduces": [
                {"kind": s["kind"], "position": s["start"], "motif": s["motif"],
                 "region": _region_at(selected, s["start"])}
                for s in introduced
            ],
            "clean": not introduced,
            "developability_delta": delta,
            "paratope": paratope_context(paratope, chain_label, position),
        })

    ranked = _rank(candidates, "user_position", preserve)
    contact = paratope_context(paratope, chain_label, position)
    liabilities = [{"kind": s["kind"], "motif": s["motif"]}
                   for s in baseline_sites
                   if s["start"] <= position <= s["start"] + len(s["motif"]) - 1]

    return {
        "chain": chain_label,
        "position": position,
        "original_aa": original,
        "imgt_number": imgt_key and f"{chain_type}{imgt_key}",
        "region": region,
        "exposure": exposure,
        "paratope": contact,
        "existing_liabilities": liabilities,
        "candidates": ranked[:max_candidates],
        "candidate_count": len(ranked),
        "basis": _RANKING_NOTE,
        "top_model_candidates": [c["variant_spec"] for c in ranked[:5]],
        "model_selection_method": _RANKING_NOTE,
        "caveat": ("This position carries no detected liability, so there is nothing to fix; "
                   "these are candidates prioritized for structural evaluation, not a "
                   "recommendation to make one. None of this predicts retained binding, "
                   "potency, expression or stability."
                   if not liabilities else
                   "This position overlaps a detected liability, which the automatic "
                   "proposals already cover; the ranking here is over all 19 substitutions "
                   "rather than only those that remove the chemistry."),
    }
