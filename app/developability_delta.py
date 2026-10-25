"""What a single substitution does to charge, pI and aggregation proxies.

A PTM fix is chosen for chemistry, but it also moves the molecule's surface
properties. Replacing an Asn with an Asp removes a deamidation motif and adds a
negative charge; replacing it with Leu removes the motif and raises a
parent-exposure-weighted hydrophobicity proxy. Mutant surface area is not calculated here. Neither effect is visible in a BLOSUM score, so this module
reports both for every candidate.

Everything here is computed with the same functions as the parent-molecule
descriptors in ``analysis``, so a delta is comparable to the value it moves.

Nothing here predicts aggregation, viscosity or clearance. These are the same
screening proxies used elsewhere in this project, differenced.
"""
from __future__ import annotations

from functools import lru_cache
import math
from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .analysis import HYDROPATHY, _window_scores

POSITIVE = "KR"
NEGATIVE = "DE"
# Charge class per residue at near-neutral pH. His is listed separately because
# its protonation is genuinely ambiguous at pH 7.4.
CHARGE_CLASS = {**{a: "positive" for a in POSITIVE},
                **{a: "negative" for a in NEGATIVE},
                "H": "weakly positive"}

PATCH_HYDROPATHY = 1.6
PATCH_MAX_CHARGED = 1

LIMITATION = (
    "Deltas in computed screening descriptors, not predicted behaviour. A pI or "
    "charge shift does not establish a change in clearance, and a hydrophobic "
    "patch count is a sequence window rule, not an aggregation prediction. "
    "Single-residue effects on a whole IgG are small by construction: these "
    "values rank candidates against each other, and do not estimate the "
    "molecule's properties."
)


@lru_cache(maxsize=4096)
def _formal_charge(sequence: str) -> int:
    """(K+R) - (D+E), the convention this project's Fv charge descriptors use."""
    return sum(a in POSITIVE for a in sequence) - sum(a in NEGATIVE for a in sequence)


@lru_cache(maxsize=2048)
def _patch_count(sequence: str) -> int:
    return len(_patches(sequence))


def _patches(sequence: str) -> list[dict[str, Any]]:
    return [w for w in _window_scores(sequence)
            if w["mean_hydropathy"] >= PATCH_HYDROPATHY
            and w["charged_residues"] <= PATCH_MAX_CHARGED]


# Every candidate at a site shares the same parent chains, and the unedited
# partner chain is identical across every candidate in the run, so these are
# computed once per distinct sequence rather than once per candidate.
@lru_cache(maxsize=2048)
def _chain_descriptors(sequence: str, ph: float) -> dict[str, Any]:
    protein = ProteinAnalysis(sequence)
    return {
        "isoelectric_point": round(protein.isoelectric_point(), 2),
        "net_charge_at_ph": round(protein.charge_at_pH(ph), 2),
        "formal_charge": _formal_charge(sequence),
        "gravy": round(protein.gravy(), 3),
        "hydrophobic_patch_count": _patch_count(sequence),
    }


@lru_cache(maxsize=2048)
def _paired_descriptors(combined: str, ph: float) -> dict[str, float]:
    protein = ProteinAnalysis(combined)
    return {"net_charge_at_ph": round(protein.charge_at_pH(ph), 2),
            "isoelectric_point": round(protein.isoelectric_point(), 2)}


def _delta(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in before.items():
        new = after[key]
        change = new - value
        out[key] = {
            "before": value,
            "after": new,
            "delta": round(change, 3) if isinstance(change, float) else change,
        }
    return out


def mutation_effect(
    chains: dict[str, str],
    chain_label: str,
    position: int,
    replacement: str,
    *,
    ph: float = 7.4,
    relative_sidechain_sasa: float | None = None,
    in_cdr: bool = False,
) -> dict[str, Any]:
    """Descriptor deltas for one substitution, in the context of both chains.

    ``chains`` maps label to sequence and should carry every chain of the Fv, so
    the paired-chain charge descriptors can move. ``position`` is one-based in
    the chain being edited. ``relative_sidechain_sasa`` is the edited residue's
    exposure when a structure was predicted; without it the exposure-weighted
    terms are reported as unavailable rather than assumed buried.
    """
    parent = chains.get(chain_label, "")
    if not 1 <= position <= len(parent):
        raise ValueError("Edit position lies outside the chain.")
    original = parent[position - 1]

    mutated_chain = parent[:position - 1] + replacement + parent[position:]
    mutated_chains = {**chains, chain_label: mutated_chain}

    per_chain: dict[str, Any] = {}
    for label, sequence in chains.items():
        before = _chain_descriptors(sequence, ph)
        after = _chain_descriptors(mutated_chains[label], ph)
        per_chain[label] = _delta(before, after)

    # Paired-chain descriptors. fv_csp is the product of per-chain formal
    # charges, the convention used by the charge-symmetry descriptor elsewhere.
    def fv(source: dict[str, str]) -> dict[str, Any]:
        heavy = source.get("H", "")
        light = next((s for label, s in source.items() if label != "H"), "")
        heavy_charge, light_charge = _formal_charge(heavy), _formal_charge(light)
        combined = "".join(source.values())
        paired = _paired_descriptors(combined, ph) if combined else {}
        return {
            "fv_formal_charge": heavy_charge + light_charge,
            "fv_charge_symmetry_parameter": heavy_charge * light_charge,
            "fv_net_charge_at_ph": paired.get("net_charge_at_ph", 0.0),
            "fv_isoelectric_point": paired.get("isoelectric_point", 0.0),
        }

    fv_delta = _delta(fv(chains), fv(mutated_chains))

    # Hold parent exposure fixed for a screening proxy, not a mutant surface
    # calculation. Positive-part KD avoids treating polar N->S as adding a
    # hydrophobic residue simply because S has a less-negative raw KD value.
    hydropathy_change = HYDROPATHY[replacement] - HYDROPATHY[original]
    positive_hydropathy_change = max(0., HYDROPATHY[replacement]) - max(0., HYDROPATHY[original])
    if (isinstance(relative_sidechain_sasa, (int, float))
            and not isinstance(relative_sidechain_sasa, bool)
            and math.isfinite(relative_sidechain_sasa) and relative_sidechain_sasa >= 0):
        exposed = {
            "status": "available",
            "relative_sidechain_sasa": round(float(relative_sidechain_sasa), 3),
            "hydropathy_delta": round(hydropathy_change, 2),
            # Retained as a raw residue-scale descriptor for API compatibility.
            "exposure_weighted_hydropathy_delta":
                round(hydropathy_change * float(relative_sidechain_sasa), 3),
            "parent_exposure_weighted_positive_hydropathy_delta":
                round(positive_hydropathy_change * float(relative_sidechain_sasa), 3),
            "basis": "parent_exposure_positive_part_kyte_doolittle_proxy",
            "definition": "Parent relative side-chain SASA × [max(0, KD(mutant)) − max(0, KD(parent))].",
            "limitation": "A residue-scale screening proxy at fixed parent exposure; not mutant solvent-accessible area or a spatial patch calculation. Recalculate the mutant structure to assess surface changes.",
        }
    else:
        exposed = {
            "status": "unavailable",
            "hydropathy_delta": round(hydropathy_change, 2),
            "reason": "No mapped side-chain exposure for this residue; "
                      "exposure-weighted terms were not computed.",
        }

    charge_before = CHARGE_CLASS.get(original, "neutral")
    charge_after = CHARGE_CLASS.get(replacement, "neutral")

    flags = _flags(original, replacement, charge_before, charge_after,
                   per_chain[chain_label], fv_delta, exposed, in_cdr)

    return {
        "edited_chain": chain_label,
        "position": position,
        "original_aa": original,
        "replacement_aa": replacement,
        "ph": ph,
        "charge_class": {"before": charge_before, "after": charge_after,
                         "changed": charge_before != charge_after},
        "edited_chain_delta": per_chain[chain_label],
        "all_chain_delta": per_chain,
        "fv_delta": fv_delta,
        "exposed_hydrophobicity": exposed,
        "flags": flags,
        "penalty": _penalty(per_chain[chain_label], exposed, flags),
        "limitation": LIMITATION,
    }


def _flags(original: str, replacement: str, charge_before: str, charge_after: str,
           chain_delta: dict[str, Any], fv_delta: dict[str, Any],
           exposed: dict[str, Any], in_cdr: bool) -> list[dict[str, str]]:
    """Plain-language notes a reviewer should see next to the substitution."""
    flags: list[dict[str, str]] = []

    pi_shift = chain_delta["isoelectric_point"]["delta"]
    if abs(pi_shift) >= 0.3:
        flags.append({
            "kind": "pi_shift",
            "severity": "note" if abs(pi_shift) < 0.6 else "attention",
            "text": f"Shifts this chain's pI by {pi_shift:+.2f} "
                    f"({chain_delta['isoelectric_point']['before']} to "
                    f"{chain_delta['isoelectric_point']['after']}).",
        })

    if charge_before != charge_after:
        if charge_after == "neutral":
            text = f"Removes a {charge_before} charge."
        elif charge_before == "neutral":
            text = f"Adds a {charge_after} charge."
        else:
            text = f"Reverses charge, {charge_before} to {charge_after}."
        flags.append({"kind": "charge_class", "severity": "attention", "text": text})

    csp = fv_delta["fv_charge_symmetry_parameter"]
    if csp["delta"]:
        flags.append({
            "kind": "charge_symmetry",
            "severity": "note",
            "text": f"Fv charge-symmetry parameter moves {csp['delta']:+g} "
                    f"({csp['before']} to {csp['after']}).",
        })

    patches = chain_delta["hydrophobic_patch_count"]["delta"]
    if patches > 0:
        flags.append({
            "kind": "hydrophobic_patch",
            "severity": "attention",
            "text": f"Creates {patches} new hydrophobic sequence window"
                    f"{'s' if patches > 1 else ''} by this screen's rule.",
        })
    elif patches < 0:
        flags.append({
            "kind": "hydrophobic_patch",
            "severity": "good",
            "text": f"Removes {abs(patches)} hydrophobic sequence window"
                    f"{'s' if abs(patches) > 1 else ''}.",
        })

    if exposed.get("status") == "available":
        weighted = exposed["parent_exposure_weighted_positive_hydropathy_delta"]
        if weighted >= 0.5:
            flags.append({
                "kind": "exposed_hydrophobicity",
                "severity": "attention",
                "text": f"Increases the parent-exposure-weighted hydrophobicity proxy "
                        f"({weighted:+.2f}; positive-part Kyte-Doolittle)"
                        + (" in a CDR. " if in_cdr else ". ")
                        + "This is not a measured or modeled spatial patch change; assess mutant coordinates.",
            })
        elif weighted <= -0.5:
            flags.append({
                "kind": "exposed_hydrophobicity",
                "severity": "good",
                "text": f"Decreases the parent-exposure-weighted hydrophobicity proxy "
                        f"({weighted:+.2f}; positive-part Kyte-Doolittle). "
                        "This is not a measured or modeled spatial patch change; assess mutant coordinates.",
            })
        if charge_before == "neutral" and charge_after != "neutral" \
                and exposed["relative_sidechain_sasa"] < 0.15:
            flags.append({
                "kind": "buried_charge",
                "severity": "attention",
                "text": "Places a charged residue at a largely buried position, "
                        "which can be structurally destabilising.",
            })

    if replacement == "P":
        flags.append({"kind": "proline", "severity": "attention",
                      "text": "Proline substitution; backbone conformation and loop "
                              "geometry are likely to change."})
    if replacement == "C":
        flags.append({"kind": "cysteine", "severity": "attention",
                      "text": "Introduces a free cysteine; mispairing and "
                              "conjugation risk require experiment."})
    if replacement == "G" and in_cdr:
        flags.append({"kind": "glycine", "severity": "note",
                      "text": "Glycine in a CDR increases local flexibility."})
    return flags


def _penalty(chain_delta: dict[str, Any], exposed: dict[str, Any],
             flags: list[dict[str, str]]) -> dict[str, Any]:
    """A small, bounded tie-break score. Never outranks conservativeness.

    This exists so that two substitutions which are equally conservative and
    have equal positional precedent can be ordered by surface consequence. It is
    applied after both of those in the ranking, by design: a surface heuristic
    should not displace a conventional, well-precedented fix.
    """
    components = {
        "pi_shift": min(1.0, abs(chain_delta["isoelectric_point"]["delta"]) / 1.0),
        "new_hydrophobic_patch": min(1.0, max(0, chain_delta["hydrophobic_patch_count"]["delta"]) / 2.0),
        "exposed_hydrophobicity": (
            min(1.0, max(0.0, exposed["parent_exposure_weighted_positive_hydropathy_delta"]) / 2.0)
            if exposed.get("status") == "available" else 0.0
        ),
        "structural_risk": min(1.0, 0.5 * sum(
            f["kind"] in {"proline", "cysteine", "buried_charge"} for f in flags)),
    }
    weights = {"pi_shift": 0.25, "new_hydrophobic_patch": 0.3,
               "exposed_hydrophobicity": 0.3, "structural_risk": 0.15}
    score = sum(components[k] * weights[k] for k in weights)
    return {
        "score": round(score, 4),
        "components": {k: round(v, 4) for k, v in components.items()},
        "weights": weights,
        "definition": "Bounded 0-1 tie-break over charge and aggregation-proxy "
                      "deltas. Applied only after conservativeness and positional "
                      "precedent, so it reorders equivalent candidates rather than "
                      "promoting an unconventional one. Not a risk probability.",
    }
