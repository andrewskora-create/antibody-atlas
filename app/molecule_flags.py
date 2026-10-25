"""Molecule-level charge and hydrophobicity flags, mapped to the outcomes they bear on.

pI, net charge and hydrophobicity are not per-residue properties, so they do not
belong in the residue table. They are reported here instead, and each one is
tied only to the outcomes the literature actually supports:

  colloidal stability / aggregation   charge near pI, hydrophobic patches
  non-specific binding               positive charge, exposed hydrophobicity
  clearance                          charge extremes, CDR hydrophobicity
  viscosity                          low net charge, charge asymmetry, hydrophobicity

Immunogenicity is deliberately absent. Charge and hydrophobicity reach it only
through aggregation, and treating a pI value as an immunogenicity signal would
overstate what these descriptors show. The direct sequence-level measure is
non-human peptide content, which is a separate axis.

Two limits apply to everything here and are repeated on every flag:

1. Every cited dataset used full-length IgG. This project computes descriptors
   from variable domains only, so a value is a rank-order proxy and not the
   molecule-level quantity the published analyses thresholded on.
2. Nothing here predicts a clearance rate, a viscosity, or an aggregation rate.
   The flags say where a molecule sits and what to measure.

The band boundaries are this project's screening policy, not published
thresholds. Each flag states that in its own `policy` field.
"""
from __future__ import annotations

from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .analysis import HYDROPATHY

# Typical mAb drug-product formulations sit acidic, around pH 5.5-6.5, so the
# colloidal margin is evaluated there by default rather than at the pH 7.4 used
# for the net-charge descriptor.
DEFAULT_FORMULATION_PH = 6.0
DEFAULT_ASSAY_PH = 7.4

VERDICT_ORDER = {"green": 0, "amber": 1, "red": 2, "unknown": 3}


def _flag(**kw) -> dict[str, Any]:
    kw.setdefault("fv_only", True)
    return kw


def _scale_charge_bands(full_length: bool) -> tuple[float, float]:
    """Charge bands differ by basis.

    An H2L2 has roughly four times the ionizable content of a single Fv, so the
    variable-domain bands would fire on almost every real molecule if reused.
    These remain this project's policy either way.
    """
    return (20.0, 32.0) if full_length else (5.0, 8.0)


def _pi(sequence: str) -> float:
    return round(ProteinAnalysis(sequence).isoelectric_point(), 2)


def _charge(sequence: str, ph: float) -> float:
    return round(ProteinAnalysis(sequence).charge_at_pH(ph), 2)


def _molecular_charge(components: list[tuple[str, int]], ph: float) -> float:
    """Sum physical chains, retaining each chain's N and C terminal groups.

    This retains Biopython's sequence-only ionization model, including its
    unmodified termini assumption. Concatenating chains creates fictitious
    peptide bonds and removes terminal groups from Fv and H2L2 calculations.
    """
    return sum(copies * ProteinAnalysis(sequence).charge_at_pH(ph)
               for sequence, copies in components)


def _molecular_pi(components: list[tuple[str, int]]) -> float:
    """Find the zero of the summed charge, not an average of chain pIs."""
    low, high = 0.0, 14.0
    for _ in range(48):
        middle = (low + high) / 2
        if _molecular_charge(components, middle) > 0:
            low = middle
        else:
            high = middle
    return round((low + high) / 2, 2)


def _cdr_hydropathy(chains: dict[str, dict]) -> dict[str, Any]:
    """Mean hydropathy over IMGT CDR residues only.

    Sharma et al. 2014 associated fast clearance with high hydrophobicity of
    particular CDRs rather than of the whole chain, so the whole-chain GRAVY
    already reported elsewhere is the wrong denominator for this flag.
    """
    residues: list[str] = []
    per_cdr: dict[str, float] = {}
    for label, chain in chains.items():
        numbering = chain.get("numbering_imgt") or chain.get("numbering") or {}
        if numbering.get("scheme") != "imgt":
            continue
        for name in ("CDR1", "CDR2", "CDR3"):
            aas = [r["aa"] for r in numbering.get("residues", []) if r["region"] == name]
            if not aas:
                continue
            residues.extend(aas)
            per_cdr[f"{label}-{name}"] = round(sum(HYDROPATHY[a] for a in aas) / len(aas), 3)
    if not residues:
        return {"mean": None, "per_cdr": {}, "residue_count": 0}
    return {"mean": round(sum(HYDROPATHY[a] for a in residues) / len(residues), 3),
            "per_cdr": per_cdr, "residue_count": len(residues)}


def _titration(sequence: str, low: float = 4.0, high: float = 8.5,
               step: float = 0.25) -> list[dict[str, float]]:
    curve, ph = [], low
    while ph <= high + 1e-9:
        curve.append({"ph": round(ph, 2), "net_charge": _charge(sequence, ph)})
        ph += step
    return curve


def molecule_flags(chains: dict[str, dict], formulation_ph: float = DEFAULT_FORMULATION_PH,
                   assay_ph: float = DEFAULT_ASSAY_PH,
                   assembly: dict[str, Any] | None = None) -> dict[str, Any]:
    """Molecule-level flags.

    ``assembly`` is the output of ``constant_regions.assemble``. When supplied,
    pI and net charge are computed on the full-length H2L2 rather than on the
    variable domains, which is the basis every published charge-and-clearance
    analysis used. Without it the values are variable-domain only and are
    labelled as rank-order proxies throughout.
    """
    if not chains:
        raise ValueError("At least one chain is required.")
    if not 0 <= formulation_ph <= 14 or not 0 <= assay_ph <= 14:
        raise ValueError("pH values must be between 0 and 14.")

    per_chain = {}
    for label in sorted(chains):
        sequence = chains[label].get("sequence")
        if not sequence:
            continue
        per_chain[label] = {"pi": _pi(sequence),
                            "net_charge_assay_ph": _charge(sequence, assay_ph),
                            "net_charge_formulation_ph": _charge(sequence, formulation_ph)}
    if not per_chain:
        raise ValueError("No chain sequences were available.")

    # Which sequence the molecule-level descriptors are computed on.
    full_length = bool(assembly and assembly.get("assembled_sequence"))
    if full_length:
        subject = assembly["assembled_sequence"]
        components = [(assembly[key]["sequence"], assembly[key]["copies"])
                      for key in ("heavy", "light") if assembly.get(key)]
        basis_note = (f"Full-length {assembly['format']}: variable domains with "
                      f"{assembly['heavy']['constant_label']} heavy constant"
                      + (f" and {assembly['light']['constant_label']} light constant"
                         if assembly.get("light") else "")
                      + (", C-terminal lysine clipped" if assembly.get("c_terminal_lysine_clipped")
                         else ", C-terminal lysine retained")
                      + ". Sequence-only charge sums the physical chains at their stated stoichiometry, "
                        "with unmodified N and C termini on each chain. pI is the zero of that summed charge.")
    else:
        subject = "".join(chains[label]["sequence"] for label in sorted(per_chain))
        components = [(chains[label]["sequence"], 1) for label in sorted(per_chain)]
        basis_note = ("Variable-domain charge sums one copy of each supplied chain, with unmodified "
                      "N and C termini on each domain. pI is the zero of that summed charge. "
                      "Constant regions and structural effects on ionization are excluded; "
                      "these are sequence-only proxies, not measured molecule properties.")
    fv = {"pi": _molecular_pi(components),
          "net_charge_assay_ph": round(_molecular_charge(components, assay_ph), 2),
          "net_charge_formulation_ph": round(_molecular_charge(components, formulation_ph), 2),
          "full_length": full_length,
          "residues": len(subject),
          "basis": basis_note}

    flags: list[dict[str, Any]] = []

    # 1. Colloidal margin: net charge approaches zero near pI, so repulsion between
    #    molecules is weakest there. The action is a formulation change, not a mutation.
    subject_label = "Molecule" if full_length else "Fv"
    margin = round(abs(fv["pi"] - formulation_ph), 2)
    verdict = "red" if margin < 0.5 else "amber" if margin < 1.5 else "green"
    flags.append(_flag(
        id="colloidal_ph_margin", label="Formulation pH margin from pI",
        value=margin, units="pH units", verdict=verdict,
        bands={"measure": f"|pI \u2212 formulation pH|, here |{fv['pi']} \u2212 {formulation_ph}|",
               "cuts": [("red", "< 0.5"), ("amber", "0.5 to 1.5"), ("green", "\u2265 1.5")],
               "direction": "higher is better"},
        detail=f"{subject_label} pI {fv['pi']} against a formulation pH of {formulation_ph}.",
        conclusion={"red": "The formulation pH sits essentially at the pI, where net charge and "
                           "electrostatic repulsion are at a minimum. Highest colloidal risk.",
                    "amber": "The formulation pH is close to the pI. Colloidal repulsion is reduced.",
                    "green": "The formulation pH is well separated from the pI."}[verdict],
        outcomes=["aggregation"],
        action="This is a formulation lever before it is an engineering one: moving the buffer pH "
               "away from the pI costs nothing in sequence. See the suggested window below.",
        source_ids=["colloidal_ph_margin", "biopython"],
        policy="Bands (<0.5 red, <1.5 amber) are this project's screening policy. Real colloidal "
               "behaviour also depends on ionic strength, buffer species, excipients, protein "
               "concentration and conformational stability, none of which is modelled."))

    # 2. Net charge. Non-monotonic for clearance: Sharma et al. 2014 reported faster
    #    clearance at both high positive and high negative net charge, so a low pI is
    #    not automatically favourable.
    net = fv["net_charge_assay_ph"]
    magnitude = abs(net)
    amber_at, red_at = _scale_charge_bands(full_length)
    verdict = "red" if magnitude >= red_at else "amber" if magnitude >= amber_at else "green"
    direction = "positive" if net > 0 else "negative" if net < 0 else "neutral"
    outcomes = ["clearance", "viscosity"] + (["nonspecific_binding"] if net > 0 else [])
    flags.append(_flag(
        id="net_charge", label=f"{subject_label} net charge", value=net,
        units=f"e at pH {assay_ph}",
        verdict=verdict,
        bands={"measure": f"|net charge| at pH {assay_ph}, here |{net}|",
               "cuts": [("green", f"< {amber_at:g}"), ("amber", f"{amber_at:g} to {red_at:g}"),
                        ("red", f"\u2265 {red_at:g}")],
               "direction": "nearer zero is milder, but the clearance relationship is "
                            "non-monotonic in either direction"},
        detail=f"Net {direction} charge of {net} at pH {assay_ph}; pI {fv['pi']}.",
        conclusion={"red": f"Strongly {direction} for {'a full-length IgG' if full_length else 'an Fv'}. "
                           "Charge extremes in either direction "
                           "have been associated with faster clearance.",
                    "amber": f"Moderately {direction}. Worth tracking alongside a poly-specificity "
                             "measurement.",
                    "green": "Net charge sits in the middle of the range."}[verdict],
        outcomes=outcomes,
        action="Measure poly-specificity (PSR, BVP ELISA or heparin retention) before treating "
               "charge as the dominant clearance risk; the published analysis found the assay the "
               "stronger discriminator.",
        source_ids=["charge_clearance", "clinical_pk_charge", "biopython"],
        policy=f"Bands (|charge| >={amber_at:g} amber, >={red_at:g} red on this basis) are this "
               "project's policy, not published cut-offs. The charge-clearance relationship is non-monotonic, so a low pI is not "
               "automatically favourable. Clearance is additionally governed by FcRn binding, "
               "target-mediated disposition and anti-drug antibodies, none of which is visible "
               "from a variable domain."))

    # 3. CDR hydrophobicity, which is the CDR-specific term the clearance work used.
    cdr = _cdr_hydropathy(chains)
    if cdr["mean"] is None:
        flags.append(_flag(id="cdr_hydrophobicity", label="CDR hydrophobicity", value=None,
                           units="mean Kyte-Doolittle", verdict="unknown",
                           detail="IMGT numbering was unavailable, so CDR residues could not be selected.",
                           conclusion="Not assessed.", outcomes=["clearance", "aggregation"],
                           action="Supply a numbered chain to assess this.",
                           source_ids=["charge_clearance"], policy="Not assessed."))
    else:
        mean = cdr["mean"]
        verdict = "red" if mean >= 0.5 else "amber" if mean >= 0.0 else "green"
        flags.append(_flag(
            id="cdr_hydrophobicity", label="CDR hydrophobicity",
            value=mean, units="mean Kyte-Doolittle over CDR residues", verdict=verdict,
            bands={"measure": "mean Kyte-Doolittle hydropathy over IMGT CDR residues",
                   "cuts": [("green", "< 0.0"), ("amber", "0.0 to 0.5"), ("red", "\u2265 0.5")],
                   "direction": "lower is better"},
            detail=f"Mean hydropathy {mean} across {cdr['residue_count']} IMGT CDR residues.",
            conclusion={"red": "CDRs are hydrophobic on average. This is the term most directly "
                               "associated with fast clearance and with self-association.",
                        "amber": "CDR hydrophobicity is mid-range.",
                        "green": "CDRs are net hydrophilic."}[verdict],
            outcomes=["clearance", "aggregation", "nonspecific_binding", "viscosity"],
            action="Confirm with HIC retention and a self-interaction measurement rather than "
                   "acting on the descriptor.",
            source_ids=["charge_clearance", "biopython"],
            policy="Bands are this project's policy. A mean over CDR residues is a cruder measure "
                   "than the per-CDR structural hydrophobicity used in the source work; per-CDR "
                   "values are reported alongside.",
            per_cdr=cdr["per_cdr"]))

    # 4. Charge asymmetry between chains, a sequence-level stand-in for the Fv charge
    #    dipole that the viscosity work implicates. The structural term is TAP's SFvCSP.
    if len(per_chain) >= 2:
        labels = sorted(per_chain)
        if full_length and assembly.get("light"):
            chain_charges = {"H": _charge(assembly["heavy"]["sequence"], assay_ph),
                             "L": _charge(assembly["light"]["sequence"], assay_ph)}
        else:
            chain_charges = {l: per_chain[l]["net_charge_assay_ph"] for l in labels}
        values = list(chain_charges.values())
        asymmetry = round(max(values) - min(values), 2)
        verdict = "amber" if asymmetry >= 6 else "green"
        # A large difference does not mean the two chains are oppositely
        # charged. +10 and +3 differ by 7 with both strongly positive, and
        # calling that "opposite charge" misdescribes the molecule and points
        # at the wrong formulation lever. The signs are checked, not inferred
        # from the magnitude of the gap.
        opposite = min(values) < 0 < max(values)
        if verdict == "green":
            conclusion = "Chain charges are comparable; no large sequence-level dipole."
        elif opposite:
            conclusion = (f"The two chains carry opposite charge ({max(values):+g} and "
                          f"{min(values):+g} e). A charge dipole across the Fv has been "
                          "associated with high-concentration viscosity.")
        else:
            direction = "positive" if min(values) >= 0 else "negative"
            conclusion = (f"The chains differ by {asymmetry:g} e but are both {direction} "
                          f"({max(values):+g} and {min(values):+g} e), so this is an imbalance "
                          "rather than a dipole. The viscosity association is reported for a "
                          "charge dipole across the Fv; a same-sign imbalance is weaker "
                          "evidence for it.")
        flags.append(_flag(
            id="charge_asymmetry",
            label="Inter-chain charge asymmetry" + (" (full-length)" if full_length else ""),
            value=asymmetry, units=f"e difference at pH {assay_ph}", verdict=verdict,
            bands={"measure": f"difference between the two chains' net charge at pH {assay_ph}, "
                              "larger minus smaller; the signs are reported separately",
                   "cuts": [("green", "< 6"), ("amber", "\u2265 6")],
                   "direction": "lower is better"},
            detail="; ".join(f"{l} {charge:+g}" for l, charge in chain_charges.items()),
            opposite_sign=opposite,
            per_chain_charge=chain_charges,
            conclusion=conclusion,
            outcomes=["viscosity"],
            action="If high concentration matters, measure viscosity directly; compute TAP's "
                   "SFvCSP from a structure for the structural version of this term.",
            source_ids=["charge_clearance", "biopython"],
            policy="A per-chain charge difference is a weak proxy for a spatial charge dipole. "
                   "The >=6 band is this project's policy."))

    # 5. Hydrophobic patch proxy, carried through rather than recomputed.
    patches = sum((chain.get("metrics") or {}).get("aggregation_proxy", {}).get("patch_count", 0)
                  for chain in chains.values())
    if any((chain.get("metrics") or {}).get("aggregation_proxy") for chain in chains.values()):
        verdict = "amber" if patches >= 3 else "green"
        flags.append(_flag(
            id="hydrophobic_patches", label="Sequence hydrophobic windows",
            value=patches, units="7-residue windows", verdict=verdict,
            bands={"measure": "count of 7-residue windows with mean hydropathy \u2265 1.6 "
                              "and at most one charged residue",
                   "cuts": [("green", "< 3"), ("amber", "\u2265 3")],
                   "direction": "lower is better"},
            detail=f"{patches} window(s) with mean hydropathy >=1.6 and <=1 charged residue.",
            conclusion=("Several contiguous hydrophobic stretches are present."
                        if verdict == "amber" else
                        "Few or no contiguous hydrophobic stretches."),
            outcomes=["aggregation", "nonspecific_binding"],
            action="A structure-based patch analysis supersedes this; run the structure step.",
            source_ids=["biopython"],
            policy="A sequence-window count is a screening proxy and not an aggregation "
                   "prediction. The >=3 band is this project's policy."))

    # Trade-offs. These descriptors pull against each other, and collapsing them into
    # one score would hide that, so a detected tension is reported as its own finding.
    tradeoffs: list[dict[str, str]] = []
    if magnitude >= amber_at:
        tradeoffs.append({
            "between": "net charge: clearance vs viscosity",
            "statement": "Raising net charge further would be expected to help high-concentration "
                         "viscosity and to work against clearance. These cannot both be optimised "
                         "by moving charge in one direction.",
            "source_ids": ["charge_clearance"]})
    if cdr["mean"] is not None and cdr["mean"] >= 0.0 and margin < 1.5:
        tradeoffs.append({
            "between": "hydrophobicity and colloidal margin",
            "statement": "Mid-to-high CDR hydrophobicity together with a formulation pH close to "
                         "the pI compounds aggregation risk: hydrophobic contact is not opposed by "
                         "charge repulsion. Address the pH margin first, since it needs no sequence change.",
            "source_ids": ["colloidal_ph_margin", "charge_clearance"]})

    # Suggested formulation window: where |net charge| is appreciable and the pH is at
    # least 1.5 units off the pI. Arithmetic on the titration curve, not a recommendation.
    curve = [{"ph": ph / 4, "net_charge": round(_molecular_charge(components, ph / 4), 2)}
             for ph in range(16, 35)]
    window = [point for point in curve
              if abs(point["ph"] - fv["pi"]) >= 1.5 and 5.0 <= point["ph"] <= 7.0]
    suggestion = None
    if window:
        best = max(window, key=lambda p: abs(p["net_charge"]))
        suggestion = {
            "ph_range": [min(p["ph"] for p in window), max(p["ph"] for p in window)],
            "highest_charge_ph": best["ph"], "net_charge_there": best["net_charge"],
            "basis": f"Points between pH 5.0 and 7.0 that are at least 1.5 units from the {subject_label} pI, "
                     "ranked by net charge magnitude. This is arithmetic on a sequence titration "
                     "curve, not a formulation recommendation, and it ignores conformational "
                     "stability, excipients and ionic strength."}

    worst = min((f["verdict"] for f in flags), key=lambda v: -VERDICT_ORDER[v])
    by_outcome: dict[str, list[str]] = {}
    for flag in flags:
        for outcome in flag["outcomes"]:
            by_outcome.setdefault(outcome, []).append(flag["id"])

    for flag in flags:
        flag["fv_only"] = not full_length

    return {
        "formulation_ph": formulation_ph, "assay_ph": assay_ph,
        "per_chain": per_chain, "fv": fv, "flags": flags,
        "basis": "full_length" if full_length else "variable_domain",
        "assembly": ({k: v for k, v in assembly.items()
                      if k not in {"assembled_sequence"}} if full_length else None),
        "worst_verdict": worst,
        "verdict_counts": {name: sum(f["verdict"] == name for f in flags)
                           for name in ("green", "amber", "red", "unknown")},
        "outcomes": by_outcome,
        "tradeoffs": tradeoffs,
        "titration": curve,
        "formulation_suggestion": suggestion,
        "excluded": {
            "immunogenicity": "Charge and hydrophobicity reach immunogenicity only through "
                              "aggregation, so no immunogenicity verdict is derived from them. "
                              "The direct sequence-level measure is non-human peptide content, "
                              "reported separately as framework identity."},
        "limitation": (
            ("Values are computed on the assembled full-length molecule, which is the basis the "
             "cited analyses used. Constant regions are one allotype, carry no glycan, and "
             "include no Fc engineering. No clearance rate, viscosity or aggregation rate is "
             "predicted, and FcRn affinity, target-mediated disposition and anti-drug antibodies "
             "are not modelled.")
            if full_length else
            ("Every cited dataset used full-length IgG; these values come from variable domains "
             "only and are rank-order proxies. Attach a constant region for values comparable "
             "with those analyses. No clearance rate, viscosity or aggregation rate is predicted.")),
    }
