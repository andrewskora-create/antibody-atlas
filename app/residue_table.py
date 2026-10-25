"""One row per residue, on one IMGT frame, carrying every annotation we have.

Every other module in this project reports in its own coordinate system. The
PTM scan speaks in one-based input index. ANARCII speaks in IMGT or Kabat
positions, prefixed by chain *type* so that 'K56' means kappa position 56 and
not lysine. The germline library is keyed by bare IMGT position. Structure
output is keyed by chain id and residue number. Paragraph scores only
IMGT CDR+/-2.

Reconciling those by hand is where mistakes happen, and it is the one piece of
work no upstream tool does for you. This module does it once: it joins every
layer onto the canonical IMGT numbering, one row per residue, and attaches a
single engineering verdict per row so the output is a list of positions to act
on rather than a set of flags to cross-reference.

Nothing here computes a new biophysical quantity. Every value is carried
through from the module that produced it, and a layer that did not run leaves
its columns null rather than a default.
"""
from __future__ import annotations

import json
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

REFERENCE_PATH = Path(__file__).resolve().parents[1] / "data" / "human_v_ogrdb_numbered.json"

# The V gene ends at the conserved Cys104; CDR3 is V(D)J junctional, so a
# V-only library cannot supply a reference residue for it. That is a real
# boundary, reported rather than papered over.
GERMLINE_MAX_IMGT = 104

# The call describes the CHEMISTRY at this position. Whether the position is
# also a predicted antigen contact is a separate axis, carried in its own
# column and folded into the next step -- a residue can be both important to
# preserve and a real degradation concern, and collapsing the two into one
# word hid whichever came second.
#
# "evaluate" rather than "fix": these are screening heuristics over a
# predicted structure, and the output of this screen is an experiment to run,
# not a mutation to make.
VERDICTS = {
    "evaluate": "Chemical liability at high priority, or exposure-assessed and above the "
                "monitoring threshold. The first candidate for characterisation and, if "
                "confirmed, substitution.",
    "watch": "Liability present, but not yet exposure-assessed or below high priority. "
             "Characterise before acting.",
    "clear": "No chemical liability detected at this position.",
    "preserve": "No chemical liability, but a predicted antigen contact. Avoid substituting "
                "here unless binding is re-tested.",
}

# What to do, given the chemistry AND the binding involvement together.
NEXT_STEPS = {
    ("evaluate", True): "Characterise degradation and binding before selecting a substitution: "
                        "every edit here also edits the paratope.",
    ("evaluate", False): "Preferred target for substitution once the liability is confirmed "
                         "analytically.",
    ("watch", True): "Monitor. A substitution here would edit the paratope, so confirm the "
                     "liability matters before considering one.",
    ("watch", False): "Characterise before acting; not yet a substitution candidate.",
    ("preserve", True): "No chemistry to act on. Avoid substituting unless binding is re-tested.",
    ("clear", False): "",
}

# What a verdict rests on. A sequence-only verdict is a real verdict -- the PTM
# motif scan either found a motif here or it did not -- but it cannot downgrade
# a buried residue, so sequence-only tables over-report 'watch'. The basis is
# carried per row so a reader can tell the difference without guessing.
BASIS_NOTES = {
    "sequence": "Motif scan only. Buried residues cannot be downgraded without a structure.",
    "exposure": "Relative side-chain exposure from the predicted model informed the priority.",
    "paratope": "Predicted antigen contact was available at this position.",
    "ensemble": "Exposure is the worst case across several structure models, not one model.",
}

GERMLINE_STATUS = {
    "match": "Matches the assigned human V allele at this position.",
    "mismatch_framework": "Differs from the assigned human V allele in framework. Framework "
                          "deviations are where humanisation back-mutations and vernier residues "
                          "sit; worth reviewing.",
    "mismatch_cdr": "Differs from the assigned human V allele inside CDR1/2. Expected for a "
                    "matured or grafted antibody and not notable on its own.",
}


@lru_cache(maxsize=1)
def _germline_index() -> dict:
    data = json.loads(REFERENCE_PATH.read_text())
    return {
        "by_name": {ref["name"]: ref for ref in data["references"]},
        "source": data["source"],
        "source_url": data["source_url"],
        "releases": data["release_versions"],
    }


def _bare(number: Any) -> str:
    """Strip the ANARCII chain-type prefix, leaving the position with insertion."""
    value = str(number or "")
    return value[1:] if value[:1] in {"H", "K", "L"} and value[1:2].isdigit() else value


def _base(number: Any) -> int | None:
    match = re.match(r"^(\d+)", _bare(number))
    return int(match[1]) if match else None


def _liability_rows(chain: dict, structure_sites: list[dict], label: str,
                    precedents: dict[tuple[str, int, str], dict] | None = None,
                    germline_by_position: dict[int, bool | None] | None = None,
                    ) -> dict[int, list[dict]]:
    """Liabilities keyed by one-based input position.

    Structure-analysed sites carry priority, exposure and comparators, so they
    are preferred. Sequence-only sites are carried through where the structure
    step did not run or did not reach them, marked so the reader can tell.
    """
    by_position: dict[int, list[dict]] = {}
    seen: set[tuple[int, str]] = set()
    for site in structure_sites:
        if site.get("chain") != label:
            continue
        position = site.get("sequence_position", site.get("start"))
        if position is None:
            continue
        for flag in site.get("flags") or [site]:
            kind = flag.get("kind") or site.get("kind")
            if (position, kind) in seen:
                continue
            seen.add((position, kind))
            by_position.setdefault(position, []).append({
                "kind": kind,
                "motif": flag.get("motif") or site.get("motif"),
                "priority": flag.get("priority", site.get("priority")),
                "exposure_assessed": flag.get("exposure_assessed", (
                    isinstance(site.get("relative_sidechain_sasa"), (int, float))
                    and math.isfinite(site["relative_sidechain_sasa"])
                    and site["relative_sidechain_sasa"] >= 0)),
                "comparators": site.get("comparators"),
                "source": "structure",
                # Present only on an ensemble run. Carried here because the
                # verdict in this table is where a reader decides whether to
                # act, so that is where "the models disagree about this" has
                # to appear -- not three sections away.
                "exposure_band_stable": site.get("exposure_band_stable"),
                "exposure_models": site.get("exposure_models"),
                "ensemble_note": site.get("ensemble_note"),
            })
    for site in (chain.get("metrics") or {}).get("ptm_sites", []):
        position, kind = site.get("start"), site.get("kind")
        if position is None or (position, kind) in seen:
            continue
        seen.add((position, kind))
        by_position.setdefault(position, []).append({
            "kind": kind, "motif": site.get("motif"), "priority": None,
            "exposure_assessed": False, "comparators": None, "source": "sequence",
        })

    # Both evidence layers are attached in one place so a liability found by
    # the structure route and one found by the motif scan cannot end up
    # annotated differently.
    for position, entries in by_position.items():
        for entry in entries:
            entry["precedent"] = (precedents or {}).get((label, position, entry["kind"]))
            entry["germline_encoded"] = _motif_is_germline(
                position, entry.get("motif"), germline_by_position)
    return by_position


def _motif_is_germline(position: int, motif: str | None,
                       germline_by_position: dict[int, bool | None] | None) -> bool | None:
    """Does every residue of the motif match the assigned human V allele?

    A motif spans more than one residue -- NG is an asparagine and the glycine
    after it -- and it is only germline-encoded if the allele supplies both.
    Unknown anywhere in the span (CDR3, an unassigned allele, a position the
    allele does not have) makes the whole answer unknown rather than False:
    "we could not tell" and "it differs" lead to different decisions.
    """
    if not motif or germline_by_position is None:
        return None
    states = [germline_by_position.get(p) for p in range(position, position + len(motif))]
    if any(state is None for state in states):
        return None
    return all(states)


def precedent_index(precedents: dict | None) -> dict[tuple[str, int, str], dict]:
    """Approved and Phase-III drugs carrying each motif, keyed to join here.

    The search already ran -- /api/fixes attaches it to the substitution
    proposals -- but it was attached to a different object than this table
    reads, so the per-residue view never showed it. It is the most directly
    relevant evidence for a framework liability, so it belongs on the row
    where the engineering call is made, not only beside the fix proposals
    three sections away.
    """
    index: dict[tuple[str, int, str], dict] = {}
    for site in ((precedents or {}).get("sites") or []):
        chain, position, kind = site.get("query_chain"), site.get("position"), site.get("kind")
        if chain is None or position is None:
            continue
        matches = site.get("matches") or []
        approved = {m["drug"] for m in matches if m.get("stage") == "Approved"}
        phase_iii = {m["drug"] for m in matches if m.get("stage") == "Phase-III"}
        # One drug can match on several reference chains; the counts are of
        # drugs, so the names shown have to be deduplicated the same way or
        # the list contradicts the number beside it.
        named, seen = [], set()
        for match in matches:
            if match["drug"] in seen:
                continue
            seen.add(match["drug"])
            named.append(f"{match['drug']} ({match['stage']})")
        index[(chain, position, kind)] = {
            "searched": True,
            "approved": len(approved),
            "phase_iii": len(phase_iii),
            "drugs": named[:6],
            "region": site.get("query_region"),
        }
    return index


def _exposure_index(structure_sites: list[dict], label: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for site in structure_sites:
        if site.get("chain") != label:
            continue
        position = site.get("sequence_position", site.get("start"))
        value = site.get("relative_sidechain_sasa")
        if position is not None and isinstance(value, (int, float)) and math.isfinite(value) and value >= 0:
            out[position] = round(float(value), 3)
    return out


def _paratope_index(paratope: dict | None, label: str) -> dict[int, dict]:
    """Paragraph scores keyed by one-based input position for this chain.

    ``predict`` returns ``by_chain[label]`` keyed by one-based sequence position
    where the IMGT number could be mapped back, and by the string ``IMGT<n>``
    where it could not. Only the mappable ones can join this table, and the
    count of the others is reported so a reader can see what was dropped.
    """
    if not paratope or paratope.get("status") != "available":
        return {}
    out: dict[int, dict] = {}
    for key, record in (paratope.get("by_chain") or {}).get(label, {}).items():
        position = _paratope_position(key)
        if position is None:
            continue
        probability = record.get("probability")
        if not isinstance(probability, (int, float)) or not math.isfinite(probability) or not 0 <= probability <= 1:
            continue
        out[position] = {"score": round(float(probability), 3),
                    "contact": bool(record.get("likely_contact"))}
    return out


def _paratope_position(key: Any) -> int | None:
    """JSON object keys become strings; explicit IMGT keys remain unmapped."""
    if isinstance(key, bool):
        return None
    if isinstance(key, int):
        return key if key > 0 else None
    if isinstance(key, str) and re.fullmatch(r"[1-9][0-9]*", key):
        return int(key)
    return None


def _unmapped_paratope(paratope: dict | None) -> int:
    if not paratope or paratope.get("status") != "available":
        return 0
    return sum(1 for scores in (paratope.get("by_chain") or {}).values()
               for key in scores if _paratope_position(key) is None)


def _verdict(liabilities: list[dict], paratope: dict | None,
             exposure: float | None, region: str = "") -> tuple[str, list[str], str]:
    """The engineering call at this position, what it rests on, and why.

    Two pieces of evidence argue against engineering a framework liability,
    and neither was reaching this function before:

    * The motif is germline-encoded. Every human antibody built on this V
      allele carries it, so it is not a defect introduced by this molecule,
      and back-mutating it is not available as a fix -- germline IS the
      back-mutation. That is evidence of tolerance, not proof: a germline
      framework motif can still deamidate if it is exposed. So it demotes
      "evaluate" to "watch" rather than clearing the position.
    * An approved or Phase-III drug carries the same motif at the same
      position. Same logic, same strength, same demotion.

    Both apply in framework only. Inside a CDR, germline encoding says little
    -- the CDR is what affinity maturation changed -- and the consequence of
    a PTM there is loss of binding rather than a CMC nuisance.
    """
    basis = ["sequence"]
    if exposure is not None:
        basis.append("exposure")
    # An ensemble run says so in the basis, because "exposure informed this"
    # means something different when it is one model's guess than when it is
    # the worst case across three.
    if any((l.get("exposure_models") or 0) > 1 for l in liabilities):
        basis.append("ensemble")
    if paratope is not None:
        basis.append("paratope")
    # A predicted contact used to return here, before any chemistry was read.
    # That made a high-priority deamidation motif on a contact residue show
    # as "preserve" with its liability invisible in the call -- two different
    # questions answered with one word, and the more urgent one lost.
    contact = bool(paratope and paratope.get("contact"))
    if not liabilities:
        return ("preserve" if contact else "clear"), basis, VERDICTS[
            "preserve" if contact else "clear"]

    high = any(l.get("priority") == "high" for l in liabilities)
    assessed = any(l.get("exposure_assessed") for l in liabilities)
    call = "evaluate" if (high or (assessed and any(l.get("priority") == "review"
                                                    for l in liabilities))) else "watch"

    # Softening is decided per LIABILITY, not per position. A position can
    # carry two chemistries -- a germline-encoded deamidation motif and an
    # introduced oxidation site, say -- and letting the first soften the
    # position's whole call would hide the second. The call is softened only
    # if EVERY liability driving it is itself supported.
    in_framework = not str(region).startswith("CDR")
    drivers = [l for l in liabilities
               if l.get("priority") == "high"
               or (l.get("exposure_assessed") and l.get("priority") == "review")]
    if call == "evaluate" and in_framework and drivers:
        per_liability = [_liability_support(l) for l in drivers]
        if all(per_liability):
            kinds = {k for support in per_liability for k in support}
            basis.extend(sorted(kinds))
            return "watch", basis, _softening_reason(drivers, per_liability)
    return call, basis, VERDICTS[call]


def _liability_support(liability: dict) -> set[str]:
    """Which arguments, if any, say this one liability need not be engineered."""
    support: set[str] = set()
    if liability.get("germline_encoded"):
        support.add("germline")
    precedent = liability.get("precedent") or {}
    if (precedent.get("approved", 0) + precedent.get("phase_iii", 0)) > 0:
        support.add("precedent")
    return support


def _softening_reason(drivers: list[dict], support: list[set[str]]) -> str:
    """Name the motif, and say what supports it, one liability at a time."""
    parts = []
    for liability, kinds in zip(drivers, support):
        label = liability.get("motif") or liability.get("kind") or "the motif"
        said = []
        if "germline" in kinds:
            said.append("germline-encoded in the assigned human V allele")
        if "precedent" in kinds:
            precedent = liability.get("precedent") or {}
            n = precedent.get("approved", 0) + precedent.get("phase_iii", 0)
            said.append(f"carried by {n} approved or Phase-III drug"
                        f"{'s' if n != 1 else ''} at this position")
        parts.append(f"{label} is " + " and ".join(said))
    return ("Framework liability demoted from evaluate to watch \u2014 " + "; ".join(parts)
            + ". Evidence of tolerance, not proof: a motif that an approved antibody "
              "tolerates can still matter in a different molecule, and an exposed germline "
              "motif can still degrade. Confirm analytically before treating it as benign.")


def _germline_status(match: bool | None, region: str) -> str | None:
    if match is None:
        return None
    if match:
        return "match"
    return "mismatch_cdr" if str(region).startswith("CDR") else "mismatch_framework"


def build_chain_rows(label: str, chain: dict, structure_sites: list[dict],
                     paratope: dict | None,
                     precedents: dict[tuple[str, int, str], dict] | None = None,
                     structure_residues: list[dict] | None = None,
                     ) -> list[dict[str, Any]]:
    canonical = chain.get("numbering_imgt") or chain.get("numbering")
    if not canonical or canonical.get("scheme") != "imgt":
        raise ValueError(f"Chain {label}: canonical IMGT numbering is required to build the residue table.")
    selected = chain.get("numbering") or canonical
    selected_by_position = {r["sequence_position"]: r for r in selected.get("residues", [])}

    humanness = chain.get("humanness") or {}
    allele_name = humanness.get("nearest_human_v_gene")
    index = _germline_index()
    allele = index["by_name"].get(allele_name) if allele_name else None
    germline_positions = (allele or {}).get("positions", {})

    # Germline is resolved for every position BEFORE the liabilities are
    # assembled: whether a motif is germline-encoded is a property of its
    # whole span, so it cannot be decided one row at a time.
    germline_detail: dict[int, tuple[str | None, bool | None, str | None]] = {}
    for residue in canonical.get("residues", []):
        bare_n = _bare(residue["number"])
        base_n = _base(residue["number"])
        reference = germline_positions.get(bare_n)
        if base_n is not None and base_n > GERMLINE_MAX_IMGT:
            match, note = None, "CDR3/FR4 is V(D)J junctional; a V-gene library has no reference residue."
        elif allele is None:
            match, note = None, "No nearest human V allele was assigned."
        elif reference is None:
            match, note = None, "This position is absent from the assigned allele."
        else:
            match, note = reference == residue["aa"], None
        germline_detail[residue["sequence_position"]] = (reference, match, note)

    liabilities = _liability_rows(
        chain, structure_sites, label, precedents,
        {position: detail[1] for position, detail in germline_detail.items()})
    # All measured/modelled residues have exposure, not only PTM-bearing ones.
    # PTM sites may carry an ensemble maximum; preserve it over the primary
    # structure's value so the table agrees with its chemistry priority.
    exposure = _exposure_index(structure_residues or [], label)
    exposure.update(_exposure_index(structure_sites, label))
    contacts = _paratope_index(paratope, label)

    rows: list[dict[str, Any]] = []
    for residue in canonical.get("residues", []):
        position = residue["sequence_position"]
        bare = _bare(residue["number"])
        alt = selected_by_position.get(position)
        site_liabilities = liabilities.get(position, [])
        site_exposure = exposure.get(position)
        site_contact = contacts.get(position)
        reference_aa, germline_match, germline_note = germline_detail[position]

        verdict, basis, reason = _verdict(
            site_liabilities, site_contact, site_exposure, residue["region"])
        is_contact = bool((site_contact or {}).get("contact"))
        next_step = NEXT_STEPS.get((verdict, is_contact), NEXT_STEPS.get((verdict, False), ""))
        status = _germline_status(germline_match, residue["region"])
        rows.append({
            "chain": label,
            "chain_type": canonical.get("chain_type"),
            "sequence_position": position,
            "aa": residue["aa"],
            "imgt": bare,
            "imgt_label": residue["number"],
            "imgt_region": residue["region"],
            "alt_scheme": selected.get("definition") if selected is not canonical else None,
            "alt_number": _bare(alt["number"]) if alt and selected is not canonical else None,
            "alt_region": alt["region"] if alt and selected is not canonical else None,
            "germline_aa": reference_aa,
            "germline_match": germline_match,
            "germline_status": status,
            "germline_note": germline_note,
            "germline_allele": allele_name,
            "liabilities": site_liabilities,
            "exposure": site_exposure,
            "paratope_score": (site_contact or {}).get("score"),
            "paratope_contact": (site_contact or {}).get("contact"),
            "verdict": verdict,
            "verdict_basis": basis,
            "verdict_reason": reason,
            # Chemistry and binding are separate axes; the step to take
            # depends on both, so it is resolved once here rather than left
            # for a reader to combine two columns in their head.
            "next_step": next_step,
            # Evidence that argues against engineering a framework liability.
            # Carried per row so the table can show why a call was softened.
            "germline_encoded": next((l.get("germline_encoded") for l in site_liabilities
                                      if l.get("germline_encoded") is not None), None),
            "precedent": next((l.get("precedent") for l in site_liabilities
                               if l.get("precedent")), None),
            # A row worth a reader's attention: a liability, a predicted contact,
            # or a framework deviation from the assigned human allele.
            "notable": bool(site_liabilities) or bool((site_contact or {}).get("contact"))
                       or status == "mismatch_framework",
        })
    return rows


def build_residue_table(chains: dict[str, dict], structure: dict | None = None,
                        paratope: dict | None = None,
                        precedents: dict | None = None) -> dict[str, Any]:
    """The consolidated per-residue table across every supplied chain."""
    structure_sites = (structure or {}).get("ptm_sites", []) or []
    precedents_by_site = precedent_index(precedents)
    rows: list[dict[str, Any]] = []
    for label in sorted(chains, key=lambda x: (x != "H", x)):
        rows.extend(build_chain_rows(label, chains[label], structure_sites, paratope,
                                     precedents_by_site, (structure or {}).get("residues", [])))

    index = _germline_index()
    counts: dict[str, int] = {name: 0 for name in VERDICTS}
    for row in rows:
        counts[row["verdict"]] += 1
    germline_counts: dict[str, int] = {name: 0 for name in GERMLINE_STATUS}
    for row in rows:
        if row["germline_status"]:
            germline_counts[row["germline_status"]] += 1
    # A layer is present if it RAN, not if it found something. The motif scan
    # running and returning nothing is a result, not a missing layer, and
    # conflating the two would tell a reader their molecule was unscreened.
    layers = {
        "numbering": True,
        "germline": any((chain.get("humanness") or {}).get("nearest_human_v_gene")
                        for chain in chains.values()),
        "liabilities": any("ptm_sites" in (chain.get("metrics") or {})
                           for chain in chains.values()) or bool(structure_sites),
        "exposure": any(r["exposure"] is not None for r in rows),
        "paratope": bool(paratope) and paratope.get("status") == "available",
    }
    return {
        "rows": rows,
        "verdict_counts": counts,
        "verdict_legend": VERDICTS,
        "basis_legend": BASIS_NOTES,
        "germline_counts": germline_counts,
        "germline_legend": GERMLINE_STATUS,
        "notable_rows": sum(r["notable"] for r in rows),
        # "Notable" is a union of three unrelated reasons, and reporting only
        # its total is what makes the table look full of positions with no
        # chemistry. The reasons are counted separately so the reader can see
        # that a row is here because it is a predicted binding-site contact,
        # not because anything is wrong with it.
        "notable_reasons": {
            "liability": sum(bool(r["liabilities"]) for r in rows),
            "contact": sum(bool(r["paratope_contact"]) for r in rows),
            "framework_mismatch": sum(r["germline_status"] == "mismatch_framework"
                                      for r in rows),
            "liability_only": sum(bool(r["liabilities"]) and not r["paratope_contact"]
                                  for r in rows),
            "contact_only": sum(bool(r["paratope_contact"]) and not r["liabilities"]
                                for r in rows),
            "both": sum(bool(r["liabilities"]) and bool(r["paratope_contact"])
                        for r in rows),
        },
        "evidence_counts": {
            "germline_encoded": sum(bool(r.get("germline_encoded")) for r in rows),
            "with_precedent": sum(bool((r.get("precedent") or {}).get("approved", 0)
                                       or (r.get("precedent") or {}).get("phase_iii", 0))
                                  for r in rows),
            "precedent_searched": sum(bool(r.get("precedent")) for r in rows),
            "softened": sum(r["verdict"] == "watch" and "germline" in (r.get("verdict_basis") or [])
                            or r["verdict"] == "watch" and "precedent" in (r.get("verdict_basis") or [])
                            for r in rows),
        },
        "notable_legend": (
            "A row is notable for any of three independent reasons: it carries a chemical "
            "liability motif, it is a predicted antigen contact, or it differs from the "
            "assigned human germline in framework. Most notable rows carry no liability \u2014 "
            "they are binding-site positions flagged so that a fix proposed elsewhere is not "
            "quietly applied here."),
        "layers": layers,
        "missing_layers": sorted(name for name, present in layers.items() if not present),
        "germline_reference": {"source": index["source"], "url": index["source_url"],
                               "releases": index["releases"]},
        "paratope_positions_unmapped": _unmapped_paratope(paratope),
        "frame": "IMGT, via ANARCII. Every row is keyed to the canonical IMGT position; "
                 "a second numbering scheme is carried alongside when one was selected.",
        "caveat": "Rows carry values through from the modules that produced them. A layer that "
                  "did not run leaves its columns empty rather than defaulting, and each row "
                  "records which layers informed its verdict. Without a structure every liability "
                  "stays at 'watch': the motif scan cannot tell a buried residue from an exposed "
                  "one, so a sequence-only table deliberately over-reports.",
    }
