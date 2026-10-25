"""Attach human constant regions to a variable domain, for molecule-level descriptors.

Every published analysis relating charge to clearance, viscosity or aggregation
measured full-length IgG. Computing pI from variable domains alone therefore
answers a different question than the one those papers answered: the constant
regions are most of the molecule by mass and carry most of its ionizable
groups, so an Fv-only pI rank-orders candidates but cannot be compared with a
published threshold.

Attaching a constant region closes that gap, because the constant regions are
invariant: within a subclass and allotype every therapeutic IgG carries the
same one. So the descriptor becomes a property of the molecule rather than of
the fragment, at the cost of the assumptions recorded in ASSUMPTIONS below.

What this does NOT do is turn a descriptor into a prediction. A full-length pI
is still a descriptor. Clearance is also governed by FcRn affinity,
target-mediated disposition and anti-drug antibodies, none of which follows
from a sequence.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

REFERENCE_PATH = Path(__file__).resolve().parents[1] / "data" / "ig_constant_regions.json"

# An IgG is two heavy and two light chains.
HEAVY_COPIES = 2
LIGHT_COPIES = 2

ASSUMPTIONS = [
    "One allotype per subclass. Allotypes differ at charged positions and shift the computed pI; "
    "the shipped sequence is the UniProt canonical, not a consensus.",
    "No glycan. Fc N-glycans, and sialylated glycoforms in particular, carry negative charge that "
    "lowers the real pI. Charge-variant profiles seen by cIEF reflect glycoforms this does not model.",
    "No engineering. IgG4 S228P, Fc silencing or half-life mutations, and charge variants are not "
    "applied; attach the result to the construct you actually have.",
    "Assembly is sequence concatenation. It is not a structural model, so no interface, no "
    "disulfide bonding and no conformation is represented.",
]


class ConstantRegionError(ValueError):
    """The requested constant region is not available."""


@lru_cache(maxsize=1)
def _reference() -> dict:
    if not REFERENCE_PATH.is_file():
        raise ConstantRegionError(
            "Constant-region reference data is missing. Run "
            "'python setup/fetch_constant_regions.py' to fetch it from UniProt.")
    return json.loads(REFERENCE_PATH.read_text())


def available() -> dict[str, Any]:
    """What can be attached, for the interface to offer."""
    data = _reference()
    heavy, light = [], []
    for key, region in sorted(data["regions"].items()):
        entry = {"key": key, "label": region["label"], "length": region["length"],
                 "accession": region["accession"], "note": region["note"]}
        (heavy if region["kind"] == "heavy" else light).append(entry)
    return {"heavy": heavy, "light": light,
            "source": data["source"], "license": data["license"],
            "release": data.get("uniprot_release"),
            "copyright": data["copyright"],
            "assumptions": ASSUMPTIONS,
            "limitations": data.get("limitations", [])}


def default_light_for(chain_type: str | None) -> str:
    """Kappa and lambda take different constant regions; follow what ANARCII called."""
    return "IGKC" if str(chain_type or "K").upper().startswith("K") else "IGLC2"


def _region(key: str, kind: str) -> dict:
    regions = _reference()["regions"]
    region = regions.get(key)
    if region is None:
        choices = ", ".join(sorted(k for k, v in regions.items() if v["kind"] == kind))
        raise ConstantRegionError(f"Unknown {kind} constant region '{key}'. Available: {choices}.")
    if region["kind"] != kind:
        raise ConstantRegionError(f"'{key}' is a {region['kind']} constant region, not {kind}.")
    return region


def assemble(chains: dict[str, dict], heavy_isotype: str = "IGHG1",
             light_constant: str | None = None,
             clip_c_terminal_lysine: bool = True) -> dict[str, Any]:
    """Build full-length chains by appending constant regions to the variable domains.

    ``clip_c_terminal_lysine`` removes the heavy-chain C-terminal lysine, which
    is clipped by carboxypeptidase in most expressed product. Clipping is
    partial and heterogeneous in reality; the flag makes the choice explicit
    rather than silent, and removing it costs the molecule two positive charges.
    """
    if "H" not in chains:
        raise ConstantRegionError("A heavy variable domain is required to assemble a full-length IgG.")
    heavy_region = _region(heavy_isotype, "heavy")

    heavy_v = chains["H"]["sequence"]
    heavy_c = heavy_region["sequence"]
    clipped = False
    if clip_c_terminal_lysine and heavy_c.endswith("K"):
        heavy_c, clipped = heavy_c[:-1], True
    heavy_full = heavy_v + heavy_c

    light_full = None
    light_region = None
    if "L" in chains:
        chain_type = (chains["L"].get("numbering_imgt") or chains["L"].get("numbering") or {}
                      ).get("chain_type")
        light_region = _region(light_constant or default_light_for(chain_type), "light")
        light_full = chains["L"]["sequence"] + light_region["sequence"]

    # Molecule-level descriptors need the real stoichiometry, not one copy of each.
    assembled = heavy_full * HEAVY_COPIES + (light_full * LIGHT_COPIES if light_full else "")

    data = _reference()
    return {
        "assembled_sequence": assembled,
        "format": "H2L2" if light_full else "H2 (no light chain supplied)",
        "complete": light_full is not None,
        "heavy": {"variable_length": len(heavy_v), "constant": heavy_region["key"],
                  "constant_label": heavy_region["label"],
                  "constant_accession": heavy_region["accession"],
                  "constant_length": len(heavy_c), "full_length": len(heavy_full),
                  "sequence": heavy_full, "copies": HEAVY_COPIES,
                  "note": heavy_region["note"]},
        "light": ({"variable_length": len(chains["L"]["sequence"]),
                   "constant": light_region["key"], "constant_label": light_region["label"],
                   "constant_accession": light_region["accession"],
                   "constant_length": light_region["length"],
                   "full_length": len(light_full), "sequence": light_full,
                   "copies": LIGHT_COPIES, "note": light_region["note"]}
                  if light_full else None),
        "c_terminal_lysine_clipped": clipped,
        "total_residues": len(assembled),
        "source": data["source"], "source_url": data["source_url"],
        "license": data["license"], "copyright": data["copyright"],
        "release": data.get("uniprot_release"), "retrieved": data.get("retrieved_utc"),
        "assumptions": ASSUMPTIONS,
        "caveat": ("Full-length values are comparable with published full-length analyses in a way "
                   "that variable-domain values are not. They remain descriptors: no clearance "
                   "rate, viscosity or aggregation rate is predicted, and FcRn affinity, "
                   "target-mediated disposition and anti-drug antibodies are not modelled."),
    }
