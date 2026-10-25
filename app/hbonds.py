"""Hydrogen-bond geometry for a variable-domain structure.

Why this exists. The substitution ranking scores conservativeness, germline
and clinical precedent, and whether a position is a predicted contact. None of
that is geometry, so a buried asparagine holding a framework hydrogen bond and
a solvent-facing one on the same loop score identically. That is the most
common way a proposal that looks clean turns out to be destabilising.

What it computes. Donors and acceptors from the coordinates, paired by the
standard distance-and-angle criteria, then annotated with how buried the bond
is and how tightly packed its surroundings are. Burial is what makes the
output actionable rather than merely descriptive: losing a solvent-exposed
hydrogen bond is usually tolerable and losing a buried one often is not.

What it does not compute. An energy. These are counts and geometry, so two
substitutions that each break one buried bond cannot be ranked against each
other here. ``app/repack.py`` exists for that question.

Source matters more than usual. Hydrogen bonds live on side chains, and side
chains are the least reliable part of a predicted backbone. Every bond this
module returns carries the ``source`` it was computed from, and nothing
downstream should report one from a predicted model with the same confidence
as one from a crystal or a map.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

#: Standard geometric criteria. Heavy-atom donor-acceptor distance, the
#: hydrogen-acceptor distance, and the donor-hydrogen-acceptor angle. These are
#: the conventional cut-offs rather than this project's invention; they are
#: permissive enough to include weak bonds, which suits a screen that is meant
#: to flag what a substitution would disturb.
MAX_DONOR_ACCEPTOR_A = 3.5
MAX_HYDROGEN_ACCEPTOR_A = 2.5
MIN_ANGLE_DEG = 120.0

#: Relative side-chain SASA below which a bond counts as buried. Matches the
#: banding structure_analysis uses, so the two sections agree about what
#: "buried" means.
BURIED_BELOW = 0.10
EXPOSED_ABOVE = 0.30

#: Radius for the packing-density count, in angstrom. Heavy atoms only.
PACKING_RADIUS_A = 8.0

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

#: Donors require explicit bonded hydrogens. Acceptor eligibility additionally
#: depends on residue/atom chemistry and protonation, never on element alone.
DONOR_ELEMENTS = {"N", "O", "S"}
CHEMICAL_ROLE_NOTE = (
    "Donors require explicit bonded hydrogens; amide, Lys and Arg nitrogens are not acceptors. "
    "Histidine acceptors use explicit ring hydrogens or HID/HIE/HIP names; ambiguous "
    "unprotonated HIS is excluded. Protonation is modeled, not experimentally measured. "
    "Unknown residue chemistries and disulfide/neutral thiol sulfur acceptors are excluded."
)


def _is_acceptor(row, attached, rows_by_residue) -> bool:
    name, atom = row["resname"], row["name"]
    has_h = bool(attached.get(row["index"]))
    if name not in {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                    "ASH", "GLH", "HID", "HIE", "HIP", "CYM", "CYX", "LYN"}:
        return False
    if atom in {"O", "OXT"}:
        return not has_h
    if atom in {"OD1", "OD2"} and name in {"ASP", "ASH"}:
        return not has_h
    if atom in {"OE1", "OE2"} and name in {"GLU", "GLH"}:
        return not has_h
    if (name, atom) in {("ASN", "OD1"), ("GLN", "OE1"), ("SER", "OG"),
                        ("THR", "OG1"), ("TYR", "OH"), ("MET", "SD")}:
        return True
    if name == "CYM" and atom == "SG":
        return not has_h
    if name == "LYN" and atom == "NZ":
        return len(attached.get(row["index"], [])) <= 2
    if name in {"HIS", "HID", "HIE", "HIP"} and atom in {"ND1", "NE2"}:
        if has_h or name == "HIP":
            return False
        if name == "HID":
            return atom == "NE2"
        if name == "HIE":
            return atom == "ND1"
        ring = rows_by_residue[(row["chain"], row["resid"])]
        return any(r["name"] in {"ND1", "NE2"} and attached.get(r["index"]) for r in ring)
    return False


@dataclass
class HydrogenBond:
    # Residues are identified by their full label including any insertion
    # code -- "112" and "112A" are different residues. Using the integer alone
    # merged them, and IMGT puts its insertion codes in CDR-H3, which is both
    # the loop most likely to be engineered and the one where a mis-assigned
    # hydrogen bond would do the most damage. Even the test fixture has two
    # such pairs.
    donor_chain: str
    donor_resid: str
    donor_resname: str
    donor_atom: str
    acceptor_chain: str
    acceptor_resid: str
    acceptor_resname: str
    acceptor_atom: str
    distance_a: float
    angle_deg: float
    kind: str              # backbone-backbone | backbone-sidechain | sidechain-sidechain
    donor_is_sidechain: bool
    acceptor_is_sidechain: bool

    def involves(self, chain: str, resid: str) -> bool:
        resid = str(resid)
        return ((self.donor_chain == chain and self.donor_resid == resid)
                or (self.acceptor_chain == chain and self.acceptor_resid == resid))

    def sidechain_of(self, chain: str, resid: str) -> bool:
        """True when this bond uses the SIDE CHAIN of the given residue.

        The distinction decides whether a substitution can break the bond at
        all: a backbone amide survives every substitution except to proline,
        while a side-chain bond is lost whenever the replacement lacks the
        donor or acceptor.
        """
        resid = str(resid)
        if self.donor_chain == chain and self.donor_resid == resid and self.donor_is_sidechain:
            return True
        return (self.acceptor_chain == chain and self.acceptor_resid == resid
                and self.acceptor_is_sidechain)


class ProtonationError(RuntimeError):
    """Raised when hydrogens could not be placed, so geometry is not available."""


def protonate(pdb_path: str | Path, ph: float = 7.4):
    """Add missing atoms and hydrogens. Returns (topology, positions in A).

    Predicted Fv models and most deposited structures carry no hydrogens, and
    hydrogen-bond geometry needs them: the donor-hydrogen-acceptor angle is
    what separates a hydrogen bond from two atoms that happen to be close.
    """
    try:
        from openmm import unit
        from pdbfixer import PDBFixer
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise ProtonationError(f"PDBFixer and OpenMM are required: {exc}") from exc
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    # Deliberately NOT building missing loops. Inventing backbone that was not
    # resolved and then reporting hydrogen bonds from it would be fabrication;
    # an unresolved region should come back as absent, not as modelled.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    positions = np.array(fixer.positions.value_in_unit(unit.nanometer)) * 10.0
    return fixer.topology, positions


def _atom_rows(topology, positions):
    rows = []
    for atom in topology.atoms():
        residue = atom.residue
        rows.append({
            "index": atom.index,
            "name": atom.name,
            "element": atom.element.symbol if atom.element is not None else "X",
            "chain": residue.chain.id,
            "resid": residue_label(residue),
            "seq_number": _numeric_part(residue.id),
            "resname": residue.name,
            "is_sidechain": atom.name not in BACKBONE_ATOMS,
            "xyz": positions[atom.index],
        })
    return rows


def find_hydrogen_bonds(topology, positions) -> list[HydrogenBond]:
    """Every donor-hydrogen...acceptor triple meeting the criteria above."""
    rows = _atom_rows(topology, positions)
    by_index = {r["index"]: r for r in rows}
    # hydrogens attached to each heavy atom, from the topology's own bonds
    attached: dict[int, list[int]] = {}
    for a, b in topology.bonds():
        for heavy, light in ((a, b), (b, a)):
            if (light.element is not None and light.element.symbol == "H"
                    and heavy.element is not None and heavy.element.symbol != "H"):
                attached.setdefault(heavy.index, []).append(light.index)

    rows_by_residue = {}
    for row in rows:
        rows_by_residue.setdefault((row["chain"], row["resid"]), []).append(row)
    acceptors = [r for r in rows if _is_acceptor(r, attached, rows_by_residue)]
    donors = [r for r in rows if r["element"] in DONOR_ELEMENTS and attached.get(r["index"])]
    if not acceptors or not donors:
        return []

    acceptor_xyz = np.array([r["xyz"] for r in acceptors])
    bonds: list[HydrogenBond] = []
    for donor in donors:
        delta = acceptor_xyz - donor["xyz"]
        distances = np.linalg.norm(delta, axis=1)
        for idx in np.where((distances <= MAX_DONOR_ACCEPTOR_A) & (distances > 0.1))[0]:
            acceptor = acceptors[idx]
            if acceptor["index"] == donor["index"]:
                continue
            # Same residue, or directly bonded: not a hydrogen bond.
            if (acceptor["chain"] == donor["chain"]
                    and acceptor["resid"] == donor["resid"]):
                continue
            # i to i+1 backbone N...O is a geometric artefact of the chain.
            if (acceptor["chain"] == donor["chain"]
                    and abs(acceptor["seq_number"] - donor["seq_number"]) <= 1
                    and not acceptor["is_sidechain"] and not donor["is_sidechain"]):
                continue
            best = None
            for h_index in attached[donor["index"]]:
                h_xyz = by_index[h_index]["xyz"]
                h_to_a = float(np.linalg.norm(acceptor["xyz"] - h_xyz))
                if h_to_a > MAX_HYDROGEN_ACCEPTOR_A:
                    continue
                v1 = donor["xyz"] - h_xyz
                v2 = acceptor["xyz"] - h_xyz
                norm = np.linalg.norm(v1) * np.linalg.norm(v2)
                if norm == 0:
                    continue
                angle = math.degrees(math.acos(
                    max(-1.0, min(1.0, float(np.dot(v1, v2) / norm)))))
                if angle < MIN_ANGLE_DEG:
                    continue
                if best is None or angle > best[1]:
                    best = (h_to_a, angle)
            if best is None:
                continue
            kind = ("sidechain-sidechain" if donor["is_sidechain"] and acceptor["is_sidechain"]
                    else "backbone-backbone" if not donor["is_sidechain"] and not acceptor["is_sidechain"]
                    else "backbone-sidechain")
            bonds.append(HydrogenBond(
                donor_chain=donor["chain"], donor_resid=donor["resid"],
                donor_resname=donor["resname"], donor_atom=donor["name"],
                acceptor_chain=acceptor["chain"], acceptor_resid=acceptor["resid"],
                acceptor_resname=acceptor["resname"], acceptor_atom=acceptor["name"],
                distance_a=round(float(distances[idx]), 2),
                angle_deg=round(best[1], 1), kind=kind,
                donor_is_sidechain=donor["is_sidechain"],
                acceptor_is_sidechain=acceptor["is_sidechain"]))
    return bonds


def residue_label(residue) -> str:
    """'112' or '112A'. The insertion code is part of the identity."""
    return f"{residue.id}{(residue.insertionCode or '').strip()}".strip()


def _numeric_part(value) -> int:
    match = re.match(r"-?\d+", str(value))
    return int(match.group(0)) if match else 0


def packing_density(topology, positions, chain: str, resid: str,
                    radius: float = PACKING_RADIUS_A) -> int | None:
    """Heavy atoms within `radius` of this residue's side-chain centroid.

    A crude but honest measure of how much room a replacement side chain
    would have. It is a count, not an energy, and glycine has no side chain
    so it returns None rather than a misleading zero.
    """
    rows = _atom_rows(topology, positions)
    resid = str(resid)
    side = [r["xyz"] for r in rows
            if r["chain"] == chain and r["resid"] == resid
            and r["is_sidechain"] and r["element"] != "H"]
    if not side:
        return None
    centroid = np.mean(np.array(side), axis=0)
    heavy = np.array([r["xyz"] for r in rows if r["element"] != "H"])
    return int(np.sum(np.linalg.norm(heavy - centroid, axis=1) <= radius))


def burial_band(relative_sasa: float | None) -> str:
    if relative_sasa is None:
        return "unknown"
    if relative_sasa < BURIED_BELOW:
        return "buried"
    if relative_sasa < EXPOSED_ABOVE:
        return "intermediate"
    return "exposed"


def summarize_position(bonds: Iterable[HydrogenBond], chain: str, resid: str,
                       relative_sasa: float | None = None,
                       packing: int | None = None,
                       source: str = "predicted") -> dict[str, Any]:
    """What a substitution at this position would have to preserve.

    `source` travels with the answer because the same geometry means different
    things from a map and from a model, and the reader has to be able to tell.
    """
    resid = str(resid)
    here = [b for b in bonds if b.involves(chain, resid)]
    sidechain = [b for b in here if b.sidechain_of(chain, resid)]
    return {
        "chain": chain,
        "position": resid,
        "source": source,
        "chemical_role_assumptions": CHEMICAL_ROLE_NOTE,
        "hydrogen_bonds_total": len(here),
        "sidechain_hydrogen_bonds": len(sidechain),
        "backbone_hydrogen_bonds": len(here) - len(sidechain),
        "burial": burial_band(relative_sasa),
        "relative_sidechain_sasa": relative_sasa,
        "packing_heavy_atoms_8a": packing,
        "bonds": [asdict(b) for b in sidechain],
        # The actionable line. A backbone bond survives any substitution but
        # proline; a side-chain bond is lost whenever the replacement lacks
        # the donor or acceptor, and losing a buried one is the case worth
        # stopping for.
        "at_risk": len(sidechain),
        "note": _position_note(len(sidechain), burial_band(relative_sasa), source),
    }


def _position_note(sidechain_bonds: int, burial: str, source: str) -> str:
    origin = ("from an experimental structure" if source == "experimental"
              else "from a predicted model, where side-chain placement is the "
                   "least reliable part")
    if sidechain_bonds == 0:
        return (f"No side-chain hydrogen bond at this position {origin}. A "
                f"substitution here does not break one, though it may still "
                f"introduce a clash or bury a charge.")
    if burial == "buried":
        return (f"{sidechain_bonds} buried side-chain hydrogen bond(s) {origin}. "
                f"Losing a buried hydrogen bond is the case most likely to cost "
                f"stability; prefer a replacement that keeps the donor or acceptor.")
    if burial == "exposed":
        return (f"{sidechain_bonds} solvent-exposed side-chain hydrogen bond(s) "
                f"{origin}. Exposed bonds are usually tolerable to lose, because "
                f"water competes for them anyway.")
    return (f"{sidechain_bonds} side-chain hydrogen bond(s) at intermediate "
            f"burial {origin}.")
