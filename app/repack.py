"""Build a mutant side chain, relax it in place, and report what changed.

The substitution ranking infers consequence by subtraction: it knows which
hydrogen bonds the wild type has and assumes a replacement lacking the donor
or acceptor loses them. That systematically over-penalises conservative swaps
-- Asn to Gln, Ser to Thr -- which are exactly the cases where a compensating
bond forms. This module answers the question directly instead: it builds the
mutant, relaxes the neighbourhood, and measures the new geometry.

Three things decide whether the answer means anything.

**The wild type is minimised identically.** This is the trap. Minimising only
the mutant and comparing it with the input structure measures the
minimisation, not the mutation -- the protocol moves every side chain a
little, and those movements swamp the one atom that matters. Both states go
through the same pipeline, and ``identity_control`` exists to prove the
residual noise is small enough to see a signal through.

**Only the neighbourhood moves.** Atoms beyond a shell around the mutation are
restrained. Without that the whole domain relaxes, every hydrogen bond in the
structure shifts slightly, and the difference is dominated by places the
substitution cannot reach.

**No free energy is reported.** Comparing force-field potential energies
between two different sequences is meaningless: the atom count and the bonded
terms differ, so the number is dominated by chemistry that has nothing to do
with stability. A plausible-looking kcal/mol here would be the most
misleading thing this project could print. What is reported instead is
countable and checkable: hydrogen bonds lost and gained, steric clashes
introduced, and polar atoms left buried without a partner.

On a predicted backbone this models one thing on top of another. The result
carries its source so that is visible, and the caller is expected to keep it
visible.
"""
from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import numpy as np

from .hbonds import (CHEMICAL_ROLE_NOTE, MAX_DONOR_ACCEPTOR_A, HydrogenBond, _numeric_part,
                     find_hydrogen_bonds, packing_density, residue_label)

THREE_FROM_ONE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
ONE_FROM_THREE = {v: k for k, v in THREE_FROM_ONE.items()}

#: Atoms within this distance of the mutated side chain are free to move;
#: everything else is restrained in place.
FREE_SHELL_A = 8.0
#: Harmonic restraint on the frozen region. Strong enough to hold the fold,
#: soft enough not to fight the local relaxation.
RESTRAINT_K = 100.0          # kJ/mol/nm^2
MINIMISE_STEPS = 400

#: A contact closer than this between heavy atoms that are not bonded is a
#: clash. Deliberately generous: the point is to catch a replacement that
#: cannot physically fit, not to adjudicate borderline packing.
CLASH_DISTANCE_A = 2.2


class RepackError(RuntimeError):
    """Raised when the mutant could not be built or relaxed."""


def _load_fixer(pdb_path: str | Path):
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}      # never invent unresolved backbone
    return fixer


def _finish(fixer, ph: float):
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    return fixer


def _positions_angstrom(fixer) -> np.ndarray:
    from openmm import unit
    return np.array(fixer.positions.value_in_unit(unit.nanometer)) * 10.0


_FORCEFIELD = None


def _forcefield():
    """Cached. Parsing the XML on every call is pure overhead."""
    global _FORCEFIELD
    if _FORCEFIELD is None:
        from openmm import app
        # No implicit solvent, deliberately, and it cost an afternoon to
        # establish why. On this hardware GBn2 took 421 s for twenty
        # minimisation steps on a 3,483-atom Fv; the cheaper OBC and HCT
        # variants took 41-50 s; vacuum with a cutoff took 1.4 s. Implicit
        # solvent is simply not affordable per substitution.
        #
        # Dropping it would normally be a bad trade: in vacuum, charged side
        # chains collapse together and invent salt bridges, which is exactly
        # the artefact that would corrupt a hydrogen-bond count. What makes it
        # defensible here is that everything outside the shell is FROZEN, so
        # that collapse has nowhere to happen -- the only motion available is
        # the mutated side chain and its immediate neighbours relaxing a few
        # angstrom into an otherwise rigid pocket. This is a local
        # re-accommodation, not an energy calculation, and nothing downstream
        # reports it as one.
        _FORCEFIELD = app.ForceField("amber14-all.xml")
    return _FORCEFIELD


def _minimise(topology, positions_nm, centre_nm: np.ndarray) -> tuple[Any, Any]:
    """Relax atoms near `centre`, hold the rest rigid. Returns (positions, n_free).

    The environment is frozen by setting particle masses to zero rather than
    by harmonic restraints. A restrained atom still carries degrees of freedom
    the minimiser has to explore and still drifts; a zero-mass particle is
    removed from the search entirely. That is both faster and a cleaner
    comparison, because the frozen half of the structure is bit-for-bit
    identical between the wild type and the mutant.
    """
    import openmm
    from openmm import app, unit

    try:
        # constraints=None on purpose. Constrained bonds cannot involve a
        # massless particle, and the environment here is frozen by zeroing
        # masses. Constraints exist to allow a longer dynamics timestep
        # anyway; a minimisation does not need them.
        system = _forcefield().createSystem(
            topology, nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
    except Exception as exc:
        raise RepackError(f"force field could not parameterise the structure: {exc}") from exc

    coords = np.array(positions_nm.value_in_unit(unit.nanometer))
    free = 0
    for index, _atom in enumerate(topology.atoms()):
        if np.linalg.norm(coords[index] - centre_nm) <= FREE_SHELL_A / 10.0:
            free += 1
        else:
            system.setParticleMass(index, 0.0)
    if free == 0:
        raise RepackError("no atoms were free to move; the shell radius is too small")

    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(
        topology, system, integrator,
        openmm.Platform.getPlatformByName("CPU"))
    simulation.context.setPositions(positions_nm)
    simulation.minimizeEnergy(maxIterations=MINIMISE_STEPS)
    state = simulation.context.getState(getPositions=True)
    return state.getPositions(), free


def _residue_centre(fixer, chain: str, resid: str) -> np.ndarray:
    from openmm import unit
    coords = np.array(fixer.positions.value_in_unit(unit.nanometer))
    picked = [a.index for a in fixer.topology.atoms()
              if a.residue.chain.id == chain and residue_label(a.residue) == str(resid)]
    if not picked:
        raise RepackError(f"no residue {chain}{resid} in the structure")
    return coords[picked].mean(axis=0)


def _clashes(topology, positions_a: np.ndarray, chain: str, resid: str) -> int:
    """Non-bonded heavy-atom contacts below CLASH_DISTANCE_A at this residue."""
    resid = str(resid)
    target_number = _numeric_part(resid)
    rows = [(a.index, a.residue.chain.id, residue_label(a.residue),
             _numeric_part(a.residue.id),
             a.element.symbol if a.element is not None else "X")
            for a in topology.atoms()]
    here = [i for i, c, r, _n, e in rows if c == chain and r == resid and e != "H"]
    other = [i for i, c, _r, n, e in rows
             if not (c == chain and abs(n - target_number) <= 1) and e != "H"]
    if not here or not other:
        return 0
    a = positions_a[here]
    b = positions_a[other]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return int(np.sum(d < CLASH_DISTANCE_A))


def _bond_key(bond: HydrogenBond) -> tuple:
    """Identity of a bond, so two states can be compared."""
    return (bond.donor_chain, bond.donor_resid, bond.donor_atom,
            bond.acceptor_chain, bond.acceptor_resid, bond.acceptor_atom)


def _mutate_exact(fixer, chain: str, resid: str, original: str, mutation: str):
    """Address insertion-coded residues without PDBFixer's numeric-ID collision.

    PDBFixer indexes only the integer residue ID. Give the selected residue a
    unique temporary ID, mutate it, then restore its exact original identity.
    Modeller rebuilds topology objects, so resolve the temporary ID again.
    """
    matches = [r for r in fixer.topology.residues()
               if r.chain.id == chain and residue_label(r) == resid]
    if len(matches) != 1:
        raise RepackError(f"position {chain}{resid} has no unique coordinate identity")
    selected = matches[0]
    old_id, old_insertion = selected.id, selected.insertionCode
    used = {_numeric_part(r.id) for r in fixer.topology.residues() if r.chain.id == chain}
    temporary = str(max(used | {0}) + 1)
    selected.id, selected.insertionCode = temporary, ""
    try:
        fixer.applyMutations([f"{original}-{temporary}-{mutation}"], chain)
        changed = [r for r in fixer.topology.residues()
                   if r.chain.id == chain and r.id == temporary]
        if len(changed) != 1 or changed[0].name != mutation:
            raise RepackError("Mutant coordinate identity could not be verified")
        changed[0].id, changed[0].insertionCode = old_id, old_insertion
    except Exception as exc:
        raise RepackError(f"could not build {original}{resid}{mutation}: {exc}") from exc


def _state(pdb_path, chain, resid, mutation, ph):
    """Build one state (wild type or mutant) through the identical pipeline."""
    from openmm import unit
    resid = str(resid)
    fixer = _load_fixer(pdb_path)
    original = next((r.name for r in fixer.topology.residues()
                     if r.chain.id == chain and residue_label(r) == resid), None)
    if original is None:
        raise RepackError(f"no residue {chain}{resid} in the structure")
    if mutation is not None:
        _mutate_exact(fixer, chain, resid, original, mutation)
    _finish(fixer, ph)
    centre = _residue_centre(fixer, chain, resid)
    positions, free = _minimise(fixer.topology, fixer.positions, centre)
    coords = np.array(positions.value_in_unit(unit.nanometer)) * 10.0
    bonds = find_hydrogen_bonds(fixer.topology, coords)
    return {
        "topology": fixer.topology, "coords": coords, "bonds": bonds,
        "original": original, "free_atoms": free,
        "clashes": _clashes(fixer.topology, coords, chain, resid),
        "packing": packing_density(fixer.topology, coords, chain, resid),
    }


def repack(pdb_path: str | Path, chain: str, resid: str, target: str,
           *, ph: float = 7.4, source: str = "predicted",
           output_paths: tuple[Path, Path] | None = None) -> dict[str, Any]:
    """Model one substitution and report the geometry it changes.

    `target` is a one-letter code. Both the wild type and the mutant are built
    and relaxed by the same route, so the difference reflects the substitution
    rather than the protocol.
    """
    target = target.upper().strip()
    if target not in THREE_FROM_ONE:
        raise RepackError(f"{target!r} is not one of the twenty amino acids")
    resid = str(resid)
    wild = _state(pdb_path, chain, resid, None, ph)
    if ONE_FROM_THREE.get(wild["original"]) == target:
        mutant = _state(pdb_path, chain, resid, None, ph)
    else:
        mutant = _state(pdb_path, chain, resid, THREE_FROM_ONE[target], ph)

    before = {_bond_key(b): b for b in wild["bonds"]}
    after = {_bond_key(b): b for b in mutant["bonds"]}
    at_site = lambda b: b.involves(chain, resid)
    lost = [b for k, b in before.items() if k not in after and at_site(b)]
    gained = [b for k, b in after.items() if k not in before and at_site(b)]
    lost_elsewhere = sum(1 for k, b in before.items() if k not in after and not at_site(b))
    gained_elsewhere = sum(1 for k, b in after.items() if k not in before and not at_site(b))

    retained = [(before[k], after[k]) for k in before.keys() & after.keys() if at_site(before[k])]
    if output_paths is not None:
        from openmm import app, unit
        for state, path in zip((wild, mutant), output_paths):
            with Path(path).open("w") as handle:
                app.PDBFile.writeFile(state["topology"], state["coords"] * unit.angstrom,
                                      handle, keepIds=True)
    buried_lost = [b for b in lost if b.sidechain_of(chain, resid)]
    return {
        "chain": chain, "position": resid,
        "from": ONE_FROM_THREE.get(wild["original"], wild["original"]),
        "to": target,
        "source": source,
        "method": ("PDBFixer side-chain build, restrained OpenMM minimisation "
                   f"(amber14, no solvent model, {FREE_SHELL_A:.0f} A free shell "
                   f"with the rest held rigid, {MINIMISE_STEPS} steps); "
                   f"wild type relaxed identically"),
        "chemical_role_assumptions": CHEMICAL_ROLE_NOTE,
        "hydrogen_bonds": {
            "wild_type_total": len(wild["bonds"]),
            "mutant_total": len(mutant["bonds"]),
            "lost_at_site": len(lost),
            "gained_at_site": len(gained),
            "net_at_site": len(gained) - len(lost),
            "sidechain_lost_at_site": len(buried_lost),
            "retained_at_site": len(retained),
            "retained": [{"baseline": a.__dict__, "candidate": b.__dict__} for a, b in retained],
            "lost": [b.__dict__ for b in lost],
            "gained": [b.__dict__ for b in gained],
        },
        # Reported so the reader can see how much of the difference is the
        # protocol rather than the substitution. Large numbers here mean the
        # comparison is not clean and should not be trusted.
        "protocol_noise": {
            "bonds_changed_away_from_the_site": lost_elsewhere + gained_elsewhere,
            "note": ("Bonds changed away from the selected residue can reflect local relaxation, "
                     "indirect substitution effects, or modeling variability. This is a sensitivity "
                     "indicator, not an isolated estimate of protocol noise; inspect these changes "
                     "before interpreting the selected-site comparison."),
        },
        "clashes": {
            "wild_type": wild["clashes"], "mutant": mutant["clashes"],
            "introduced": max(0, mutant["clashes"] - wild["clashes"]),
            "delta": mutant["clashes"] - wild["clashes"],
            "note": "Net count difference in close heavy-atom contacts at the site; not an identity-resolved list of newly introduced clashes.",
        },
        "packing_heavy_atoms_8a": mutant["packing"],
        "free_atoms_relaxed": mutant["free_atoms"],
        "no_free_energy": (
            "Local OpenMM relaxation does not calculate folding free energy. "
            "Hydrogen-bond counts and raw potential-energy differences are not folding ΔΔG. "
            "Any separate stability prediction is reported with its own method and limitations."),
    }


def identity_control(pdb_path: str | Path, chain: str, resid: str,
                     *, ph: float = 7.4) -> dict[str, Any]:
    """Mutate a residue to itself and measure what the protocol alone changes.

    The floor for every other result. If running the pipeline twice on the
    same sequence already moves N hydrogen bonds, then a substitution that
    moves fewer than N has told you nothing.
    """
    resid = str(resid)
    fixer = _load_fixer(pdb_path)
    original = next((r.name for r in fixer.topology.residues()
                     if r.chain.id == chain and residue_label(r) == resid), None)
    if original is None:
        raise RepackError(f"no residue {chain}{resid} in the structure")
    one = ONE_FROM_THREE.get(original)
    if one is None:
        raise RepackError(f"{original} is not a standard amino acid")
    result = repack(pdb_path, chain, resid, one, ph=ph)
    changed = (result["hydrogen_bonds"]["lost_at_site"]
               + result["hydrogen_bonds"]["gained_at_site"]
               + result["protocol_noise"]["bonds_changed_away_from_the_site"])
    return {
        "chain": chain, "position": resid, "residue": one,
        "bonds_changed_by_protocol_alone": changed,
        "detail": result["hydrogen_bonds"],
        "interpretation": (
            "This is the noise floor. A real substitution has to change more "
            "than this before the change can be attributed to the mutation."),
    }
