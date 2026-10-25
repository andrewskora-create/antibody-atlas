"""Hydrogen-bond geometry and side-chain repacking.

The repacking result is only meaningful if the wild type goes through the
identical pipeline: minimising just the mutant and comparing it with the input
measures the minimisation, not the mutation. ``identity_control`` is the proof
of that, and it is the most important test here -- if the protocol alone moves
hydrogen bonds, no substitution result from this module can be believed.
"""
from __future__ import annotations

import unittest
from pathlib import Path

from app import hbonds

ROOT = Path(__file__).resolve().parents[1]
FV = ROOT / "tests" / "fixtures" / "example_fv.pdb"

try:
    import openmm  # noqa: F401
    import pdbfixer  # noqa: F401
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

needs_openmm = unittest.skipUnless(
    HAVE_OPENMM, "needs OpenMM and PDBFixer, which the structure path installs")


@needs_openmm
class Geometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.topology, cls.positions = hbonds.protonate(FV)
        cls.bonds = hbonds.find_hydrogen_bonds(cls.topology, cls.positions)

    def test_an_fv_has_a_beta_sheet_worth_of_backbone_bonds(self):
        """An Fv is two beta sandwiches. Far fewer than this means the angle
        criteria or the hydrogen placement is wrong."""
        backbone = [b for b in self.bonds if b.kind == "backbone-backbone"]
        self.assertGreater(len(backbone), 80)

    def test_every_bond_satisfies_the_stated_criteria(self):
        for bond in self.bonds:
            self.assertLessEqual(bond.distance_a, hbonds.MAX_DONOR_ACCEPTOR_A)
            self.assertGreaterEqual(bond.angle_deg, hbonds.MIN_ANGLE_DEG)
            self.assertLessEqual(bond.angle_deg, 180.0)

    def test_a_residue_is_not_bonded_to_itself(self):
        for bond in self.bonds:
            same = (bond.donor_chain == bond.acceptor_chain
                    and bond.donor_resid == bond.acceptor_resid)
            self.assertFalse(same)

    def test_insertion_codes_are_part_of_residue_identity(self):
        """IMGT 112 and 112A are different residues, and IMGT puts its
        insertion codes in CDR-H3 -- the loop most likely to be engineered.
        Keying on the integer alone merges them and mis-assigns bonds."""
        labels = {(b.donor_chain, b.donor_resid) for b in self.bonds}
        labels |= {(b.acceptor_chain, b.acceptor_resid) for b in self.bonds}
        lettered = {x for x in labels if not str(x[1]).isdigit()}
        self.assertTrue(lettered, "fixture should contain insertion-coded residues")
        for chain, resid in lettered:
            bare = "".join(c for c in str(resid) if c.isdigit())
            if (chain, bare) in labels:
                a = [b for b in self.bonds if b.involves(chain, resid)]
                b_ = [b for b in self.bonds if b.involves(chain, bare)]
                self.assertNotEqual(
                    {id(x) for x in a}, {id(x) for x in b_},
                    f"{chain}{resid} and {chain}{bare} resolved to the same bonds")

    def test_sidechain_classification_is_consistent(self):
        for bond in self.bonds:
            if bond.kind == "backbone-backbone":
                self.assertFalse(bond.donor_is_sidechain or bond.acceptor_is_sidechain)
            if bond.kind == "sidechain-sidechain":
                self.assertTrue(bond.donor_is_sidechain and bond.acceptor_is_sidechain)

    def test_burial_bands_match_the_structure_analysis_thresholds(self):
        self.assertEqual(hbonds.burial_band(0.05), "buried")
        self.assertEqual(hbonds.burial_band(0.20), "intermediate")
        self.assertEqual(hbonds.burial_band(0.50), "exposed")
        self.assertEqual(hbonds.burial_band(None), "unknown")

    def test_summary_carries_its_source(self):
        summary = hbonds.summarize_position(self.bonds, "H", "35", 0.05, 50, "predicted")
        self.assertEqual(summary["source"], "predicted")
        self.assertIn("predicted model", summary["note"])
        experimental = hbonds.summarize_position(self.bonds, "H", "35", 0.05, 50,
                                                 "experimental")
        self.assertIn("experimental", experimental["note"])


@needs_openmm
class Repacking(unittest.TestCase):
    def test_identity_control_is_quiet(self):
        """Mutating a residue to itself must change nothing.

        This is the noise floor for every other result in the module. If the
        protocol alone moves hydrogen bonds, a substitution that moves fewer
        than that has told you nothing.
        """
        from app.repack import identity_control
        result = identity_control(FV, "H", "35")
        self.assertLessEqual(result["bonds_changed_by_protocol_alone"], 1)

    def test_truncation_loses_side_chain_bonds(self):
        from app import repack as rp
        from app.hbonds import find_hydrogen_bonds, protonate
        topology, positions = protonate(FV)
        bonds = find_hydrogen_bonds(topology, positions)
        target = None
        for residue in topology.residues():
            label = hbonds.residue_label(residue)
            count = sum(1 for b in bonds
                        if b.sidechain_of(residue.chain.id, label))
            if count >= 2 and residue.name not in {"GLY", "ALA", "PRO"}:
                target = (residue.chain.id, label)
                break
        self.assertIsNotNone(target, "no suitable residue in the fixture")
        result = rp.repack(FV, target[0], target[1], "G")
        self.assertGreater(result["hydrogen_bonds"]["lost_at_site"], 0)
        self.assertEqual(result["hydrogen_bonds"]["gained_at_site"], 0,
                         "glycine has no side chain and cannot gain a bond")

    def test_no_free_energy_is_claimed(self):
        """A plausible kcal/mol between two different sequences would be the
        most misleading number this project could print."""
        from app import repack as rp
        result = rp.repack(FV, "H", "35", "A")
        self.assertNotIn("ddg", str(result).lower())
        self.assertIn("Local OpenMM relaxation does not calculate folding free energy", result["no_free_energy"])

    def test_bad_target_is_refused(self):
        from app import repack as rp
        with self.assertRaises(rp.RepackError):
            rp.repack(FV, "H", "35", "Z")

    def test_missing_residue_is_refused(self):
        from app import repack as rp
        with self.assertRaises(rp.RepackError):
            rp.repack(FV, "H", "9999", "A")


if __name__ == "__main__":
    unittest.main()
