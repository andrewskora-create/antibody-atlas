"""Ingesting an experimental structure.

The machinery downstream of folding always took a PDB path, but only the
folding step ever wrote one, so a crystal or cryo-EM structure could not be
used at all. These tests cover the parts that are easy to get wrong and whose
failure is quiet: which chain is which, where the variable domain ends, and
the fixed-column PDB format.

The co-crystal tests need a file that is not redistributed, so they skip on a
fresh clone rather than fail.
"""
from __future__ import annotations

import unittest
from pathlib import Path

from app import structure_input

ROOT = Path(__file__).resolve().parents[1]
COCRYSTAL = ROOT / "validation" / ".cocrystal_work" / "9pi9.cif"
FIXTURE_FV = ROOT / "tests" / "fixtures" / "example_fv.pdb"

needs_cocrystal = unittest.skipUnless(
    COCRYSTAL.is_file(),
    "needs the co-crystal; fetch it with: python validation/paratope_cocrystal.py")


class Format(unittest.TestCase):
    def test_unknown_file_is_refused_with_a_reason(self):
        tmp = ROOT / "tests" / "_not_a_structure.txt"
        tmp.write_text("hello\n")
        try:
            with self.assertRaises(structure_input.StructureInputError) as caught:
                structure_input.inspect(tmp)
            self.assertIn("PDB or mmCIF", str(caught.exception))
        finally:
            tmp.unlink()

    def test_missing_file_is_refused(self):
        with self.assertRaises(structure_input.StructureInputError):
            structure_input.inspect(ROOT / "does-not-exist.pdb")


class PredictedFvRoundTrip(unittest.TestCase):
    """The ingester must handle an ordinary Fv, not only a Fab complex."""

    def test_predicted_model_is_recognised(self):
        report = structure_input.inspect(FIXTURE_FV)
        roles = {c.chain_id: c.role for c in report.chains}
        self.assertEqual(roles.get("H"), "heavy")
        self.assertEqual(roles.get("L"), "light")
        self.assertEqual(report.fv_pairs, [("H", "L")])
        self.assertFalse(report.antigen_chains)

    def test_prepared_file_parses_back(self):
        """Fixed-column PDB output is unforgiving; a one-column slip writes an
        insertion code into the residue-number field and every parser rejects
        the file."""
        from Bio.PDB import PDBParser
        out = ROOT / "tests" / "_prepared_roundtrip.pdb"
        try:
            structure_input.prepare(FIXTURE_FV, out)
            model = next(PDBParser(QUIET=True).get_structure("x", str(out)).get_models())
            self.assertEqual({c.id for c in model}, {"H", "L"})
        finally:
            out.unlink(missing_ok=True)


class Reconciliation(unittest.TestCase):
    def test_full_match_is_full_coverage(self):
        seq = "EVQLVESGGGLVQPGGSLRLSCAAS"
        result = structure_input.reconcile(seq, seq)
        self.assertEqual(result["coverage"], 1.0)
        self.assertEqual(result["unresolved_runs"], [])
        self.assertTrue(result["usable"])

    def test_a_gap_is_reported_not_skipped(self):
        """Unresolved density must surface as a gap. Silently closing it
        shifts every downstream position by the length of the gap."""
        full = "AAAACCCCDDDDEEEEFFFFGGGG"
        observed = full[:8] + full[12:]          # four residues unresolved
        result = structure_input.reconcile(observed, full)
        self.assertLess(result["coverage"], 1.0)
        self.assertTrue(result["unresolved_runs"])
        self.assertEqual(sum(r["length"] for r in result["unresolved_runs"]), 4)

    def test_low_coverage_is_marked_unusable(self):
        result = structure_input.reconcile("AAAA", "AAAA" + "C" * 40)
        self.assertFalse(result["usable"])
        self.assertIn("guesswork", result["note"])


@needs_cocrystal
class ComplexWithAntigen(unittest.TestCase):
    def test_chains_are_classified_without_being_told(self):
        report = structure_input.inspect(COCRYSTAL)
        roles = {c.chain_id: c.role for c in report.chains}
        self.assertEqual(roles["A"], "heavy")
        self.assertEqual(roles["B"], "light")
        self.assertEqual(sorted(report.antigen_chains), ["E", "F"])

    def test_copies_are_paired_by_proximity(self):
        report = structure_input.inspect(COCRYSTAL)
        self.assertEqual(sorted(report.fv_pairs), [("A", "B"), ("C", "D")])

    def test_multiple_copies_are_reported_not_averaged(self):
        report = structure_input.inspect(COCRYSTAL)
        self.assertTrue(any("copies" in w for w in report.warnings))

    def test_constant_domain_is_dropped_from_the_fv(self):
        out = ROOT / "tests" / "_prepared_complex.pdb"
        try:
            report = structure_input.prepare(COCRYSTAL, out, copy=0)
            heavy = report["chains"]["H"]
            self.assertGreater(heavy["residues_in_file"], heavy["variable_domain_residues"])
            self.assertGreater(heavy["constant_or_untyped_dropped"], 50)
            self.assertEqual(report["source"], "experimental")
        finally:
            out.unlink(missing_ok=True)

    def test_epitope_distance_is_measured_and_cdr_dominated(self):
        out = ROOT / "tests" / "_prepared_epitope.pdb"
        try:
            report = structure_input.prepare(COCRYSTAL, out, copy=0)
            epitope = structure_input.epitope_distances(out, report["antigen_chains"])
            self.assertEqual(epitope["status"], "measured")
            self.assertGreater(epitope["contact_count"], 10)
            # An epitope footprint sits on the CDRs. If most contacts were in
            # framework, the chain mapping or the numbering would be wrong.
            cdr = 0
            for contact in epitope["contacts"]:
                number = int("".join(ch for ch in contact[1:] if ch.isdigit()))
                if 27 <= number <= 38 or 56 <= number <= 65 or 105 <= number <= 117:
                    cdr += 1
            self.assertGreater(cdr / len(epitope["contacts"]), 0.5,
                               "epitope contacts should be CDR-dominated")
        finally:
            out.unlink(missing_ok=True)

    def test_no_antigen_reports_unavailable_rather_than_zero(self):
        out = ROOT / "tests" / "_prepared_noantigen.pdb"
        try:
            structure_input.prepare(COCRYSTAL, out, copy=0, keep_antigen=False)
            epitope = structure_input.epitope_distances(out, [])
            self.assertEqual(epitope["status"], "unavailable")
            self.assertEqual(epitope["residues"], [])
        finally:
            out.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
