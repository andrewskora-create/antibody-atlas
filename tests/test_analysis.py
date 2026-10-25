import unittest

from app.analysis import clean_sequence, ptm_sites, sequence_metrics


class SequenceAnalysisTests(unittest.TestCase):
    def test_fasta_and_ambiguous_residue_rejection(self):
        self.assertEqual(clean_sequence(">H\naCdE\nFG\n"), "ACDEFG")
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            clean_sequence("ACDX")
        with self.assertRaisesRegex(ValueError, "multiple FASTA"):
            clean_sequence(">H\nACD\n>L\nEFG")

    def test_sequon_requires_non_proline_middle_residue(self):
        sequons = [site for site in ptm_sites("NPTNATS") if site["kind"] == "N-linked glycosylation sequon"]
        self.assertEqual([(site["start"], site["motif"]) for site in sequons], [(4, "NAT")])

    def test_charge_changes_with_ph(self):
        low = sequence_metrics("ACDEFGHIKLMNPQRSTVWY", 5)
        high = sequence_metrics("ACDEFGHIKLMNPQRSTVWY", 9)
        self.assertGreater(low["net_charge_at_ph"], high["net_charge_at_ph"])


if __name__ == "__main__":
    unittest.main()
