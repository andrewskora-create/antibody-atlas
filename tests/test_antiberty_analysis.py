import math
import unittest

from app.antiberty_analysis import _summarize_residues, analyze_antiberty


from tests.support import antiberty_checkpoint, requires_antiberty
# Public ImmuneBuilder example; no user sequence is used in tests.
PUBLIC_HEAVY = "EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS"


class ResidueSummaryTests(unittest.TestCase):
    def test_region_aggregation_and_ties_preserve_evidence(self):
        numbering = {"residues": [{"sequence_position": i, "number": f"H{i}",
                                    "region": "CDR1" if i <= 2 else "Framework"}
                                   for i in range(1, 10)]}
        values = [-4.0, -4.0] + [-1.0] * 8
        residues, regions = _summarize_residues("A" * 10, values, numbering)
        self.assertEqual(len(residues), 10)
        self.assertEqual(sum(r["lowest_decile"] for r in residues), 2)
        self.assertEqual(residues[0]["rank_within_chain"], residues[1]["rank_within_chain"])
        self.assertEqual(residues[-1]["region"], "Unassigned")
        self.assertEqual(sum(r["residue_count"] for r in regions), 10)
        self.assertEqual(sum(r["sum_log_likelihood"] for r in regions), sum(values))
        self.assertAlmostEqual(residues[0]["observed_probability"], math.exp(-4))

    def test_missing_or_invalid_residue_scores_fail(self):
        for values in ([-1], [-1, float("nan")], [-1, 0.2]):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    _summarize_residues("AA", values, {})


@requires_antiberty()
class LocalAntiBertyIntegrationTests(unittest.TestCase):
    def test_detailed_outputs_reproduce_official_runner(self):
        import torch
        from antiberty import AntiBERTyRunner
        torch.set_num_threads(1)
        runner = AntiBERTyRunner(device="cpu", checkpoint_path=str(antiberty_checkpoint()))
        detail = analyze_antiberty(runner, {"H": PUBLIC_HEAVY}, scheme="kabat", batch_size=16)["H"]
        official = float(runner.pseudo_log_likelihood([PUBLIC_HEAVY], batch_size=16)[0])
        self.assertAlmostEqual(detail["pseudo_log_likelihood"], official, places=6)
        self.assertEqual(len(detail["residue_scores"]), len(PUBLIC_HEAVY))
        self.assertEqual(sum(r["residue_count"] for r in detail["region_summaries"]), len(PUBLIC_HEAVY))
        self.assertEqual(detail["scheme"], "kabat")
        self.assertTrue(all(math.isfinite(r["masked_log_probability"]) for r in detail["residue_scores"]))
        for probabilities in detail["class_probabilities"].values():
            self.assertAlmostEqual(sum(probabilities.values()), 1, places=6)
        self.assertEqual(len(detail["provenance"]["checkpoint_artifacts"]), 2)
        self.assertTrue(all(len(a["sha256"]) == 64 for a in detail["provenance"]["checkpoint_artifacts"]))


if __name__ == "__main__":
    unittest.main()
