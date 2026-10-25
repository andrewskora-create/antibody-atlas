"""Percentiles over the therapeutic panel, and the guard that withholds them.

The panel used to be 36 hand-picked molecules, where "observed min to max" was
the only honest summary: a percentile over a curated subset describes the
curation. It is now every approved human or humanised paired whole mAb in the
snapshot, which is what makes a position within it interpretable.

Two things have to hold. The percentile must be arithmetically right, and it
must be absent when the panel is too small -- a silently-computed percentile
over a handful of references is worse than none, because it reads as precision.
"""
from __future__ import annotations

import unittest

from app.profiler_interpretation import (MIN_REFERENCES_FOR_PERCENTILE, _benchmark,
                                         _percentile_of, _quantile)


def refs(values):
    return [{"value": v} for v in values]


class PercentileArithmetic(unittest.TestCase):
    def test_empirical_definition_is_fraction_at_or_below(self):
        ordered = [1.0, 2.0, 3.0, 4.0]
        self.assertEqual(_percentile_of(0.5, ordered), 0.0)
        self.assertEqual(_percentile_of(1.0, ordered), 25.0)   # at or below, inclusive
        self.assertEqual(_percentile_of(3.5, ordered), 75.0)
        self.assertEqual(_percentile_of(9.0, ordered), 100.0)

    def test_quantile_interpolates(self):
        ordered = [0.0, 10.0]
        self.assertAlmostEqual(_quantile(ordered, 0.5), 5.0)
        self.assertAlmostEqual(_quantile(ordered, 0.0), 0.0)
        self.assertAlmostEqual(_quantile(ordered, 1.0), 10.0)

    def test_single_reference_quantile_does_not_divide_by_zero(self):
        self.assertEqual(_quantile([7.0], 0.25), 7.0)


class PercentileGuard(unittest.TestCase):
    def test_small_panel_gets_no_percentile_and_says_so(self):
        small = refs(range(MIN_REFERENCES_FOR_PERCENTILE - 1))
        benchmark = _benchmark("k", "L", 10.0, "u", small, "up", "scale")
        self.assertIsNone(benchmark["percentile"])
        self.assertIn("Withheld", benchmark["percentile_basis"])
        # The range must survive: withholding the percentile must not remove
        # the comparison that was there before.
        self.assertIsNotNone(benchmark["min"])
        self.assertIsNotNone(benchmark["max"])
        self.assertIn(benchmark["query_relation"],
                      {"within_range", "above_range", "below_range"})

    def test_large_panel_gets_a_percentile(self):
        large = refs(range(MIN_REFERENCES_FOR_PERCENTILE + 1))
        benchmark = _benchmark("k", "L", 10.0, "u", large, "up", "scale")
        self.assertIsNotNone(benchmark["percentile"])
        self.assertIn("quartiles", benchmark)
        self.assertLessEqual(benchmark["quartiles"]["p5"], benchmark["quartiles"]["p25"])
        self.assertLessEqual(benchmark["quartiles"]["p25"], benchmark["quartiles"]["p75"])
        self.assertLessEqual(benchmark["quartiles"]["p75"], benchmark["quartiles"]["p95"])

    def test_non_finite_query_gets_no_percentile(self):
        large = refs(range(MIN_REFERENCES_FOR_PERCENTILE + 1))
        benchmark = _benchmark("k", "L", None, "u", large, "up", "scale")
        self.assertIsNone(benchmark["percentile"])


class WhatThePercentileAdds(unittest.TestCase):
    def test_within_range_can_still_be_an_extreme_position(self):
        """The reason this work was done, expressed as a test.

        A molecule inside the observed range but above almost every reference
        reads as unremarkable under min/max and as extreme under a percentile.
        """
        values = [float(i) for i in range(100)] + [1000.0]   # one far outlier
        benchmark = _benchmark("k", "L", 99.0, "u", refs(values), "up", "scale")
        self.assertEqual(benchmark["query_relation"], "within_range")
        self.assertGreaterEqual(benchmark["percentile"], 99.0)

    def test_percentile_is_consistent_with_relation_at_the_edges(self):
        values = [float(i) for i in range(MIN_REFERENCES_FOR_PERCENTILE + 5)]
        top = _benchmark("k", "L", max(values), "u", refs(values), "up", "scale")
        self.assertEqual(top["percentile"], 100.0)
        self.assertEqual(top["query_relation"], "within_range")
        bottom = _benchmark("k", "L", min(values) - 1, "u", refs(values), "up", "scale")
        self.assertEqual(bottom["percentile"], 0.0)
        self.assertEqual(bottom["query_relation"], "below_range")


class SequencePanelAgreement(unittest.TestCase):
    def test_both_percentile_thresholds_match(self):
        """Two modules compute percentiles; a reader must not get one in the
        sequence table and none in the descriptor table for the same panel."""
        from app.sequence_panel import MIN_REFERENCES_FOR_PERCENTILE as other
        self.assertEqual(MIN_REFERENCES_FOR_PERCENTILE, other)


if __name__ == "__main__":
    unittest.main()


class DistributedStubContents(unittest.TestCase):
    """What ships in the panel stub, and what must never be dropped from it.

    The stub is what a fresh clone rehydrates from, so a field missing here is
    a feature silently switched off for every user. Two have gone wrong:
    `sequence_descriptors` was produced by a one-off script that was never
    committed, and a rebuild without it disables every sequence-descriptor
    comparison with no error anywhere.
    """

    @classmethod
    def setUpClass(cls):
        import json
        from pathlib import Path
        cls.path = Path(__file__).resolve().parents[1] / "data" / "profiler_comparators.stripped.json"
        if not cls.path.is_file():
            raise unittest.SkipTest("panel stub not present")
        cls.panel = json.loads(cls.path.read_text())

    def test_no_reference_carries_sequence_text(self):
        for reference in self.panel["references"]:
            self.assertNotIn("sequences", reference,
                             f"{reference.get('id')} still carries Thera-SAbDab sequence text")
        self.assertTrue(self.panel.get("sequences_removed"))

    def test_every_reference_keeps_its_sequence_descriptors(self):
        missing = [r.get("id") for r in self.panel["references"]
                   if not r.get("sequence_descriptors")]
        self.assertEqual(missing, [],
                         "sequence_panel.panel_ranges() silently returns {} without these, "
                         "which disables every sequence descriptor comparison")

    def test_sequence_descriptors_have_the_fields_the_bands_need(self):
        from app.sequence_panel import METRICS
        for reference in self.panel["references"][:5]:
            for chain, values in reference["sequence_descriptors"].items():
                for metric in METRICS:
                    self.assertIn(metric, values,
                                  f"{reference.get('id')} {chain} is missing {metric}")

    def test_antiberty_keeps_what_the_comparison_reads(self):
        """The trim must not remove the likelihood or the provenance signature."""
        for reference in self.panel["references"][:5]:
            antiberty = reference.get("antiberty") or {}
            if antiberty.get("status") != "complete":
                continue
            for chain, record in antiberty["chains"].items():
                self.assertIn("pseudo_log_likelihood", record)
                self.assertIn("provenance", record)

    def test_antiberty_per_residue_detail_is_not_shipped(self):
        """It is ~90% of the file and nothing reads it."""
        for reference in self.panel["references"]:
            for record in ((reference.get("antiberty") or {}).get("chains") or {}).values():
                self.assertNotIn("residue_scores", record)
