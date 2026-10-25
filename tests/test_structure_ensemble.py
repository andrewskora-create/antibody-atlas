"""Ensemble exposure: the flag-if-any rule, and the checkpoint guard.

Two behaviours matter here and neither is about numerics.

The first is that a corrupt checkpoint must be detected rather than skipped.
This is not hypothetical: antibody_model_3 in a working installation here was
a truncated 135 MB download where models 2 and 4 are 214 MB. Nothing noticed,
because nothing used anything but model 1. An ensemble that quietly dropped it
would run at three-quarters strength and say it ran at full.

The second is that exposure is reconciled by maximum, not mean. A site exposed
in one model of three is a site worth an assay, and averaging it away is the
failure this feature exists to prevent.
"""
from __future__ import annotations

import unittest

from app import structure_ensemble as ens


def residue(chain, position, exposure, aa="N", region="CDR1"):
    return {"chain": chain, "sequence_position": position, "aa": aa,
            "region": region, "number": f"{chain}{position}",
            "relative_sidechain_sasa": exposure}


def analysis(*residues):
    return {"residues": list(residues), "sites": []}


class Banding(unittest.TestCase):
    def test_bands_match_structure_analysis(self):
        self.assertEqual(ens._band(None), "unknown")
        self.assertEqual(ens._band(0.05), "low")
        self.assertEqual(ens._band(0.10), "intermediate")
        self.assertEqual(ens._band(0.29), "intermediate")
        self.assertEqual(ens._band(0.30), "high")

    def test_band_order_covers_every_band_band_can_return(self):
        produced = {ens._band(v) for v in (None, 0.0, 0.2, 0.9)}
        self.assertTrue(produced <= set(ens.BAND_ORDER))


class Consensus(unittest.TestCase):
    def test_exposure_is_the_maximum_not_the_mean(self):
        result = ens.consensus([
            (1, analysis(residue("H", 31, 0.05))),
            (2, analysis(residue("H", 31, 0.50))),
            (4, analysis(residue("H", 31, 0.05))),
        ])
        row = result["residues"][0]
        self.assertEqual(row["exposure_max"], 0.5)
        self.assertEqual(row["band"], "high")      # mean would be 0.2 -> intermediate
        self.assertFalse(row["band_stable"])
        self.assertEqual(row["models"], 3)

    def test_agreement_is_reported_as_stable(self):
        result = ens.consensus([
            (1, analysis(residue("H", 31, 0.42))),
            (2, analysis(residue("H", 31, 0.45))),
        ])
        row = result["residues"][0]
        self.assertTrue(row["band_stable"])
        self.assertEqual(row["bands_seen"], ["high"])
        self.assertEqual(result["band_disagreements"], 0)

    def test_disagreements_are_counted_and_sorted_first(self):
        result = ens.consensus([
            (1, analysis(residue("H", 31, 0.05), residue("H", 32, 0.40))),
            (2, analysis(residue("H", 31, 0.60), residue("H", 32, 0.41))),
        ])
        self.assertEqual(result["band_disagreements"], 1)
        self.assertEqual(result["residues"][0]["sequence_position"], 31)

    def test_residue_absent_from_one_model_still_reported(self):
        result = ens.consensus([
            (1, analysis(residue("H", 31, 0.4))),
            (2, analysis(residue("H", 31, 0.4), residue("L", 50, 0.2))),
        ])
        light = [r for r in result["residues"] if r["chain"] == "L"][0]
        self.assertEqual(light["models"], 1)

    def test_unmeasured_residue_does_not_become_buried(self):
        """A missing measurement must stay unknown, never read as low exposure."""
        result = ens.consensus([(1, analysis(residue("H", 31, None)))])
        row = result["residues"][0]
        self.assertIsNone(row["exposure_max"])
        self.assertEqual(row["band"], "unknown")

    def test_empty_input_is_unavailable_not_a_crash(self):
        self.assertEqual(ens.consensus([])["status"], "unavailable")


class SiteReconciliation(unittest.TestCase):
    def test_single_model_value_is_preserved_for_audit(self):
        rows = ens.consensus([
            (1, analysis(residue("H", 31, 0.05))),
            (2, analysis(residue("H", 31, 0.55))),
        ])["residues"]
        sites = [{"chain": "H", "sequence_position": 31, "kind": "deamidation",
                  "relative_sidechain_sasa": 0.05, "exposure_band": "low"}]
        updated = ens.reconcile_sites(sites, rows)[0]
        self.assertEqual(updated["exposure_band"], "high")
        self.assertEqual(updated["single_model_exposure_band"], "low")
        self.assertEqual(updated["single_model_exposure"], 0.05)
        self.assertIn("disagree", updated["ensemble_note"])

    def test_stable_band_change_is_explained(self):
        rows = ens.consensus([
            (1, analysis(residue("H", 31, 0.44))),
            (2, analysis(residue("H", 31, 0.46))),
        ])["residues"]
        sites = [{"chain": "H", "sequence_position": 31,
                  "relative_sidechain_sasa": 0.44, "exposure_band": "low"}]
        updated = ens.reconcile_sites(sites, rows)[0]
        self.assertIn("Band changed", updated["ensemble_note"])

    def test_site_with_no_matching_residue_is_untouched(self):
        sites = [{"chain": "L", "sequence_position": 99, "exposure_band": "low"}]
        updated = ens.reconcile_sites(sites, [])[0]
        self.assertEqual(updated, sites[0])


class KeyContract(unittest.TestCase):
    """The field name reconcile_sites is fed from, pinned by a test.

    This exists because of a bug that every other test in this file missed.
    The fold worker read `structure_analysis["sites"]`; analyze_structure
    returns `ptm_sites`. reconcile_sites was handed an empty list every time,
    returned an empty list, and the ensemble's worst-case exposure reached no
    PTM site. Nothing failed -- the unit tests above call reconcile_sites
    directly with synthetic input, so they were testing a function the product
    never actually used on real data.
    """

    def test_analyze_structure_names_its_site_list_ptm_sites(self):
        import inspect
        from app import structure_analysis
        source = inspect.getsource(structure_analysis.analyze_structure)
        self.assertIn('"ptm_sites": sites', source,
                      "analyze_structure's site key changed; the fold worker and "
                      "residue_table both read 'ptm_sites' and must be updated together")

    def test_fold_worker_reconciles_the_key_that_exists(self):
        import inspect
        from app import main
        source = inspect.getsource(main._fold)
        self.assertIn('structure_analysis["ptm_sites"] = structure_ensemble.reconcile_sites',
                      source, "the ensemble must reconcile ptm_sites, not a key that is never set")

    def test_residue_table_reads_the_same_key(self):
        import inspect
        from app import residue_table
        source = inspect.getsource(residue_table.build_residue_table)
        self.assertIn('.get("ptm_sites"', source)


class CheckpointGuard(unittest.TestCase):
    def test_plan_never_returns_fewer_than_one_model(self):
        result = ens.plan(0)
        self.assertGreaterEqual(len(result["model_ids"]), min(1, len(result["available"])))

    def test_plan_reports_what_it_could_not_use(self):
        """Whatever this installation has, the plan must account for all four."""
        result = ens.plan(4)
        accounted = set(result["available"]) | {int(k) for k in result["unusable"]}
        self.assertEqual(accounted, set(ens.ALL_MODEL_IDS))

    def test_plan_flags_when_it_came_up_short(self):
        result = ens.plan(4)
        if len(result["model_ids"]) < 4:
            self.assertTrue(result["short"])
            self.assertTrue(result["note"])


if __name__ == "__main__":
    unittest.main()
