import unittest
from unittest.mock import patch

from app.profiler_interpretation import annotate_antiberty, annotate_propermab


def anti(value=-1, model_hash="weights"):
    return {"pseudo_log_likelihood": value, "provenance": {"version": "1",
            "checkpoint_artifacts": [{"filename": "model.safetensors", "sha256": model_hash},
                                     {"filename": "config.json", "sha256": "config"}]}}


def protocol():
    return {"software_versions": {name: "1" for name in ("ImmuneBuilder", "anarcii", "freesasa", "openmm", "numpy", "scipy", "numba")},
            "structure_model": "Test model", "preparation": "Test preparation",
            "structure_checkpoint_artifacts": [{"filename": "model", "sha256": "test_hash"}]}


def pm(value=10, ph=7.4, adapter="1"):
    return {"ph": ph, "source_commit": "source", "adapter_version": adapter,
            "model_checkpoint": None, "comparison_protocol": protocol(),
            "features": [{"key": "hyd_asa", "label": "Apolar accessible area", "value": value, "unit": "Å²"}]}


def panel():
    return {"description": "Two synthetic test controls", "protocol": protocol(), "references": [
        {"id": "a", "name": "Control A", "chain_types": {"H": "H", "L": "K"},
         "antiberty": {"status": "complete", "chains": {"H": anti(-2), "L": anti(-0.2)}},
         "propermab": {"status": "complete", "results": pm(10)}},
        {"id": "b", "name": "Control B", "chain_types": {"H": "H", "L": "L"},
         "antiberty": {"status": "complete", "chains": {"H": anti(-0.5), "L": anti(-4)}},
         "propermab": {"status": "complete", "results": pm(20)}},
    ]}


class ProfilerInterpretationTests(unittest.TestCase):
    @patch("app.profiler_interpretation._load_panel", side_effect=panel)
    def test_anti_comparison_matches_actual_subclass_and_preserves_input(self, _):
        raw = {"L": anti(-1)}
        result = annotate_antiberty(raw, {"L": {"numbering": {"chain_type": "K"}}})
        benchmark = result["L"]["summary"]["benchmarks"][0]
        self.assertEqual(benchmark["reference_n"], 1)
        self.assertEqual(benchmark["reference_values"][0]["name"], "Control A")
        self.assertEqual(benchmark["query_relation"], "below_range")
        self.assertNotIn("summary", raw["L"])

    @patch("app.profiler_interpretation._load_panel", side_effect=panel)
    def test_unknown_subclass_or_changed_checkpoint_withholds_comparison(self, _):
        for raw, chains in (({"L": anti()}, {}),
                            ({"H": anti(model_hash="changed")}, {"H": {"numbering": {"chain_type": "H"}}})):
            result = next(iter(annotate_antiberty(raw, chains).values()))
            benchmark = result["summary"]["benchmarks"][0]
            self.assertEqual(benchmark["reference_n"], 0)
            self.assertIsNone(benchmark["min"])
            self.assertEqual(benchmark["query_relation"], "unavailable")

    @patch("app.profiler_interpretation._load_panel", side_effect=panel)
    def test_pm_actual_range_and_protocol_rejection(self, _):
        result = annotate_propermab(pm(25))
        benchmark = result["features"][0]["benchmark"]
        self.assertEqual((benchmark["min"], benchmark["median"], benchmark["max"]), (10, 15, 20))
        self.assertEqual(benchmark["query_relation"], "above_range")
        self.assertIn("not validated failure", result["conclusion"])
        for query in (pm(ph=6), pm(adapter="changed")):
            withheld = annotate_propermab(query)["features"][0]["benchmark"]
            self.assertEqual(withheld["query_relation"], "unavailable")
            self.assertEqual(withheld["reference_n"], 0)

    @patch("app.profiler_interpretation._load_panel", return_value={})
    def test_missing_panel_retains_intrinsic_scale_without_inventing_range(self, _):
        result = annotate_propermab(pm())
        benchmark = result["features"][0]["benchmark"]
        self.assertIsNone(benchmark["max"])
        self.assertIn("Å²", benchmark["scale_description"])
        self.assertEqual(benchmark["reference_values"], [])

    @patch("app.profiler_interpretation._load_panel", side_effect=panel)
    def test_equality_is_within_range_not_outside(self, _):
        self.assertEqual(annotate_propermab(pm(10))["features"][0]["benchmark"]["query_relation"], "within_range")

    @patch("app.profiler_interpretation._load_panel", side_effect=panel)
    def test_missing_metadata_or_changed_structure_checkpoint_withholds(self, _):
        missing = pm()
        missing.pop("comparison_protocol")
        changed = pm()
        changed["comparison_protocol"]["structure_checkpoint_artifacts"][0]["sha256"] = "different"
        changed_version = pm()
        changed_version["comparison_protocol"]["software_versions"]["openmm"] = "different"
        for query in (missing, changed, changed_version):
            with self.subTest(query=query):
                result = annotate_propermab(query)
                self.assertEqual(result["features"][0]["benchmark"]["query_relation"], "unavailable")


    def test_individual_reference_protocol_is_checked(self):
        refs=panel()
        refs['references'][1]['propermab']['results']['comparison_protocol']['preparation']='different'
        with patch('app.profiler_interpretation._load_panel',return_value=refs):
            b=annotate_propermab(pm(15))['features'][0]['benchmark']
        self.assertEqual(b['reference_n'],1)
        self.assertEqual(b['max'],10)

    def test_native_override_requires_matching_fingerprint_for_patches(self):
        from app.profiler_interpretation import _protocol_compatible
        q=protocol();q['software_versions'].update({'pdb2pqr':'1','apbs-binary':'1'})
        override={**q,'nanoshaper_override_sha256':'other-executable'}
        self.assertFalse(_protocol_compatible(override,q,'hyd_patch_area'))
        self.assertTrue(_protocol_compatible(override,override,'hyd_patch_area'))
        self.assertTrue(_protocol_compatible(override,q,'hyd_asa'))


if __name__ == "__main__":
    unittest.main()
