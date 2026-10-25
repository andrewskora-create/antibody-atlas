"""Synthetic-only checks for evidence-bounded variant comparisons."""
import copy
import hashlib
import unittest

from app.variant_comparison import compare_records


def record(sequence="ANSTM", analysis_id="a"):
    digest = hashlib.sha256(f"H:{sequence}".encode()).hexdigest()
    sites = [
        {"kind": "N-linked glycosylation sequon", "start": 2, "motif": sequence[1:4]},
        {"kind": "Oxidation-susceptible residue", "start": 5, "motif": sequence[4]},
    ]
    return {
        "analysis_id": analysis_id, "scheme": "imgt", "ph": 7.4,
        "chains": {"H": {
            "sequence": sequence,
            "numbering": {"residues": [{"sequence_position": i, "number": f"H{i}", "region": "Framework"}
                                        for i in range(1, len(sequence) + 1)]},
            "metrics": {"isoelectric_point": 7.0, "net_charge_at_ph": 1.0,
                        "gravy": -0.1, "ptm_sites": sites},
            "humanness": {"name": "Human V-gene framework identity", "locus": "IGH",
                          "reference_source": "synthetic", "reference_release": "1",
                          "definition": "synthetic identity", "percent": 70.0},
        }},
        "provenance": {"input_sha256": digest, "reference_sha256": "reference"},
    }


def add_structure(item, *, complete=True):
    digest = item["provenance"]["input_sha256"]
    item["provenance"].update(structure_input_sha256=digest,
                              structure_pdb_sha256=f"pdb-{item['analysis_id']}",
                              structure_model="model 1", structure_model_version="1")
    length = len(item["chains"]["H"]["sequence"])
    mapped = length if complete else length - 1
    item["structure_analysis"] = {
        "scheme": "imgt", "ph": 7.4, "method": "surface",
        "patch_definition": {"hydrophobic": "synthetic"},
        "chain_mapping": {"H": {"input_residues": length, "model_residues": mapped,
                                "mapped_residues": mapped, "coverage_percent": 100 * mapped / length}},
        "mapping_coverage_percent": {"H": 100 * mapped / length},
        "patch_counts": {"hydrophobic": 2},
        "patches": [{"kind": "hydrophobic", "exposed_sidechain_sasa_sum_a2": 25.0}],
        "ptm_sites": [{"chain": "H", "sequence_position": 2,
                       "kind": "N-linked glycosylation sequon", "priority": "review",
                       "exposure_band": "low", "exposure_percent": 5.0}],
        "unassessed_ptm_sites": [],
    }


class VariantComparisonTests(unittest.TestCase):
    def test_mutation_ptm_and_metrics_without_structures(self):
        a, b = record(), record("ANATM", "b")
        b["chains"]["H"]["metrics"].update(isoelectric_point=7.3, gravy=-0.2,
                                                   ptm_sites=[{"kind": "Oxidation-susceptible residue",
                                                               "start": 5, "motif": "M"}])
        result = compare_records(a, b)
        self.assertEqual(result["mutation_count"], 1)
        self.assertEqual((result["mutations"][0]["position"], result["mutations"][0]["from"],
                          result["mutations"][0]["to"]), (3, "S", "A"))
        self.assertEqual(result["chains"]["H"]["metrics"]["isoelectric_point"]["delta"], 0.3)
        self.assertEqual(len(result["ptm"]["removed"]), 1)
        self.assertEqual(len(result["ptm"]["retained"]), 1)
        self.assertEqual(result["ptm"]["retained"][0]["priority_transition"]["status"], "unavailable")
        self.assertEqual(result["structure"]["status"], "unavailable")

    def test_indel_and_context_mismatch_rejected(self):
        a, b = record(), record("ANSTM", "b")
        b["chains"]["H"]["sequence"] += "A"
        with self.assertRaisesRegex(ValueError, "insertions and deletions"):
            compare_records(a, b)
        b = record("ANSTM", "b")
        b["ph"] = 6.5
        with self.assertRaisesRegex(ValueError, "same analysis pH"):
            compare_records(a, b)
        b = record("ANSTM", "b")
        b["scheme"] = "kabat"
        with self.assertRaisesRegex(ValueError, "same IMGT or Kabat"):
            compare_records(a, b)
        b = record("ANSTM", "b")
        b["chains"]["L"] = copy.deepcopy(b["chains"]["H"])
        with self.assertRaisesRegex(ValueError, "same chain set"):
            compare_records(a, b)

    def test_owned_structures_enable_site_transitions_and_full_mapping_patches(self):
        a, b = record(), record("ANSTM", "b")
        add_structure(a)
        add_structure(b)
        b["structure_analysis"]["ptm_sites"][0].update(priority="high", exposure_band="high")
        b["structure_analysis"]["patch_counts"]["hydrophobic"] = 3
        b["structure_analysis"]["patches"][0]["exposed_sidechain_sasa_sum_a2"] = 35.0
        result = compare_records(a, b)
        self.assertEqual(result["structure"]["status"], "available")
        self.assertEqual(result["structure"]["patch_counts"]["hydrophobic"]["delta"], 1)
        self.assertEqual(result["structure"]["patch_areas"]["hydrophobic"]["delta"], 10.0)
        site = next(s for s in result["ptm"]["retained"] if s["kind"] == "N-linked glycosylation sequon")
        self.assertEqual(site["priority_transition"]["status"], "changed")
        self.assertEqual(site["exposure_transition"]["candidate"], "high")

    def test_stale_structure_and_partial_mapping_cannot_imply_improvement(self):
        a, b = record(), record("ANSTM", "b")
        add_structure(a)
        add_structure(b, complete=False)
        result = compare_records(a, b)
        self.assertEqual(result["structure"]["patch_counts"], {})
        self.assertIn("withheld", result["structure"]["aggregate_reason"])
        b["provenance"]["structure_input_sha256"] = "wrong"
        result = compare_records(a, b)
        self.assertEqual(result["structure"]["status"], "unavailable")
        self.assertEqual(result["ptm"]["retained"][0]["exposure_transition"]["status"], "unavailable")

    def test_antiberty_requires_matching_checkpoint(self):
        a, b = record(), record("ANSTM", "b")
        for item, value in ((a, -2.0), (b, -1.8)):
            item["chains"]["H"]["antiberty"] = {
                "mean_log_likelihood": value, "scheme": "imgt",
                "provenance": {"package": "antiberty", "version": "1", "checkpoint_artifacts": [{"sha256": "same"}]},
            }
        metric = compare_records(a, b)["chains"]["H"]["metrics"]["antiberty_mean_log_likelihood"]
        self.assertAlmostEqual(metric["delta"], 0.2)
        b["chains"]["H"]["antiberty"]["provenance"]["checkpoint_artifacts"][0]["sha256"] = "other"
        metric = compare_records(a, b)["chains"]["H"]["metrics"]["antiberty_mean_log_likelihood"]
        self.assertEqual(metric["status"], "unavailable")

    def test_propermab_sequence_and_structure_stages(self):
        a, b = record(), record("ANSTM", "b")
        add_structure(a)
        add_structure(b)
        for item, charge, area in ((a, 1.0, 20.0), (b, 2.0, 30.0)):
            provenance = item["provenance"]
            provenance["propermab_input_sha256"] = provenance["input_sha256"]
            provenance["propermab_structure_pdb_sha256"] = provenance["structure_pdb_sha256"]
            item["propermab"] = {
                "source_commit": "source", "adapter_version": "adapter", "comparison_protocol": {"ph": 7.4},
                "numbering_scheme": "IMGT", "ph": 7.4, "scope": "Fv descriptors",
                "input_pdb_sha256": provenance["structure_pdb_sha256"],
                "patch_pipeline": "pipeline", "patch_conditions": "same conditions",
                "features": [{"key": "vh_charge", "label": "VH charge", "value": charge, "unit": "e"},
                             {"key": "hyd_patch_area", "label": "Hydrophobic patch area", "value": area, "unit": "Å²"}],
            }
        pm = compare_records(a, b)["propermab"]
        self.assertEqual(pm["status"], "available")
        self.assertEqual({f["stage"] for f in pm["features"]}, {"sequence", "surface_patch"})
        self.assertEqual(next(f for f in pm["features"] if f["key"] == "hyd_patch_area")["delta"], 10.0)
        b["provenance"]["propermab_structure_pdb_sha256"] = "wrong"
        pm = compare_records(a, b)["propermab"]
        self.assertEqual([f["key"] for f in pm["features"]], ["vh_charge"])
        self.assertIn("matched input PDBs", pm["unavailable"][0]["reason"])
        b["propermab"]["comparison_protocol"] = {"ph": 7.4, "changed": True}
        self.assertEqual(compare_records(a, b)["propermab"]["status"], "unavailable")


if __name__ == "__main__":
    unittest.main()
