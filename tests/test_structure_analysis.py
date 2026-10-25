import json
import unittest
from pathlib import Path

from app.structure_analysis import _group_physical_sites, _sidechain_charge, _site_priority, analyze_structure, extract_trace

HEAVY = "EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS"
LIGHT = "VIWMTQSPSSLSASVGDRVTITCQASQDIRFYLNWYQQKPGKAPKLLISDASNMETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLPFTFGPGTKVDFK"
PDB = Path(__file__).parent / "fixtures" / "example_fv.pdb"


class StructureAnalysisTests(unittest.TestCase):
    def score(self, kind, motif, exposure, region="Framework"):
        return _site_priority({"kind": kind, "motif": motif, "start": 31, "end": 31 + len(motif) - 1},
                              {"chain": "H", "aa": motif[0], "number": "H31", "region": region,
                               "relative_sidechain_sasa": exposure, "sidechain_sasa_a2": 40.0}, "imgt")

    def test_strong_motifs_retain_review_when_buried(self):
        for motif in ("NG", "NS", "DG", "DS"):
            kind = "Asparagine deamidation motif" if motif[0] == "N" else "Aspartate isomerization motif"
            with self.subTest(motif=motif):
                self.assertEqual(self.score(kind, motif, 0, region="CDR1")["priority"], "review")
                self.assertEqual(self.score(kind, motif, 0.299, region="CDR1")["priority"], "review")
                self.assertEqual(self.score(kind, motif, 0.30, region="CDR1")["priority"], "high")

    def test_framework_sites_are_never_high_priority(self):
        """A framework liability is a CMC question, not a high-priority finding.

        The chemistry is identical to a CDR site, so the screen still reports it,
        but it is capped at review and routed to the approved/Phase-III
        comparator search instead of being escalated on exposure alone.
        """
        for motif in ("NG", "NS", "DG", "DS"):
            kind = "Asparagine deamidation motif" if motif[0] == "N" else "Aspartate isomerization motif"
            for region in ("Framework", "Unassigned"):
                with self.subTest(motif=motif, region=region):
                    site = self.score(kind, motif, 0.95, region=region)
                    self.assertEqual(site["priority"], "review")
                    self.assertTrue(site["framework_priority_capped"])
                    self.assertIn("not rated high", " ".join(site["reasons"]))
            # The same motif in a CDR is still allowed to reach high.
            cdr = self.score(kind, motif, 0.95, region="CDR3")
            self.assertEqual(cdr["priority"], "high")
            self.assertFalse(cdr["framework_priority_capped"])

    def test_exposure_alone_never_makes_oxidation_high(self):
        for motif in ("M", "W"):
            for exposure, priority in ((0.099, "lower"), (0.10, "review"), (0.9, "review")):
                with self.subTest(motif=motif, exposure=exposure):
                    site = self.score("Oxidation-susceptible residue", motif, exposure)
                    self.assertEqual(site["priority"], priority)
                    self.assertEqual(site["exposure_percent"], round(exposure * 100, 1))

    def test_glycosylation_is_not_downgraded_by_native_burial(self):
        for exposure in (0, 0.1, 0.9, None):
            site = self.score("N-linked glycosylation sequon", "NAT", exposure)
            self.assertEqual(site["priority"], "review")

    def test_other_asn_asp_motifs_need_exposure_for_review(self):
        for kind, motif in (("Asparagine deamidation motif", "NT"),
                            ("Aspartate isomerization motif", "DT")):
            self.assertEqual(self.score(kind, motif, 0.099)["priority"], "lower")
            self.assertEqual(self.score(kind, motif, 0.10)["priority"], "review")
            self.assertEqual(self.score(kind, motif, 0.90)["priority"], "review")

    def test_missing_exposure_is_unknown_not_zero_or_low(self):
        for exposure in (None, float("nan"), -1):
            site = self.score("Oxidation-susceptible residue", "M", exposure)
            self.assertIsNone(site["relative_sidechain_sasa"])
            self.assertIsNone(site["exposure_percent"])
            self.assertEqual(site["exposure_band"], "unknown")
            self.assertEqual(site["priority"], "review")
            self.assertIn("Unknown", site["risk_label"])
            json.dumps(site, allow_nan=False)

    def test_cdr_impact_does_not_inflate_chemical_priority(self):
        framework = self.score("Oxidation-susceptible residue", "W", 0.15)
        cdr = self.score("Oxidation-susceptible residue", "W", 0.15, "CDR1")
        self.assertEqual(framework["priority"], cdr["priority"])
        self.assertNotEqual(framework["functional_impact"], cdr["functional_impact"])

    def test_charge_changes_with_ph(self):
        self.assertGreater(_sidechain_charge("H", 5), _sidechain_charge("H", 8))
        self.assertLess(_sidechain_charge("D", 8), _sidechain_charge("D", 2))

    def test_surface_and_ptm_mapping(self):
        result = analyze_structure(PDB, {"H": HEAVY, "L": LIGHT})
        self.assertEqual(len(result["residues"]), len(HEAVY) + len(LIGHT))
        self.assertGreater(result["patch_counts"]["hydrophobic"], 0)
        self.assertIn(result["overall_ptm_risk"]["call"], {"Follow-up recommended", "Review sites", "No priority site identified"})
        self.assertTrue(all(0 <= site["relative_sidechain_sasa"] for site in result["ptm_sites"]))
        self.assertTrue(any(site["chain"] == "H" and site["number"] == "H59"
                            and site["priority"] == "high" for site in result["ptm_sites"]))
        heats = [r["surface_heat"] for r in result["residues"]]
        self.assertTrue(all(0 <= value <= 1 for heat in heats for value in heat.values()))
        self.assertGreater(max(h["hydrophobic"] for h in heats), 0.7)
        self.assertEqual(min(h["hydrophobic"] for h in heats), 0)
        json.dumps(result, allow_nan=False)

    def test_kabat_regions_follow_selected_scheme(self):
        imgt = analyze_structure(PDB, {"H": HEAVY, "L": LIGHT}, scheme="imgt")
        kabat = analyze_structure(PDB, {"H": HEAVY, "L": LIGHT}, scheme="kabat")
        def region(data, position):
            return next(r["region"] for r in data["residues"]
                        if r["chain"] == "H" and r["sequence_position"] == position)
        self.assertEqual([region(imgt, p) for p in (30, 34, 62)], ["CDR1", "Framework", "Framework"])
        self.assertEqual([region(kabat, p) for p in (30, 34, 62)], ["Framework", "CDR1", "CDR2"])
        self.assertTrue(all(s["scheme"] == "kabat" for s in kabat["ptm_sites"]))

    def test_overlapping_chemistry_counts_one_physical_site(self):
        flags = [{"chain": "H", "start": 31, "kind": "N-linked glycosylation sequon",
                  "motif": "NST", "priority_points": 4, "priority": "high"},
                 {"chain": "H", "start": 31, "kind": "Asparagine deamidation motif",
                  "motif": "NS", "priority_points": 3, "priority": "review"}]
        grouped = _group_physical_sites(flags)
        self.assertEqual(len(grouped), 1)
        self.assertEqual(len(grouped[0]["flags"]), 2)
        self.assertEqual(grouped[0]["priority"], "high")

    def test_single_residue_difference_preserves_mapped_structure(self):
        altered_light = LIGHT.replace("I", "M", 1)
        result = analyze_structure(PDB, {"H": HEAVY, "L": altered_light})
        self.assertTrue(result["overall_ptm_risk"]["call"].startswith("Incomplete evidence"))
        self.assertEqual(result["chain_mapping"]["L"]["mapped_residues"], len(LIGHT) - 1)
        self.assertTrue(any(s["chain"] == "L" and s["motif"] == "M"
                            for s in result["unassessed_ptm_sites"]))
        self.assertTrue(all(s["exposure_percent"] is None and s["exposure_band"] == "unknown"
                            for s in result["unassessed_ptm_sites"]))
        self.assertEqual(len(result["residues"]), len(HEAVY) + len(LIGHT))

    def test_unrelated_light_chain_is_withheld_but_heavy_is_shown(self):
        result = analyze_structure(PDB, {"H": HEAVY, "L": "A" * len(LIGHT)})
        self.assertEqual(result["chain_mapping"]["L"]["mapped_residues"] < 0.9 * len(LIGHT), True)
        self.assertTrue(result["mapping_warnings"])
        self.assertGreater(len(result["residues"]), 0)

    def test_trace_remains_available_without_sequence_mapping(self):
        trace = extract_trace(PDB)
        self.assertEqual(len(trace["residues"]), len(HEAVY) + len(LIGHT))


if __name__ == "__main__":
    unittest.main()
