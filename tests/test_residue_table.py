"""The consolidated per-residue table and the molecule-level flags.

Numbering dictionaries are built by hand here rather than by calling ANARCII,
so the joins and the verdict policy are tested in isolation and the suite stays
fast.
"""
import unittest

from app.attribution import BY_ID
from app.molecule_flags import molecule_flags
from app.residue_table import build_residue_table


def numbering(chain_type, rows, scheme="imgt", definition="IMGT"):
    return {"chain_type": chain_type, "scheme": scheme, "definition": definition,
            "residues": [{"sequence_position": i, "number": f"{chain_type}{n}",
                          "number_base": int(str(n).rstrip("ABC") or 0), "aa": aa,
                          "region": region}
                         for i, (n, aa, region) in enumerate(rows, 1)]}


# IMGT 1-3 framework, 36 inside CDR1, 110 inside CDR3.
IMGT_ROWS = [(1, "E", "Framework"), (2, "V", "Framework"), (3, "Q", "Framework"),
             (36, "D", "CDR1"), (37, "T", "CDR1"), (110, "N", "CDR3"), (111, "G", "CDR3")]
KABAT_ROWS = [(1, "E", "Framework"), (2, "V", "Framework"), (3, "Q", "Framework"),
              (31, "D", "CDR1"), (32, "T", "CDR1"), (98, "N", "CDR3"), (99, "G", "CDR3")]


def chain(liabilities=(), scheme_rows=None):
    imgt = numbering("H", IMGT_ROWS)
    selected = numbering("H", scheme_rows or KABAT_ROWS, "kabat", "KABAT") if scheme_rows is not None \
        else numbering("H", KABAT_ROWS, "kabat", "KABAT")
    return {"sequence": "EVQDTNG", "numbering": selected, "numbering_imgt": imgt,
            "humanness": {"nearest_human_v_gene": None},
            "metrics": {"ptm_sites": list(liabilities)}}


class ResidueTableTests(unittest.TestCase):
    def test_rows_are_keyed_to_imgt_and_carry_the_second_scheme(self):
        table = build_residue_table({"H": chain()})
        row = next(r for r in table["rows"] if r["sequence_position"] == 4)
        self.assertEqual(row["imgt"], "36")
        self.assertEqual(row["imgt_region"], "CDR1")
        self.assertEqual(row["alt_number"], "31")
        self.assertEqual(row["alt_scheme"], "KABAT")

    def test_region_disagreement_between_schemes_is_preserved_not_reconciled(self):
        """A position may be framework in IMGT and CDR in Kabat. Both are shown."""
        kabat = [(1, "E", "Framework"), (2, "V", "Framework"), (3, "Q", "Framework"),
                 (31, "D", "CDR1"), (32, "T", "CDR1"), (98, "N", "CDR3"), (99, "G", "CDR3")]
        kabat[0] = (1, "E", "CDR1")
        table = build_residue_table({"H": chain(scheme_rows=kabat)})
        row = table["rows"][0]
        self.assertEqual(row["imgt_region"], "Framework")
        self.assertEqual(row["alt_region"], "CDR1")

    def test_cdr3_has_no_germline_reference_and_says_why(self):
        table = build_residue_table({"H": chain()})
        row = next(r for r in table["rows"] if r["imgt"] == "110")
        self.assertIsNone(row["germline_match"])
        self.assertIsNone(row["germline_status"])
        self.assertIn("junctional", row["germline_note"])

    def test_sequence_only_liability_is_watch_and_declares_its_basis(self):
        table = build_residue_table({"H": chain([
            {"start": 6, "kind": "Asparagine deamidation motif", "motif": "NG"}])})
        row = next(r for r in table["rows"] if r["sequence_position"] == 6)
        self.assertEqual(row["verdict"], "watch")
        self.assertEqual(row["verdict_basis"], ["sequence"])
        self.assertTrue(row["notable"])

    def test_positions_without_a_liability_are_clear_not_unassessed(self):
        table = build_residue_table({"H": chain()})
        self.assertTrue(all(r["verdict"] == "clear" for r in table["rows"]))
        self.assertEqual(table["notable_rows"], 0)
        self.assertNotIn("unassessed", table["verdict_legend"])

    def test_exposed_high_priority_liability_becomes_fix(self):
        structure = {"ptm_sites": [{"chain": "H", "sequence_position": 6, "priority": "high",
                                    "relative_sidechain_sasa": 0.44, "region": "CDR3",
                                    "kind": "Asparagine deamidation motif", "motif": "NG"}]}
        table = build_residue_table({"H": chain()}, structure=structure)
        row = next(r for r in table["rows"] if r["sequence_position"] == 6)
        self.assertEqual(row["verdict"], "evaluate")
        self.assertIn("exposure", row["verdict_basis"])
        self.assertEqual(row["exposure"], 0.44)

    def test_a_contact_does_not_suppress_the_chemistry_at_the_same_position(self):
        """Chemistry and binding are separate axes, reported separately.

        This used to assert the opposite: a predicted contact returned
        "preserve" before any liability was read, so a high-priority
        deamidation motif sitting on a contact residue showed no chemical
        call at all. A residue can be both important to preserve and a
        serious degradation concern, and one word cannot carry both. The
        call now describes the chemistry; the contact is its own column,
        and the step to take is resolved from the two together.
        """
        structure = {"ptm_sites": [{"chain": "H", "sequence_position": 6, "priority": "high",
                                    "relative_sidechain_sasa": 0.44, "region": "CDR3",
                                    "kind": "Asparagine deamidation motif", "motif": "NG"}]}
        paratope = {"status": "available",
                    "by_chain": {"H": {6: {"imgt": "110", "aa": "N", "probability": 0.91,
                                           "likely_contact": True}}}}
        table = build_residue_table({"H": chain()}, structure=structure, paratope=paratope)
        row = next(r for r in table["rows"] if r["sequence_position"] == 6)
        self.assertEqual(row["verdict"], "evaluate")
        self.assertTrue(row["paratope_contact"])
        self.assertEqual(row["paratope_score"], 0.91)
        self.assertIn("edits the paratope", row["next_step"])

    def test_a_contact_with_no_liability_is_still_preserve(self):
        paratope = {"status": "available",
                    "by_chain": {"H": {6: {"imgt": "110", "aa": "N", "probability": 0.91,
                                           "likely_contact": True}}}}
        table = build_residue_table({"H": chain()}, paratope=paratope)
        row = next(r for r in table["rows"] if r["sequence_position"] == 6)
        self.assertEqual(row["verdict"], "preserve")

    def test_unmappable_paratope_entries_are_counted_not_silently_dropped(self):
        paratope = {"status": "available",
                    "by_chain": {"H": {6: {"probability": 0.2, "likely_contact": False},
                                       "IMGT112": {"probability": 0.8, "likely_contact": True}}}}
        table = build_residue_table({"H": chain()}, paratope=paratope)
        self.assertEqual(table["paratope_positions_unmapped"], 1)

    def test_absent_layers_are_reported(self):
        table = build_residue_table({"H": chain()})
        self.assertEqual(table["missing_layers"], ["exposure", "germline", "paratope"])

    def test_canonical_imgt_numbering_is_required(self):
        bad = chain()
        bad["numbering_imgt"] = numbering("H", KABAT_ROWS, "kabat", "KABAT")
        with self.assertRaises(ValueError):
            build_residue_table({"H": bad})


def flag_chain(sequence, cdr_aas=("D", "T"), label="H"):
    rows = [(1, sequence[0], "Framework")] + \
           [(36 + i, aa, "CDR1") for i, aa in enumerate(cdr_aas)]
    return {"sequence": sequence, "numbering": numbering(label, rows),
            "numbering_imgt": numbering(label, rows),
            "metrics": {"aggregation_proxy": {"patch_count": 0}}}


class MoleculeFlagTests(unittest.TestCase):
    def test_formulation_ph_sitting_on_the_pi_is_red(self):
        # A sequence with pI near 6 evaluated at a formulation pH of 6.
        result = molecule_flags({"H": flag_chain("DEHAAAK")}, formulation_ph=6.0)
        flag = next(f for f in result["flags"] if f["id"] == "colloidal_ph_margin")
        self.assertLess(flag["value"], 1.5)
        self.assertIn(flag["verdict"], {"amber", "red"})
        self.assertEqual(flag["outcomes"], ["aggregation"])

    def test_charge_extreme_fires_and_emits_the_clearance_viscosity_tradeoff(self):
        result = molecule_flags({"H": flag_chain("KKKKKKKKKKRRRR")})
        flag = next(f for f in result["flags"] if f["id"] == "net_charge")
        self.assertEqual(flag["verdict"], "red")
        self.assertIn("clearance", flag["outcomes"])
        self.assertIn("nonspecific_binding", flag["outcomes"])
        self.assertTrue(any("viscosity" in t["between"] for t in result["tradeoffs"]))

    def test_negative_charge_is_flagged_without_claiming_nonspecific_binding(self):
        """Clearance is non-monotonic, so a low pI is flagged too -- but a
        negative Fv is not asserted to bind non-specifically."""
        result = molecule_flags({"H": flag_chain("DDDDDDDDDDEEEE")})
        flag = next(f for f in result["flags"] if f["id"] == "net_charge")
        self.assertEqual(flag["verdict"], "red")
        self.assertIn("clearance", flag["outcomes"])
        self.assertNotIn("nonspecific_binding", flag["outcomes"])

    def test_hydrophobic_cdrs_are_flagged_with_per_cdr_values(self):
        result = molecule_flags({"H": flag_chain("AIIIVVV", cdr_aas=("I", "V", "L", "F"))})
        flag = next(f for f in result["flags"] if f["id"] == "cdr_hydrophobicity")
        self.assertEqual(flag["verdict"], "red")
        self.assertIn("H-CDR1", flag["per_cdr"])

    def test_immunogenicity_is_never_an_outcome(self):
        result = molecule_flags({"H": flag_chain("KKKKKKKKKKRRRR")})
        self.assertNotIn("immunogenicity", result["outcomes"])
        for flag in result["flags"]:
            self.assertNotIn("immunogenicity", flag["outcomes"])
        self.assertIn("immunogenicity", result["excluded"])

    def test_every_flag_is_attributed_and_marked_variable_domain_only(self):
        result = molecule_flags({"H": flag_chain("EVQDTNG"), "L": flag_chain("DIQMTQS")})
        self.assertTrue(result["flags"])
        for flag in result["flags"]:
            self.assertTrue(flag["source_ids"], flag["id"])
            for source_id in flag["source_ids"]:
                self.assertIn(source_id, BY_ID, f"{flag['id']} cites unregistered {source_id}")
            self.assertTrue(flag["fv_only"])
            self.assertIn("policy", flag)

    def test_suggested_formulation_window_keeps_its_distance_from_the_pi(self):
        result = molecule_flags({"H": flag_chain("EVQDTNGKK")})
        suggestion = result["formulation_suggestion"]
        if suggestion:
            pi = result["fv"]["pi"]
            for bound in suggestion["ph_range"]:
                self.assertGreaterEqual(abs(bound - pi), 1.5)

    def test_rejects_impossible_ph(self):
        with self.assertRaises(ValueError):
            molecule_flags({"H": flag_chain("EVQDTNG")}, formulation_ph=20)
        with self.assertRaises(ValueError):
            molecule_flags({})


if __name__ == "__main__":
    unittest.main()


class AnalyzeEndpointTests(unittest.TestCase):
    """The table and the flags must reach the API response, not just exist."""

    HEAVY = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
             "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
    LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
             "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")

    @classmethod
    def setUpClass(cls):
        from fastapi.testclient import TestClient
        from app.main import app
        cls.client = TestClient(app)
        cls.body = cls.client.post("/api/analyze", json={
            "heavy": cls.HEAVY, "light": cls.LIGHT, "scheme": "kabat",
            "ph": 7.4, "formulation_ph": 6.0}).json()

    def test_table_is_returned_with_both_schemes_per_row(self):
        rows = self.body["residue_table"]["rows"]
        self.assertTrue(rows)
        self.assertTrue(all(r["imgt"] for r in rows))
        self.assertTrue(any(r["alt_number"] for r in rows))
        self.assertTrue(any(r["alt_region"] != r["imgt_region"] for r in rows),
                        "IMGT and Kabat should disagree on region somewhere in a real Fv")

    def test_sequence_only_response_claims_no_exposure_or_paratope(self):
        table = self.body["residue_table"]
        self.assertEqual(table["missing_layers"], ["exposure", "paratope"])
        self.assertEqual(table["verdict_counts"]["evaluate"], 0)
        self.assertEqual(table["verdict_counts"]["preserve"], 0)
        self.assertTrue(all(r["verdict_basis"] == ["sequence"] for r in table["rows"]))

    def test_flags_are_returned_and_exclude_immunogenicity(self):
        flags = self.body["molecule_flags"]
        self.assertTrue(flags["flags"])
        self.assertNotIn("immunogenicity", flags["outcomes"])
        self.assertIn("formulation_ph", flags)

    def test_impossible_formulation_ph_is_rejected_by_the_schema(self):
        response = self.client.post("/api/analyze",
                                    json={"heavy": self.HEAVY, "formulation_ph": 99})
        self.assertEqual(response.status_code, 422)
