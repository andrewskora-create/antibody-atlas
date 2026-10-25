"""Per-position substitution ranking, and batch screening from the browser."""
import time
import unittest

from app.analysis import clean_sequence
from app.cdr_substitutions import propose_at_position
from app.numbering import number_chain
from tests.support import requires_dataset

HEAVY = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
         "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
         "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")


def chains():
    out = {}
    for label, raw in (("H", HEAVY), ("L", LIGHT)):
        sequence = clean_sequence(raw)
        canonical = number_chain(sequence, "imgt")
        out[label] = {"sequence": sequence, "numbering": canonical,
                      "numbering_imgt": canonical}
    return out


class PositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chains = chains()

    def test_every_substitution_is_considered(self):
        result = propose_at_position(self.chains, "H", 99, max_candidates=19)
        self.assertEqual(result["candidate_count"], 19)
        self.assertEqual(result["original_aa"], "W")
        self.assertNotIn("W", {c["replacement_aa"] for c in result["candidates"]})

    def test_a_position_with_no_liability_still_ranks(self):
        """The point of this entry point: the automatic path only covers
        detected liabilities, and the user often wants somewhere else."""
        position = next(i for i, aa in enumerate(self.chains["H"]["sequence"], 1)
                        if aa == "G")
        result = propose_at_position(self.chains, "H", position)
        self.assertTrue(result["candidates"])
        self.assertIn("no detected liability", result["caveat"])

    def test_tryptophan_site_uses_blended_selection_not_a_fixed_replacement(self):
        result = propose_at_position(self.chains, "H", 99, max_candidates=19)
        clean = [c for c in result["candidates"] if c["clean"] and not c["edits_preserved_position"]]
        scores = [c["selection_score"]["score"] for c in clean]
        self.assertEqual(scores, sorted(scores, reverse=True))
        self.assertEqual(len(result["top_model_candidates"]), 5)
        self.assertTrue(all(c["selection_score"]["basis"] in
                            {"similarity_clinical_blend", "similarity_only"}
                            for c in result["candidates"]))

    def test_each_candidate_is_rescanned_for_what_it_creates(self):
        result = propose_at_position(self.chains, "H", 99, max_candidates=19)
        for candidate in result["candidates"]:
            self.assertEqual(candidate["clean"], not candidate["introduces"])

    def test_imgt_position_and_region_are_reported(self):
        result = propose_at_position(self.chains, "H", 99)
        self.assertEqual(result["imgt_number"], "H107")
        self.assertEqual(result["region"], "CDR3")

    def test_bad_input_is_refused_clearly(self):
        with self.assertRaises(ValueError):
            propose_at_position(self.chains, "H", 9999)
        with self.assertRaises(ValueError):
            propose_at_position(self.chains, "Z", 1)
        with self.assertRaises(ValueError):
            propose_at_position({"H": {"sequence": HEAVY}}, "H", 10)


class BatchWebTests(unittest.TestCase):
    """The browser runner. It shares scoring with the command-line runner but
    has its own job plumbing, which is what these cover."""

    PANEL = (f">MAB01|H\n{HEAVY}\n>MAB01|L\n{LIGHT}\n"
             f">MAB02|H\n{HEAVY}\n>MAB02|L\n{LIGHT}\n")

    def _finish(self, job, limit=180):
        from app import batch_web
        for _ in range(limit):
            current = batch_web.get(job["id"])
            if current["status"] == "complete":
                return current
            time.sleep(0.5)
        self.fail("batch did not finish")

    def test_a_panel_runs_and_produces_one_row_per_molecule(self):
        from app import batch_web
        job = self._finish(batch_web.start(
            self.PANEL, {"scheme": "imgt", "ph": 7.4, "precedents": False}))
        self.assertEqual(job["total"], 2)
        self.assertEqual(len(job["rows"]), 2)
        self.assertEqual({r["name"] for r in job["rows"]}, {"MAB01", "MAB02"})
        self.assertTrue(all(r["status"] == "complete" for r in job["results"]))

    def test_a_sequence_that_is_not_an_antibody_fails_alone(self):
        """Parse errors and screening errors are different failures. This one
        is valid protein that ANARCII cannot number, so it must fail by itself
        and leave the rest of the panel's results intact."""
        from app import batch_web
        panel = self.PANEL + ">NOTANAB|H\n" + ("GSGSGSGSAP" * 12) + "\n"
        job = self._finish(batch_web.start(
            panel, {"scheme": "imgt", "ph": 7.4, "precedents": False}))
        statuses = {r["name"]: r["status"] for r in job["results"]}
        self.assertEqual(statuses["NOTANAB"], "failed")
        self.assertEqual(statuses["MAB01"], "complete")
        self.assertEqual(statuses["MAB02"], "complete")
        self.assertEqual(len(job["rows"]), 2)
        self.assertTrue(next(r["error"] for r in job["results"]
                             if r["name"] == "NOTANAB"))

    def test_an_unreadable_residue_rejects_the_panel_before_it_starts(self):
        """A typo is caught up front rather than after a long run, and the
        message names the construct and the line."""
        from app import batch_web
        panel = self.PANEL + ">TYPO|H\n" + HEAVY[:60] + "BOU" + HEAVY[60:] + "\n"
        with self.assertRaises(ValueError) as caught:
            batch_web.start(panel, {"scheme": "imgt", "ph": 7.4, "precedents": False})
        message = str(caught.exception)
        self.assertIn("TYPO", message)
        self.assertIn("B", message)

    def test_the_isotype_moves_the_rows_onto_a_full_length_basis(self):
        from app import batch_web
        fv = self._finish(batch_web.start(
            self.PANEL, {"scheme": "imgt", "ph": 7.4, "precedents": False}))
        full = self._finish(batch_web.start(
            self.PANEL, {"scheme": "imgt", "ph": 7.4, "precedents": False,
                         "isotype": "IGHG1"}))
        self.assertEqual(fv["rows"][0]["charge_basis"], "variable_domain")
        self.assertEqual(full["rows"][0]["charge_basis"], "full_length")
        self.assertNotEqual(fv["rows"][0]["pI"], full["rows"][0]["pI"])

    def test_csv_carries_every_column_any_row_has(self):
        from app import batch_web
        job = self._finish(batch_web.start(
            self.PANEL, {"scheme": "imgt", "ph": 7.4, "precedents": False}))
        csv_text = batch_web.to_csv(job["id"])
        header = csv_text.splitlines()[0].split(",")
        self.assertIn("name", header)
        self.assertIn("residues_notable", header)
        self.assertEqual(len(csv_text.splitlines()), 3)

    def test_an_oversized_panel_is_refused_and_points_at_the_cli(self):
        from app import batch_web
        panel = "".join(f">M{i}|H\n{HEAVY}\n" for i in range(batch_web.MAX_CONSTRUCTS + 5))
        with self.assertRaises(ValueError) as caught:
            batch_web.start(panel, {"scheme": "imgt", "ph": 7.4})
        self.assertIn("run_batch.py", str(caught.exception))

    def test_an_empty_panel_is_refused(self):
        from app import batch_web
        with self.assertRaises(ValueError):
            batch_web.start("   ", {"scheme": "imgt", "ph": 7.4})

    def test_an_unknown_job_is_not_invented(self):
        from app import batch_web
        self.assertIsNone(batch_web.get("nosuchjob"))
        self.assertEqual(batch_web.to_csv("nosuchjob"), "")


if __name__ == "__main__":
    unittest.main()


class PanelFormatTests(unittest.TestCase):
    """Formats real users actually paste.

    A panel labelled the way people write labels -- no '>', the chain word in
    the middle of the phrase, nucleotide sequence -- used to reach the CSV
    reader and fail with "need at least a name column and a heavy column",
    which described a file the user never had.
    """

    @staticmethod
    def _dna(protein, seed=5):
        import random
        from app.translation import CODONS
        back = {}
        for codon, aa in CODONS.items():
            back.setdefault(aa, []).append(codon)
        rng = random.Random(seed)
        return "".join(rng.choice(back[a]) for a in protein)

    def panel(self):
        return (f"PUB_10 variable light chain coding sequence\n{self._dna(LIGHT)}\n"
                f"PUB_11 variable heavy chain coding sequence\n{self._dna(HEAVY)}\n"
                f"PUB_11 variable light chain coding sequence:\n{self._dna(LIGHT)}\n")

    def test_a_pasted_panel_is_not_mistaken_for_csv(self):
        from app.batch_web import _looks_like_csv
        self.assertFalse(_looks_like_csv(self.panel()))
        self.assertTrue(_looks_like_csv("name,heavy,light\nMAB01,EVQL,DIQM\n"))

    def test_headers_without_a_leading_marker_are_read(self):
        from app.batch_input import parse_fasta
        parsed = parse_fasta(self.panel())
        self.assertEqual(parsed.problems, [])
        self.assertEqual(set(parsed.constructs), {"PUB_10", "PUB_11"})

    def test_the_chain_word_may_sit_mid_header(self):
        from app.batch_input import _chain_token
        self.assertEqual(_chain_token("PUB_1 variable heavy chain coding sequence"),
                         ("PUB_1", "H"))
        self.assertEqual(_chain_token("PUB_1 variable light chain coding sequence"),
                         ("PUB_1", "L"))
        # Both chains of one construct must land on the same name.
        heavy = _chain_token("X variable heavy chain coding sequence")[0]
        light = _chain_token("X variable light chain coding sequence")[0]
        self.assertEqual(heavy, light)

    def test_a_header_with_no_chain_word_is_still_reported(self):
        from app.batch_input import _chain_token
        self.assertEqual(_chain_token("MAB01 coding sequence")[1], None)

    def test_nucleotide_records_are_translated_before_being_checked(self):
        from app.batch_input import parse_fasta
        parsed = parse_fasta(self.panel())
        self.assertEqual(len(parsed.constructs["PUB_11"]["H"]), len(HEAVY))
        self.assertEqual(parsed.constructs["PUB_11"]["H"], HEAVY)

    def test_a_light_only_construct_is_screened_rather_than_rejected(self):
        from app.batch_input import parse_fasta
        parsed = parse_fasta(self.panel())
        self.assertIn("PUB_10", parsed.constructs)
        self.assertEqual(sorted(parsed.constructs["PUB_10"]), ["L"])
        self.assertTrue(any("no heavy chain" in w["detail"] for w in parsed.warnings))

    def test_a_light_only_construct_falls_back_to_the_fv_charge_basis(self):
        """Only the charge basis needs a heavy chain; losing the whole
        molecule over it would be the wrong trade."""
        from app import batch_web
        job = batch_web.start(self.panel(), {"scheme": "imgt", "ph": 7.4,
                                             "precedents": False, "isotype": "IGHG1"})
        for _ in range(240):
            current = batch_web.get(job["id"])
            if current["status"] == "complete":
                break
            time.sleep(0.5)
        self.assertTrue(all(r["status"] == "complete" for r in current["results"]))
        rows = {r["name"]: r for r in current["rows"]}
        self.assertEqual(rows["PUB_10"]["charge_basis"], "variable_domain")
        self.assertEqual(rows["PUB_11"]["charge_basis"], "full_length")
        self.assertIn("Full-length basis unavailable", rows["PUB_10"]["charge_basis_note"])


class BatchReportTests(unittest.TestCase):
    """A batch row is a way into the molecule, not a dead end."""

    PANEL = (f">MAB01|H\n{HEAVY}\n>MAB01|L\n{LIGHT}\n"
             f">MAB02|H\n{HEAVY}\n>MAB02|L\n{LIGHT}\n")

    @classmethod
    def setUpClass(cls):
        from app import batch_web
        job = batch_web.start(cls.PANEL, {"scheme": "imgt", "ph": 7.4, "precedents": False})
        for _ in range(240):
            current = batch_web.get(job["id"])
            if current["status"] == "complete":
                break
            time.sleep(0.5)
        cls.job = current

    def test_records_are_kept_so_a_row_can_be_reported_on(self):
        from app import batch_web
        self.assertIsNotNone(batch_web.record_for(self.job["id"], "MAB01"))
        self.assertIsNone(batch_web.record_for(self.job["id"], "NOPE"))

    def test_the_polled_status_does_not_carry_every_record(self):
        """The page polls this once a second; records are large."""
        self.assertNotIn("records", self.job)

    def test_input_sequences_are_recoverable_for_the_full_screen(self):
        from app import batch_web
        chains = batch_web.sequences_for(self.job["id"], "MAB01")
        self.assertEqual(sorted(chains), ["H", "L"])
        self.assertEqual(chains["H"], HEAVY)

    def test_a_full_report_can_be_built_from_a_batch_record(self):
        from app import batch_web
        from app.report import make_report
        pdf = make_report(batch_web.record_for(self.job["id"], "MAB01"), title="MAB01")
        self.assertTrue(pdf.startswith(b"%PDF"))
        import re as _re
        self.assertEqual(len(_re.findall(rb"/Type\s*/Page\b", pdf)), 5)
