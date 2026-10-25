"""The report has to say what each number is measured against.

Three things in the output were read as meaning something they do not:

  * "N notable" in the residue table is a union of three unrelated reasons.
    Most notable rows carry no chemical liability at all -- they are predicted
    antigen contacts, flagged so that a fix proposed elsewhere is not quietly
    applied to the binding site. Reported as one total, that looks like a
    screen crying wolf.
  * A green or red molecule-level flag did not say which boundary the value
    cleared, so the colour was unfalsifiable from the page.
  * Two sections quote percentiles against two different populations -- the
    natural repertoire for the sequence-only descriptors, the clinical panel
    for the panel comparison -- and nothing said so.

These tests hold the data shape that the UI needs to say all three.
"""
from __future__ import annotations

import json
import re
import unittest
from pathlib import Path

from app.molecule_flags import molecule_flags
from app.residue_table import (build_residue_table, precedent_index,
                               _motif_is_germline)
from tests.test_residue_table import chain as numbered_chain

ROOT = Path(__file__).resolve().parents[1]
START = ROOT / "static" / "start.html"

HEAVY = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
         "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGYFDYWGQGTLVTVSS")
LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSG"
         "SGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGQGTKVEIK")


def _chains() -> dict:
    return {"H": {"sequence": HEAVY, "metrics": {}},
            "L": {"sequence": LIGHT, "metrics": {}}}


class NotableReasonsTests(unittest.TestCase):
    """'Notable' is three reasons, and the reader needs them apart."""

    def setUp(self):
        # One real liability at sequence position 6, numbered by hand so the
        # join is exercised without calling ANARCII.
        self.table = build_residue_table({"H": numbered_chain(
            [{"start": 6, "kind": "Asparagine deamidation motif", "motif": "NG"}])})

    def test_a_liability_is_actually_counted(self):
        self.assertEqual(self.table["notable_reasons"]["liability"], 1)
        self.assertEqual(self.table["notable_reasons"]["contact"], 0)

    def test_the_reasons_are_reported_separately(self):
        reasons = self.table["notable_reasons"]
        self.assertEqual(
            set(reasons),
            {"liability", "contact", "framework_mismatch",
             "liability_only", "contact_only", "both"})

    def test_the_split_adds_back_up(self):
        r = self.table["notable_reasons"]
        self.assertEqual(r["liability_only"] + r["both"], r["liability"])
        self.assertEqual(r["contact_only"] + r["both"], r["contact"])

    def test_no_reason_can_exceed_the_union(self):
        r = self.table["notable_reasons"]
        for key in ("liability", "contact", "framework_mismatch"):
            self.assertLessEqual(r[key], self.table["notable_rows"], key)

    def test_the_legend_says_a_notable_row_need_not_carry_a_liability(self):
        legend = self.table["notable_legend"].lower()
        self.assertIn("contact", legend)
        self.assertIn("liability", legend)


class MoleculeFlagBandTests(unittest.TestCase):
    """A colour the page cannot justify is a colour the reader cannot check."""

    def setUp(self):
        self.flags = molecule_flags(_chains())["flags"]

    def test_every_assessed_flag_carries_its_boundaries(self):
        for flag in self.flags:
            if flag["verdict"] == "unknown":
                continue
            with self.subTest(flag=flag["id"]):
                self.assertIn("bands", flag)
                self.assertTrue(flag["bands"]["measure"])
                self.assertTrue(flag["bands"]["cuts"])

    def test_the_verdict_is_always_one_of_the_bands_shown(self):
        # The UI marks the live band by matching the verdict string. If no
        # band matched, the row would show boundaries with none highlighted
        # and the colour would again be unexplained.
        for flag in self.flags:
            if flag["verdict"] == "unknown":
                continue
            with self.subTest(flag=flag["id"]):
                self.assertIn(flag["verdict"],
                              [verdict for verdict, _rule in flag["bands"]["cuts"]])

    def test_an_unassessed_flag_claims_no_bands(self):
        for flag in self.flags:
            if flag["verdict"] == "unknown":
                self.assertIsNone(flag.get("bands"), flag["id"])

    def test_the_payload_survives_the_trip_to_the_browser(self):
        json.dumps(self.flags)


class ReferencePopulationTests(unittest.TestCase):
    """Each percentile has to name the population it is a percentile of."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_both_percentile_sections_carry_a_reference_box(self):
        self.assertGreaterEqual(self.html.count('class="refbox"'), 2)

    def test_the_sequence_descriptors_disclaim_the_clinical_panel(self):
        self.assertIn("This is not the clinical panel", self.html)

    def test_the_panel_section_disclaims_the_repertoire(self):
        self.assertIn("composition-based reference comparisons, not surface-patch measurements",
                      self.html)

    def test_every_residue_table_column_has_a_definition(self):
        # The headers are assembled from one list; a column added there
        # without an entry in RT_COLS would render as a bare word again.
        block = re.search(r"const RT_COLS = \{(.*?)\n\};", self.html, re.S)
        self.assertIsNotNone(block)
        defined = set(re.findall(r"^\s*'([^']+)':", block.group(1), re.M))
        emitted = {"Chain", "Pos", "AA", "IMGT", "Region", "Germline", "Exposure",
                   "Contact", "Predicted contact", "To epitope", "Call", "Liability",
                   "Next step"}
        self.assertEqual(emitted - defined, set())


class SplitTableTests(unittest.TestCase):
    """Chemistry and binding surface are different questions.

    One table with a filter made them look like the same finding, which is
    what produced 'why are so many residues notable with no PTM impact'.
    They are now separate tables with separate columns.
    """

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_three_tables_exist(self):
        for fn in ("renderLiabilityTable", "renderParatopeResidueTable",
                   "renderAllResiduesTable"):
            self.assertIn(f"function {fn}(", self.html)

    def test_the_single_combined_table_is_gone(self):
        self.assertNotIn("renderResidueTable", self.html)

    def test_the_liability_table_keeps_the_contact_flag(self):
        # Whether a liability sits on the binding surface is the first thing
        # that decides how to fix it, so it stays in the liability table even
        # though contacts have a table of their own.
        self.assertIn("'Predicted contact'", self.html)

    def test_the_binding_surface_table_says_nothing_is_wrong(self):
        self.assertIn("Nothing is wrong with these", self.html)


class PanelTests(unittest.TestCase):
    """The report is a search, not a scroll."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_sections_are_collapsed_into_panels(self):
        self.assertIn("function panelise(", self.html)
        self.assertIn("$('#results').innerHTML = panelise(out);", self.html)

    def test_only_a_couple_of_sections_open_by_default(self):
        # The point of the accordion is that the page opens short. If this
        # set grows, it stops being an accordion.
        block = re.search(r"const PANEL_OPEN_BY_DEFAULT = new Set\(\[(.*?)\]\);",
                          self.html, re.S)
        self.assertIsNotNone(block)
        self.assertLessEqual(len(re.findall(r"'([^']+)'", block.group(1))), 2)

    def test_printing_reveals_every_section(self):
        # A collapsed <details> hides its body by UA stylesheet, so both the
        # print rule and the beforeprint handler are load-bearing for the PDF.
        self.assertIn("details.panel>.p-body{display:block !important}", self.html)
        self.assertIn("'beforeprint'", self.html)
        self.assertIn("'afterprint'", self.html)

    def test_the_filter_wiring_covers_every_bar(self):
        # There is more than one .rt table now; a single querySelector would
        # wire the first bar and leave the rest dead.
        self.assertIn("document.querySelectorAll('.rt-bar').forEach", self.html)
        self.assertNotIn("const bar = document.querySelector('.rt-bar');", self.html)


class CascadeTests(unittest.TestCase):
    """Two rules claimed the same selector and the loser was the new layout."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_verdict_is_declared_once(self):
        # `.verdict{display:flex}` from the old badge row overrode the new
        # verdict masthead and laid its heading, call, tiles and key numbers
        # out side by side. Chrome shrank them; Safari overflowed.
        self.assertEqual(len(re.findall(r"^\s*\.verdict\{", self.html, re.M)), 1)

    def test_the_pill_rule_cannot_reach_the_key_numbers(self):
        # `.v` is the verdict pill, but `.keynum .v` and `.cell .v` are also
        # `.v` and were picking up its inline-block pill padding.
        self.assertNotIn("\n  .v{", self.html)
        self.assertIn("span.v{", self.html)

    def test_the_basis_badge_font_shorthand_is_valid(self):
        # `font: 700 9.5px/1.6 inherit` is invalid -- `inherit` is not a font
        # family -- so the whole declaration was dropped.
        self.assertNotIn("font:700 9.5px/1.6 inherit", self.html)


if __name__ == "__main__":
    unittest.main()


class GermlineAndPrecedentTests(unittest.TestCase):
    """A germline framework motif is not a defect this molecule introduced.

    The call used to rest on priority and exposure alone, so a deamidation
    motif that every antibody on the assigned V allele carries came back as
    "evaluate" -- a substitution that would move the framework AWAY from human
    germline to remove a motif human germline puts there. Two pieces of
    evidence now reach the verdict, and both argue only for framework
    positions: inside a CDR, germline encoding says little, because the CDR
    is what affinity maturation changed.
    """

    KIND = "Asparagine deamidation motif"

    def _structure(self, position, priority="high", motif="EV"):
        return {"ptm_sites": [{"chain": "H", "sequence_position": position,
                               "kind": self.KIND, "motif": motif, "priority": priority,
                               "relative_sidechain_sasa": 0.5}]}

    def _precedents(self, position, region="Framework", matches=(("DrugA", "Approved"),)):
        return {"sites": [{"query_chain": "H", "position": position, "kind": self.KIND,
                           "query_region": region,
                           "matches": [{"drug": d, "stage": s} for d, s in matches]}]}

    def _row(self, position, structure=None, precedents=None):
        table = build_residue_table({"H": numbered_chain()}, structure=structure,
                                    precedents=precedents)
        return table, next(r for r in table["rows"]
                           if r["sequence_position"] == position)

    # -- the motif-span germline test ------------------------------------
    def test_a_motif_is_germline_only_if_every_residue_is(self):
        self.assertTrue(_motif_is_germline(5, "NG", {5: True, 6: True}))
        self.assertFalse(_motif_is_germline(5, "NG", {5: True, 6: False}))

    def test_unknown_anywhere_in_the_span_makes_the_answer_unknown(self):
        # "we could not tell" and "it differs" lead to different decisions,
        # so they must not collapse into the same False.
        self.assertIsNone(_motif_is_germline(5, "NG", {5: True, 6: None}))
        self.assertIsNone(_motif_is_germline(5, "NG", {5: True}))
        self.assertIsNone(_motif_is_germline(5, "NG", None))

    # -- the precedent join ----------------------------------------------
    def test_one_drug_matching_twice_is_counted_and_named_once(self):
        index = precedent_index(self._precedents(
            1, matches=(("DrugA", "Approved"), ("DrugA", "Approved"),
                        ("DrugB", "Phase-III"))))
        entry = index[("H", 1, self.KIND)]
        self.assertEqual((entry["approved"], entry["phase_iii"]), (1, 1))
        self.assertEqual(entry["drugs"], ["DrugA (Approved)", "DrugB (Phase-III)"])

    def test_the_precedent_reaches_the_row(self):
        # It was being computed and attached to the fix proposals instead, so
        # the table where the call is made never saw it.
        _table, row = self._row(1, self._structure(1), self._precedents(1))
        self.assertEqual(row["precedent"]["approved"], 1)

    def test_an_unsearched_position_is_not_a_zero_result(self):
        _table, row = self._row(1, self._structure(1))
        self.assertIsNone(row["precedent"])

    # -- the demotion ------------------------------------------------------
    def test_an_exposed_framework_liability_is_a_fix_by_default(self):
        _table, row = self._row(1, self._structure(1))
        self.assertEqual(row["verdict"], "evaluate")

    def test_clinical_precedent_demotes_a_framework_fix_to_watch(self):
        _table, row = self._row(1, self._structure(1), self._precedents(1))
        self.assertEqual(row["verdict"], "watch")
        self.assertIn("precedent", row["verdict_basis"])
        self.assertIn("approved or Phase-III", row["verdict_reason"])

    def test_the_demotion_says_it_is_evidence_not_proof(self):
        _table, row = self._row(1, self._structure(1), self._precedents(1))
        self.assertIn("not proof", row["verdict_reason"])

    def test_a_cdr_liability_is_never_demoted(self):
        _table, row = self._row(4, self._structure(4, motif="DT"),
                                self._precedents(4, region="CDR1"))
        self.assertEqual(row["imgt_region"], "CDR1")
        self.assertEqual(row["verdict"], "evaluate")

    def test_a_searched_but_empty_precedent_softens_nothing(self):
        _table, row = self._row(1, self._structure(1), self._precedents(1, matches=()))
        self.assertEqual(row["verdict"], "evaluate")
        self.assertEqual(row["precedent"]["approved"], 0)

    def test_one_unsupported_liability_keeps_the_whole_position_at_fix(self):
        # A position can carry two chemistries. Letting a germline-encoded
        # deamidation motif soften the position would hide an introduced
        # oxidation site sitting on the same residue.
        structure = {"ptm_sites": [
            {"chain": "H", "sequence_position": 1, "kind": self.KIND, "motif": "EV",
             "priority": "high", "relative_sidechain_sasa": 0.5},
            {"chain": "H", "sequence_position": 1, "kind": "Methionine oxidation",
             "motif": "EV", "priority": "high", "relative_sidechain_sasa": 0.5}]}
        table = build_residue_table({"H": numbered_chain()}, structure=structure,
                                    precedents=self._precedents(1))
        row = next(r for r in table["rows"] if r["sequence_position"] == 1)
        self.assertEqual(row["verdict"], "evaluate")

    def test_the_reason_names_each_motif_it_softened(self):
        _table, row = self._row(1, self._structure(1), self._precedents(1))
        self.assertIn("EV is", row["verdict_reason"])

    def test_the_reason_warns_that_precedent_does_not_transfer(self):
        # A motif an approved antibody tolerates can still matter in another
        # candidate; the sentence has to say so.
        _table, row = self._row(1, self._structure(1), self._precedents(1))
        self.assertIn("can still matter in a different molecule", row["verdict_reason"])

    def test_the_counts_are_reported_for_the_interface(self):
        table, _row = self._row(1, self._structure(1), self._precedents(1))
        self.assertEqual(table["evidence_counts"]["with_precedent"], 1)


class DescriptorLabellingTests(unittest.TestCase):
    """An acronym with no expansion is not a label."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_gravy_is_spelled_out(self):
        self.assertIn("Grand Average of Hydropathy", self.html)

    def test_linear_surface_proxies_are_not_displayed(self):
        self.assertNotIn("Hydrophobicity appears twice, on purpose", self.html)
        self.assertIn("const rows = displayedMoleculeFlags(mf).map", self.html)
        self.assertIn("const BAND_ORDER = ['isoelectric_point','net_charge_at_ph'];", self.html)

    def test_the_structure_is_not_buried_in_go_deeper(self):
        self.assertIn("<h2>3D structure", self.html)
        self.assertIn("'3d-structure'", self.html)
        self.assertLess(self.html.index("<h2>3D structure"),
                        self.html.index("out.push('<h2>Go deeper</h2>')"))

    def test_the_viewer_still_draws_when_its_panel_starts_open(self):
        # It rendered lazily on the details 'toggle' event; a panel that
        # starts open never fires one.
        self.assertIn("drawOnce();", self.html)


class GlyphTests(unittest.TestCase):
    """A symbol the page does not define is an acronym it does not expand."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_demotion_arrow_has_its_own_hover_target(self):
        # The reader points at the symbol, not at the word beside it.
        self.assertIn('class="v-down" title="Demoted from evaluate', self.html)

    def test_the_arrow_is_only_drawn_where_it_is_also_defined(self):
        # Both come out of rtBadgeKey, which is given the rows so it can
        # count them; the key appears exactly when some row carries one.
        self.assertIn("Call ↓</b> on ${softened.length}", self.html)
        self.assertIn("const softened = (rows||[]).filter", self.html)

    def test_the_key_distinguishes_a_plain_watch(self):
        self.assertIn("is the ordinary case", self.html)


class ViewerLegendTests(unittest.TestCase):
    """A coloured structure with no key is a picture, not a result."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_scales_are_declared_once(self):
        # The legend and the atom colours must come from the same numbers;
        # a key that drifts from the screen is a wrong answer the reader
        # cannot check.
        self.assertIn("const HYDRO_SCALE = {", self.html)
        self.assertIn("const CHARGE_SCALE = {", self.html)
        self.assertIn("bandColour(HYDRO_SCALE, hydroScore(e.aa, e.rsa))", self.html)
        self.assertIn("bandColour(CHARGE_SCALE, chargeScore(e.aa, e.rsa))", self.html)

    def test_no_threshold_is_written_twice(self):
        # Each cut-off should appear only in its scale definition.
        for cut in ("0.45", "0.28", "0.14"):
            self.assertEqual(self.html.count(f"above: {cut}"), 1, cut)

    def test_every_mode_has_a_key(self):
        for mode in ("ptm", "cdr", "plain"):
            self.assertIn(f"{mode}: {{title:", self.html)

    def test_the_key_says_the_scale_is_exposure_weighted(self):
        # The one thing that explains why a cartoon can look alarming while
        # the solvent-accessible surface is mild.
        self.assertIn("multiplied by relative side-chain exposure", self.html)
        self.assertIn("buried hydrophobe stays cool", self.html)

    def test_the_key_follows_the_mode(self):
        self.assertIn("box.innerHTML = viewerLegend(mode,", self.html)


class SpatialEvidenceTests(unittest.TestCase):
    def test_surface_followup_comes_from_structural_categories(self):
        html = START.read_text(encoding="utf-8")
        shortlist = html.split('function renderShortlist(', 1)[1].split('// ---- liability fix', 1)[0]
        self.assertIn('surfaceCategories', shortlist)
        self.assertNotIn('flag.action', shortlist)
        self.assertIn('Structural coverage is incomplete', shortlist)
        self.assertIn('Suggested assay and decision it informs', shortlist)


class ChargeSignTests(unittest.TestCase):
    """A gap between two charges does not tell you their signs.

    The flag described any difference of 6 e or more as "opposite charge on
    the two chains". +10 and +3 differ by 7 with both strongly positive, so
    the sentence could misdescribe the molecule and point at the wrong
    formulation lever.
    """

    def _flag(self, heavy, light):
        flags = molecule_flags({"H": {"sequence": heavy, "metrics": {}},
                                "L": {"sequence": light, "metrics": {}}})["flags"]
        return next(f for f in flags if f["id"] == "charge_asymmetry")

    POSITIVE = "K" * 14 + "GSGSGSGSGS"
    MILD = "K" * 4 + "GSGSGSGSGS"
    NEGATIVE = "E" * 14 + "GSGSGSGSGS"

    def test_a_same_sign_gap_is_not_called_opposite(self):
        flag = self._flag(self.POSITIVE, self.MILD)
        self.assertGreaterEqual(flag["value"], 6)
        self.assertFalse(flag["opposite_sign"])
        self.assertNotIn("opposite", flag["conclusion"].lower())

    def test_a_same_sign_gap_is_described_as_an_imbalance(self):
        flag = self._flag(self.POSITIVE, self.MILD)
        self.assertIn("imbalance rather than a dipole", flag["conclusion"])

    def test_the_dipole_claim_is_weakened_when_it_is_not_a_dipole(self):
        # The viscosity association in the literature is for a dipole.
        flag = self._flag(self.POSITIVE, self.MILD)
        self.assertIn("weaker", flag["conclusion"])

    def test_genuinely_opposite_signs_are_still_called_opposite(self):
        flag = self._flag(self.POSITIVE, self.NEGATIVE)
        self.assertTrue(flag["opposite_sign"])
        self.assertIn("opposite charge", flag["conclusion"])

    def test_the_per_chain_charges_are_carried_not_just_the_gap(self):
        flag = self._flag(self.POSITIVE, self.NEGATIVE)
        self.assertEqual(set(flag["per_chain_charge"]), {"H", "L"})


class UnassessedColourTests(unittest.TestCase):
    """A placeholder must not look like a measurement.

    Missing exposure used to fall back to rsa 0.4, which painted an
    unassessed residue in the middle of the ramp -- indistinguishable by eye
    from a real, mild result.
    """

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_there_is_no_numeric_fallback_for_missing_exposure(self):
        self.assertNotIn("typeof e.rsa === 'number' ? e.rsa : 0.4", self.html)

    def test_unreached_residues_get_their_own_colour(self):
        self.assertIn("const NOT_ASSESSED = 0x9aa0a6;", self.html)
        self.assertIn("if(!e) return NOT_ASSESSED;", self.html)
        self.assertIn("countUnassessedAtoms(component.structure, index, mode)", self.html)

    def test_that_colour_is_outside_every_ramp(self):
        # Otherwise "not assessed" would still read as a position on the scale.
        for scale in ("HYDRO_SCALE", "CHARGE_SCALE"):
            block = re.search(scale + r" = \{(.*?)\n\};", self.html, re.S)
            self.assertIsNotNone(block, scale)
            self.assertNotIn("9aa0a6", block.group(1))

    def test_glycine_is_not_counted_as_unassessed(self):
        # No side chain means no side-chain SASA; that is a fact about
        # glycine, not a gap in the analysis.
        self.assertIn("if(e.aa !== 'G')", self.html)

    def test_the_key_says_grey_is_not_a_low_value(self):
        self.assertIn("not a low value", self.html)


class ChemistryVersusBindingTests(unittest.TestCase):
    """Two questions, two columns, one step resolved from both."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_call_is_documented_as_chemistry_only(self):
        self.assertIn("What the CHEMISTRY at this position calls for", self.html)

    def test_the_liability_table_carries_the_next_step(self):
        self.assertIn("'Next step'", self.html)
        self.assertIn("case 'Next step':", self.html)

    def test_fix_is_gone_as_a_verdict(self):
        # "evaluate": the output of a screening heuristic over a predicted
        # structure is an experiment to run, not a mutation to make.
        self.assertNotIn(".v-fix{", self.html)
        self.assertNotIn("'fix'", self.html)

    def test_the_next_step_table_covers_every_combination(self):
        from app.residue_table import NEXT_STEPS, VERDICTS
        for verdict in VERDICTS:
            self.assertTrue(
                any(key[0] == verdict for key in NEXT_STEPS),
                f"no next step defined for {verdict}")


class ShortlistTests(unittest.TestCase):
    """The output of a screen should be a list of measurements."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_shortlist_exists_and_runs_after_the_evidence(self):
        self.assertIn("function renderShortlist(", self.html)
        self.assertLess(self.html.index("out.push(renderShortlist("),
                        self.html.index("out.push('<h2>Go deeper</h2>')"))

    def test_it_invents_no_measurement(self):
        # A concern whose module states no experiment is omitted, not given
        # a plausible-sounding one.
        self.assertIn("category.action ||", self.html)

    def test_it_says_absence_is_not_a_pass(self):
        self.assertIn("Absence from this list is not a pass", self.html)

    def test_a_superseded_sequence_proxy_is_demoted_in_it(self):
        self.assertNotIn("Lowest priority of the concerns here", self.html)


class PackagingTests(unittest.TestCase):
    """A fresh installer hits these; the maintainer never does.

    run_local.sh takes a different path from run.py, so neither the missing
    checkpoint gate nor the port disagreement could fire for anyone already
    running the server.
    """

    def setUp(self):
        self.run_py = (ROOT / "run.py").read_text(encoding="utf-8")

    def test_no_engine_checkpoint_gates_startup(self):
        # AntiBERTy carries its own checkpoint inside the installed pip
        # package; ImmuneBuilder and DeepSP are fetched after install.
        # Listing any of them turns a missing optional engine into
        # "the package is incomplete" and refuses to start at all.
        block = re.search(r"REQUIRED_FILES = \((.*?)\n\)", self.run_py, re.S)
        self.assertIsNotNone(block)
        self.assertNotIn("models/", block.group(1))

    def test_every_required_file_is_actually_shipped(self):
        for rel in re.findall(r'ROOT / "([^"]+)"',
                              re.search(r"REQUIRED_FILES = \((.*?)\n\)",
                                        self.run_py, re.S).group(1)):
            self.assertTrue((ROOT / rel).is_file(), f"REQUIRED_FILES lists a missing file: {rel}")

    def test_no_deliberately_unshipped_dataset_gates_startup(self):
        """The licence-restricted datasets are absent from every package.

        The project does not redistribute the Thera-SAbDab-derived corpora,
        so they are missing by design. ReferenceDataMissing explains that at
        the point of use and names the build command; demanding them at
        startup pre-empts a message written for exactly this case, and made
        the published package refuse to boot.
        """
        from app.reference_data import OPTIONAL_DATASETS
        optional = {filename for filename, _d, _c in OPTIONAL_DATASETS.values()}
        listed = set(re.findall(r'ROOT / "data/([^"]+)"',
                                re.search(r"REQUIRED_FILES = \((.*?)\n\)",
                                          self.run_py, re.S).group(1)))
        self.assertEqual(listed & optional, set())

    def test_the_launchers_agree_on_the_port(self):
        shell = (ROOT / "run_local.sh").read_text(encoding="utf-8")
        self.assertIn('PORT="${PORT:-8000}"', shell)
        self.assertIn('os.environ.get("PORT", "8000")', self.run_py)

    def test_no_document_still_advertises_the_old_port(self):
        for name in ("README.md", "START_HERE.md", "QUICKSTART.md"):
            path = ROOT / name
            if path.is_file():
                self.assertNotIn("8765", path.read_text(encoding="utf-8"), name)

    def test_the_antiberty_checkpoint_is_resolved_not_assumed(self):
        provenance = (ROOT / "app" / "source_provenance.py").read_text(encoding="utf-8")
        self.assertNotIn("models/antiberty", provenance)
        self.assertIn("get_weights(download=False)", provenance)

    def test_no_module_resolves_the_checkpoint_by_a_path_alone(self):
        """models/antiberty was wrong in six places, each failing silently.

        The startup gate refused to boot. The provenance record always read
        "Not installed". The test gate skipped tests that had therefore never
        run anywhere. The batch CLI reported the engine unavailable on every
        machine. Two build scripts pinned it as an explicit checkpoint_path,
        harmless only because get_weights() ignores a path that does not
        exist. A release script packed the tree.

        The rule is not "never name the path" -- app/main.py names it as a
        first guess and then falls through, which is correct and documented.
        The rule is that naming it is never the ONLY way a module finds the
        checkpoint: anything that mentions it must also ask the package.
        """
        offenders = []
        for path in sorted(ROOT.rglob("*.py")):
            if {".venv", "__pycache__", "tests"} & set(path.parts):
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            mentions = [line for line in text.splitlines()
                        if ("models/antiberty" in line or 'models" / "antiberty' in line)
                        and not line.lstrip().startswith("#")]
            if mentions and "get_weights" not in text:
                offenders.append(f"{path.relative_to(ROOT)}: {mentions[0].strip()[:70]}")
        self.assertEqual(offenders, [])

    def test_the_batch_path_resolves_the_checkpoint(self):
        panel = (ROOT / "run_panel.py").read_text(encoding="utf-8")
        self.assertIn("get_weights(download=False)", panel)

    def test_the_test_gate_resolves_it_too(self):
        # A gate that is false everywhere does not protect a test, it
        # deletes it.
        support = (ROOT / "tests" / "support.py").read_text(encoding="utf-8")
        self.assertIn("get_weights(download=False)", support)

    def test_the_ui_test_fixture_is_not_in_the_tree(self):
        self.assertFalse((ROOT / "static" / "_uitest_9pi9.cif").exists())


class ViewerInteractionTests(unittest.TestCase):
    """A hydrophobic map drawn on a ribbon cannot show a patch."""

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_surface_properties_open_on_a_surface(self):
        self.assertIn("SURFACE_MODES = new Set(['hydrophobic', 'charge'])", self.html)
        self.assertIn("let style = styleFor(mode);", self.html)

    def test_the_markup_agrees_with_the_default(self):
        self.assertIn('data-style="surface" aria-pressed="true"', self.html)
        self.assertIn('data-style="cartoon" aria-pressed="false"', self.html)

    def test_a_chosen_style_is_not_overridden(self):
        self.assertIn("styleChosen = true;", self.html)
        self.assertIn("if(!styleChosen){", self.html)

    def test_clicking_a_residue_reads_it_out(self):
        self.assertIn("stage.signals.clicked.add", self.html)
        self.assertIn("stage.signals.clicked.removeAll()", self.html)

    def test_the_readout_cannot_disagree_with_the_colours(self):
        # Both come from the same index and the same score functions.
        self.assertIn("hydroScore(e.aa, e.rsa).toFixed(2)", self.html)
        self.assertIn("chargeScore(e.aa, e.rsa).toFixed(2)", self.html)

    def test_a_table_row_focuses_its_residue(self):
        self.assertIn("function wireResidueFocus()", self.html)
        self.assertIn("wireResidueFocus();", self.html)
        self.assertIn("if(!tr) return;", self.html)
        self.assertIn("AtlasWorkspace.selectPosition(tr.dataset.chain", self.html)


class FixWorkbenchTests(unittest.TestCase):
    """Ten sites used to mean ten stacked cards and ten tables.

    Comparing two sites meant scrolling between them with nothing held
    still. One list and one detail pane: the list is the comparison, the
    pane is the depth.
    """

    def setUp(self):
        self.html = START.read_text(encoding="utf-8")

    def test_the_stacked_cards_are_gone(self):
        self.assertIn("function renderFixWorkbench(", self.html)
        self.assertNotIn('out.push(`<div class="liab"><div class="liab-head">', self.html)

    def test_every_pane_is_rendered_up_front(self):
        # The "test" buttons are bound once by querySelectorAll after render,
        # and printing has to carry every site. Building panes on demand
        # would break both, silently.
        self.assertIn("wb-pane", self.html)
        self.assertNotIn("innerHTML = renderFixPane", self.html)

    def test_printing_reveals_every_site(self):
        self.assertIn(".wb-pane[hidden]{display:block !important}", self.html)

    def test_the_list_is_keyboard_navigable(self):
        # A roving tabindex is meaningless without the keys that move it.
        self.assertIn("ArrowDown: current + 1, ArrowUp: current - 1", self.html)
        self.assertIn('role="tablist"', self.html)
        self.assertIn('role="tabpanel"', self.html)

    def test_it_is_wired_after_render(self):
        self.assertIn("wireFixWorkbench();", self.html)

    def test_one_signal_builder_feeds_list_and_pane(self):
        # Two derivations of "is this on a contact" is how a list row and its
        # own detail pane end up disagreeing.
        self.assertIn("function fixSiteSignals(site)", self.html)
        self.assertEqual(self.html.count("fixSiteSignals(site)"), 3)
