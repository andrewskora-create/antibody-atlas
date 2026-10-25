"""End-to-end check that the numpy DeepSP re-implementation matches upstream.

DeepSP ships Keras HDF5 weights and a TensorFlow predictor. This project
re-implements inference in numpy so it does not take a TensorFlow dependency
for three small convolutional networks, and uses ANARCII where DeepSP uses
legacy ANARCI. Both substitutions are places where the answer could silently
drift, so both are checked against the sixteen reference molecules DeepSP
publishes with the model:

  * alignment -- ANARCII must place every residue in the same column ANARCI did
  * inference -- the numpy forward pass must reproduce the published descriptors

The fixture carries upstream's own inputs, alignments and outputs, so the
comparison is against DeepSP, not against a value this project recorded once.
"""
from __future__ import annotations

import json
import unittest
from pathlib import Path

from app import deepsp, numbering
from tests.support import requires_deepsp

FIXTURE = Path(__file__).parent / "fixtures" / "deepsp_reference_molecules.json"

#: float32 accumulation order differs between Keras and numpy; anything above
#: this is a real disagreement rather than rounding.
TOLERANCE = 1e-3


def _fixture() -> dict:
    with FIXTURE.open() as handle:
        return json.load(handle)


class AlignmentGrid(unittest.TestCase):
    def test_grid_is_the_shape_deepsp_expects(self):
        self.assertEqual(deepsp.ALIGNED_LENGTH, 272)
        self.assertEqual(len(deepsp._H_GRID), 145)
        self.assertEqual(len(deepsp._L_GRID), 127)

    def test_cdrh3_insertion_ladder_order(self):
        # DeepSP numbers CDR-H3 insertions 111A..111H then 112I..112A, i.e. the
        # second half runs backwards. Getting this wrong mirrors CDR-H3 without
        # changing the length, so nothing downstream would notice.
        ladder = deepsp._H_GRID[110:129]
        self.assertEqual(
            list(ladder),
            ["111", "111A", "111B", "111C", "111D", "111E", "111F", "111G", "111H",
             "112I", "112H", "112G", "112F", "112E", "112D", "112C", "112B", "112A", "112"])

    def test_light_grid_has_no_insertion_columns(self):
        self.assertTrue(all(label.isdigit() for label in deepsp._L_GRID))


class AlignmentAgainstUpstream(unittest.TestCase):
    """ANARCII must reproduce the alignment DeepSP's ANARCI produced."""

    def test_every_reference_molecule_aligns_identically(self):
        differing = []
        for molecule in _fixture()["molecules"]:
            heavy = numbering.number_chain(molecule["heavy"], "imgt")
            light = numbering.number_chain(molecule["light"], "imgt")
            result = deepsp.align(heavy, light)
            if result["aligned"] != molecule["aligned"]:
                columns = [i for i, (a, b) in
                           enumerate(zip(result["aligned"], molecule["aligned"])) if a != b]
                differing.append((molecule["name"], columns[:10]))
            self.assertEqual(result["dropped"], [],
                             f"{molecule['name']}: residues fell off the grid")
        self.assertEqual(differing, [], "ANARCII and ANARCI disagree on column placement")


class InferenceAgainstUpstream(unittest.TestCase):
    @requires_deepsp()
    def test_published_descriptors_are_reproduced(self):
        fixture = _fixture()
        molecules = fixture["molecules"]
        predicted = deepsp.predict_aligned([m["aligned"] for m in molecules])
        worst, worst_key = 0.0, None
        compared = 0
        for molecule, values in zip(molecules, predicted):
            for key, expected in molecule["descriptors"].items():
                self.assertIn(key, values)
                relative = abs(values[key] - expected) / max(abs(expected), 1e-6)
                compared += 1
                if relative > worst:
                    worst, worst_key = relative, f"{molecule['name']}/{key}"
        self.assertEqual(compared, len(molecules) * 30)
        self.assertLess(worst, TOLERANCE,
                        f"worst deviation {worst:.2e} at {worst_key}")

    @requires_deepsp()
    def test_full_pass_from_raw_sequence(self):
        """Sequence in, descriptors out, with no upstream alignment handed to us."""
        molecule = _fixture()["molecules"][0]
        heavy = numbering.number_chain(molecule["heavy"], "imgt")
        light = numbering.number_chain(molecule["light"], "imgt")
        result = deepsp.analyze(heavy, light)
        self.assertEqual(result["status"], "ok")
        for key, expected in molecule["descriptors"].items():
            self.assertAlmostEqual(result["descriptors"][key], expected,
                                   delta=max(abs(expected), 1.0) * TOLERANCE)

    @requires_deepsp()
    def test_highlights_carry_percentiles_when_reference_present(self):
        molecule = _fixture()["molecules"][0]
        heavy = numbering.number_chain(molecule["heavy"], "imgt")
        light = numbering.number_chain(molecule["light"], "imgt")
        result = deepsp.analyze(heavy, light)
        self.assertTrue(result["highlights"])
        for row in result["highlights"]:
            self.assertIn("risk", row)
            if result.get("reference"):
                self.assertIsNotNone(row["percentile"])
                self.assertGreaterEqual(row["percentile"], 0.0)
                self.assertLessEqual(row["percentile"], 100.0)


class Guards(unittest.TestCase):
    def test_heavy_only_is_skipped_not_guessed(self):
        result = deepsp.analyze({"scheme": "imgt", "residues": []}, None)
        self.assertEqual(result["status"], "skipped")

    def test_kabat_numbering_is_refused(self):
        with self.assertRaises(ValueError):
            deepsp.align({"scheme": "kabat", "residues": []},
                         {"scheme": "kabat", "residues": []})

    def test_overflowing_cdrh3_is_reported_as_unreliable(self):
        """A CDR-H3 longer than the grid must not quietly return numbers."""
        molecule = _fixture()["molecules"][0]
        heavy = numbering.number_chain(molecule["heavy"], "imgt")
        light = numbering.number_chain(molecule["light"], "imgt")
        # Synthesise insertions past 111H/112I, which the grid cannot hold.
        heavy = dict(heavy, residues=list(heavy["residues"]) + [
            {"number": "H111J", "number_base": 111, "aa": "W",
             "region": "CDR3", "sequence_position": 999}])
        alignment = deepsp.align(heavy, light)
        self.assertTrue(alignment["dropped_from_cdr"])
        if deepsp.available():
            result = deepsp.analyze(heavy, light)
            self.assertEqual(result["status"], "unreliable")
            self.assertIn("truncated", result["warning"])


if __name__ == "__main__":
    unittest.main()
