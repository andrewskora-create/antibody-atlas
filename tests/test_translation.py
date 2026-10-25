"""Nucleotide input, and the trap it sets."""
import random
import unittest

from app.analysis import clean_sequence
from app.translation import (CODONS, describe, looks_like_nucleotide, resolve,
                             reverse_complement, to_protein, translate)

VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
      "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
      "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
SIGNAL = "MGWSCIILFLVATATGVHS"

_BACK = {}
for _codon, _aa in CODONS.items():
    _BACK.setdefault(_aa, []).append(_codon)


def back_translate(protein, seed=11):
    rng = random.Random(seed)
    return "".join(rng.choice(_BACK[a]) for a in protein)


class DetectionTests(unittest.TestCase):
    def test_a_protein_is_not_mistaken_for_dna(self):
        self.assertFalse(looks_like_nucleotide(VH))
        self.assertFalse(looks_like_nucleotide(VL))

    def test_an_all_acgt_peptide_is_too_short_to_be_called_dna(self):
        """A, C, G and T are amino acids. Short all-ACGT strings stay protein."""
        self.assertFalse(looks_like_nucleotide("ACGTACGTACGTACGT"))

    def test_a_real_gene_is_detected(self):
        self.assertTrue(looks_like_nucleotide(back_translate(VH)))


class TranslationTests(unittest.TestCase):
    def test_bare_variable_domain_round_trips(self):
        result = to_protein(back_translate(VH))
        self.assertEqual(result["protein"], VH)
        self.assertEqual(result["frame"], "+1")
        self.assertEqual(result["chain_type"], "H")

    def test_signal_peptide_and_vector_flanks_are_removed_not_just_tolerated(self):
        dna = "GCTAGCGCCACC" + back_translate(SIGNAL + VH) + "TGA" + "GGATCCGGTACC"
        result = to_protein(dna)
        self.assertEqual(result["protein"], VH)
        self.assertIn(SIGNAL, result["leader"])

    def test_the_leader_is_dropped_before_any_descriptor_sees_it(self):
        """A secretion leader is strongly hydrophobic. Leaving it in would
        inflate GRAVY and invent hydrophobic windows that are not in the
        folded protein, which is exactly the kind of silent wrong answer this
        whole module exists to avoid."""
        from app.analysis import sequence_metrics
        clean = sequence_metrics(VH)
        with_leader = sequence_metrics(SIGNAL + VH)
        self.assertGreater(with_leader["gravy"], clean["gravy"])
        self.assertGreater(with_leader["aggregation_proxy"]["patch_count"],
                           clean["aggregation_proxy"]["patch_count"])
        translated = sequence_metrics(to_protein(back_translate(SIGNAL + VH))["protein"])
        self.assertEqual(translated["gravy"], clean["gravy"])
        self.assertEqual(translated["aggregation_proxy"]["patch_count"],
                         clean["aggregation_proxy"]["patch_count"])

    def test_reverse_complement_is_found(self):
        dna = reverse_complement(back_translate(SIGNAL + VH))
        result = to_protein(dna)
        self.assertEqual(result["protein"], VH)
        self.assertTrue(result["reverse_complement"])
        self.assertTrue(result["frame"].startswith("-"))

    def test_the_frame_is_chosen_by_numbering_not_by_stop_count(self):
        """Shifting the start by one base must not shift the answer."""
        dna = "G" + back_translate(VH)
        self.assertEqual(to_protein(dna)["protein"], VH)
        self.assertEqual(to_protein(dna)["frame"], "+2")

    def test_light_chain_is_recognised_as_light(self):
        self.assertEqual(to_protein(back_translate(VL))["chain_type"], "K")

    def test_whitespace_numbering_and_case_are_tolerated(self):
        dna = back_translate(VH)
        formatted = "  1 " + dna[:60].lower() + "\n 61 " + dna[60:].lower()
        self.assertEqual(to_protein(formatted)["protein"], VH)

    def test_non_antibody_dna_fails_with_the_frames_it_tried(self):
        with self.assertRaises(ValueError) as caught:
            to_protein("ATGCATGCAT" * 12)
        message = str(caught.exception)
        self.assertIn("reading frame", message)
        self.assertIn("+1", message)

    def test_a_sequence_that_is_not_nucleotide_is_rejected_clearly(self):
        with self.assertRaises(ValueError) as caught:
            to_protein("ACGT" * 20 + "ZZZ")
        self.assertIn("unexpected symbols", str(caught.exception))

    def test_stops_and_unknown_codons_are_represented(self):
        self.assertEqual(translate("TAATAGTGA"), "***")
        self.assertEqual(translate("NNN"), "X")


class DiagnosticTests(unittest.TestCase):
    """The old failure said only that nothing was recognised, which named no
    cause and sent the user looking in the wrong place."""

    def test_a_short_fragment_is_called_truncated(self):
        self.assertIn("truncated", describe(VH[:40]))

    def test_short_dna_says_so_rather_than_being_called_protein(self):
        message = describe("ATGGGATGGAGCTGTATC")
        self.assertIn("DNA", message)
        self.assertIn("60", message)

    def test_an_overlong_input_suggests_splitting_the_chains(self):
        self.assertIn("separately", describe(VH + VL + VH + VL))

    def test_a_few_ambiguity_codes_do_not_stop_dna_being_recognised(self):
        dna = list(back_translate(VH))
        for i in (10, 50, 120):
            dna[i] = "R"
        self.assertTrue(looks_like_nucleotide("".join(dna)))

    def test_asparagines_do_not_make_a_protein_look_like_dna(self):
        """N is both an ambiguity code and an amino acid; only unambiguous
        bases count towards the nucleotide score."""
        self.assertFalse(looks_like_nucleotide("ACGTN" * 30))
        self.assertTrue(looks_like_nucleotide("ACGT" * 40))


class EntryPointTests(unittest.TestCase):
    def test_resolve_passes_protein_through_untouched(self):
        sequence, record = resolve(VH)
        self.assertEqual(sequence, VH)
        self.assertIsNone(record)

    def test_clean_sequence_translates_dna_by_default(self):
        self.assertEqual(clean_sequence(back_translate(VH)), VH)

    def test_clean_sequence_can_be_told_not_to(self):
        """The API translates up front so it can report it, then cleans with
        translation off; that path must not translate twice."""
        dna = back_translate(VH)
        self.assertEqual(clean_sequence(dna, allow_nucleotide=False), dna)


if __name__ == "__main__":
    unittest.main()
