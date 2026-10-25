"""The modeling shortlist combines observations without inventing missing evidence."""
import unittest
from unittest.mock import patch

from app.cdr_substitutions import _rank, _selection_score, propose_at_position


def candidate(replacement, *, similarity=.8, clinical=0, position=10, clean=True,
              germline=0):
    return {
        'original_aa': 'A', 'replacement_aa': replacement, 'edit_position': position,
        'blosum62': None, 'physicochemical': {'distance': 1 - similarity},
        'clinical_precedent': {'status': 'unavailable' if clinical is None else 'available',
                               'percent': clinical},
        'germline_precedent': {'percent': germline}, 'clean': clean,
        'variant_spec': f'H:A{position}{replacement}',
    }


class SelectionRankingTests(unittest.TestCase):
    def test_clinical_observation_can_change_similarity_only_order(self):
        high_similarity = candidate('V', similarity=.95, clinical=0)
        high_frequency = candidate('S', similarity=.65, clinical=80)
        ranked = _rank([high_similarity, high_frequency], 'user_position')
        self.assertEqual(ranked[0]['replacement_aa'], 'S')
        self.assertEqual(ranked[0]['selection_score']['score'], .725)
        self.assertEqual(ranked[0]['selection_score']['weights'],
                         {'residue_similarity': .5, 'clinical_frequency': .5})

    def test_observed_zero_and_unavailable_are_distinct(self):
        ranked = _rank([candidate('V', clinical=0), candidate('S', clinical=None)], 'user_position')
        by_aa = {c['replacement_aa']: c['selection_score'] for c in ranked}
        self.assertEqual(by_aa['V']['score'], .4)
        self.assertEqual(by_aa['V']['clinical_frequency'], 0)
        self.assertEqual(by_aa['S']['score'], .8)
        self.assertIsNone(by_aa['S']['clinical_frequency'])
        self.assertEqual(by_aa['S']['basis'], 'similarity_only')

    def test_safety_precedence_is_never_displaced_by_high_blend(self):
        protected = candidate('V', similarity=1, clinical=100, position=12)
        new_liability = candidate('C', similarity=1, clinical=100, clean=False)
        clean = candidate('S', similarity=.1, clinical=0)
        ranked = _rank([protected, new_liability, clean], 'user_position', {12})
        self.assertEqual([c['replacement_aa'] for c in ranked], ['S', 'C', 'V'])

    def test_germline_only_breaks_ties_in_blend(self):
        a = candidate('S', clinical=10, germline=95)
        b = candidate('V', clinical=10, germline=5)
        self.assertEqual(_rank([b, a], 'user_position')[0]['replacement_aa'], 'S')

    def test_invalid_clinical_values_use_explicit_fallback(self):
        for value in (float('nan'), float('inf'), -1, 101, True, '80'):
            with self.subTest(value=value):
                c = candidate('S', clinical=value)
                c['conservativeness'] = {'score': .8}
                result = _selection_score(c)
                self.assertEqual(result['score'], .8)
                self.assertEqual(result['basis'], 'similarity_only')
                self.assertIsNone(result['clinical_frequency'])

    def test_all_nineteen_remain_available_and_top_five_follow_ranking(self):
        numbering = {'scheme': 'imgt', 'chain_type': 'H', 'residues': [
            {'sequence_position': 1, 'number': 'H27', 'insertion': '', 'region': 'CDR1', 'aa': 'A'}]}
        chains = {'H': {'sequence': 'A', 'numbering': numbering, 'numbering_imgt': numbering}}
        clinical = {'status': 'available', 'percent': 20, 'observed': 2, 'total': 10}
        with patch('app.cdr_substitutions.ptm_sites', return_value=[]), \
                patch('app.cdr_substitutions._germline_frequency', return_value={'percent': 10}), \
                patch('app.cdr_substitutions._clinical_frequency', return_value=clinical), \
                patch('app.cdr_substitutions.mutation_effect', return_value={}), \
                patch('app.cdr_substitutions._blosum62', return_value=None):
            result = propose_at_position(chains, 'H', 1, max_candidates=19)
        self.assertEqual(len(result['candidates']), 19)
        self.assertEqual(result['top_model_candidates'],
                         [c['variant_spec'] for c in result['candidates'][:5]])
        self.assertTrue(all('selection_score' in c for c in result['candidates']))
        self.assertIn('equal-weight blend', result['basis'])
        self.assertIn('not a prediction of folding stability', result['basis'])


if __name__ == '__main__':
    unittest.main()
