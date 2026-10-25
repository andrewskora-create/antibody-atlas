import json
import unittest

from tests.support import requires_dataset
from unittest.mock import patch

from app.ptm_precedents import PRECEDENTS, REFERENCE_PATH, matching_precedents, precedent_catalog, sequence_precedents


@requires_dataset("therapeutic_ptm_sequences")
class PrecedentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = json.loads(REFERENCE_PATH.read_text())['antibodies']

    def test_verified_reference_sites_match_their_input_chain(self):
        for reference in self.references.values():
            chains = {label: {'sequence': sequence} for label, sequence in reference['chains'].items()}
            result = sequence_precedents(chains)
            for site in reference['validated_sites']:
                with self.subTest(site=site['precedent_id']):
                    case = next(c for c in result['cases'] if c['id'] == site['precedent_id'])
                    match = case['matches'][0]
                    self.assertEqual(match['identity_percent'], 100)
                    self.assertEqual(match['query_sequence_position'], site['sequence_position'])
                    self.assertEqual(match['matched_motif'], site['expected_motif'])
                    self.assertEqual(match['match_type'], 'Exact reference chain sequence')
                    self.assertEqual(match['identity_matches'], match['identity_denominator'])
                    self.assertEqual(len(match['query_context']), match['local_positions'])
            self.assertEqual(result['searched_cases'], len(PRECEDENTS))

    def test_same_chemistry_at_different_framework_position_is_not_match(self):
        sequence = self.references['infliximab']['chains']['H']
        changed = list(sequence)
        changed[17] = 'L'  # Remove documented Met18; other Met flags remain.
        changed[5] = 'M'
        cases = sequence_precedents({'H': {'sequence': ''.join(changed)}})['cases']
        self.assertNotIn('infliximab-hc-m18-oxidation', {c['id'] for c in cases})

    def test_heavy_sequence_in_light_slot_cannot_match_heavy_precedent(self):
        sequence = self.references['infliximab']['chains']['H']
        self.assertFalse(sequence_precedents({'L': {'sequence': sequence}})['cases'])

    def test_same_anchor_with_different_chemistry_context_does_not_match(self):
        sequence = self.references['bevacizumab']['chains']['H']
        changed = sequence[:84] + 'T' + sequence[85:]  # N84 remains, NS becomes NT.
        cases = sequence_precedents({'H': {'sequence': changed}})['cases']
        self.assertNotIn('bevacizumab-hc-n84-deamidation', {c['id'] for c in cases})
        glyco = self.references['cetuximab']['chains']['H']
        changed = glyco[:88] + 'A' + glyco[89:]  # NDT becomes NAT: still a sequon, different context.
        cases = sequence_precedents({'H': {'sequence': changed}})['cases']
        self.assertNotIn('cetuximab-hc-n88-glycosylation', {c['id'] for c in cases})

    def test_imgt_and_kabat_preserve_site_match_but_label_selected_convention(self):
        sequence = self.references['cetuximab']['chains']['H']
        outputs = [sequence_precedents({'H': {'sequence': sequence}}, scheme) for scheme in ('imgt','kabat')]
        matches = [next(c for c in output['cases'] if c['id']=='cetuximab-hc-n88-glycosylation')['matches'][0] for output in outputs]
        self.assertEqual(matches[0]['reference_imgt_number'], matches[1]['reference_imgt_number'])
        self.assertEqual(matches[0]['query_sequence_position'], matches[1]['query_sequence_position'])
        self.assertEqual(matches[1]['query_number'], 'H85')
        self.assertEqual(matches[1]['query_region'], 'Framework')
        self.assertEqual(matches[1]['query_scheme'], 'kabat')

    def test_lambda_cannot_match_kappa_reference(self):
        from app.numbering import number_chain
        sequence = self.references['trastuzumab']['chains']['L']
        numbered = number_chain(sequence, 'imgt')
        self.assertEqual(numbered['chain_type'], 'K')
        numbered['chain_type'] = 'L'
        for row in numbered['residues']:
            row['number'] = 'L' + row['number'][1:]
        with patch('app.ptm_precedents.number_chain', return_value=numbered):
            cases = sequence_precedents({'L': {'sequence': sequence}})['cases']
        self.assertNotIn('trastuzumab-lc-n30-deamidation', {c['id'] for c in cases})

    def test_unnumberable_optional_chain_does_not_hide_valid_chain_matches(self):
        from app.numbering import number_chain
        sequence = self.references['bevacizumab']['chains']['H']
        def number_or_reject(value, scheme):
            if value == 'AAAA':
                raise ValueError('No variable domain')
            return number_chain(value, scheme)
        with patch('app.ptm_precedents.number_chain', side_effect=number_or_reject):
            result = sequence_precedents({'H': {'sequence': sequence}, 'L': {'sequence':'AAAA'}})
        self.assertIn('bevacizumab-hc-n84-deamidation', {c['id'] for c in result['cases']})
        self.assertEqual(result['skipped_chains'][0]['chain'], 'L')

    def test_pyroglutamate_requires_actual_input_n_terminus(self):
        sequence = self.references['rituximab']['chains']['H']
        exact = sequence_precedents({'H': {'sequence': sequence}})['cases']
        self.assertIn('rituximab-hc-q1-pyroglutamate', {c['id'] for c in exact})
        extended = sequence_precedents({'H': {'sequence': 'AA' + sequence}})['cases']
        self.assertNotIn('rituximab-hc-q1-pyroglutamate', {c['id'] for c in extended})

    def test_broad_region_matching_is_disabled(self):
        self.assertEqual(matching_precedents('Asparagine deamidation motif','Framework','N'), [])
        for case in precedent_catalog()['cases']:
            self.assertTrue(case['numbering'])
            self.assertTrue(case['conditions'])
            self.assertTrue(case['url'].startswith('https://'))

    def _synthetic(self, query_positions):
        def numbered(positions):
            return {'chain_type':'H','scheme':'imgt','residues':[
                {'sequence_position':i+1,'number':f'H{p}','aa':aa,'region':'Framework'}
                for i,(p,aa) in enumerate(zip(positions,'NSA'))]}
        reference = {'antibodies':{'synthetic':{'chains':{'H':'NSA'},'validated_sites':[
            {'precedent_id':'bevacizumab-hc-n84-deamidation','validated':True,'chain':'H',
             'sequence_position':1,'observed_motif':'NS','expected_motif':'NS'}]}}}
        with patch('app.ptm_precedents._reference_data',return_value=reference), \
             patch('app.ptm_precedents._number_reference',return_value=numbered([1,2,3])), \
             patch('app.ptm_precedents.number_chain',return_value=numbered(query_positions)):
            return sequence_precedents({'H':{'sequence':'NSA'}})

    def test_all_motif_positions_must_align(self):
        # Same NS letters and N anchor, but S belongs to a different canonical position.
        self.assertEqual(self._synthetic([1,3,4])['cases'], [])

    def test_identity_denominator_includes_union_gaps(self):
        match = self._synthetic([1,2,4])['cases'][0]['matches'][0]
        self.assertEqual((match['identity_matches'],match['identity_denominator']), (2,4))
        self.assertEqual(match['identity_percent'],50)
        self.assertEqual((match['query_context'],match['reference_context']), ('NS-A','NSA-'))
        self.assertEqual(match['local_matches'],2)
        self.assertEqual(match['match_type'],'Aligned site and motif; differing neighborhood')


if __name__ == '__main__':
    unittest.main()
