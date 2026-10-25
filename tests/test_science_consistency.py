"""Cross-view regressions using public fixture sequences and synthetic evidence."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from app.constant_regions import assemble
from app.molecule_flags import molecule_flags, _molecular_pi
from app.residue_table import _paratope_index, _unmapped_paratope, _liability_rows, _verdict
from app.structure_analysis import _group_physical_sites, _site_priority
from app.structure_ensemble import reconcile_sites
from tests.test_constant_regions import chains


class MolecularChargeTests(unittest.TestCase):
    def check_charge(self, chain_data, assembly=None):
        result = molecule_flags(chain_data, assembly=assembly)
        components = ([(assembly[k]['sequence'], assembly[k]['copies'])
                       for k in ('heavy', 'light') if assembly.get(k)] if assembly else
                      [(c['sequence'], 1) for c in chain_data.values()])
        expected = lambda ph: sum(n * ProteinAnalysis(s).charge_at_pH(ph) for s, n in components)
        self.assertEqual(result['fv']['net_charge_assay_ph'], round(expected(7.4), 2))
        self.assertEqual(result['fv']['net_charge_formulation_ph'], round(expected(6.0), 2))
        for point in result['titration']:
            self.assertEqual(point['net_charge'], round(expected(point['ph']), 2))
        pi = result['fv']['pi']
        self.assertGreaterEqual(expected(pi - 0.006), 0)
        self.assertLessEqual(expected(pi + 0.006), 0)
        return result

    def test_fv_h2l2_and_h2_preserve_all_physical_chain_termini(self):
        c = chains()
        self.check_charge(c)
        for clipped in (True, False):
            a = assemble(c, clip_c_terminal_lysine=clipped)
            self.check_charge(c, a)
        h = {'H': c['H']}
        self.check_charge(h, assemble(h))
        self.check_charge(h)
        self.check_charge({'L': c['L']})

    def test_pi_unchanged_when_all_chain_copy_counts_are_doubled(self):
        c = chains()
        components = [(item['sequence'], 1) for item in c.values()]
        self.assertEqual(_molecular_pi(components), _molecular_pi([(s, 2) for s, _ in components]))

    def test_asymmetry_value_and_shown_charges_use_same_basis(self):
        c = chains()
        a = assemble(c)
        result = self.check_charge(c, a)
        flag = next(f for f in result['flags'] if f['id'] == 'charge_asymmetry')
        expected = {label: round(ProteinAnalysis(a[key]['sequence']).charge_at_pH(7.4), 2)
                    for label, key in (('H', 'heavy'), ('L', 'light'))}
        self.assertEqual(flag['per_chain_charge'], expected)
        self.assertEqual(flag['value'], round(max(expected.values()) - min(expected.values()), 2))


class ResidueEvidenceTests(unittest.TestCase):
    def test_contact_without_liability_keeps_exposure_and_ensemble_ptm_overrides_primary(self):
        from app.residue_table import build_residue_table
        from tests.test_residue_table import chain
        structure = {'residues': [
            {'chain': 'H', 'sequence_position': 2, 'relative_sidechain_sasa': .6},
            {'chain': 'H', 'sequence_position': 6, 'relative_sidechain_sasa': .05}],
            'ptm_sites': [{'chain': 'H', 'sequence_position': 6, 'priority': 'high',
                           'relative_sidechain_sasa': .8, 'kind': 'Asparagine deamidation motif',
                           'motif': 'NG', 'exposure_models': 2}]}
        paratope = {'status': 'available', 'by_chain': {'H': {2: {'probability': .9, 'likely_contact': True}}}}
        rows = build_residue_table({'H': chain()}, structure, paratope)['rows']
        self.assertEqual(rows[1]['verdict'], 'preserve')
        self.assertEqual(rows[1]['exposure'], .6)
        self.assertEqual(rows[5]['exposure'], .8)

    def test_json_export_roundtrip_retains_contact_positions(self):
        p = {'status': 'available', 'by_chain': {'H': {
            27: {'probability': .9, 'likely_contact': True},
            'IMGT28': {'probability': .8, 'likely_contact': True}}}}
        loaded = json.loads(json.dumps(p))
        self.assertEqual(_paratope_index(p, 'H'), _paratope_index(loaded, 'H'))
        self.assertEqual(_unmapped_paratope(loaded), 1)

    def test_grouped_flags_keep_individual_priority_for_softening(self):
        # NTS has a glycosylation sequon plus a weaker deamidation motif.
        residue = {'chain': 'H', 'aa': 'N', 'number': 'H10', 'region': 'Framework',
                   'relative_sidechain_sasa': .05, 'sidechain_sasa_a2': 5}
        glyco = {'kind': 'N-linked glycosylation sequon', 'motif': 'NTS', 'start': 10}
        deamidation = {'kind': 'Asparagine deamidation motif', 'motif': 'NT', 'start': 10}
        grouped = _group_physical_sites([_site_priority(s, residue, 'imgt') for s in (glyco, deamidation)])
        precedent = {('H', 10, glyco['kind']): {'approved': 1}}
        liabilities = _liability_rows({}, grouped, 'H', precedent)[10]
        self.assertEqual([l['priority'] for l in liabilities], ['review', 'lower'])
        self.assertEqual(_verdict(liabilities, None, .05, 'Framework')[0], 'watch')

    def test_ensemble_max_exposure_recomputes_chemistry_and_keeps_audit(self):
        residue = {'chain': 'H', 'aa': 'N', 'number': 'H31', 'region': 'CDR1',
                   'relative_sidechain_sasa': .05, 'sidechain_sasa_a2': 5}
        sites = _group_physical_sites([_site_priority(
            {'kind': 'Asparagine deamidation motif', 'motif': 'NG', 'start': 31}, residue, 'imgt')])
        row = {'chain': 'H', 'sequence_position': 31, 'exposure_max': .8, 'exposure_min': .05,
               'band': 'high', 'models': 2, 'spread': .75, 'band_stable': False,
               'bands_seen': ['low', 'high']}
        result = reconcile_sites(sites, [row])[0]
        self.assertEqual(result['single_model_priority'], 'review')
        self.assertEqual(result['priority'], 'high')
        self.assertEqual(result['flags'][0]['priority'], 'high')
        self.assertEqual(result['exposure_percent'], 80)
        self.assertIn('80.0%', result['reasons'][0])


class CrossViewTests(unittest.TestCase):
    def test_batch_csv_keeps_evaluate_count_in_legacy_header(self):
        from run_batch import _row
        row = _row('example', {'residue_table': {'verdict_counts': {'evaluate': 3, 'watch': 2}}})
        self.assertEqual(row['residues_fix'], 3)

    def test_batch_worker_retains_computed_precedents(self):
        from run_batch import _screen_one
        numbering = {'scheme': 'imgt', 'chain_type': 'H', 'residues': [
            {'sequence_position': i, 'number': f'H{i}', 'number_base': i,
             'region': 'Framework', 'aa': aa} for i, aa in enumerate('NTS', 1)]}
        corpus = {'sites': [{'query_chain': 'H', 'position': 1,
                            'kind': 'N-linked glycosylation sequon',
                            'matches': [{'drug': 'Public fixture', 'stage': 'Approved'}]}]}
        with tempfile.TemporaryDirectory() as directory, \
             patch('app.numbering.number_chain', return_value=numbering), \
             patch('app.humanness.framework_identity', return_value={}), \
             patch('app.sequence_panel.annotate_chain', return_value={}), \
             patch('app.sequence_panel.conclusion', return_value={}), \
             patch('app.clinical_frameworks.clinical_framework_precedents', return_value=corpus), \
             patch('app.ptm_precedents.sequence_precedents', return_value={'cases': []}):
            result = _screen_one(('example', {'H': 'NTS'}, {
                'out': directory, 'threads': 1, 'scheme': 'imgt', 'ph': 7.4,
                'fold': False, 'substitutions': False, 'pdf': False}))
            self.assertEqual(result['status'], 'complete', result)
            record = json.loads((Path(directory) / 'example.json').read_text())
            self.assertEqual(record['residue_table']['rows'][0]['precedent']['approved'], 1)

    def test_unassessed_chain_cannot_receive_germline_match_tile(self):
        from app.report_brief import signals
        assessment = {'sites': [], 'patch_risk': {'categories': []}}
        data = {'chains': {'H': {}, 'L': {}}, 'residue_table': {
            'germline_counts': {'match': 1}, 'rows': [
                {'chain': 'H', 'imgt': '10', 'imgt_region': 'Framework', 'germline_match': True},
                {'chain': 'L', 'imgt': '10', 'imgt_region': 'Framework', 'germline_match': None}]}}
        tiles = signals(data, assessment)[0]
        self.assertEqual(next(state for name, state, detail in tiles if name == 'Germline FR'), 'unknown')
        data['residue_table']['rows'][1]['germline_match'] = True
        tiles = signals(data, assessment)[0]
        self.assertEqual(next(state for name, state, detail in tiles if name == 'Germline FR'), 'good')

    def test_sequence_only_is_not_a_clear_structural_cdr_assessment(self):
        from app.report_brief import signals
        tiles = signals({'chains': {'H': {}}}, {'sites': [], 'patch_risk': {'categories': []}})[0]
        self.assertEqual(next(state for name, state, detail in tiles if name == 'CDR liabilities'), 'unknown')

    def test_unknown_negative_patch_does_not_hide_high_positive_patch(self):
        from app.report_brief import signals
        assessment = {'sites': [], 'patch_risk': {'categories': [
            {'kind': 'negative', 'priority': 'unknown', 'count': None, 'largest_area_a2': None},
            {'kind': 'positive', 'priority': 'high', 'count': 1, 'largest_area_a2': 200}]}}
        tiles = signals({}, assessment)[0]
        self.assertEqual(next(state for name, state, detail in tiles if name == 'Charge surface'), 'high')


if __name__ == '__main__':
    unittest.main()
