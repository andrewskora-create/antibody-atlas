"""Regression checks for structure ownership and common imported panel formats."""
import unittest
import tempfile
import json
from pathlib import Path
from concurrent.futures import Future
from unittest.mock import patch

from fastapi.testclient import TestClient

from app import main as api
from app.batch_input import parse_csv, parse_fasta
from tests.test_position_and_batch import HEAVY, LIGHT


class StructureOwnership(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(api.app)
        self.jid = 'review-structure'
        self.aid = 'review-analysis'
        self.sequences = {'H': HEAVY}
        api.jobs[self.jid] = {
            'status': 'complete', 'scheme': 'imgt', 'ph': 7.4,
            'input_sha256': api._sequence_digest(self.sequences),
            'structure_analysis': {'owner': 'baseline'},
        }
        api.analysis_records[self.aid] = {
            'scheme': 'imgt', 'ph': 7.4,
            'chains': {'H': {'sequence': HEAVY}},
            'provenance': {'input_sha256': api._sequence_digest(self.sequences)},
        }
        self.numbering = patch.object(api, 'number_chain', return_value={
            'scheme': 'imgt', 'chain_type': 'H', 'residues': []})
        self.numbering.start()
        self.framework = patch.object(api, 'framework_identity', return_value={})
        self.framework.start()

    def tearDown(self):
        self.numbering.stop()
        self.framework.stop()
        api.jobs.pop(self.jid, None)
        api.analysis_records.pop(self.aid, None)

    def test_other_sequence_or_conditions_rejected_by_both_substitution_routes(self):
        for route in ('/api/substitutions', '/api/substitutions/position'):
            for changed in ({'heavy': 'A' + HEAVY[1:]}, {'ph': 6.0}, {'scheme': 'kabat'}):
                payload = {'heavy': HEAVY, 'fold_job_id': self.jid, 'position': 99,
                           'repack_top': 0, **changed}
                with self.subTest(route=route, changed=changed):
                    response = self.client.post(route, json=payload)
                    self.assertEqual(response.status_code, 422, response.text)
                    self.assertIn('must match', response.json()['detail'])

    def test_expired_structure_does_not_silently_fall_back_to_sequence_ranking(self):
        response = self.client.post('/api/substitutions/position', json={
            'heavy': HEAVY, 'position': 99, 'fold_job_id': 'expired', 'repack_top': 0})
        self.assertEqual(response.status_code, 422)

    def test_matching_structure_and_user_preserved_positions_reach_proposer(self):
        with patch.object(api, 'propose_substitutions', return_value={'chains': {}}) as propose, \
                patch.object(api, 'build_residue_table', return_value={}), \
                patch('app.clinical_frameworks.clinical_framework_precedents', return_value={'sites': []}):
            response = self.client.post('/api/substitutions', json={
                'heavy': HEAVY, 'fold_job_id': self.jid, 'preserve': {'H': [99]}})
        self.assertEqual(response.status_code, 200, response.text)
        self.assertEqual(propose.call_args.kwargs['preserve'], {'H': {99}})
        self.assertEqual(propose.call_args.kwargs['structure'], {'owner': 'baseline'})

    def test_upload_cannot_attach_another_molecule_to_existing_analysis(self):
        with patch('app.structure_input.prepare') as prepare:
            response = self.client.post('/api/structure',
                files={'file': ('example.pdb', b'placeholder')},
                data={'heavy': LIGHT, 'analysis_id': self.aid})
        self.assertEqual(response.status_code, 422, response.text)
        prepare.assert_not_called()

    def test_upload_rejects_invalid_ph_before_structure_work(self):
        with patch('app.structure_input.prepare') as prepare:
            response = self.client.post('/api/structure',
                files={'file': ('example.pdb', b'placeholder')}, data={'ph': '15'})
        self.assertEqual(response.status_code, 422, response.text)
        prepare.assert_not_called()


class PanelImport(unittest.TestCase):
    def test_batch_charge_value_remains_numeric_beside_its_flag(self):
        from app.batch_web import _row
        row = _row('Example', {'molecule_flags': {
            'fv': {'net_charge_assay_ph': 4.25},
            'flags': [{'id': 'net_charge', 'verdict': 'watch'}]}})
        self.assertEqual(row['net_charge'], 4.25)
        self.assertEqual(row['flag_net_charge'], 'watch')

    def test_batch_chain_labels_are_validated_before_scoring(self):
        from app.batch_web import _screen_one
        with patch('app.numbering.number_chain', return_value={'chain_type': 'K'}):
            with self.assertRaisesRegex(ValueError, 'chain label'):
                _screen_one('Example', {'H': LIGHT}, {'scheme': 'imgt', 'ph': 7.4})

    def test_wrapped_fasta_preserves_short_terminal_lines(self):
        for tail_length in (1, 5, 11):
            with self.subTest(tail_length=tail_length):
                panel = parse_fasta(f'>Example|H\n{HEAVY[:-tail_length]}\n{HEAVY[-tail_length:]}\n')
                self.assertTrue(panel.ok, panel.report())
                self.assertEqual(panel.constructs['Example']['H'], HEAVY)

    def test_plain_labels_still_allow_short_wrapped_terminal_lines(self):
        panel = parse_fasta(f'Example heavy chain\n{HEAVY[:-5]}\n{HEAVY[-5:]}\n')
        self.assertTrue(panel.ok, panel.report())
        self.assertEqual(panel.constructs['Example']['H'], HEAVY)

    def test_tab_delimited_and_bom_csv_inputs_are_read(self):
        from app.batch_web import _looks_like_csv
        for delimiter in (',', '\t'):
            with self.subTest(delimiter=delimiter):
                text = f'\ufeffname{delimiter}heavy{delimiter}light\nExample{delimiter}{HEAVY}{delimiter}{LIGHT}\n'
                self.assertTrue(_looks_like_csv(text))
                panel = parse_csv(text)
                self.assertTrue(panel.ok, panel.report())
                self.assertEqual(panel.constructs['Example'], {'H': HEAVY, 'L': LIGHT})

    def test_light_only_csv_matches_light_only_fasta_support(self):
        panel = parse_csv(f'name,light\nExample,{LIGHT}\n')
        self.assertTrue(panel.ok, panel.report())
        self.assertEqual(panel.constructs['Example'], {'L': LIGHT})
        self.assertTrue(panel.warnings)

    def test_malformed_csv_is_rejected_before_screening(self):
        panel = parse_csv(f'name,heavy\nExample,{HEAVY},extra\n')
        self.assertFalse(panel.ok)
        self.assertFalse(panel.constructs)

    def test_csv_names_cannot_escape_output_directory(self):
        for name in ('../outside', '/absolute', '..\\outside'):
            with self.subTest(name=name):
                panel = parse_csv(f'name,heavy\n{name},{HEAVY}\n')
                self.assertFalse(panel.ok)
                self.assertFalse(panel.constructs)

    def test_unicode_report_download_uses_rfc5987_filename(self):
        from app import batch_web
        name = 'Antibody α'
        with patch.object(batch_web, 'record_for', return_value={}), \
                patch('app.report.make_report', return_value=b'%PDF-placeholder'):
            response = TestClient(api.app).get('/api/batch/example/report/' + name)
        self.assertEqual(response.status_code, 200)
        self.assertIn("filename*=UTF-8''Antibody%20%CE%B1", response.headers['content-disposition'])


class BatchResume(unittest.TestCase):
    def _exercise(self, old_chains, new_chains, old_options=None, new_options=None):
        from app.batch_input import ParsedPanel
        from run_batch import run_batch, input_options_fingerprint, save_state
        options = {'resume': True, 'workers': 1, 'threads': 1, 'ph': 7.4,
                   'scheme': 'imgt', 'fold': False, 'pdf': False}
        previous = {**options, **(old_options or {})}
        current = {**options, **(new_options or {})}
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            (out / 'Example.json').write_text('{}')
            save_state(out, {'molecules': {
                'Example': {'status': 'complete', 'row': {'name': 'Example'},
                            'input_options_sha256': input_options_fingerprint(old_chains, previous)},
                'Removed': {'status': 'complete', 'row': {'name': 'Removed'}}}})
            future = Future()
            future.set_result({'name': 'Example', 'status': 'complete', 'row': {'name': 'Example'}})
            with patch('run_batch.ProcessPoolExecutor') as executor, patch('run_batch.signal.signal'):
                executor.return_value.__enter__.return_value.submit.return_value = future
                state = run_batch(ParsedPanel({'Example': new_chains}, [], []), out, current)
                called = executor.called
            saved = json.loads((out / 'batch_state.json').read_text())
        self.assertEqual(set(saved['molecules']), {'Example'})
        self.assertEqual(state['molecules']['Example']['input_options_sha256'],
                         input_options_fingerprint(new_chains, current))
        return called

    def test_same_name_with_different_sequence_reruns(self):
        self.assertTrue(self._exercise({'H': HEAVY}, {'H': 'A' + HEAVY[1:]}))

    def test_changed_ph_or_output_options_rerun(self):
        for options in ({'ph': 6.0}, {'pdf': True}, {'fold': True}):
            with self.subTest(options=options):
                self.assertTrue(self._exercise({'H': HEAVY}, {'H': HEAVY}, new_options=options))

    def test_identical_analysis_reuses_output_despite_worker_count_change(self):
        self.assertFalse(self._exercise({'H': HEAVY}, {'H': HEAVY}, new_options={'workers': 3}))


class AnalysisContextHandoffs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from tests.test_translation import back_translate
        cls.dna = {'heavy': back_translate(HEAVY), 'light': back_translate(LIGHT)}
        cls.analysis = api.analyze(api.AnalyzeRequest(**cls.dna))
        cls.expected = {'H': HEAVY, 'L': LIGHT}

    @classmethod
    def tearDownClass(cls):
        api.analysis_records.pop(cls.analysis['analysis_id'], None)

    def test_sequence_only_substitution_table_keeps_analysis_germline_and_liabilities(self):
        # Exercise the actual table join. Only candidate enumeration is stubbed;
        # numbering, metrics, germline assignment and table construction are real.
        with patch.object(api, 'propose_substitutions', return_value={'chains': {}}):
            result = api.substitutions(api.SubstitutionRequest(**self.dna))
        left, right = self.analysis['residue_table'], result['residue_table']
        self.assertEqual(right['germline_counts'], left['germline_counts'])
        self.assertNotIn('germline', right['missing_layers'])
        self.assertEqual([(r['chain'], r['sequence_position'], r['germline_aa'], r['germline_match'])
                          for r in right['rows']],
                         [(r['chain'], r['sequence_position'], r['germline_aa'], r['germline_match'])
                          for r in left['rows']])
        self.assertGreater(sum(bool(r['liabilities']) for r in right['rows']), 0)
        self.assertEqual(sum(bool(r['liabilities']) for r in right['rows']),
                         sum(bool(r['liabilities']) for r in left['rows']))

    def test_fold_dispatch_receives_translated_sequences_matching_analysis(self):
        with patch.object(api.importlib.util, 'find_spec', return_value=object()), \
                patch.object(api.threading, 'Thread') as thread:
            job = api.fold(api.FoldRequest(**self.dna, analysis_id=self.analysis['analysis_id']))
        try:
            self.assertEqual(thread.call_args.kwargs['args'][1], self.expected)
            self.assertEqual(job['input_sha256'], self.analysis['provenance']['input_sha256'])
        finally:
            api.jobs.pop(job['id'], None)

    def test_structure_upload_translates_before_protein_validation(self):
        def prepare(source, prepared, **kwargs):
            self.assertEqual(kwargs['sequences'], self.expected)
            prepared.write_text('END\n')
            return {'source_format': 'PDB', 'copy_used': 0, 'copies_available': 1,
                    'antigen_chains': []}
        with patch('app.structure_input.prepare', side_effect=prepare), \
                patch.object(api, 'analyze_structure', return_value={}):
            response = TestClient(api.app).post('/api/structure',
                files={'file': ('fixture.pdb', b'END\n')},
                data={**self.dna, 'analysis_id': self.analysis['analysis_id']})
        self.assertEqual(response.status_code, 200, response.text)
        job = response.json()
        try:
            self.assertEqual(job['sequences_from_structure'], self.expected)
            self.assertEqual(job['input_sha256'], self.analysis['provenance']['input_sha256'])
        finally:
            api.jobs.pop(job['id'], None)
            (api.RUNS / f"{job['id']}.pdb").unlink(missing_ok=True)


if __name__ == '__main__':
    unittest.main()
