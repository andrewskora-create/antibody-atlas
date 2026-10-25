import unittest

from tests.support import requires_propermab
from unittest.mock import patch
from tests.test_structure_analysis import HEAVY, LIGHT, PDB
from app.main import analyze, AnalyzeRequest
from app.propermab import run, _terminal_identity

@requires_propermab()
class PropermabTests(unittest.TestCase):
    def test_real_offline_sequence_descriptors(self):
        record = analyze(AnalyzeRequest(heavy=HEAVY, light=LIGHT))
        result = run(record)
        features = {f['key']: f['value'] for f in result['features']}
        self.assertEqual(features['fv_charge'], features['vh_charge'] + features['vl_charge'])
        self.assertEqual(features['n_charged_res_fv'], 40)
        self.assertIsNone(result['model_checkpoint'])
        self.assertTrue(any(x['feature']=='Structure descriptors' for x in result['unavailable']))

    def test_mismatched_structure_retains_sequence_results(self):
        record = analyze(AnalyzeRequest(heavy=HEAVY, light=LIGHT))
        with patch('app.propermab._prepare_imgt_pdb', side_effect=ValueError('Model differs from input.')):
            result = run(record, PDB)
        self.assertEqual(len(result['features']), 5)
        self.assertTrue(any(x['reason']=='Model differs from input.' for x in result['unavailable']))

    def test_failed_repair_retains_original_measurements(self):
        from types import SimpleNamespace
        baseline={'structure_partial':True,'features':[{'key':'hyd_patch_area','value':100}],
                  'structure_coverage':{'L':{'omitted_positions':[108]}},'unavailable':[]}
        with patch('app.propermab._run_single',return_value=baseline), patch('app.propermab.subprocess.run',return_value=SimpleNamespace(returncode=1)):
            result=run({'chains':{'H':{'sequence':HEAVY},'L':{'sequence':LIGHT+'S'}},'ph':7.4},PDB)
        self.assertEqual(result['features'][0]['value'],100)
        self.assertEqual(result['terminal_repair']['status'],'failed')

class TerminalIdentityTests(unittest.TestCase):
    def test_exact_and_terminal_omissions(self):
        self.assertEqual(_terminal_identity(LIGHT,LIGHT,'L')['omitted_positions'],[])
        c=_terminal_identity(LIGHT+'S',LIGHT,'L')
        self.assertEqual(c['omitted_positions'],[len(LIGHT)+1])
        self.assertEqual(c['input_positions'],list(range(1,len(LIGHT)+1)))
        c=_terminal_identity('A'+LIGHT,LIGHT,'L')
        self.assertEqual(c['omitted_positions'],[1])
        self.assertEqual(c['input_positions'][0],2)

    def test_reject_substitutions_internal_gaps_and_excessive_trimming(self):
        for observed in (LIGHT[:30]+'A'+LIGHT[31:], LIGHT[:30]+LIGHT[31:], LIGHT[:-3], ''):
            with self.subTest(length=len(observed)), self.assertRaises(ValueError):
                _terminal_identity(LIGHT,observed,'L')

    def test_partial_structure_cannot_use_complete_structure_reference(self):
        from app.profiler_interpretation import _protocol_compatible
        from app.propermab import comparison_protocol
        protocol=comparison_protocol(7.4,True)
        partial={**protocol,'partial_structure':True}
        self.assertFalse(_protocol_compatible(partial,protocol,'hyd_asa'))
        self.assertTrue(_protocol_compatible(partial,protocol,'fv_charge'))

if __name__ == '__main__':
    unittest.main()
