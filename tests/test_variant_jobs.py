import copy
import unittest
from unittest.mock import patch
from fastapi.testclient import TestClient
from app import main as api
from app.variant_jobs import Variant, VariantRequest, candidate_sequences, _run, variant_jobs


def record(seq, aid):
    return {'analysis_id':aid, 'scheme':'imgt','ph':7.4,
        'chains':{'H':{'sequence':seq,'metrics':{'ptm_sites':[]},'numbering':{'residues':[]}}},
        'provenance':{'input_sha256':api._sequence_digest({'H':seq})}}


class VariantJobsTests(unittest.TestCase):
    def test_substitutions_reject_wrong_source_repeated_positions_and_indels(self):
        base=record('ACDE','base')
        good=Variant(name='good',mutations=[{'chain':'H','position':2,'from_aa':'C','to_aa':'S'}])
        self.assertEqual(candidate_sequences(base,good),({'H':'ASDE'},1))
        for v in [Variant(name='wrong',mutations=[{'chain':'H','position':2,'from_aa':'N','to_aa':'S'}]),
                  Variant(name='indel',heavy='ACD'),
                  Variant(name='repeat',mutations=[{'chain':'H','position':2,'from_aa':'C','to_aa':'S'}]*2)]:
            with self.assertRaises(ValueError): candidate_sequences(base,v)

    def test_model_failure_never_reuses_baseline_structure(self):
        baseline=record('A'*100,'variant-test-base')
        candidate=record('V'+'A'*99,'variant-test-candidate')
        api.analysis_records[candidate['analysis_id']]=copy.deepcopy(candidate)
        calls=[]
        def fold(jid,seqs,ph,scheme):
            calls.append(seqs['H'])
            if len(calls)==1:
                api.jobs[jid].update(status='complete',method='test-model',model_version='1',model_sha256='base-pdb',pdb_url='/test',
                    structure_analysis={'scheme':scheme,'ph':ph,'method':'test', 'patch_definition':'test',
                     'chain_mapping':{'H':{'input_residues':100,'mapped_residues':100,'model_residues':100}},'ptm_sites':[]})
            else: api.jobs[jid].update(status='failed',error='Deliberate candidate fold failure')
        v=Variant(name='candidate',heavy=candidate['chains']['H']['sequence'])
        request=VariantRequest(analysis_id=baseline['analysis_id'],variants=[v],run_structure=True)
        jid='variant-test-job';variant_jobs[jid]={'status':'queued'}
        try:
            with patch.object(api,'_fold',side_effect=fold),patch.object(api,'analyze',return_value=copy.deepcopy(candidate)):
                _run(jid,baseline,[(v,{'H':candidate['chains']['H']['sequence']},1)],request)
            result=variant_jobs[jid]['results'][0]
            self.assertEqual(len(calls),2)
            self.assertNotEqual(calls[0],calls[1])
            self.assertEqual(result['status'],'complete')
            self.assertEqual(result['model_status']['structure']['status'],'unavailable')
            self.assertEqual(result['comparison']['structure']['status'],'unavailable')
            self.assertEqual(result['comparison']['structure']['patch_counts'],{})
        finally:
            api.analysis_records.pop(candidate['analysis_id'],None);variant_jobs.pop(jid,None)

    def test_batch_orders_singles_first_and_rejects_invalid_before_dispatch(self):
        base=record('ACDE','variant-order-test');api.analysis_records[base['analysis_id']]=base
        client=TestClient(api.app)
        try:
            with patch('app.variant_jobs.threading') as thread_module:
                thread=thread_module.Thread
                response=client.post('/api/variants',json={'analysis_id':base['analysis_id'],'variants':[
                    {'name':'double','heavy':'VSDE'},{'name':'single','heavy':'VCDE'}]})
                self.assertEqual(response.status_code,200)
                prepared=thread.call_args.kwargs['args'][2]
                self.assertEqual([item[0].name for item in prepared],['single','double'])
                thread.reset_mock()
                bad=client.post('/api/variants',json={'analysis_id':base['analysis_id'],'variants':[{'name':'bad','heavy':'ACD'}]})
                self.assertEqual(bad.status_code,422);thread.assert_not_called()
        finally:api.analysis_records.pop(base['analysis_id'],None)
