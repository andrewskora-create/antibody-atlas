import copy
import unittest

from tests.support import requires_dataset
from unittest.mock import patch
from app import clinical_frameworks as cf

@requires_dataset("clinical_frameworks")
class ClinicalFrameworkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference=next(r for r in cf._reference_data()['references'] if r['drug'].lower()=='bevacizumab' and r['family']=='H')
    def run_match(self,query=None,reference=None,scheme='imgt',motifs=None):
        ref=copy.deepcopy(reference or self.reference)
        q=copy.deepcopy(query or self.reference)
        data={'version':'test','coverage':{'reference_chains':1},'references':[ref]}
        with patch.object(cf,'_reference_data',return_value=data),patch.object(cf,'number_chain',side_effect=lambda seq,s:q['numbering'][s]):
            if motifs is None:return cf.clinical_framework_precedents({'H':{'sequence':q['sequence']}},scheme)
            with patch.object(cf,'ptm_sites',return_value=motifs):return cf.clinical_framework_precedents({'H':{'sequence':q['sequence']}},scheme)
    def target(self,result):return next(s for s in result['sites'] if s['position']==84 and s['motif']=='NS')
    def test_exact_public_framework(self):
        result=self.run_match(); match=self.target(result)['matches'][0]
        self.assertEqual(match['framework_identity_percent'],100)
        self.assertEqual(match['match_type'],'Exact framework')
        self.assertEqual(match['query_context'],match['reference_context'])
    def test_kabat_consistency(self):
        site=self.target(self.run_match(scheme='kabat'))
        self.assertEqual(site['matches'][0]['framework_identity_percent'],100)
        self.assertNotEqual(site['number'],site['query_imgt_number'])
    def test_subtype_mismatch(self):
        ref=copy.deepcopy(self.reference)
        ref['numbering']['imgt']['chain_type']='K'
        self.assertEqual(self.target(self.run_match(reference=ref))['matches'],[])
    def test_motif_mismatch(self):
        ref=copy.deepcopy(self.reference)
        next(r for r in ref['numbering']['imgt']['residues'] if r['sequence_position']==85)['aa']='A'
        self.assertEqual(self.target(self.run_match(reference=ref))['matches'],[])
    def test_reference_cdr_not_framework(self):
        ref=copy.deepcopy(self.reference)
        next(r for r in ref['numbering']['imgt']['residues'] if r['sequence_position']==84)['region']='CDR2'
        self.assertEqual(self.target(self.run_match(reference=ref))['matches'],[])
    def test_union_gap_and_threshold(self):
        ref=copy.deepcopy(self.reference)
        rows=ref['numbering']['imgt']['residues']
        removable=[r for r in rows if r['region']=='Framework' and r['sequence_position']<70]
        rows.remove(removable[0])
        match=self.target(self.run_match(reference=ref))['matches'][0]
        self.assertEqual(match['framework_matches'],match['framework_denominator']-1)
        self.assertEqual(match['match_type'],'Similar framework and aligned site')
        for r in removable[1:12]: rows.remove(r)
        self.assertEqual(self.target(self.run_match(reference=ref))['matches'],[])
    def test_query_region_change_rejects_reference(self):
        q=copy.deepcopy(self.reference)
        next(r for r in q['numbering']['imgt']['residues'] if r['sequence_position']==84)['region']='CDR2'
        self.assertEqual(self.target(self.run_match(query=q))['matches'],[])
    def test_mislabeled_query_chain_is_skipped(self):
        with patch.object(cf,'number_chain',return_value=self.reference['numbering']['imgt']):
            result=cf.clinical_framework_precedents({'L':{'sequence':self.reference['sequence']}})
        self.assertEqual(result['sites'],[])
        self.assertEqual(len(result['skipped_chains']),1)
    def test_exact_cdr_matches_despite_different_framework(self):
        ref=copy.deepcopy(self.reference)
        for r in ref['numbering']['imgt']['residues']:
            if r['region']=='Framework':r['aa']='A' if r['aa']!='A' else 'V'
        site=next(s for s in self.run_match(reference=ref)['sites'] if s['position']==52)
        match=site['matches'][0]
        self.assertEqual(site['query_region'],'CDR2')
        self.assertEqual(match['reference_region'],'CDR2')
        self.assertEqual(match['match_type'],'Exact CDR2 sequence')
        self.assertEqual(match['region_identity_percent'],100)
        self.assertLess(match['framework_identity_percent'],90)
    def test_changed_cdr_outside_motif_rejected(self):
        ref=copy.deepcopy(self.reference)
        r=next(r for r in ref['numbering']['imgt']['residues'] if r['region']=='CDR2' and r['sequence_position']>53)
        r['aa']='A' if r['aa']!='A' else 'V'
        site=next(s for s in self.run_match(reference=ref)['sites'] if s['position']==52)
        self.assertEqual(site['matches'],[])
    def test_cdr_region_mismatch_rejected(self):
        ref=copy.deepcopy(self.reference)
        next(r for r in ref['numbering']['imgt']['residues'] if r['sequence_position']==52)['region']='CDR1'
        site=next(s for s in self.run_match(reference=ref)['sites'] if s['position']==52)
        self.assertEqual(site['matches'],[])
    def test_cdr_gap_rejected(self):
        ref=copy.deepcopy(self.reference)
        rows=ref['numbering']['imgt']['residues']
        rows.remove(next(r for r in rows if r['region']=='CDR2' and r['sequence_position']>53))
        site=next(s for s in self.run_match(reference=ref)['sites'] if s['position']==52)
        self.assertEqual(site['matches'],[])
    def test_all_eligible_rows_and_precomputed_schemes(self):
        data=cf._reference_data()
        self.assertEqual(data['coverage']['eligible_drugs'],365)
        self.assertEqual(data['coverage']['represented_drugs'],365)
        self.assertEqual(data['coverage']['reference_chains'],len(data['references']))
        for r in data['references']:
            self.assertIn(r['stage'],{'Approved','Phase-III'})
            for scheme in ('imgt','kabat'):
                self.assertEqual(''.join(x['aa'] for x in r['numbering'][scheme]['residues']),r['sequence'])

if __name__=='__main__':unittest.main()
