import unittest
from app.patch_assessment import assess_patches

class PatchRiskTests(unittest.TestCase):
    def data(self):
        return {'chains':{'H':{},'L':{}},'provenance':{'structure_pdb_sha256':'same'},
                'structure_analysis':{'mapping_coverage_percent':{'H':100,'L':100},'patches':[
                    {'kind':'hydrophobic','near_cdr':True,'exposed_sidechain_sasa_sum_a2':180}]}}
    def calls(self,d):return {x['kind']:x for x in assess_patches(d)['categories']}
    def test_presence_medium_absence_lower(self):
        c=self.calls(self.data());self.assertEqual(c['hydrophobic']['priority'],'review');self.assertEqual(c['negative']['priority'],'lower')
    def test_missing_mapping_is_unknown_not_lower(self):
        d=self.data();d['structure_analysis']['mapping_coverage_percent']['L']=99
        c=self.calls(d);self.assertEqual(c['hydrophobic']['priority'],'review');self.assertTrue(c['hydrophobic']['provisional']);self.assertEqual(c['hydrophobic']['observed_priority'],'review')
        self.assertEqual(c['negative']['priority'],'unknown')
    def test_incomplete_chain_mapping_preserves_patch_category_labels(self):
        d=self.data()
        d['structure_analysis']['mapping_coverage_percent']['L']=99
        d['structure_analysis']['chain_mapping']={'H':{'input_residues':2},'L':{'input_residues':2}}
        c=self.calls(d)
        self.assertEqual([c[k]['label'] for k in ('hydrophobic','positive','negative')],
                         ['Hydrophobic patches','Positive-charge patches','Negative-charge patches'])
        self.assertEqual(c['hydrophobic']['call'],'Provisional · medium follow-up priority')
        self.assertTrue(any('L input positions 1, 2' in r for r in c['hydrophobic']['reasons']))

    def test_compatible_above_reference_high_and_wrong_model_not_used(self):
        d=self.data();d['propermab']={'input_pdb_sha256':'same','summary':{'benchmarks':[{'key':'pos_patch_area','reference_n':4,'query_relation':'above_range','query_value':300}]}}
        self.assertEqual(self.calls(d)['positive']['priority'],'high')
        d['propermab']['input_pdb_sha256']='other';self.assertEqual(self.calls(d)['positive']['priority'],'lower')
    def test_within_reference_nonzero_patch_remains_medium(self):
        d=self.data();d['propermab']={'input_pdb_sha256':'same','summary':{'benchmarks':[{'key':'neg_patch_area','reference_n':4,'query_relation':'within_range','query_value':50}]}}
        self.assertEqual(self.calls(d)['negative']['priority'],'review')
    def test_partial_propermab_and_missing_sasa_do_not_escalate(self):
        d=self.data();d['propermab']={'structure_partial':True,'input_pdb_sha256':'same','summary':{'benchmarks':[{'key':'hyd_patch_area','reference_n':4,'query_relation':'above_range'}]}}
        self.assertEqual(self.calls(d)['hydrophobic']['priority'],'review')
        d['structure_analysis']['residues']=[{'aa':'M','relative_sidechain_sasa':None}]
        self.assertEqual(self.calls(d)['hydrophobic']['priority'],'review')

    def anchored(self):
        d=self.data()
        d['chains']={'H':{'sequence':'A'},'L':{'sequence':'A'}}
        d['structure_analysis']['residues']=[{'chain':c,'sequence_position':1,'aa':'A','relative_sidechain_sasa':0.5} for c in ('H','L')]
        d['propermab']={'input_pdb_sha256':'same','summary':{'benchmarks':[
            {'key':key,'reference_n':35,'query_relation':'within_range','query_value':value,'reference_anchor':{'value':100}}
            for key,value in [('hyd_patch_area',110),('hyd_patch_area_cdr',90)]]}}
        return d

    def test_herceptin_both_inclusive_boundaries_lower_comparative(self):
        c=self.calls(self.anchored())['hydrophobic']
        self.assertTrue(c['herceptin_like']);self.assertEqual(c['priority'],'lower')
        self.assertIn('Low comparative risk',c['call'])

    def test_herceptin_one_outside_missing_or_zero_anchor_not_lower(self):
        for mode in ('outside','missing','zero','nan','above'):
            d=self.anchored();b=d['propermab']['summary']['benchmarks'][1]
            if mode=='outside':b['query_value']=89.99
            if mode=='missing':del b['reference_anchor']
            if mode=='zero':b['reference_anchor']['value']=0
            if mode=='nan':b['query_value']=float('nan')
            if mode=='above':b['query_relation']='above_range'
            self.assertFalse(self.calls(d)['hydrophobic']['herceptin_like'],mode)

    def test_herceptin_coverage_model_protocol_and_disable_gates(self):
        for mode in ('partial','coverage','model','disabled','five'):
            d=self.anchored()
            if mode=='partial':d['propermab']['structure_partial']=True
            if mode=='coverage':d['structure_analysis']['mapping_coverage_percent']['L']=99
            if mode=='model':d['propermab']['input_pdb_sha256']='other'
            if mode=='disabled':d['assessment_policy']={'herceptin_tolerance_percent':0}
            if mode=='five':d['assessment_policy']={'herceptin_tolerance_percent':5}
            self.assertFalse(self.calls(d)['hydrophobic']['herceptin_like'],mode)

    def test_herceptin_missing_or_duplicate_exposure_rows_fail_closed(self):
        for mode in ('empty','missing','duplicate','identity'):
            d=self.anchored();rows=d['structure_analysis']['residues']
            if mode=='empty':rows.clear()
            if mode=='missing':rows.pop()
            if mode=='duplicate':rows.append(dict(rows[0]))
            if mode=='identity':rows[0]['aa']='W'
            self.assertFalse(self.calls(d)['hydrophobic']['herceptin_like'],mode)
