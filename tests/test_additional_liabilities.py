import unittest
from app.analysis import ptm_sites
from app.structure_analysis import _site_priority, _group_physical_sites


class AdditionalLiabilitiesTests(unittest.TestCase):
    def test_only_actual_input_n_terminus_gets_conditional_pyro_flag(self):
        self.assertTrue(any(s['kind']=='N-terminal pyroglutamate formation' for s in ptm_sites('QAQE')))
        self.assertFalse(any(s['kind']=='N-terminal pyroglutamate formation' for s in ptm_sites('AQAQE')))

    def test_cysteine_identity_not_count_or_parity(self):
        numbering={'scheme':'imgt','chain_type':'H','residues':[
            {'sequence_position':1,'number':'H23','aa':'C'},
            {'sequence_position':2,'number':'H104','aa':'C'},
            {'sequence_position':3,'number':'H104A','aa':'C'},
            {'sequence_position':4,'number':'H50','aa':'C'}]}
        flags=[s for s in ptm_sites('CCCC',numbering) if s['kind']=='Noncanonical cysteine']
        self.assertEqual([s['start'] for s in flags],[3,4])
        self.assertFalse(any(s['kind']=='Noncanonical cysteine' for s in ptm_sites('CCCC')))
        self.assertTrue(all('May form' in s['interpretation'] for s in flags))

    def test_additional_chemistries_do_not_become_oxidation(self):
        flags=ptm_sites('QDP')
        scored=[]
        for flag in flags:
            result=_site_priority(flag,{'chain':'H','aa':flag['motif'][0],'number':'H1',
                'region':'Framework','relative_sidechain_sasa':0,'sidechain_sasa_a2':0},'imgt')
            self.assertEqual(result['priority'],'lower' if flag['start']==1 else 'review')
            scored.append(result)
        chemistries={f['chemistry'] for site in _group_physical_sites(scored) for f in site['flags']}
        self.assertEqual(chemistries,{'pyroglutamate','cleavage'})

    def test_cdr_trp_oxidation_keeps_exposure_review(self):
        flag=next(s for s in ptm_sites('AWA') if 'Oxidation' in s['kind'])
        site=_site_priority(flag,{'chain':'H','aa':'W','number':'H58','region':'CDR2',
            'relative_sidechain_sasa':0.6,'sidechain_sasa_a2':80},'imgt')
        self.assertEqual(site['priority'],'review')
        self.assertIn('CDR',site['functional_impact'])
