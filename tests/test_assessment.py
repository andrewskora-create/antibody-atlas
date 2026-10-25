import unittest
from app.assessment import build_assessment


class AssessmentTests(unittest.TestCase):
    def site(self, position, kind, priority='review', region='Framework', exposure=.5):
        return dict(chain='H', sequence_position=position, number=str(position), start=position,
                    kind=kind, motif='DS', priority=priority, region=region, relative_sidechain_sasa=exposure)

    def test_focus_and_unknown_are_not_safety_calls(self):
        ox='Oxidation-susceptible residue'
        sites=[self.site(1,ox),self.site(2,ox,region='CDR3'),self.site(3,'Aspartate isomerization motif','high'),self.site(4,ox,exposure=None)]
        result=build_assessment({'structure_analysis':{'ptm_sites':sites}})
        self.assertEqual(result['counts'],dict(high=1,medium=1,unknown=1,deferred=1))
        self.assertEqual(result['status'],'Partial')
        self.assertEqual(len(result['paragraphs']),2)

    def test_clinical_context_matches_chain_position_and_chemistry_without_downgrade(self):
        kind='Aspartate isomerization motif'
        s=self.site(62,kind,'high')
        clinical=dict(query_chain='H',position=62,kind=kind,number='69',motif='DS',matches=[dict(drug='Public reference',stage='Approved')])
        case=dict(chemistry=kind,antibody='Reference',site='HC62',matches=[dict(query_chain='H',query_sequence_position=62)])
        result=build_assessment({'structure_analysis':{'ptm_sites':[s]},'clinical_framework_precedents':{'sites':[clinical,{**clinical,'query_chain':'L'}]},'ptm_precedents':{'cases':[case,{**case,'chemistry':'Other'}]}})
        self.assertEqual(result['counts']['high'],1)
        self.assertEqual(len(result['clinical_sites']),1)
        self.assertEqual(len(result['documented_matches']),1)
        self.assertIn('retained-motif match in Public reference',result['paragraphs'][1])
        self.assertIn('does not prove PTM tolerance',result['paragraphs'][1])

    def test_without_structure_all_sites_unknown(self):
        result=build_assessment({'chains':{'H':{'metrics':{'ptm_sites':[{'start':4},{'start':4},{'start':8}]}}}})
        self.assertEqual(result['counts']['unknown'],2)
        self.assertEqual(result['counts']['high'],0)
        self.assertTrue(result['status'].startswith('Preliminary'))
