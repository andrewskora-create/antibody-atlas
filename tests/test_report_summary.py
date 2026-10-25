import re
import unittest
from app.report import make_report, priority_sites

class ShortReportTests(unittest.TestCase):
    def test_report_page_count_is_fixed_regardless_of_input_size(self):
        sites=[{'chain':'H','number':f'H{i}','sequence_position':i,'region':'Framework',
                'kind':'Oxidation-susceptible residue','motif':'M','priority':'review',
                'relative_sidechain_sasa':0.2} for i in range(1,1201)]
        cases=[{'antibody':'Extremely long therapeutic name '*100,'site':'HC Met1',
                'observation':'Detailed evidence '*500,'url':'https://example.org/study',
                'matches':[{'query_chain':'H','query_sequence_position':1,'query_number':'H1','query_region':'Framework'}]}]*1000
        data={'chains':{},'structure_analysis':{'ptm_sites':sites},'ptm_precedents':{'cases':cases}}
        pdf=make_report(data,'Long report title '*1000)
        # Five fixed pages: overview, priority sites, biophysical properties,
        # precedents. The count must not grow with input size.
        self.assertEqual(len(re.findall(rb'/Type\s*/Page\b',pdf)),5)
        self.assertLess(len(pdf),50000)

    def test_missing_models_still_produce_bounded_report(self):
        pdf=make_report({'chains':{}})
        self.assertEqual(len(re.findall(rb'/Type\s*/Page\b',pdf)),5)

    def test_priority_selection_preserves_high_sites_and_deduplicates(self):
        base={'chain':'H','region':'Framework','relative_sidechain_sasa':.5}
        sites=[{**base,'sequence_position':1,'priority':'lower'},
               {**base,'sequence_position':2,'priority':'high'},
               {**base,'sequence_position':2,'priority':'review'},
               {**base,'sequence_position':3,'priority':'unassessed'}]
        chosen=priority_sites({'structure_analysis':{'ptm_sites':sites}})
        self.assertEqual([s['sequence_position'] for s in chosen],[2,3,1])
