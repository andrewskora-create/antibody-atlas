import unittest
from copy import deepcopy
from app.model_selection import select_propermab_structure


def fixture():
    report={'chains':{'H':{},'L':{}},'scheme':'imgt','ph':7.4,
            'provenance':{'structure_pdb_sha256':'original'},
            'structure_analysis':{'original':True}}
    result={'structure_selection':'terminal_repaired','original_input_pdb_sha256':'original',
            'input_pdb_sha256':'repaired','terminal_repair':{'status':'complete','repaired_pdb_sha256':'repaired'},
            'selected_structure_analysis':{'mapping_coverage_percent':{'H':100,'L':100},'scheme':'imgt','ph':7.4}}
    return report,result

class ModelSelectionTests(unittest.TestCase):
    def test_selects_full_model_and_preserves_original_across_repeated_selection(self):
        report,result=fixture()
        self.assertTrue(select_propermab_structure(report,result))
        self.assertEqual(report['provenance']['structure_pdb_sha256'],'repaired')
        self.assertEqual(report['original_structure_analysis'],{'original':True})
        self.assertTrue(select_propermab_structure(report,result))
        self.assertEqual(report['original_structure_analysis'],{'original':True})
        result['selected_structure_analysis']['ph']=6
        self.assertEqual(report['structure_analysis']['ph'],7.4)

    def test_failed_lineage_or_settings_never_promotes(self):
        for mutate in (
            lambda r:r.update(original_input_pdb_sha256='other'),
            lambda r:r.update(input_pdb_sha256='other'),
            lambda r:r['terminal_repair'].update(status='failed'),
            lambda r:r['selected_structure_analysis']['mapping_coverage_percent'].update(L=99),
            lambda r:r['selected_structure_analysis']['mapping_coverage_percent'].pop('L'),
            lambda r:r['selected_structure_analysis'].update(scheme='kabat'),
            lambda r:r['selected_structure_analysis'].update(ph=6),
        ):
            report,result=fixture();before=deepcopy(report);mutate(result)
            with self.assertRaises(ValueError):select_propermab_structure(report,result)
            self.assertEqual(report,before)

    def test_original_result_not_promoted(self):
        report,result=fixture();result.pop('structure_selection')
        self.assertFalse(select_propermab_structure(report,result))
        self.assertEqual(report['provenance']['structure_pdb_sha256'],'original')
