import unittest

from app.humanness import framework_identity
from app.numbering import number_chain
from app.report import make_report
from app.analysis import sequence_metrics

HEAVY = "EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS"


class HumannessReportTests(unittest.TestCase):
    def test_framework_identity_is_position_auditable(self):
        numbering = number_chain(HEAVY, 'imgt')
        result = framework_identity(HEAVY, numbering)
        self.assertEqual((result['nearest_human_v_gene'], result['matches'], result['framework_positions']),
                         ('IGHV3-30*02', 78, 80))
        self.assertEqual(result['percent'], 97.5)
        self.assertTrue(any(c['name'].startswith('IGHV3-48*') for c in result['candidate_alleles']))
        self.assertEqual(len(result['input_sha256_prefix']), 16)

    def test_sequence_only_report_is_a_pdf(self):
        numbering = number_chain(HEAVY, 'imgt')
        report = {'scheme': 'imgt', 'chains': {'H': {'sequence': HEAVY, 'numbering': numbering,
                  'humanness': framework_identity(HEAVY, numbering),
                  'metrics': sequence_metrics(HEAVY)}}}
        pdf = make_report(report, 'Test antibody')
        self.assertTrue(pdf.startswith(b'%PDF-'))
        self.assertGreater(len(pdf), 1000)


if __name__ == '__main__':
    unittest.main()
