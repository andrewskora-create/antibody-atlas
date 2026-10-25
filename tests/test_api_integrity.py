import unittest

from fastapi.testclient import TestClient

from app.main import app, jobs, propermab_jobs
from tests.test_structure_analysis import HEAVY


class ApiIntegrityTests(unittest.TestCase):
    def test_report_is_built_from_server_record(self):
        client = TestClient(app)
        analysis = client.post('/api/analyze', json={'heavy': HEAVY, 'scheme': 'kabat'})
        self.assertEqual(analysis.status_code, 200)
        record = analysis.json()
        self.assertEqual(len(record['provenance']['input_sha256']), 64)
        self.assertEqual(len(record['provenance']['reference_sha256']), 64)
        pdf = client.post('/api/report', json={'analysis_id': record['analysis_id']})
        self.assertEqual(pdf.status_code, 200)
        self.assertTrue(pdf.content.startswith(b'%PDF-'))
        forged = client.post('/api/report', json={'analysis_id': 'missing', 'report': record})
        self.assertEqual(forged.status_code, 404)

    def test_propermab_report_rejects_another_structure_for_same_analysis(self):
        client = TestClient(app)
        record = client.post('/api/analyze', json={'heavy': HEAVY}).json()
        aid = record['analysis_id']
        jobs['test-new-fold'] = {'status':'complete','analysis_id':aid}
        propermab_jobs['test-old-pm'] = {'status':'complete','analysis_id':aid,
                                        'fold_job_id':'test-old-fold', 'results':{}}
        try:
            response = client.post('/api/report', json={'analysis_id':aid,
                'fold_job_id':'test-new-fold','propermab_job_id':'test-old-pm'})
            self.assertEqual(response.status_code,422)
            self.assertIn('different structure',response.json()['detail'])
        finally:
            jobs.pop('test-new-fold',None)
            propermab_jobs.pop('test-old-pm',None)


if __name__ == '__main__':
    unittest.main()
