"""Stability estimates stay attached to the matching modeled substitution."""
from unittest.mock import patch
from fastapi.testclient import TestClient
from app import main as api


def test_local_preview_caches_coordinates_and_carries_separate_prediction(tmp_path):
    jid = 'stability-cache-test'
    (tmp_path / f'{jid}.pdb').write_text('test coordinate identity')
    api.jobs[jid] = {'status': 'complete', 'ph': 7.4, 'structure_analysis': {
        'residues': [{'chain': 'H', 'sequence_position': 9, 'pdb_number': '12A', 'aa': 'N'}]}}
    args = dict(fold_job_id=jid, chain='H', position=9, replacement='Q', ph=7.4)
    def model(path, chain, resid, target, **kwargs):
        assert resid == '12A'
        for output in kwargs['output_paths']:
            output.write_text('coordinates')
        return {'from': 'N', 'to': target, 'method': 'local model',
                'hydrogen_bonds': {}, 'no_free_energy': 'Geometry only'}
    energy = {'status': 'available', 'value_kcal_mol': -0.3, 'method': 'Test predictor'}
    client = TestClient(api.app)
    try:
        with patch.object(api, 'RUNS', tmp_path), patch('app.repack.repack', side_effect=model) as repack, \
             patch('app.stability_energy.predict_stability', return_value=energy) as predict:
            first = client.post('/api/repack', json=args)
            second = client.post('/api/repack', json=args)
            assert first.status_code == second.status_code == 200
            assert first.json()['free_energy'] == energy
            assert second.json()['mutant_pdb_url'] == first.json()['mutant_pdb_url']
            assert repack.call_count == predict.call_count == 1
            predict.assert_called_with(tmp_path/f'{jid}.pdb', 'H', '12A', 'N', 'Q')
            # Different sequence position cannot reuse a prediction.
            assert client.post('/api/repack', json={**args, 'position': 10}).status_code == 422
            # New structure content invalidates coordinate and energy caches.
            (tmp_path/f'{jid}.pdb').write_text('changed coordinates')
            assert client.post('/api/repack', json=args).status_code == 200
            assert repack.call_count == predict.call_count == 2
    finally:
        api.jobs.pop(jid, None)
        api.repack_results.clear()
        api.repack_previews.clear()


def test_default_topfive_and_on_demand_bounds():
    assert api.PositionRequest(heavy='EVQL', position=1).repack_top == 5
    assert api.PositionRequest(heavy='EVQL', position=1, repack_top=19).repack_top == 19
