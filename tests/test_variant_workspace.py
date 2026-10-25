"""Scientific and API regressions for structure-linked engineering comparisons."""
import copy
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB import PDBParser
from fastapi.testclient import TestClient

from app import main as api
from app.variant_comparison import compare_records
from app.variant_jobs import VariantRequest


def record(h='ADK', l='AEK'):
    chains = {label: {'sequence': seq, 'metrics': {'ptm_sites': []},
                     'numbering': {'chain_type': 'K' if label == 'L' else 'H', 'residues': []}}
              for label, seq in [('H', h), ('L', l)]}
    return {'chains': chains, 'scheme': 'imgt', 'ph': 7.4, 'provenance': {}}


def test_light_chain_change_reaches_whole_fv_charge():
    a, b = record(), record(l='AKK')
    result = compare_records(a, b)['molecule']
    metric = result['metrics']['net_charge_at_ph']
    expected = sum(ProteinAnalysis(b['chains'][key]['sequence']).charge_at_pH(7.4) for key in ('H', 'L'))
    assert metric['candidate'] == round(expected, 2)
    assert 1.9 < metric['delta'] < 2.1
    assert result['per_chain']['H']['metrics']['net_charge_at_ph']['delta'] == 0
    assert result['per_chain']['L']['metrics']['net_charge_at_ph']['delta'] > 1.9
    assert result['per_chain']['L']['copies'] == 1


def test_igg_charge_counts_two_mutated_chains_and_same_constants():
    from app.constant_regions import assemble
    a, b = record(), record(l='AKK')
    result = compare_records(a, b, isotype='IGHG1')['molecule']
    assembled = assemble(b['chains'], 'IGHG1')
    expected = sum(assembled[key]['copies'] * ProteinAnalysis(assembled[key]['sequence']).charge_at_pH(7.4)
                   for key in ('heavy', 'light'))
    assert result['metrics']['net_charge_at_ph']['candidate'] == round(expected, 2)
    assert 3.9 < result['metrics']['net_charge_at_ph']['delta'] < 4.1
    assert result['per_chain']['L']['copies'] == 2
    assert result['format'] == 'H2L2'
    assert 'full-length structure' in result['structure_scope']
    # Same molecular pI is computed from the summed charge, not average chain pIs.
    pi = result['metrics']['isoelectric_point']['candidate']
    charge_at_pi = sum(assembled[key]['copies'] * ProteinAnalysis(assembled[key]['sequence']).charge_at_pH(pi)
                      for key in ('heavy', 'light'))
    assert abs(charge_at_pi) < 0.3


def test_isotype_selection_validates_and_missing_reference_stays_unknown():
    from app.constant_regions import ConstantRegionError
    assert VariantRequest(analysis_id='x', variants=[{'name':'v','heavy':'ADK'}]).isotype == 'fv'
    with pytest.raises(ValueError):
        VariantRequest(analysis_id='x', variants=[{'name':'v','heavy':'ADK'}], isotype='invented')
    with patch('app.constant_regions.assemble', side_effect=ConstantRegionError('Missing reference')):
        result = compare_records(record(), record(l='AKK'), isotype='IGHG1')['molecule']
    assert result['status'] == 'unavailable'
    assert result['metrics'] == {}


def test_structure_selection_uses_pdb_insertion_code_not_sequence_index():
    a, b = record(), record(l='AKK')
    for item in (a, b):
        digest = api._sequence_digest({key: chain['sequence'] for key, chain in item['chains'].items()})
        item['provenance'] = dict(input_sha256=digest, structure_input_sha256=digest,
            structure_model='test', structure_model_version='1', structure_pdb_sha256='pdb')
        item['structure_analysis'] = {'scheme':'imgt', 'ph':7.4, 'chain_mapping':{'H':{},'L':{}},
            'residues':[{'chain':'L','sequence_position':2,'pdb_number':'112A'}]}
    mutation = compare_records(a,b)['mutations'][0]
    assert mutation['candidate_structure'] == {'chain':'L','residue':'112A'}
    b['provenance']['structure_input_sha256'] = 'wrong'
    assert compare_records(a,b)['mutations'][0]['candidate_structure'] is None


def test_repack_api_exports_owned_files_resolves_position_and_rejects_ph_mismatch(tmp_path):
    jid = 'workspace-repack-test'
    path = tmp_path / f'{jid}.pdb'; path.write_text('MODEL\nENDMDL\n')
    api.jobs[jid] = {'status':'complete', 'structure_source':'experimental',
                    'structure_analysis': {'ph':7.4, 'residues':[
                        {'chain':'H', 'sequence_position':12, 'pdb_number':'15'}]}}
    args = dict(fold_job_id=jid, chain='H', position=12, replacement='A', ph=7.4)
    def fake_repack(pdb_path, chain, resid, target, **kw):
        assert resid == '15'  # sequence 12 is not PDB 12
        for destination in kw['output_paths']:
            assert destination.parent == tmp_path
            destination.write_text('MODEL\nENDMDL\n')
        return {'method':'same relaxation', 'no_free_energy':'No ΔΔG.', 'hydrogen_bonds':{},'position':resid}
    client = TestClient(api.app)
    try:
        with patch.object(api, 'RUNS', tmp_path), patch('app.repack.repack', side_effect=fake_repack) as run:
            response = client.post('/api/repack', json=args)
            assert response.status_code == 200, response.text
            body = response.json()
            assert body['sequence_position'] == 12
            assert body['structures']['selection'] == {'chain':'H','residue':'15'}
            for key in ('parent_pdb_url','mutant_pdb_url'):
                assert client.get(body[key]).status_code == 200
                assert str(tmp_path) not in body[key]
            assert body['free_energy']['status'] == 'unavailable'
            assert client.get('/api/repack/not-registered/baseline/pdb').status_code == 404
            assert client.post('/api/repack', json={**args,'ph':6}).status_code == 422
            assert client.post('/api/repack', json={**args,'position':13}).status_code == 422
            assert run.call_count == 1
    finally:
        api.jobs.pop(jid, None)
        api.repack_previews.clear()


def test_actual_repack_pdbs_are_relaxed_parent_and_mutant(tmp_path):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    from app.repack import repack
    fixture = Path(__file__).parent / 'fixtures' / 'example_fv.pdb'
    if not fixture.exists():
        pytest.skip('Public example structure not present')
    paths = (tmp_path/'parent.pdb',tmp_path/'mutant.pdb')
    result = repack(fixture, 'H', '35', 'A', output_paths=paths)
    parser = PDBParser(QUIET=True)
    parent, mutant = [parser.get_structure(str(i), str(path))[0] for i,path in enumerate(paths)]
    assert mutant['H'][35].resname == 'ALA'
    assert parent['H'][35].resname != 'ALA'
    assert 'retained' in result['hydrogen_bonds']
    assert result['hydrogen_bonds']['retained_at_site'] == len(result['hydrogen_bonds']['retained'])
    # Atoms far outside the selected neighbourhood remain in a common frame.
    assert (parent['L'][1]['CA'].coord == mutant['L'][1]['CA'].coord).all()
    assert result['clashes']['delta'] == result['clashes']['mutant'] - result['clashes']['wild_type']


def test_repack_refuses_uploaded_structure_of_a_different_sequence():
    from tests.test_api_endpoints import FV, HEAVY, LIGHT
    client = TestClient(api.app)
    with FV.open('rb') as handle:
        uploaded = client.post('/api/structure',
            files={'file':('example_fv.pdb', handle, 'chemical/x-pdb')},
            data={'heavy': HEAVY, 'light': LIGHT, 'ph':'7.4','scheme':'imgt'})
    assert uploaded.status_code == 200
    job = uploaded.json()
    assert not job['ingestion']['usable']  # Fixture is a different public antibody.
    with patch('app.repack.repack') as run:
        response = client.post('/api/repack',json={
            'fold_job_id':job['id'], 'chain':'H','position':35,'replacement':'A','ph':7.4})
        assert response.status_code == 422
        assert 'coordinate mapping' in response.json()['detail']
        run.assert_not_called()
