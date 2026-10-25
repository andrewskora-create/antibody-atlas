"""Exposure remains available for a selected residue without a PTM motif."""
import math
from unittest.mock import patch

import pytest

from app.cdr_substitutions import _exposure_at, propose_at_position
from tests.test_position_and_batch import chains


def test_nonmotif_position_uses_its_own_chain_sequence_mapping():
    structure = {'ptm_sites': [], 'residues': [
        {'chain':'L','sequence_position':52,'pdb_number':'60','relative_sidechain_sasa':.01},
        {'chain':'H','sequence_position':52,'pdb_number':'61','relative_sidechain_sasa':.55},
        {'chain':'H','sequence_position':61,'pdb_number':'70','relative_sidechain_sasa':.02}]}
    exposure = _exposure_at(structure, 'H', 52)
    assert exposure['status'] == 'available'
    assert exposure['relative_sidechain_sasa'] == .55
    assert exposure['exposure_band'] == 'high'
    assert exposure['pdb_number'] == '61'
    assert exposure['source'] == 'mapped_residue'
    assert _exposure_at(structure, 'H', 53)['status'] == 'unavailable'


def test_valid_ensemble_ptm_value_takes_precedence():
    structure = {'ptm_sites': [{'chain':'H','sequence_position':52,'relative_sidechain_sasa':.55}],
                 'residues': [{'chain':'H','sequence_position':52,'relative_sidechain_sasa':.05}]}
    assert _exposure_at(structure,'H',52)['relative_sidechain_sasa'] == .55
    structure['ptm_sites'][0]['relative_sidechain_sasa'] = None
    assert _exposure_at(structure,'H',52)['relative_sidechain_sasa'] == .05


@pytest.mark.parametrize('value', [None, math.nan, math.inf, -1, True, '0.2'])
def test_invalid_exposure_stays_unknown(value):
    structure = {'residues':[{'chain':'H','sequence_position':52,'relative_sidechain_sasa':value}]}
    assert _exposure_at(structure,'H',52)['status'] == 'unavailable'


@pytest.mark.parametrize('value,band', [(0,'low'),(.1,'intermediate'),(.3,'high')])
def test_bands_use_same_boundaries_as_structural_screen(value, band):
    structure = {'residues':[{'chain':'H','sequence_position':52,'relative_sidechain_sasa':value}]}
    assert _exposure_at(structure,'H',52)['exposure_band'] == band


def test_user_position_passes_mapped_exposure_into_all_candidate_calculations():
    data = chains()
    # Select a serine that the sequence PTM scan does not itself flag.
    from app.analysis import ptm_sites
    blocked = {s['start'] for s in ptm_sites(data['H']['sequence'], data['H']['numbering_imgt'])}
    position = next(i for i,aa in enumerate(data['H']['sequence'],1) if aa == 'S' and i not in blocked)
    structure = {'ptm_sites': [], 'residues': [
        {'chain':'H','sequence_position':position,'pdb_number':'61','relative_sidechain_sasa':.55}]}
    with patch('app.cdr_substitutions.mutation_effect', return_value={'status':'unavailable'}) as calculate:
        result = propose_at_position(data,'H',position,structure=structure)
    assert result['exposure']['status'] == 'available'
    assert calculate.call_count == 19
    assert all(call.kwargs['relative_sidechain_sasa'] == .55 for call in calculate.call_args_list)
