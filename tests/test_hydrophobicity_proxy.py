"""Polar substitutions must not be presented as creating spatial patches."""
import pytest
from app.developability_delta import mutation_effect


def effect(replacement, exposure=.5, original='N'):
    return mutation_effect({'H': 'EVQLVESG' + original + 'GG'}, 'H', 9, replacement,
                           relative_sidechain_sasa=exposure, in_cdr=True)


def test_asparagine_to_serine_does_not_add_hydrophobic_burden():
    result = effect('S')
    descriptor = result['exposed_hydrophobicity']
    assert descriptor['hydropathy_delta'] > 0  # Less hydrophilic is not hydrophobic.
    assert descriptor['parent_exposure_weighted_positive_hydropathy_delta'] == 0
    assert result['penalty']['components']['exposed_hydrophobicity'] == 0
    assert not any(f['kind'] == 'exposed_hydrophobicity' for f in result['flags'])


def test_asparagine_to_isoleucine_increases_labeled_proxy():
    result = effect('I')
    descriptor = result['exposed_hydrophobicity']
    assert descriptor['parent_exposure_weighted_positive_hydropathy_delta'] == 2.25
    assert result['penalty']['components']['exposed_hydrophobicity'] > 0
    flag = next(f for f in result['flags'] if f['kind'] == 'exposed_hydrophobicity')
    assert 'proxy' in flag['text']
    assert 'not a measured or modeled spatial patch change' in flag['text']
    assert 'Adds exposed hydrophobic surface' not in flag['text']
    assert descriptor['basis'] == 'parent_exposure_positive_part_kyte_doolittle_proxy'


def test_reverse_edit_uses_same_positive_part_rule():
    result = effect('N', original='I')
    assert result['exposed_hydrophobicity']['parent_exposure_weighted_positive_hydropathy_delta'] == -2.25
    assert result['penalty']['components']['exposed_hydrophobicity'] == 0
    assert next(f for f in result['flags'] if f['kind'] == 'exposed_hydrophobicity')['severity'] == 'good'


@pytest.mark.parametrize('exposure', [None, float('nan'), float('inf'), -1, True])
def test_absent_or_invalid_parent_exposure_never_invents_surface_evidence(exposure):
    result = effect('I', exposure)
    assert result['exposed_hydrophobicity']['status'] == 'unavailable'
    assert result['penalty']['components']['exposed_hydrophobicity'] == 0
    assert not any(f['kind'] == 'exposed_hydrophobicity' for f in result['flags'])
