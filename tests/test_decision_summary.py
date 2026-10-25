"""Overview must distinguish observed concerns, lower-priority motifs, and gaps."""
from app.report_brief import signals


def tiles(data=None, categories=None):
    result = signals(data or {}, {'sites': [], 'patch_risk': {'categories': categories or []}})
    return {name: (state, detail) for name, state, detail in result[0]}


def patch(kind, priority):
    return {'kind': kind, 'priority': priority, 'count': 1,
            'largest_area_a2': None if priority == 'unknown' else 100}


def test_lower_priority_cdr_is_visible_without_escalating_it():
    data = {'structure_analysis': {'ptm_sites': [
        {'region': 'CDR1', 'priority': 'lower', 'relative_sidechain_sasa': 0.1}]}}
    state, detail = tiles(data)['CDR liabilities']
    assert state == 'good'
    assert 'No high/review-priority' in detail
    assert '1 lower-priority motif sites' in detail
    assert 'none flagged' not in detail


def test_missing_scm_does_not_create_concentration_concern():
    state, detail = tiles(categories=[patch('negative', 'lower'), patch('positive', 'lower')])['Concentration']
    assert state == 'good'
    assert 'No charge-patch concern flagged' in detail
    assert 'optional SCM not assessed' in detail


def test_concentration_review_names_the_observed_concern():
    state, detail = tiles(categories=[patch('negative', 'lower'), patch('positive', 'high')])['Concentration']
    assert state == 'review'
    assert 'Positive patch concern' in detail
    assert 'optional SCM not assessed' in detail


def test_incomplete_charge_coverage_cannot_read_as_a_clear_concentration_proxy():
    state, detail = tiles(categories=[patch('negative', 'lower'), patch('positive', 'unknown')])['Concentration']
    assert state == 'unknown'
    assert 'evidence incomplete' in detail


def test_unknown_lower_priority_exposure_is_not_called_assessed():
    data = {'structure_analysis': {'ptm_sites': [
        {'region': 'CDR1', 'priority': 'lower', 'relative_sidechain_sasa': None}]}}
    state, detail = tiles(data)['CDR liabilities']
    assert state == 'unknown'
    assert 'lower-priority motif sites' not in detail
    assert '1 unassessed' in detail
