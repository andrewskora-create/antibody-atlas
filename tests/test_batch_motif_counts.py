"""A sequence-only batch must not erase detected motif positions."""
from app.batch_web import _row

def test_sequence_batch_counts_motif_positions_without_exposure():
    row = _row('example', {'assessment': {'sites': []}, 'residue_table': {'rows': [
        {'imgt_region': 'CDR1', 'liabilities': [{'kind': 'one'}, {'kind': 'two'}]},
        {'imgt_region': 'Framework', 'liabilities': [{'kind': 'one'}]},
        {'imgt_region': 'Unassigned', 'liabilities': [{'kind': 'one'}]},
        {'imgt_region': 'CDR2', 'liabilities': []},
    ]}})
    assert row['cdr_liabilities'] == 1
    assert row['framework_liabilities'] == 1
    assert row['unassigned_liabilities'] == 1
