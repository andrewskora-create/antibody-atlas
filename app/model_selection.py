"""Select a fully validated terminal-repaired model without losing lineage."""
from copy import deepcopy


def select_propermab_structure(report, result):
    if result.get('structure_selection')!='terminal_repaired':
        return False
    repair=result.get('terminal_repair') or {}
    structure=result.get('selected_structure_analysis') or {}
    original=result.get('original_input_pdb_sha256')
    selected=result.get('input_pdb_sha256')
    provenance=report['provenance']
    expected=provenance.get('original_structure_pdb_sha256') or provenance.get('structure_pdb_sha256')
    coverage=structure.get('mapping_coverage_percent') or {}
    if (repair.get('status')!='complete' or not original or original!=expected or not selected
        or selected!=repair.get('repaired_pdb_sha256') or set(coverage)!=set(report['chains'])
        or any(value!=100 for value in coverage.values())
        or structure.get('scheme')!=report.get('scheme') or structure.get('ph')!=report.get('ph')):
        raise ValueError('Repaired structure provenance, full-input coverage or analysis settings do not match this record.')
    if provenance.get('structure_selection')!='terminal_repaired':
        report['original_structure_analysis']=deepcopy(report.get('structure_analysis'))
    report['structure_analysis']=deepcopy(structure)
    provenance.update(original_structure_pdb_sha256=original,structure_pdb_sha256=selected,
                      structure_selection='terminal_repaired',
                      structure_preparation='Full-input terminal reconstruction; identity and peptide geometry checked. Rebuilt positions remain model hypotheses.')
    return True
