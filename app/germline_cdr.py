"""Exact full IMGT CDR1/2 identity to local functional human V alleles.

Amino-acid identity cannot establish ancestry or absence of historical mutation.
The functional-only OGRDB library is the same versioned reference used for humanness.
"""
from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import re

REFERENCE_PATH = Path(__file__).resolve().parents[1] / 'data' / 'human_v_ogrdb_numbered.json'
LIMITATION = ('Exact amino-acid identity does not establish germline ancestry, prove that no mutations occurred, '
              'or demonstrate PTM occurrence, tolerance, or clinical safety. CDR identity is evaluated using IMGT boundaries.')


@lru_cache(maxsize=1)
def _references():
    return json.loads(REFERENCE_PATH.read_text())


def _key(row):
    return re.sub(r'^[HKL]', '', row['number'])


def _base(key):
    match = re.match(r'^([0-9]+)', key)
    return int(match[1]) if match else -1


def annotate_germline_cdr(site, querycanonical, queryselected):
    """Return annotation; site.position is one-based input index, site.motif optional.

    Numbering dictionaries have the existing number_chain schema. querycanonical
    must be IMGT; queryselected is the user's selected numbering convention.
    Complete keyed CDR dictionaries must equal (union positions, no missing/extra
    reference residues), rather than only matching a short PTM motif.
    """
    data = _references()
    locus = {'H': 'IGH', 'K': 'IGK', 'L': 'IGL'}.get(querycanonical.get('chain_type'))
    out = {'status': 'unassessed', 'summary': '', 'region': None, 'alleles': [],
           'reference_source': data['source'], 'reference_url': data['source_url'],
           'release': data['release_versions'].get(locus), 'sequence': '',
           'positions': [], 'motif_matched': False, 'limitation': LIMITATION}
    if not locus or querycanonical.get('scheme', 'imgt').lower() != 'imgt':
        out['summary'] = 'Human germline CDR identity unassessed: canonical IMGT numbering or chain subtype unavailable.'
        return out
    position = site.get('position', site.get('sequence_position'))
    anchor = next((r for r in querycanonical['residues'] if r['sequence_position'] == position), None)
    selected = next((r for r in queryselected['residues'] if r['sequence_position'] == position), None)
    if anchor is None or selected is None:
        out['summary'] = 'Human germline CDR identity unassessed: site could not be mapped.'
        return out
    b = _base(_key(anchor))
    region = 'CDR1' if 27 <= b <= 38 else 'CDR2' if 56 <= b <= 65 else 'CDR3' if 105 <= b <= 117 else 'Framework'
    out['region'] = region
    if region == 'CDR3' or selected.get('region') == 'CDR3':
        out['summary'] = 'Full CDR3 germline identity unassessed: a V-only reference cannot establish V/(D)/J and junctional contributions.'
        return out
    if region not in {'CDR1', 'CDR2'} or selected.get('region') not in {'CDR1', 'CDR2'}:
        out['summary'] = ('No germline CDR claim: the site is outside IMGT CDR1/2 or outside the selected CDR convention.')
        return out
    low, high = (27, 38) if region == 'CDR1' else (56, 65)
    rows = [r for r in querycanonical['residues'] if low <= _base(_key(r)) <= high]
    query = {_key(r): r['aa'] for r in rows}
    out['sequence'] = ''.join(r['aa'] for r in rows)
    out['positions'] = [r['number'] for r in rows]
    motif = site.get('motif', '')
    bypos = {r['sequence_position']: r for r in querycanonical['residues']}
    motif_rows = [bypos.get(position + offset) for offset in range(len(motif))]
    motif_in_cdr = bool(motif) and all(r is not None and _key(r) in query for r in motif_rows)
    motif_valid = motif_in_cdr and ''.join(r['aa'] for r in motif_rows) == motif
    for ref in data['references']:
        if ref['locus'] != locus or ref.get('functional') is False:
            continue
        candidate = {k: aa for k, aa in ref['positions'].items() if low <= _base(k) <= high}
        if query and query == candidate:
            out['alleles'].append(ref['name'])
    out['alleles'] = sorted(set(out['alleles']))
    out['status'] = 'exact_match' if out['alleles'] else 'no_exact_match'
    out['motif_matched'] = bool(out['alleles']) and bool(motif_valid)
    out['summary'] = (f'Exact IMGT {region} amino-acid match to functional human V allele(s).'
                      if out['alleles'] else f'No exact full IMGT {region} amino-acid match in the local functional human V-allele library.')
    if motif and not motif_in_cdr:
        out['limitation'] += ' The complete PTM motif crosses the IMGT CDR boundary; the CDR match does not establish a complete germline motif match.'
    return out
