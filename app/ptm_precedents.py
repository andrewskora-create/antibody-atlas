"""Curated observed PTMs in approved antibody molecules; no sequence is sent out.

Query associations require a verified public reference site and aligned chemistry motif.
Literature numbering is preserved. Approval is not approval of a modified species.
This small catalog is deliberately not a prevalence or acceptable-limit dataset.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import json
import hashlib
from pathlib import Path
import re

from .numbering import number_chain
from typing import Any

CATALOG_VERSION = "2026-10-04-expanded-v3"
CATALOG_NOTE = (
    "Selected experimentally documented cases in approved antibody molecules, with framework cases first. "
    "Some studies used biosimilars or stressed material. A precedent does not establish an acceptable PTM level, "
    "prove the same liability in another antibody, or change its screening concern. This general library "
    "is not a query result; sequence-specific matching is reported separately with canonical site and motif checks. "
    "Fc/constant-domain sites are not variable-framework precedents. This is a curated subset, not a complete "
    "Phase III/approved-drug survey; absence of a match is not evidence of absence."
)

EVIDENCE_PATH = Path(__file__).resolve().parent.parent / "data" / "therapeutic_ptm_evidence.json"
PRECEDENTS: tuple[dict[str, Any], ...] = tuple(json.loads(EVIDENCE_PATH.read_text())["cases"])



REFERENCE_PATH = Path(__file__).resolve().parent.parent / "data" / "therapeutic_ptm_sequences.json"
SEQUENCE_NOTE = (
    "Only input-specific matches are shown: the query and documented therapeutic site must share heavy/light "
    "chain family and light-chain subtype, canonical IMGT position, observed residue and the complete aligned immediate chemistry motif. "
    "Displayed query and reference regions use the selected CDR convention; publication regions retain their original definition. "
    "Chain identity includes every position in the union of the two numbered variable domains; gaps count as mismatches, and unnumbered terminal extensions are excluded. "
    "The aligned local window extends five numbered positions either side of the site. Differences remain explicit. "
    "A homologous-site match does not establish equal conformation, PTM rate, clinical consequence, or acceptable PTM level, "
    "and never changes screening priority. No match means no match in this curated evidence catalog. All matching is local."
)


def matching_precedents(kind: str, region: str, residue: str | None = None) -> list[dict[str, Any]]:
    """Deprecated broad matching deliberately returns no query associations."""
    return []


def precedent_catalog() -> dict[str, Any]:
    """General library only; never use as a query-specific result."""
    return {"version": CATALOG_VERSION, "note": CATALOG_NOTE, "cases": deepcopy(list(PRECEDENTS))}


@lru_cache(maxsize=1)
def _reference_data() -> dict[str, Any]:
    from .reference_data import load
    return load('therapeutic_ptm_sequences')


@lru_cache(maxsize=32)
def _number_reference(sequence: str, scheme: str) -> dict[str, Any]:
    return number_chain(sequence, scheme)


def _canonical(residue: dict[str, Any]) -> str:
    # H/K/L denote chain type, not distinct numeric alignment coordinates.
    return re.sub(r"^[HKL]", "", residue["number"])


#: IMGT positions whose insertion letters run BACKWARDS towards the base.
#:
#: IMGT grows each CDR outward from both anchors towards its middle, so the
#: lower of the two central positions counts up and the higher counts down.
#: In N-to-C order a long CDR2 reads 60, 60A, 60B, 61C, 61B, 61A, 61, and a
#: long CDR3 reads 111, 111A, 111B, 112B, 112A, 112. The intuition that 112
#: precedes 112A is the Kabat convention, where insertions always follow their
#: base; IMGT inverts that on the descending side only.
#:
#: Verified against 600 numbered therapeutic chains: position 61 was observed
#: emitting 61C, 61B, 61A in that order, exactly as 112 does, while 60, 81 and
#: 111 ascend normally. 33 is the CDR1 counterpart of 61 and 112 and is
#: included by symmetry; it did not occur in that sample, so it is the one
#: entry here resting on the rule rather than on an observation.
DESCENDING_INSERTION_POSITIONS = frozenset({33, 61, 112})


def _position_order(position: str) -> tuple[int, int]:
    """Sort key putting IMGT positions in N-to-C order.

    Plain string sorting is wrong twice over: it orders 112 before 112A, and
    it puts 36 after 114 because it compares digits as text. Both errors are
    silent -- the set of positions is right and only the order is wrong -- so
    they surface as a reader misreading which residues were involved.
    """
    match = re.fullmatch(r"(\d+)([A-Z]*)", position)
    if not match:
        raise ValueError("Unexpected canonical antibody position.")
    base, suffix = int(match[1]), match[2]
    offset = 0
    for letter in suffix:
        offset = offset * 26 + ord(letter) - ord('A') + 1
    return base, -offset if base in DESCENDING_INSERTION_POSITIONS else offset


def position_order(position: str) -> tuple[str, int, int]:
    """Sort key for a position that may carry its chain prefix, e.g. H112A.

    The same ordering as _position_order, lifted to the labels the structure
    and paratope paths actually carry. Heavy sorts before light; K and L are
    folded together because they are the same coordinate system.
    """
    match = re.fullmatch(r"([HKL]?)(\d+[A-Z]*)", str(position).strip())
    if not match:
        raise ValueError(f"Unexpected antibody position {position!r}.")
    chain = "H" if match[1] == "H" else "L" if match[1] else ""
    base, offset = _position_order(match[2])
    return chain, base, offset


def _family(numbering: dict[str, Any]) -> str:
    return 'H' if numbering['chain_type'] == 'H' else 'L'


def _motif_length(chemistry: str) -> int:
    return 3 if chemistry == 'N-linked glycosylation sequon' else 1 if chemistry in {'Oxidation-susceptible residue', 'N-terminal pyroglutamate formation'} else 2


def sequence_precedents(chains: dict[str, Any], scheme: str = 'imgt') -> dict[str, Any]:
    """Match observed literature sites through canonical IMGT, without risk transfer."""
    if scheme not in {'imgt', 'kabat'}:
        raise ValueError('Numbering scheme must be IMGT or Kabat.')
    references = _reference_data()
    metadata = references if "antibodies" in references else {}
    # Also accept a metadata wrapper used by future catalog releases.
    references = references.get('antibodies', references)
    indexed = {site['precedent_id']: (ref, site)
               for ref in references.values() if isinstance(ref, dict)
               for site in ref.get('validated_sites', []) if site.get('validated')}
    queries = {}
    skipped_chains = []
    for label, value in chains.items():
        sequence = value['sequence'] if isinstance(value, dict) else value
        try:
            canonical = number_chain(sequence, 'imgt')
            selected = value.get('numbering') if isinstance(value, dict) else None
            if not selected or selected.get('scheme') != scheme:
                selected = canonical if scheme == 'imgt' else number_chain(sequence, scheme)
        except ValueError:
            skipped_chains.append({'chain': label, 'reason': 'Antibody numbering unavailable; this chain was not compared.'})
            continue
        if label not in {'H', 'L'} or _family(canonical) != label:
            continue
        queries[label] = (sequence, canonical, selected)
    cases, searched = [], 0
    for precedent in PRECEDENTS:
        entry = indexed.get(precedent['id'])
        if not entry:
            continue
        reference, site = entry
        chain = site['chain']
        ref_sequence = reference['chains'][chain]
        position = site['sequence_position']
        length = _motif_length(precedent['chemistry'])
        motif = ref_sequence[position - 1:position - 1 + length]
        # Fail closed if the public reference no longer agrees with documented chemistry.
        if len(motif) != length or motif[0] != precedent['residue'] or site.get('observed_motif', motif) != motif:
            continue
        if site.get('expected_motif') and site['expected_motif'] != motif:
            continue
        ref_numbered = _number_reference(ref_sequence, 'imgt')
        if _family(ref_numbered) != chain:
            continue
        ref_by_seq = {r['sequence_position']:r for r in ref_numbered['residues']}
        motif_rows = [ref_by_seq.get(p) for p in range(position, position + length)]
        if any(row is None for row in motif_rows):
            continue
        searched += 1
        if chain not in queries:
            continue
        query_sequence, query_numbered, query_selected = queries[chain]
        if query_numbered['chain_type'] != ref_numbered['chain_type']:
            continue
        ref_map = {_canonical(r):r for r in ref_numbered['residues']}
        query_map = {_canonical(r):r for r in query_numbered['residues']}
        motif_keys = [_canonical(r) for r in motif_rows]
        query_rows = [query_map.get(key) for key in motif_keys]
        if any(row is None for row in query_rows):
            continue
        query_position = query_rows[0]['sequence_position']
        if precedent['chemistry'] == 'N-terminal pyroglutamate formation' and (position != 1 or query_position != 1):
            continue
        if [r['sequence_position'] for r in query_rows] != list(range(query_position, query_position + length)):
            continue
        if ''.join(r['aa'] for r in query_rows) != motif:
            continue
        positions = sorted(set(ref_map) | set(query_map), key=_position_order)
        identity_matches = sum(key in ref_map and key in query_map and ref_map[key]['aa'] == query_map[key]['aa'] for key in positions)
        anchor = positions.index(motif_keys[0])
        local_keys = positions[max(0, anchor - 5):anchor + 6]
        query_context = ''.join(query_map[key]['aa'] if key in query_map else '-' for key in local_keys)
        reference_context = ''.join(ref_map[key]['aa'] if key in ref_map else '-' for key in local_keys)
        local_matches = sum(a == b and a != '-' for a,b in zip(query_context, reference_context))
        selected_query = {r['sequence_position']:r for r in query_selected['residues']}[query_position]
        ref_selected = ref_numbered if scheme == 'imgt' else _number_reference(ref_sequence, scheme)
        selected_reference = {r['sequence_position']:r for r in ref_selected['residues']}[position]
        match_type = ('Exact reference chain sequence' if query_sequence == ref_sequence and identity_matches == len(positions)
                      else 'Exact numbered variable domain; terminal sequence differs' if identity_matches == len(positions)
                      else 'Aligned site and motif; exact local window' if local_matches == len(local_keys)
                      else 'Aligned site and motif; differing neighborhood')
        match = {
            'query_chain': chain, 'query_sequence_position': query_position,
            'query_number': selected_query['number'], 'query_region': selected_query['region'],
            'query_imgt_number': query_rows[0]['number'], 'query_scheme': scheme,
            'reference_imgt_number': motif_rows[0]['number'], 'reference_position': position,
            'reference_number': selected_reference['number'], 'reference_region': selected_reference['region'],
            'reference_scheme': scheme, 'reference_chain_type': ref_numbered['chain_type'],
            'query_chain_type': query_numbered['chain_type'],
            'identity_percent': round(identity_matches / len(positions) * 100, 1),
            'identity_matches': identity_matches, 'identity_denominator': len(positions),
            'query_context': query_context, 'reference_context': reference_context,
            'local_matches': local_matches, 'local_positions': len(local_keys),
            'local_numbered_positions': local_keys, 'context_anchor_index': local_keys.index(motif_keys[0]),
            'matched_motif': motif, 'motif_imgt_positions': motif_keys, 'match_type': match_type,
            'reference_source_url': reference.get('source_url'),
        }
        case = deepcopy(precedent)
        case.update(matches=[match], sequence_source_url=reference.get('source_url') or metadata.get('source_url'),
                    reference_sequence_sha256=hashlib.sha256(ref_sequence.encode()).hexdigest(), association='Canonical aligned site and identical chemistry motif; no transfer of PTM risk')
        cases.append(case)
    return {'version': CATALOG_VERSION, 'note': SEQUENCE_NOTE, 'cases': cases,
            'searched_cases': searched, 'unmatched_count': searched - len(cases),
            'unavailable_cases': len(PRECEDENTS) - searched, 'skipped_chains': skipped_chains,
            'reference_version': metadata.get('version'), 'sequence_source_url': metadata.get('source_url'),
            'reference_sha256': hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest(),
            'upstream_source_sha256': metadata.get('source_sha256')}
