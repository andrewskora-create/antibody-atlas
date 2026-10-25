"""Accept a nucleotide sequence and translate it to the protein to screen.

This exists because of a specific trap. A, C, G and T are all valid
amino-acid one-letter codes, so a DNA sequence passes every "is this a
protein" check unchallenged and arrives at the numbering step as a long
nonsense peptide. The user then gets "ANARCII could not recognize a complete
antibody variable domain", which is true but tells them nothing about what
actually went wrong.

Reading frame is not guessed. Every frame is translated and the numbering
model decides: the correct frame is the one that contains a recognisable
antibody variable domain. That also handles a signal peptide or vector
sequence at either end without needing to find the start codon, because the
numbering reports where the domain begins.
"""
from __future__ import annotations

import re
from typing import Any

CODONS = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMPLEMENT = str.maketrans('ACGTUNRYSWKMBDHV', 'TGCAANYRSWMKVHDB')
# IUPAC ambiguity codes appear routinely in sequencing output. They are not
# translatable, but a handful of them must not stop a sequence being
# recognised as DNA -- that was the old failure: one R in a read and the whole
# thing was treated as protein, giving a numbering error that named no cause.
UNAMBIGUOUS = set('ACGTU')
AMBIGUITY = set('NRYSWKMBDHV')
NUCLEOTIDES = UNAMBIGUOUS | AMBIGUITY

# Below this length, an all-ACGT string is as likely to be a short peptide as
# a gene fragment, and an antibody variable domain needs ~330 nt anyway.
MIN_NUCLEOTIDE_LENGTH = 60
# Allow a few ambiguity codes without abandoning the call.
MIN_NUCLEOTIDE_FRACTION = 0.9


def _strip(raw: str) -> str:
    lines = [l.strip() for l in str(raw).strip().splitlines() if l.strip()]
    if lines and lines[0].startswith('>'):
        lines = lines[1:]
    return re.sub(r'[\s0-9]+', '', ''.join(lines)).upper()


def looks_like_nucleotide(raw: str) -> bool:
    """Is this DNA or RNA rather than protein?

    Deliberately conservative. A, C, G and T are amino acids too, so this
    only fires on a sequence long enough that an all-ACGT protein would be
    an extraordinary coincidence.
    """
    sequence = _strip(raw)
    if len(sequence) < MIN_NUCLEOTIDE_LENGTH:
        return False
    # Judge on unambiguous bases only: 'N' and friends are also amino acids,
    # so counting them towards the nucleotide score would make a protein with
    # a few asparagines look like DNA.
    hits = sum(c in UNAMBIGUOUS for c in sequence)
    return hits / len(sequence) >= MIN_NUCLEOTIDE_FRACTION


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def translate(sequence: str, frame: int = 0) -> str:
    """Translate one frame. Unknown codons become X, stops become *."""
    coding = sequence[frame:]
    coding = coding[:len(coding) - len(coding) % 3]
    return ''.join(CODONS.get(coding[i:i + 3], 'X') for i in range(0, len(coding), 3))


def candidate_frames(sequence: str, try_reverse: bool = True) -> list[dict[str, Any]]:
    """All six frames, or three if the reverse strand is excluded."""
    out = []
    strands = [('+', sequence)] + ([('-', reverse_complement(sequence))] if try_reverse else [])
    for strand, strand_sequence in strands:
        for frame in range(3):
            protein = translate(strand_sequence, frame)
            out.append({'strand': strand, 'frame': frame + 1,
                        'label': f'{strand}{frame + 1}', 'protein': protein,
                        'stops': protein.count('*'), 'length': len(protein)})
    return out


def to_protein(raw: str, try_reverse: bool = True) -> dict[str, Any]:
    """Translate a nucleotide sequence, choosing the frame by numbering it.

    Returns the chosen frame's protein plus everything needed to show the user
    what happened, because a silent translation is worse than no translation:
    they need to be able to see which frame was used and check it.
    """
    sequence = _strip(raw)
    unknown = sorted(set(sequence) - NUCLEOTIDES)
    if unknown:
        raise ValueError(
            f"Not a nucleotide sequence: unexpected symbols {', '.join(unknown)}.")
    if len(sequence) < MIN_NUCLEOTIDE_LENGTH:
        raise ValueError(
            f"A nucleotide sequence needs at least {MIN_NUCLEOTIDE_LENGTH} bases; "
            f"this has {len(sequence)}.")

    frames = candidate_frames(sequence, try_reverse)

    # The numbering model is the frame oracle. A frame that contains a
    # recognisable antibody V domain is the right one; no other heuristic
    # (longest ORF, fewest stops, first ATG) is as reliable, and all of them
    # fail on a construct that carries a signal peptide or vector flanks.
    from .numbering import number_chain
    chosen = None
    for candidate in frames:
        # A stop inside the domain rules the frame out, but a stop after it
        # (the natural end of the construct) does not, so numbering decides.
        if not candidate['protein']:
            continue
        for segment in candidate['protein'].split('*'):
            if len(segment) < 90:
                continue
            try:
                numbering = number_chain(segment, 'imgt')
            except Exception:
                continue
            # Trim to the variable domain itself. A construct usually carries a
            # signal peptide, and a secretion leader is strongly hydrophobic --
            # keeping it would inflate GRAVY, shift pI and invent hydrophobic
            # windows that are not in the folded protein. The numbering already
            # located the domain, so use its bounds.
            start = int(numbering.get('domain_start') or 1)
            end = int(numbering.get('domain_end') or len(segment))
            chosen = {**candidate, 'domain': segment[start - 1:end],
                      'context': segment,
                      'leader': segment[:start - 1],
                      'trailer': segment[end:],
                      'chain_type': numbering.get('chain_type')}
            break
        if chosen:
            break

    if chosen is None:
        summary = '; '.join(f"{c['label']}: {c['length']} aa, {c['stops']} stop(s)"
                            for c in frames)
        raise ValueError(
            'Translated this nucleotide sequence in every reading frame, but none '
            'contained a recognisable antibody variable domain. Check the sequence is '
            f'an antibody V region and is not truncated. Frames tried -- {summary}.')

    protein = chosen['domain']
    leader, trailer = chosen['leader'], chosen['trailer']
    trimmed = []
    if leader:
        trimmed.append(f"{len(leader)} residue{'s' if len(leader) != 1 else ''} before it"
                       + (' (signal peptide or vector sequence)' if len(leader) > 8 else ''))
    if trailer:
        trimmed.append(f"{len(trailer)} after it")
    return {
        'translated': True,
        'protein': protein,
        'nucleotide_length': len(sequence),
        'frame': chosen['label'],
        'strand': chosen['strand'],
        'reverse_complement': chosen['strand'] == '-',
        'chain_type': chosen['chain_type'],
        'full_frame_protein': chosen['protein'],
        'translated_context': chosen['context'],
        'leader': leader,
        'trailer': trailer,
        'stops_in_frame': chosen['stops'],
        'note': (f"Translated {len(sequence)} bases in frame {chosen['label']}"
                 + (' on the reverse complement' if chosen['strand'] == '-' else '')
                 + f". Screening the {len(protein)}-residue variable domain"
                 + (f", discarding {' and '.join(trimmed)}." if trimmed else '.')
                 + " The reading frame was chosen by numbering each frame and keeping the"
                   " one containing an antibody variable domain."),
        'caveat': ('Only the variable domain is screened. A signal peptide is strongly '
                   'hydrophobic and would distort GRAVY, pI and the hydrophobic-window '
                   'count if it were left in, so anything outside the numbered domain is '
                   'removed. Check the translated sequence below before trusting a result '
                   'derived from it.'),
    }


def resolve(raw: str, try_reverse: bool = True) -> tuple[str, dict[str, Any] | None]:
    """Return (sequence to screen, translation record or None if already protein)."""
    if not looks_like_nucleotide(raw):
        return raw, None
    record = to_protein(raw, try_reverse)
    return record['protein'], record


def describe(raw: str) -> str:
    """Explain what a sequence looks like, for an error message.

    "ANARCII could not recognize a complete antibody variable domain" is true
    and useless. This says what was actually received so the cause is visible
    without the user having to share the sequence with anyone.
    """
    sequence = _strip(raw)
    if not sequence:
        return 'The field was empty after whitespace was removed.'

    length = len(sequence)
    unambiguous = sum(c in UNAMBIGUOUS for c in sequence)
    fraction = unambiguous / length
    protein_only = sorted(set(sequence) - NUCLEOTIDES)

    parts = [f'Received {length} characters.']

    if fraction >= MIN_NUCLEOTIDE_FRACTION and length >= MIN_NUCLEOTIDE_LENGTH:
        parts.append('It was read as a nucleotide sequence and translated, but no reading '
                     'frame contained a recognisable variable domain.')
    elif fraction >= MIN_NUCLEOTIDE_FRACTION and length < MIN_NUCLEOTIDE_LENGTH:
        parts.append(f'It looks like DNA but is only {length} bases, below the '
                     f'{MIN_NUCLEOTIDE_LENGTH} needed to attempt translation. A variable '
                     'domain needs roughly 330 bases.')
    elif fraction >= 0.7:
        odd = ''.join(protein_only[:8])
        parts.append(f'{unambiguous} of {length} characters are A/C/G/T/U, so this may be DNA '
                     f'with {length - unambiguous} non-base characters in it'
                     + (f' ({odd})' if odd else '')
                     + '. It was treated as protein. Remove the stray characters if it is DNA.')
    else:
        parts.append(f'It was read as a protein sequence of {length} residues.')
        if length < 95:
            parts.append('A complete variable domain is about 110 to 130 residues, so this is '
                         'probably truncated -- a CDR alone or a partial domain cannot be '
                         'numbered.')
        elif length > 400:
            parts.append('That is longer than a variable domain. A full-length chain or a '
                         'scFv may need splitting; paste VH and VL separately.')
        if protein_only:
            unexpected = [c for c in protein_only if c not in 'ACDEFGHIKLMNPQRSTVWY']
            if unexpected:
                parts.append('Unsupported symbols: ' + ', '.join(unexpected[:8]) + '.')

    return ' '.join(parts)
