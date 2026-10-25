"""Read a panel of antibodies from one multi-FASTA file, and refuse bad input early.

A batch of 100 takes hours. Discovering on molecule 80 that row 40 had a stray
'X' wastes most of a day, so everything that can be checked from the file is
checked before any work starts, and every problem is reported at once with the
line it came from rather than one per run.

Pairing is by name. These header conventions are all accepted:

    >MAB01|H            >MAB01|L
    >MAB01_VH           >MAB01_VL
    >MAB01_HC           >MAB01_LC
    >MAB01 heavy        >MAB01 light
    >MAB01:H            >MAB01:K
    >MAB01 H            >MAB01 L

The chain word may be anywhere in the header, not only at the end, so a label
written the way people actually write them also works:

    MAB01 variable heavy chain coding sequence
    MAB01 variable light chain coding sequence

The leading '>' is optional. A plain name line followed by its sequence is
read the same way, because that is how sequences arrive when pasted out of a
document rather than exported as FASTA.

Nucleotide input is translated before anything is checked, so the residue and
length checks apply to the protein that will actually be screened rather than
to the bases.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .analysis import AA

HEAVY_TOKENS = {'h', 'hc', 'vh', 'heavy', 'heavychain', 'heavy_chain'}
LIGHT_TOKENS = {'l', 'lc', 'vl', 'light', 'lightchain', 'light_chain', 'k', 'kappa', 'lambda'}
SEPARATORS = re.compile(r'[|:/\s_,;-]+')
FV_MIN, FV_MAX = 90, 150


def validate_name(name: str) -> None:
    """Construct names become local output filenames; keep them in one directory."""
    if (not name.strip() or name in {'.', '..'} or
            any(c in name for c in '/\\') or
            any(ord(c) < 32 or ord(c) == 127 for c in name) or
            len(name.encode('utf-8')) > 180):
        raise ValueError('use a molecule name under 180 UTF-8 bytes without path separators or control characters.')


@dataclass
class Construct:
    name: str
    chains: dict[str, str] = field(default_factory=dict)
    lines: dict[str, int] = field(default_factory=dict)


@dataclass
class ParsedPanel:
    constructs: dict[str, dict[str, str]]
    problems: list[dict]
    warnings: list[dict]

    @property
    def ok(self) -> bool:
        return not self.problems

    def report(self) -> str:
        out = [f"{len(self.constructs)} construct(s) read."]
        for kind, items in (('PROBLEM', self.problems), ('warning', self.warnings)):
            for item in items:
                where = f" (line {item['line']})" if item.get('line') else ''
                out.append(f"  {kind}: {item['name']}{where}: {item['detail']}")
        return "\n".join(out)


# Words that carry no information about which chain this is, so they can be
# dropped from the name without making two constructs collide.
NOISE_TOKENS = {'variable', 'chain', 'coding', 'sequence', 'domain', 'region',
                'protein', 'aa', 'nt', 'dna', 'cds', 'seq'}


def _chain_token(header: str) -> tuple[str, str | None]:
    """Split a header into a construct name and a chain assignment.

    The chain word is looked for anywhere in the header, not only in the last
    field. A header like "MAB01 variable heavy chain coding sequence" ends in
    "sequence", and insisting the chain be last rejected it -- which is how
    people actually label sequences pasted out of a document.
    """
    parts = [p for p in SEPARATORS.split(header.strip()) if p]
    if not parts:
        return header.strip(), None

    chain = None
    kept: list[str] = []
    for part in parts:
        lowered = part.lower()
        # First chain word wins; a later one is treated as part of the name so
        # that "heavy_vs_light_comparison" cannot silently become two chains.
        if chain is None and lowered in HEAVY_TOKENS:
            chain = 'H'
            continue
        if chain is None and lowered in LIGHT_TOKENS:
            chain = 'L'
            continue
        if lowered in NOISE_TOKENS:
            continue
        kept.append(part)

    if chain is None:
        return header.strip(), None
    name = SEPARATORS.sub('_', ' '.join(kept)).strip('_')
    return (name or header.strip()), chain


SEQUENCE_LINE = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYUXBZacdefghiklmnpqrstvwyuxbz\s0-9*-]+$')


def _is_sequence_line(line: str) -> bool:
    """Does this look like sequence rather than a label?

    A label almost always contains a character no sequence does, or a word
    with a vowel pattern no sequence has. Requiring a minimum length and a
    high residue-letter fraction separates them without needing a '>'.
    """
    stripped = re.sub(r'[\s0-9*-]+', '', line)
    if len(stripped) < 12:
        return False
    return bool(SEQUENCE_LINE.match(line))


def parse_fasta(text: str) -> ParsedPanel:
    records: list[tuple[str, str, int]] = []
    name, buffer, start = None, [], 0
    explicit_headers = any(line.lstrip().startswith('>') for line in text.splitlines())
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip().rstrip(':')
        if not line:
            continue
        # A leading '>' marks a header, but it is optional: a label line
        # followed by its sequence is the shape sequences take when they are
        # copied out of a document, and refusing it helps nobody.
        # Once a FASTA header is present, every following non-header line is
        # sequence, including a short final wrapped line. The heuristic is
        # only for the optional plain-label format.
        short_continuation = (name is not None and bool(SEQUENCE_LINE.fullmatch(line))
                              and not re.search(r'\s', line))
        plain_chain_label = bool(re.search(
            r'\b(?:heavy|light|vh|vl|hc|lc|kappa|lambda)\b', line, re.I))
        header = line.startswith('>') or (not explicit_headers and
                                          (plain_chain_label or
                                           (not _is_sequence_line(line) and
                                            not short_continuation)))
        if header:
            if name is not None:
                records.append((name, ''.join(buffer), start))
            name = line[1:].strip() if line.startswith('>') else line
            buffer, start = [], number
        else:
            buffer.append(re.sub(r'\s+', '', line).upper())
    if name is not None:
        records.append((name, ''.join(buffer), start))

    problems: list[dict] = []
    warnings: list[dict] = []
    found: dict[str, Construct] = {}

    for header, sequence, line in records:
        construct, chain = _chain_token(header)
        if chain is None:
            problems.append({'name': header, 'line': line,
                             'detail': "no chain could be identified. End the header with a chain "
                                       "token, for example '>Name|H' or '>Name_VL'."})
            continue
        try:
            validate_name(construct)
        except ValueError as exc:
            problems.append({'name': header, 'line': line, 'detail': str(exc)})
            continue
        if not sequence:
            problems.append({'name': header, 'line': line, 'detail': 'record has no sequence.'})
            continue
        # Translate first. A, C, G and T are amino acids, so DNA passes the
        # residue check untouched and then fails the length check as a
        # 'protein' three times too long -- a confusing way to learn that the
        # input was nucleotide.
        from .translation import looks_like_nucleotide, to_protein
        if looks_like_nucleotide(sequence):
            try:
                sequence = to_protein(sequence)['protein']
            except ValueError as exc:
                problems.append({'name': header, 'line': line,
                                 'detail': f'nucleotide sequence could not be translated: {exc}'})
                continue
        bad = sorted(set(sequence) - AA)
        if bad:
            problems.append({'name': header, 'line': line,
                             'detail': f"unsupported residue(s) {', '.join(bad)}; use the 20 standard amino acids."})
            continue
        entry = found.setdefault(construct, Construct(construct))
        if chain in entry.chains:
            problems.append({'name': header, 'line': line,
                             'detail': f"a {chain} chain for '{construct}' was already read at line "
                                       f"{entry.lines[chain]}. Construct names must be unique per chain."})
            continue
        entry.chains[chain] = sequence
        entry.lines[chain] = line

    constructs: dict[str, dict[str, str]] = {}
    for name_, entry in found.items():
        if 'H' not in entry.chains:
            # A light chain alone is screenable -- the single-molecule page
            # accepts one, and so should this. Paired-chain descriptors are
            # simply unavailable, which the warning says.
            warnings.append({'name': name_, 'line': entry.lines.get('L'),
                             'detail': 'no heavy chain; it will be screened as a light-only '
                                       'domain and paired-chain descriptors will be unavailable.'})
        if 'L' not in entry.chains:
            warnings.append({'name': name_, 'line': entry.lines.get('H'),
                             'detail': 'no light chain; it will be screened as a heavy-only domain '
                                       'and paired-chain descriptors will be unavailable.'})
        for chain, sequence in entry.chains.items():
            if not FV_MIN <= len(sequence) <= FV_MAX:
                warnings.append({'name': name_, 'line': entry.lines[chain],
                                 'detail': f"{chain} chain is {len(sequence)} residues; a variable domain is "
                                           f"usually {FV_MIN}-{FV_MAX}. Structure prediction will be skipped "
                                           f"if this is a full-length chain."})
        constructs[name_] = dict(entry.chains)

    if not records:
        problems.append({'name': '(file)', 'line': None, 'detail': 'no sequence records found.'})
    return ParsedPanel(constructs=constructs, problems=problems, warnings=warnings)


def parse_csv(text: str) -> ParsedPanel:
    """name,heavy,light - the other format people actually have."""
    import csv
    import io
    text = text.lstrip('\ufeff\r\n')
    first = next((line for line in text.splitlines() if line.strip()), '')
    delimiter = '\t' if '\t' in first and ',' not in first else ','
    rows = list(csv.DictReader(io.StringIO(text), delimiter=delimiter))
    problems, warnings, constructs = [], [], {}
    if not rows:
        return ParsedPanel({}, [{'name': '(file)', 'line': None, 'detail': 'no rows found.'}], [])
    fields = {(f or '').strip().lower(): f for f in (rows[0].keys() or [])}
    name_key = next((fields[k] for k in ('name', 'construct', 'id', 'antibody') if k in fields), None)
    heavy_key = next((fields[k] for k in ('heavy', 'h', 'vh', 'heavy_chain') if k in fields), None)
    light_key = next((fields[k] for k in ('light', 'l', 'vl', 'light_chain') if k in fields), None)
    if not name_key or not (heavy_key or light_key):
        return ParsedPanel({}, [{'name': '(file)', 'line': 1,
                                 'detail': 'need a name column and at least one heavy or light column.'}], [])
    for index, row in enumerate(rows, start=2):
        name = (row.get(name_key) or '').strip()
        if not name:
            problems.append({'name': '(blank)', 'line': index, 'detail': 'missing name.'}); continue
        try:
            validate_name(name)
        except ValueError as exc:
            problems.append({'name': name, 'line': index, 'detail': str(exc)})
            continue
        if None in row:
            problems.append({'name': name, 'line': index,
                             'detail': 'more fields than header columns; quote names containing delimiters.'})
            continue
        chains = {}
        for chain, key in (('H', heavy_key), ('L', light_key)):
            if not key:
                continue
            sequence = re.sub(r'\s+', '', (row.get(key) or '')).upper()
            if not sequence:
                continue
            from .translation import looks_like_nucleotide, to_protein
            if looks_like_nucleotide(sequence):
                try:
                    sequence = to_protein(sequence)['protein']
                except ValueError as exc:
                    problems.append({'name': name, 'line': index,
                                     'detail': f'{chain}: nucleotide sequence could not be translated: {exc}'})
                    continue
            bad = sorted(set(sequence) - AA)
            if bad:
                problems.append({'name': name, 'line': index,
                                 'detail': f"{chain}: unsupported residue(s) {', '.join(bad)}."})
                continue
            chains[chain] = sequence
        if not chains:
            problems.append({'name': name, 'line': index, 'detail': 'no usable heavy or light chain.'}); continue
        if name in constructs:
            problems.append({'name': name, 'line': index, 'detail': 'duplicate name.'}); continue
        constructs[name] = chains
        for absent in {'H', 'L'} - chains.keys():
            warnings.append({'name': name, 'line': index,
                             'detail': f'no {absent} chain; paired-chain descriptors will be unavailable.'})
    return ParsedPanel(constructs, problems, warnings)


def load(path: str | Path) -> ParsedPanel:
    path = Path(path)
    text = path.read_text()
    if path.suffix.lower() in {'.csv', '.tsv'}:
        return parse_csv(text)
    return parse_fasta(text)
