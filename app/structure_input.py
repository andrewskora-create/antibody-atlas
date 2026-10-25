"""Read an experimental structure and work out what is in it.

The screen can only predict a structure today, which means every geometric
result is a model on top of a model. A crystal or cryo-EM structure removes
one of those layers, and if it contains the antigen it removes an inference
too: epitope distance stops being guessed from a paratope prediction and
becomes a measurement.

Nothing here assumes a layout. A deposited antibody structure can have its
chains named anything, carry several copies in the asymmetric unit, include
constant domains, and may or may not have the antigen present. So the first
job is inventory: parse the file, classify every chain by numbering it, and
report what was found rather than requiring the user to describe it.

Deliberate non-goals. Nothing is repaired, no missing loop is built, and no
chain is renamed on disk. An unresolved region must stay absent: inventing
backbone and then reporting hydrogen bonds from it would be fabrication.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any

from .ptm_precedents import position_order

#: Shorter than this and a polypeptide is a peptide, a tag or an artefact,
#: not something to classify as a chain of interest.
MIN_CHAIN_RESIDUES = 25

#: A variable domain is ~110-130 residues. A Fab chain is roughly double, and
#: numbering still finds the V domain inside it, so no upper bound is imposed.
SUPPORTED_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}


@dataclass
class ChainRecord:
    chain_id: str
    residues: int
    sequence: str
    role: str                      # heavy | light | antigen | unclassified
    chain_type: str | None = None  # H, K or L from numbering
    domain_start: int | None = None
    domain_end: int | None = None
    note: str = ""
    #: Author residue numbers, so a finding can be reported in the numbering
    #: the depositor used rather than one this tool invented.
    author_numbers: list[int] = field(default_factory=list)


class StructureInputError(ValueError):
    """The file could not be read, or holds nothing that can be analysed."""


def _parser_for(path: Path):
    from Bio.PDB import MMCIFParser, PDBParser
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True)
    if suffix in {".pdb", ".ent"}:
        return PDBParser(QUIET=True)
    raise StructureInputError(
        f"Unsupported file type {suffix!r}. Supply a PDB or mmCIF structure "
        f"(.pdb, .ent, .cif or .mmcif).")


def _chain_sequence(chain):
    """One-letter sequence and author numbering for the standard residues."""
    from Bio.PDB.Polypeptide import protein_letters_3to1 as three2one
    letters, numbers = [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue                      # heteroatoms, waters, ligands
        name = residue.get_resname()
        # Resolve once by two spellings. Resolving it twice by different
        # routes silently drops residues and corrupts the sequence -- this
        # bit the co-crystal validation once already.
        aa = three2one.get(name.capitalize()) or three2one.get(name)
        if aa:
            letters.append(aa)
            numbers.append(int(residue.id[1]))
    return "".join(letters), numbers


def _classify(sequence: str) -> tuple[str, str | None, int | None, int | None, str]:
    """Numbering decides the role. An antibody chain numbers; an antigen does not."""
    from .numbering import number_chain
    try:
        numbering = number_chain(sequence, "imgt")
    except Exception as exc:
        return ("antigen", None, None, None,
                f"No antibody variable domain found ({str(exc)[:80]}), so this is "
                f"treated as antigen or other protein.")
    chain_type = numbering.get("chain_type")
    role = "heavy" if chain_type == "H" else "light" if chain_type in {"K", "L"} else "unclassified"
    start, end = numbering.get("domain_start"), numbering.get("domain_end")
    note = ""
    if end is not None and len(sequence) > end + 40:
        note = (f"Variable domain spans residues {start}-{end} of {len(sequence)}; "
                f"the remainder is presumably a constant domain and is not numbered.")
    return role, chain_type, start, end, note


def _inventory(path: str | Path) -> dict[str, Any]:
    """Parse a structure and classify every chain in its first model."""
    path = Path(path)
    if not path.is_file():
        raise StructureInputError(f"No such file: {path}")
    if path.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise StructureInputError(
            f"Unsupported file type {path.suffix!r}. Supply a PDB or mmCIF "
            f"structure (.pdb, .ent, .cif or .mmcif).")
    try:
        structure = _parser_for(path).get_structure("input", str(path))
    except StructureInputError:
        raise
    except Exception as exc:
        raise StructureInputError(f"Could not parse {path.name}: {str(exc)[:160]}") from exc

    models = list(structure.get_models())
    if not models:
        raise StructureInputError(f"{path.name} contains no model.")

    chains: list[ChainRecord] = []
    for chain in models[0]:
        sequence, numbers = _chain_sequence(chain)
        if len(sequence) < MIN_CHAIN_RESIDUES:
            if sequence:
                chains.append(ChainRecord(
                    chain_id=chain.id, residues=len(sequence), sequence=sequence,
                    role="unclassified", author_numbers=numbers,
                    note=f"Only {len(sequence)} standard residues; too short to classify."))
            continue
        role, chain_type, start, end, note = _classify(sequence)
        chains.append(ChainRecord(
            chain_id=chain.id, residues=len(sequence), sequence=sequence,
            role=role, chain_type=chain_type, domain_start=start, domain_end=end,
            note=note, author_numbers=numbers))

    heavy = [c for c in chains if c.role == "heavy"]
    light = [c for c in chains if c.role == "light"]
    antigen = [c for c in chains if c.role == "antigen"]
    return {
        "file": path.name,
        "format": "mmCIF" if path.suffix.lower() in {".cif", ".mmcif"} else "PDB",
        "models_in_file": len(models),
        "chains": [asdict(c) for c in chains],
        "heavy_chains": [c.chain_id for c in heavy],
        "light_chains": [c.chain_id for c in light],
        "antigen_chains": [c.chain_id for c in antigen],
        "has_antigen": bool(antigen),
        "copies": min(len(heavy), len(light)) if heavy and light else 0,
        "multiple_models_note": (
            f"{len(models)} models in the file; only the first is used. For an NMR "
            f"ensemble or a multi-state cryo-EM deposition that is a choice, not a "
            f"default to rely on." if len(models) > 1 else ""),
    }

# --- public report ---------------------------------------------------------

#: Heavy-atom separation counted as an interface contact.
CONTACT_A = 5.0
#: Below this many contacts two chains are not a paired Fv. A real VH/VL
#: interface buries a large surface and gives hundreds; a crystal-packing
#: neighbour gives a handful.
MIN_INTERFACE_CONTACTS = 40
#: An antigen chain counts as ENGAGED by a given Fv above this many contacts.
#: Set from measurement, not taste: in a two-copy co-crystal the cognate
#: antigen gave 145 and 141 contacts while the neighbouring copy in the
#: lattice gave 24 and 26. Below the threshold a chain is reported as an
#: incidental contact and never used as the epitope, because measuring
#: distance to a lattice neighbour is worse than not measuring.
MIN_ANTIGEN_CONTACTS = 60


@dataclass
class StructureReport:
    path: Path
    file_format: str
    models_in_file: int
    chains: list[ChainRecord]
    fv_pairs: list[tuple[str, str]]
    antigen_chains: list[str]
    warnings: list[str] = field(default_factory=list)
    pair_detail: list[dict[str, Any]] = field(default_factory=list)

    @property
    def copies(self) -> int:
        return len(self.fv_pairs)

    def as_dict(self) -> dict[str, Any]:
        """JSON-serialisable form, for the inspect endpoint.

        The chain records carry the full observed sequence, which the caller
        needs: it is what they would paste into the form so the structure and
        the submitted sequence describe the same molecule.
        """
        return {
            "file": self.path.name,
            "format": self.file_format,
            "models_in_file": self.models_in_file,
            "chains": [asdict(c) for c in self.chains],
            "fv_pairs": [list(p) for p in self.fv_pairs],
            "antigen_chains": list(self.antigen_chains),
            "has_antigen": bool(self.antigen_chains),
            "copies": self.copies,
            "pair_detail": self.pair_detail,
            "warnings": list(self.warnings),
        }


def _heavy_atoms(chain, limit_to: tuple[int, int] | None = None):
    """Heavy atoms of the standard residues, optionally only the V domain."""
    out, index = [], 0
    for residue in chain:
        if residue.id[0] != " ":
            continue
        index += 1
        if limit_to and not (limit_to[0] <= index <= limit_to[1]):
            continue
        out.extend(a for a in residue if a.element != "H")
    return out


def inspect(path: str | Path) -> StructureReport:
    """Read a structure and work out what is in it, without being told.

    A deposited antibody structure can name its chains anything, carry
    several copies in the asymmetric unit, include constant domains, and may
    or may not have the antigen present. Nothing here is assumed: every chain
    is classified by numbering it, and heavy chains are matched to light ones
    by measured interface contact rather than by chain order or naming.
    """
    from Bio.PDB.NeighborSearch import NeighborSearch

    path = Path(path)
    inv = _inventory(path)
    chains = [ChainRecord(**c) for c in inv["chains"]]
    model = next(_parser_for(path).get_structure("input", str(path)).get_models())
    spans = {c.chain_id: (c.domain_start, c.domain_end)
             for c in chains if c.domain_start}

    def contacts_between(a: str, b: str) -> int:
        atoms_b = _heavy_atoms(model[b], spans.get(b))
        if not atoms_b:
            return 0
        search = NeighborSearch(atoms_b)
        return sum(1 for atom in _heavy_atoms(model[a], spans.get(a))
                   if search.search(atom.coord, CONTACT_A))

    pairs: list[tuple[str, str]] = []
    detail: list[dict[str, Any]] = []
    warnings: list[str] = []
    used: set[str] = set()
    for heavy in inv["heavy_chains"]:
        scored = sorted(((contacts_between(heavy, light), light)
                         for light in inv["light_chains"] if light not in used),
                        reverse=True)
        if not scored or scored[0][0] < MIN_INTERFACE_CONTACTS:
            warnings.append(
                f"Heavy chain {heavy} has no light chain across an interface; "
                f"treated as a heavy-only domain.")
            continue
        best, light = scored[0]
        used.add(light)
        pairs.append((heavy, light))
        antigens = sorted(((contacts_between(light, ag) + contacts_between(heavy, ag), ag)
                           for ag in inv["antigen_chains"]), reverse=True)
        engaged = [{"chain": ag, "contacts": n} for n, ag in antigens
                   if n >= MIN_ANTIGEN_CONTACTS]
        incidental = [{"chain": ag, "contacts": n} for n, ag in antigens
                      if 0 < n < MIN_ANTIGEN_CONTACTS]
        detail.append({"heavy": heavy, "light": light, "interface_contacts": best,
                       "engaged_antigen": [a["chain"] for a in engaged],
                       "antigen_contacts": engaged,
                       "incidental_contacts": incidental})
        if incidental:
            warnings.append(
                f"Chain(s) {', '.join(i['chain'] for i in incidental)} touch the "
                f"{heavy}/{light} Fv below the engagement threshold, most likely "
                f"lattice neighbours. They are excluded from epitope measurement.")

    if len(pairs) > 1:
        warnings.append(
            f"{len(pairs)} Fv copies in the asymmetric unit. They are reported "
            f"separately and never averaged: copies can differ, and a mean over "
            f"them would describe no molecule that exists.")
    orphans = [c for c in inv["light_chains"] if c not in used]
    if orphans:
        warnings.append(f"Light chain(s) {', '.join(orphans)} were not paired.")
    if inv["models_in_file"] > 1:
        warnings.append(
            f"{inv['models_in_file']} models in the file; only the first is used. "
            f"For an NMR ensemble or a multi-state cryo-EM deposition that is a "
            f"choice, not a default to rely on.")
    return StructureReport(
        path=path, file_format=inv["format"], models_in_file=inv["models_in_file"],
        chains=chains, fv_pairs=pairs, antigen_chains=list(inv["antigen_chains"]),
        warnings=warnings, pair_detail=detail)


# --- normalised output -----------------------------------------------------

def _imgt_labels(sequence: str) -> dict[int, str]:
    """1-based residue index -> IMGT label without the chain prefix."""
    from .numbering import number_chain
    numbering = number_chain(sequence, "imgt")
    return {r["sequence_position"]: re.sub(r"^[HKL]", "", r["number"])
            for r in numbering["residues"]}


def _atom_line(serial: int, atom, resname: str, chain: str,
               resseq: int, icode: str) -> str:
    """One ATOM record in fixed-column PDB format.

    The columns are not negotiable and the failure mode is quiet: slip the
    insertion code one place left and it lands in the residue-number field,
    after which every parser either rejects the file or reads the wrong
    residue. Written out explicitly rather than by string formatting that
    happens to line up.
    """
    name = atom.get_name()
    # Atom names are themselves column-sensitive: a one-character element
    # starts in column 14, a two-character one in column 13.
    padded = f" {name:<3}" if len(name) < 4 else name
    x, y, z = atom.coord
    element = (atom.element or "").strip().rjust(2)
    return (f"ATOM  {serial:>5} {padded}{' ':1}{resname:>3} {chain:1}"
            f"{resseq:>4}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {element}\n")


def prepare(source: str | Path, out_path: str | Path, *, copy: int = 0,
            keep_antigen: bool = True, ph: float = 7.4,
            model_index: int = 0,
            sequences: dict[str, str] | None = None) -> dict[str, Any]:
    """Write a normalised Fv: chains renamed H/L, renumbered to IMGT.

    This is what lets an experimental structure reach the rest of the screen.
    Everything downstream -- exposure, surface patches, the paratope
    predictor -- was written against a predicted model, so it expects chains
    called H and L carrying IMGT numbering. A deposited file satisfies
    neither, and the constant domain has to go too: it is not part of the Fv
    and would distort exposure and patch areas.

    `copy` selects which Fv when the asymmetric unit holds several. They are
    never merged. `keep_antigen` writes the engaged antigen alongside, which
    is what makes a measured epitope distance possible.
    """
    source, out_path = Path(source), Path(out_path)
    report = inspect(source)
    if not report.fv_pairs:
        raise StructureInputError(
            f"No paired Fv found in {source.name}; nothing to prepare.")
    if copy >= len(report.fv_pairs):
        raise StructureInputError(
            f"Copy {copy} requested but {source.name} holds "
            f"{len(report.fv_pairs)}.")
    heavy_id, light_id = report.fv_pairs[copy]
    detail = report.pair_detail[copy]
    # `model_index` matters for an NMR ensemble or a multi-state cryo-EM
    # deposition, where the file holds several conformers and picking the
    # first silently throws the rest away. The caller chooses; the choice is
    # recorded in the returned report.
    _models = list(_parser_for(source).get_structure("input", str(source)).get_models())
    if not _models:
        raise StructureInputError(f"{Path(source).name} contains no model.")
    if not 0 <= model_index < len(_models):
        raise StructureInputError(
            f"Model {model_index + 1} was requested but the file holds "
            f"{len(_models)}.")
    model = _models[model_index]
    by_id = {c.chain_id: c for c in report.chains}

    lines: list[str] = []
    serial = 1
    chain_report: dict[str, Any] = {}
    for label, chain_id in (("H", heavy_id), ("L", light_id)):
        record = by_id[chain_id]
        labels = _imgt_labels(record.sequence)
        kept = dropped = 0
        index = 0
        from Bio.PDB.Polypeptide import protein_letters_3to1 as three2one
        for residue in model[chain_id]:
            if residue.id[0] != " ":
                continue
            name = residue.get_resname()
            if not (three2one.get(name.capitalize()) or three2one.get(name)):
                continue
            index += 1
            imgt = labels.get(index)
            if imgt is None:
                # Outside the variable domain: constant region, or a residue
                # the numbering could not place. Dropped rather than renumbered
                # arbitrarily, because an invented number is worse than a gap.
                dropped += 1
                continue
            kept += 1
            resseq, icode = (int(imgt[:-1]), imgt[-1]) if imgt[-1].isalpha() else (int(imgt), " ")
            for atom in residue:
                lines.append(_atom_line(serial, atom, name, label, resseq, icode))
                serial += 1
        lines.append("TER\n")
        chain_report[label] = {
            "source_chain": chain_id,
            "residues_in_file": record.residues,
            "variable_domain_residues": kept,
            "constant_or_untyped_dropped": dropped,
            # Kept so reconciliation compares against what this chain actually
            # resolves, rather than against the whole deposited chain including
            # its constant domain.
            "observed_sequence": record.sequence,
        }

    written_antigen: list[str] = []
    if keep_antigen:
        for antigen_id in detail["engaged_antigen"]:
            if antigen_id not in model:
                continue
            written_antigen.append(antigen_id)
            for residue in model[antigen_id]:
                if residue.id[0] != " ":
                    continue
                icode = residue.id[2] if residue.id[2].strip() else " "
                for atom in residue:
                    lines.append(_atom_line(serial, atom, residue.get_resname(),
                                            antigen_id, residue.id[1], icode))
                    serial += 1
            lines.append("TER\n")
    lines.append("END\n")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(lines))

    # Reconcile against what the user actually typed, if they typed it. Doing
    # it here rather than leaving it to the caller means a structure that is
    # the WRONG MOLECULE is caught at the point of ingestion, before anything
    # downstream computes geometry from it and presents the result as
    # measured.
    reconciliation: dict[str, Any] = {}
    usable, problems = True, []
    for label, submitted in (sequences or {}).items():
        observed = (chain_report.get(label) or {}).get("observed_sequence")
        if not observed or not submitted:
            continue
        check = reconcile(observed, submitted)
        reconciliation[label] = check
        if not check["usable"]:
            usable = False
            problems.append(f"{label}: {check['note']}")

    return {
        "path": str(out_path),
        "source": "experimental",
        "source_format": report.file_format,
        "source_file": source.name,
        "copy": copy,
        "copy_used": copy,
        "model_used": model_index,
        "models_available": len(_models),
        "copies_available": len(report.fv_pairs),
        "chains": chain_report,
        "antigen_chains": written_antigen,
        "has_antigen": bool(written_antigen),
        "reconciliation": reconciliation,
        "usable": usable,
        "problems": problems,
        "warnings": report.warnings,
        "note": ("Chains renamed to H and L and renumbered to IMGT so the rest of "
                 "the screen can read them. Constant domains and residues the "
                 "numbering could not place are dropped, not renumbered. Nothing "
                 "is repaired and no missing loop is built: an unresolved region "
                 "stays absent, because inventing backbone and then reporting "
                 "geometry from it would be fabrication."),
    }


# --- reconciliation with the submitted sequence ----------------------------

#: Below this fraction of the submitted sequence resolved, exposure and
#: hydrogen-bond results would describe a fragment rather than the molecule.
MIN_COVERAGE = 0.85


def reconcile(observed: str, submitted: str) -> dict[str, Any]:
    """Compare what a structure resolves against the sequence that was submitted.

    A deposited chain is normally a subsequence of what the user types, with
    gaps where there was no density, so this is a gapped containment check
    rather than a global alignment.

    Unresolved runs are reported individually and never closed up. Silently
    joining the two sides of a gap would renumber everything downstream of it
    by the length of the gap, which is the kind of error that produces
    confident answers about the wrong residues.
    """
    from difflib import SequenceMatcher

    observed = (observed or "").strip().upper()
    submitted = (submitted or "").strip().upper()
    if not submitted:
        raise StructureInputError("No submitted sequence to reconcile against.")

    blocks = SequenceMatcher(None, submitted, observed, autojunk=False).get_matching_blocks()
    matched = sum(b.size for b in blocks)
    covered: set[int] = set()
    for block in blocks:
        covered.update(range(block.a, block.a + block.size))

    runs, start = [], None
    for i in range(len(submitted) + 1):
        missing = i < len(submitted) and i not in covered
        if missing and start is None:
            start = i
        elif not missing and start is not None:
            runs.append({"start": start + 1, "end": i, "length": i - start})
            start = None

    coverage = len(covered) / len(submitted)
    # Two different questions, which must not share a denominator. "Is this
    # the same molecule" is matched over what was OBSERVED; "how much did it
    # resolve" is coverage over what was SUBMITTED. Dividing both by the
    # submitted length once caused a 25-residue purification tag to be
    # refused as a different molecule when every resolved residue matched.
    identity = matched / len(observed) if observed else 0.0

    usable = identity >= 0.90 and coverage >= MIN_COVERAGE
    if identity < 0.90:
        note = (f"Only {identity:.0%} of the resolved residues match the submitted "
                f"sequence. This is a different molecule, not missing density, and "
                f"proceeding would be guesswork.")
    elif coverage < MIN_COVERAGE:
        note = (f"Only {coverage:.0%} of the submitted sequence is resolved. Too "
                f"much is missing to describe the whole domain; anything said "
                f"about the absent regions would be guesswork.")
    elif runs:
        note = (f"{len(runs)} unresolved run(s), {sum(r['length'] for r in runs)} "
                f"residue(s) total. Those positions have no coordinates and are "
                f"reported as unassessed rather than as buried.")
    else:
        note = "Every submitted residue is resolved in the structure."

    return {
        "coverage": round(coverage, 4),
        "identity_of_resolved": round(identity, 4),
        "matched_residues": matched,
        "submitted_length": len(submitted),
        "observed_length": len(observed),
        "unresolved_runs": runs,
        "usable": usable,
        "note": note,
    }


# --- measured epitope distance ---------------------------------------------

#: Heavy-atom separation at or below which a residue contacts the antigen.
EPITOPE_CONTACT_A = 4.5
#: Beyond this, report "away from the interface" rather than a number that
#: implies more precision than it carries.
EPITOPE_FAR_A = 15.0


def epitope_distances(prepared_pdb: str | Path,
                      antigen_chains: list[str]) -> dict[str, Any]:
    """Minimum distance from each antibody residue to the antigen, by IMGT label.

    This is the reason to supply a complex rather than an antibody alone. The
    predicted route has no antigen in it: Paragraph predicts the PARATOPE --
    positions likely to contact something -- so a distance to it is a distance
    to an inference. With the antigen present this becomes a measurement, and
    the prediction can be scored against it rather than trusted.

    Takes a file already through `prepare`, so chains are H/L and residue
    numbers are IMGT.
    """
    from Bio.PDB import PDBParser
    from Bio.PDB.NeighborSearch import NeighborSearch

    prepared_pdb = Path(prepared_pdb)
    model = next(PDBParser(QUIET=True).get_structure("prepared", str(prepared_pdb)).get_models())
    present = [c for c in antigen_chains if c in model]
    if not present:
        return {
            "status": "unavailable",
            "reason": ("No antigen chain is present in this structure, so there is "
                       "no epitope to measure a distance to. A paratope prediction "
                       "is an inference about where antigen would bind; it is not "
                       "an epitope and no distance is reported from it here."),
            "residues": [], "contacts": [], "contact_count": 0,
        }

    antigen_atoms = [a for cid in present for r in model[cid] if r.id[0] == " "
                     for a in r if a.element != "H"]
    search = NeighborSearch(antigen_atoms)

    rows, contacts = [], []
    for label in ("H", "L"):
        if label not in model:
            continue
        for residue in model[label]:
            if residue.id[0] != " ":
                continue
            atoms = [a for a in residue if a.element != "H"]
            if not atoms:
                continue
            nearest = None
            for atom in atoms:
                for other in search.search(atom.coord, EPITOPE_FAR_A):
                    d = float((((atom.coord - other.coord) ** 2).sum()) ** 0.5)
                    if nearest is None or d < nearest:
                        nearest = d
            number = f"{label}{residue.id[1]}{residue.id[2].strip()}"
            in_contact = nearest is not None and nearest <= EPITOPE_CONTACT_A
            rows.append({
                "position": number,
                "chain": label,
                "distance_a": None if nearest is None else round(nearest, 2),
                "contact": in_contact,
                "beyond_cutoff": nearest is None,
            })
            if in_contact:
                contacts.append(number)

    from .ptm_precedents import position_order
    contacts.sort(key=position_order)
    return {
        "status": "measured",
        "source": "experimental",
        "antigen_chains": present,
        "contact_threshold_a": EPITOPE_CONTACT_A,
        "residues": rows,
        "contacts": contacts,
        "contact_count": len(contacts),
        "note": (f"Heavy-atom minima to the antigen in the supplied coordinates. "
                 f"Residues beyond {EPITOPE_FAR_A:g} A are reported as beyond the "
                 f"cutoff rather than with a number. This is one deposited complex: "
                 f"it is not evidence about a different antigen, a different "
                 f"epitope, or the unbound state."),
    }
