#!/usr/bin/env python
"""Validate the paratope predictor against a held-out antibody-antigen cocrystal.

Runs the prediction on sacituzumab bound to Trop2 (PDB 9PI9, 1.56 A X-ray) and
scores it against the experimentally observed contacts. Everything is
downloaded at run time; no structure is redistributed with this project.

    python validation/paratope_cocrystal.py

Why this is a fair test:

  * Held out, and the script proves it rather than asserting it. Paragraph's
    own train/validation/test lists are fetched from its repository and the
    entry is checked against them. 9PI9 was released in December 2025, long
    after those splits were fixed.
  * The interface is real. Trop2 is a protein antigen, not a short peptide in
    a groove, so the observed paratope is a genuine protein-protein interface.
  * The coordinates are good. At 1.56 A the side-chain positions are reliable,
    so the contact set is a measurement rather than an estimate. Contact sets
    derived from low-resolution cryo-EM are not.
  * Two independent copies. The asymmetric unit holds two Fab-antigen
    complexes, so the same interface is measured twice and the spread between
    them bounds how much of any score is crystallographic noise.
  * Sacituzumab is a clinically approved therapeutic, which is the kind of
    molecule this project is used on.

The prediction is run three ways: on each crystal copy, and on an
ABodyBuilder2 model built from sequence alone. The last of those is what the
pipeline actually does, and it is the number to believe for ordinary use.

One antibody is a sanity check, not a benchmark.
"""
from __future__ import annotations

import csv
import io
import json
import math
import sys
import urllib.request
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
WORK = Path(__file__).resolve().parent / ".cocrystal_work"

ENTRY = "9PI9"
DESCRIPTION = "Sacituzumab Fab bound to Trop2 (1.56 A X-ray, released 2025-12-10)"
# (heavy chain id, light chain id) for each Fab copy in the asymmetric unit.
FAB_COPIES = (("A", "B"), ("C", "D"))
ANTIGEN_CHAINS = ("E", "F")
THRESHOLDS = (4.0, 4.5, 5.0)
# Paragraph scores only IMGT CDR +/- 2: 25-40, 54-67, 103-119.
WINDOWS = ((25, 40, "CDR1"), (54, 67, "CDR2"), (103, 119, "CDR3"))

SPLIT_URLS = [
    f"https://raw.githubusercontent.com/oxpig/Paragraph/main/training_data/{dataset}/{split}.csv"
    for dataset in ("Expanded", "PECAN")
    for split in ("train_set", "val_set", "test_set")
]


def region_of(imgt: str) -> str:
    n = int("".join(c for c in imgt if c.isdigit()))
    for lo, hi, name in WINDOWS:
        if lo <= n <= hi:
            return name
    return "other"


def fetch(url: str, target: Path) -> Path:
    WORK.mkdir(exist_ok=True)
    if not target.is_file():
        print(f"  downloading {url}")
        with urllib.request.urlopen(url, timeout=120) as response:
            target.write_bytes(response.read())
    return target


def verify_held_out(entry: str) -> bool:
    """Check the entry against Paragraph's own published splits.

    A hold-out claim a reader cannot check is just an assertion, so this
    downloads the lists the model was trained on and looks.
    """
    seen: set[str] = set()
    for url in SPLIT_URLS:
        name = "_".join(url.rsplit("/", 2)[-2:])
        try:
            path = fetch(url, WORK / name)
        except Exception as exc:                      # offline, or layout changed
            print(f"  could not fetch {name}: {type(exc).__name__}. "
                  "Hold-out NOT verified in this run.")
            return False
        for row in csv.reader(io.StringIO(path.read_text())):
            if row:
                seen.add(row[0].strip().lower()[:4])
    held_out = entry.lower()[:4] not in seen
    print(f"  {len(seen)} unique PDB entries across Paragraph's six split files")
    print(f"  {entry} present in them: {'YES - NOT held out' if not held_out else 'no'}")
    return held_out


def _chain_rows(model, chain_id):
    from Bio.PDB.Polypeptide import protein_letters_3to1 as three2one
    rows = []
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        name = residue.get_resname()
        # Resolve the one-letter code once. Resolving it twice by different
        # routes silently drops residues and corrupts the sequence.
        aa = three2one.get(name.capitalize()) or three2one.get(name)
        if aa:
            rows.append((residue, aa))
    return rows


def _imgt_map(rows):
    from app.numbering import number_chain
    numbering = number_chain("".join(a for _, a in rows), "imgt")
    return {x["sequence_position"]: x["number"].lstrip("HKL")
            for x in numbering["residues"]}


def contacts(cif: Path, heavy: str, light: str) -> dict[float, set[tuple[str, str]]]:
    """Antibody residues near antigen, keyed by H/L label and IMGT number."""
    from Bio.PDB import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch

    model = next(MMCIFParser(QUIET=True).get_structure("x", str(cif)).get_models())
    antigen = [a for cid in ANTIGEN_CHAINS if cid in model
               for r in model[cid] if r.id[0] == " " for a in r if a.element != "H"]
    search = NeighborSearch(antigen)

    out: dict[float, set[tuple[str, str]]] = {t: set() for t in THRESHOLDS}
    for label, cid in (("H", heavy), ("L", light)):
        rows = _chain_rows(model, cid)
        position_to_imgt = _imgt_map(rows)
        for i, (residue, _) in enumerate(rows, start=1):
            imgt = position_to_imgt.get(i)
            if not imgt:
                continue  # constant domain, outside the variable numbering
            for threshold in THRESHOLDS:
                if any(search.search(a.coord, threshold) for a in residue if a.element != "H"):
                    out[threshold].add((label, imgt))
    return out


def sequences(cif: Path, heavy: str, light: str) -> dict[str, str]:
    """Variable-domain sequences, for the model-from-sequence path."""
    from Bio.PDB import MMCIFParser
    model = next(MMCIFParser(QUIET=True).get_structure("x", str(cif)).get_models())
    out = {}
    for label, cid in (("H", heavy), ("L", light)):
        rows = _chain_rows(model, cid)
        position_to_imgt = _imgt_map(rows)
        out[label] = "".join(aa for i, (_, aa) in enumerate(rows, start=1)
                             if position_to_imgt.get(i))
    return out


def write_imgt_pdb(cif: Path, heavy: str, light: str, out_pdb: Path) -> None:
    """Re-chain to H/L and renumber to IMGT, which is what Paragraph expects."""
    from Bio.PDB import MMCIFParser

    model = next(MMCIFParser(QUIET=True).get_structure("x", str(cif)).get_models())
    lines, serial = [], 1
    for label, cid in (("H", heavy), ("L", light)):
        rows = _chain_rows(model, cid)
        position_to_imgt = _imgt_map(rows)
        for i, (residue, _) in enumerate(rows, start=1):
            imgt = position_to_imgt.get(i)
            if not imgt:
                continue
            insert = " "
            if imgt[-1].isalpha():
                imgt, insert = imgt[:-1], imgt[-1]
            for atom in residue:
                if atom.element == "H" or atom.get_altloc() not in (" ", "A"):
                    continue
                x, y, z = atom.coord
                name = atom.get_name()
                nm = f" {name:<3}" if len(name) < 4 else name
                lines.append(
                    f"ATOM  {serial:>5} {nm} {residue.get_resname():>3} {label}"
                    f"{int(imgt):>4}{insert}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                    f"  1.00{atom.get_bfactor():>6.2f}          {atom.element:>2}")
                serial += 1
    out_pdb.write_text("\n".join(lines) + "\nEND\n")


def build_model(seqs: dict[str, str], out_pdb: Path) -> Path | None:
    """ABodyBuilder2 model from sequence alone, as the pipeline does it."""
    if out_pdb.is_file():
        return out_pdb
    try:
        from ImmuneBuilder import ABodyBuilder2
    except ImportError as exc:
        print(f"  ImmuneBuilder unavailable ({exc}); skipping the sequence-built model.")
        return None
    # model_ids=[1] matches what app/main.py does. ABodyBuilder2 ships a
    # four-model ensemble, but the pipeline runs a single model, so validating
    # the ensemble would measure something the user never gets. This is also
    # the single-model variance that TAP2 and Park & Izadi both warn about.
    print("  predicting structure from sequence (ABodyBuilder2, model 1, as the pipeline does)")
    ABodyBuilder2(model_ids=[1]).predict({"H": seqs["H"], "L": seqs["L"]}).save(str(out_pdb))
    return out_pdb


def score(predicted: dict, truth: set, cutoff: float) -> dict:
    tp = sum(1 for k, p in predicted.items() if p >= cutoff and k in truth)
    fp = sum(1 for k, p in predicted.items() if p >= cutoff and k not in truth)
    fn = sum(1 for k, p in predicted.items() if p < cutoff and k in truth)
    tn = len(predicted) - tp - fp - fn
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "precision": precision,
            "recall": recall,
            "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
            "mcc": ((tp * tn) - (fp * fn)) / denominator if denominator else 0.0}


def auc(predicted: dict, truth: set) -> float | None:
    # Ascending, so rank 1 is the lowest score, as Mann-Whitney U requires.
    pairs = sorted((p, k in truth) for k, p in predicted.items())
    positives = sum(1 for _, t in pairs if t)
    negatives = len(pairs) - positives
    if not positives or not negatives:
        return None
    ranks, i = {}, 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[k] = (i + j) / 2 + 1
        i = j + 1
    total = sum(ranks[k] for k, (_, t) in enumerate(pairs) if t)
    return (total - positives * (positives + 1) / 2) / (positives * negatives)


def main() -> int:
    from app.paratope import predict, CONTACT_CUTOFF, available

    status = available()
    if not status["available"]:
        print(f"Paragraph is not installed: {status['reason']}")
        return 1

    print(f"Paratope validation against {ENTRY}: {DESCRIPTION}\n")

    print("Checking the entry is held out of Paragraph's training data")
    if not verify_held_out(ENTRY):
        print("\n  WARNING: this entry is in Paragraph's splits, or the check could not run. "
              "Scores below are not a held-out result.")

    print("\nStructure")
    cif = fetch(f"https://files.rcsb.org/download/{ENTRY}.cif", WORK / f"{ENTRY.lower()}.cif")

    inputs: list[tuple[Path, str, set]] = []
    truth_by_threshold_first: dict | None = None

    for index, (heavy, light) in enumerate(FAB_COPIES, start=1):
        label = f"crystal copy {index} (chains {heavy}/{light})"
        print(f"\nObserved contacts, {label}")
        truth_by_threshold = contacts(cif, heavy, light)
        for threshold in THRESHOLDS:
            print(f"  {threshold} A: {len(truth_by_threshold[threshold])} residues")
        if truth_by_threshold_first is None:
            truth_by_threshold_first = truth_by_threshold
        pdb = WORK / f"{ENTRY.lower()}_copy{index}_imgt.pdb"
        if not pdb.is_file():
            write_imgt_pdb(cif, heavy, light, pdb)
        inputs.append((pdb, label, truth_by_threshold[4.5]))

    # The pipeline never has a crystal structure; it has a sequence. Score that too.
    heavy, light = FAB_COPIES[0]
    seqs = sequences(cif, heavy, light)
    print(f"\nSequence-built model (VH {len(seqs['H'])} aa, VL {len(seqs['L'])} aa)")
    built = build_model(seqs, WORK / f"{ENTRY.lower()}_abb2.pdb")
    if built:
        inputs.append((built, "ABodyBuilder2 model from sequence", inputs[0][2]))

    print("\nPrediction versus observation")
    best = None
    for pdb, label, truth in inputs:
        result = predict(pdb, heavy_chain_id="H", light_chain_id="L")
        if result["status"] != "available":
            print(f"  {label}: failed, {result.get('reason')}")
            continue
        predicted = {(c, rec["imgt"]): rec["probability"]
                     for c, residues in result["by_chain"].items()
                     for rec in residues.values()}
        unreachable = sorted(truth - set(predicted))
        m = score(predicted, truth, CONTACT_CUTOFF)
        area = auc(predicted, truth)
        print(f"\n  {label}")
        print(f"    {len(predicted)} positions scored; {len(truth) - len(unreachable)}/{len(truth)} "
              f"observed contacts lie inside the scored window")
        if unreachable:
            print(f"    outside the scored window, so unreachable: "
                  f"{', '.join(c + n for c, n in unreachable)}")
        print(f"    TP {m['tp']}  FP {m['fp']}  FN {m['fn']}  TN {m['tn']}")
        print(f"    precision {m['precision']:.2f}  recall {m['recall']:.2f}  "
              f"MCC {m['mcc']:.2f}  AUC {area:.3f}" if area is not None else "    AUC undefined")
        if best is None:
            best = (predicted, truth)

    if best and truth_by_threshold_first:
        predicted, truth = best
        print("\n  Sensitivity to the contact distance (crystal copy 1)")
        for threshold in THRESHOLDS:
            t = truth_by_threshold_first[threshold]
            m = score(predicted, t, CONTACT_CUTOFF)
            area = auc(predicted, t)
            print(f"    {threshold} A ({len(t):>2} contacts): precision {m['precision']:.2f}  "
                  f"recall {m['recall']:.2f}  MCC {m['mcc']:.2f}  "
                  + (f"AUC {area:.3f}" if area is not None else "AUC -"))

        print("\n  By loop (crystal copy 1, 4.5 A)")
        print(f"    {'loop':<10}{'scored':<8}{'observed':<10}{'TP':<4}{'FP':<4}{'FN':<4}{'precision'}")
        for chain in ("H", "L"):
            for name in ("CDR1", "CDR2", "CDR3"):
                keys = [k for k in predicted if k[0] == chain and region_of(k[1]) == name]
                tp = sum(1 for k in keys if predicted[k] >= CONTACT_CUTOFF and k in truth)
                fp = sum(1 for k in keys if predicted[k] >= CONTACT_CUTOFF and k not in truth)
                fn = sum(1 for k in keys if predicted[k] < CONTACT_CUTOFF and k in truth)
                observed = sum(1 for k in keys if k in truth)
                precision = f"{tp/(tp+fp):.2f}" if tp + fp else "-"
                print(f"    {chain + ' ' + name:<10}{len(keys):<8}{observed:<10}"
                      f"{tp:<4}{fp:<4}{fn:<4}{precision}")

    print("\nOne antibody is a sanity check, not a benchmark. See "
          "PTM_FIX_ENGINEERING.md section 5b for what this does and does not establish.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
