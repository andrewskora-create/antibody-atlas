#!/usr/bin/env python3
"""Fetch human immunoglobulin constant-region sequences from UniProt.

Unlike the Thera-SAbDab corpora, this data IS redistributable: UniProt applies
CC BY 4.0 to all copyrightable parts of its databases and permits copying so
long as the copyright notice travels with it. The built file therefore ships in
the repository, and this script exists so anyone can rebuild and verify it.

    python setup/fetch_constant_regions.py            # rebuild
    python setup/fetch_constant_regions.py --verify   # check the shipped file
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "data" / "ig_constant_regions.json"
BASE = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

# Heavy-chain entries are requested as the "-1" ISOFORM, which is the secreted
# form. UniProt's displayed isoform for the gamma chains is the MEMBRANE-BOUND
# one: P01857 as served is 399 aa and carries a transmembrane helix at 351-371,
# where the secreted constant region is 330 aa ending ...SLSLSPGK. Appending the
# membrane form to an Fv would add a hydrophobic helix that is not in the drug
# and would corrupt every hydrophobicity and charge number downstream.
#
# accession -> (key, kind, label, note)
ENTRIES = {
    "P01857-1": ("IGHG1", "heavy", "IgG1",
               "Most common therapeutic isotype. Allotypes (G1m) differ at a small number of "
               "positions, some of them charged; this is the UniProt canonical sequence."),
    "P01859-1": ("IGHG2", "heavy", "IgG2",
               "Shorter hinge with four cysteines; reduced effector function."),
    "P01860-1": ("IGHG3", "heavy", "IgG3",
               "Long hinge and the most allotypic variation of the four subclasses, including "
               "hinge-length differences between allotypes. Treat this sequence as one allotype, "
               "not as representative."),
    "P01861-1": ("IGHG4", "heavy", "IgG4",
               "Therapeutic IgG4 is usually engineered with the S228P hinge substitution to "
               "prevent Fab-arm exchange; that substitution is NOT applied here."),
    "P01834": ("IGKC", "light", "Kappa", "Single human kappa constant region."),
    "P0CG04": ("IGLC1", "light", "Lambda (IGLC1)", "One of several lambda constant genes."),
    "P0DOY2": ("IGLC2", "light", "Lambda (IGLC2)", "The lambda constant most used in therapeutics."),
    "P0DOY3": ("IGLC3", "light", "Lambda (IGLC3)", "Differs from IGLC2 at a small number of positions."),
    "A0M8Q6": ("IGLC7", "light", "Lambda (IGLC7)", "A less common lambda constant gene."),
}

COPYRIGHT = ("UniProtKB data are distributed under the Creative Commons Attribution 4.0 "
             "International (CC BY 4.0) License. Copyright 2002-2026 UniProt Consortium. "
             "See https://www.uniprot.org/help/license")


def fetch(accession: str) -> tuple[str, str]:
    request = urllib.request.Request(BASE.format(acc=accession),
                                     headers={"User-Agent": "antibody-atlas/constant-regions"})
    with urllib.request.urlopen(request, timeout=45) as response:
        release = response.headers.get("x-uniprot-release", "")
        text = response.read().decode()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines or not lines[0].startswith(">"):
        raise SystemExit(f"{accession}: unexpected response from UniProt")
    sequence = "".join(lines[1:])
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise SystemExit(f"{accession}: non-standard residues in returned sequence")
    return sequence, release


def check_secreted(key: str, kind: str, sequence: str) -> None:
    """Refuse a membrane-bound gamma chain.

    The secreted heavy constant region ends at the C-terminal lysine of CH3
    (...SLSLSPGK, or ...SLSLSLGK for IgG4). The membrane form continues into a
    transmembrane helix instead. Catching that here keeps a silently wrong
    sequence out of every descriptor that uses it.
    """
    if kind != "heavy":
        if not 100 <= len(sequence) <= 115:
            raise SystemExit(f"{key}: light constant region of {len(sequence)} aa is out of range")
        return
    if not sequence.endswith("GK"):
        raise SystemExit(f"{key}: does not end at the CH3 C-terminal lysine "
                         f"(ends ...{sequence[-12:]}). This looks like the membrane-bound "
                         f"isoform; request the '-1' secreted isoform.")
    if not 300 <= len(sequence) <= 400:
        raise SystemExit(f"{key}: heavy constant region of {len(sequence)} aa is out of range")
    window = sequence[-120:]
    for start in range(len(window) - 19):
        stretch = window[start:start + 20]
        if sum(a in "AVLIMFWCY" for a in stretch) >= 17:
            raise SystemExit(f"{key}: a 20-residue hydrophobic stretch was found near the C "
                             f"terminus ({stretch}). That is a transmembrane helix, not a "
                             f"secreted constant region.")


def build() -> dict:
    regions, release = {}, ""
    for accession, (key, kind, label, note) in ENTRIES.items():
        sequence, rel = fetch(accession)
        check_secreted(key, kind, sequence)
        release = release or rel
        regions[key] = {"key": key, "kind": kind, "label": label, "accession": accession,
                        "length": len(sequence), "sequence": sequence, "note": note,
                        "sha256": hashlib.sha256(sequence.encode()).hexdigest()}
        print(f"  {key:<7} {accession}  {len(sequence):>4} aa", file=sys.stderr)
    return {
        "source": "UniProtKB (reviewed/Swiss-Prot)",
        "source_url": "https://www.uniprot.org/",
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "copyright": COPYRIGHT,
        "redistribution": "Permitted with the copyright notice reproduced; this file ships in-tree.",
        "uniprot_release": release,
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "limitations": [
            "Constant-region allotypes differ at charged positions and shift the computed pI; "
            "each entry is one allotype, not a consensus.",
            "No glycan is modelled. Sialylated glycoforms carry negative charge that lowers the "
            "real pI, and the C-terminal lysine of the heavy chain is usually clipped in product.",
            "Engineering common in therapeutics (IgG4 S228P, Fc silencing, half-life or charge "
            "variants) is not applied.",
        ],
        "regions": regions,
    }


def canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="Rebuild and compare against the shipped file instead of writing.")
    args = parser.parse_args()

    print("Fetching human Ig constant regions from UniProt…", file=sys.stderr)
    payload = build()

    if args.verify:
        if not OUT.is_file():
            print(f"FAIL: {OUT} is missing.", file=sys.stderr)
            return 1
        shipped = json.loads(OUT.read_text())
        drift = [k for k, v in payload["regions"].items()
                 if shipped.get("regions", {}).get(k, {}).get("sha256") != v["sha256"]]
        if drift:
            print(f"FAIL: sequence drift in {', '.join(drift)}", file=sys.stderr)
            return 1
        print(f"OK: all {len(payload['regions'])} constant regions match the shipped file.",
              file=sys.stderr)
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(canonical(payload))
    print(f"Wrote {OUT} ({len(payload['regions'])} regions, UniProt release "
          f"{payload['uniprot_release']}).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
