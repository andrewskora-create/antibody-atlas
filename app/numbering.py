"""ANARCII numbering and explicit CDR definitions."""
from __future__ import annotations

from functools import lru_cache
from threading import RLock
from typing import Any


CDR_RANGES = {
    "imgt": {"H": [(27, 38), (56, 65), (105, 117)],
             "L": [(27, 38), (56, 65), (105, 117)]},
    "kabat": {"H": [(31, 35), (50, 65), (95, 102)],
              "L": [(24, 34), (50, 56), (89, 97)]},
}
_numbering_lock = RLock()


@lru_cache(maxsize=1)
def _model():
    from anarcii import Anarcii
    return Anarcii(seq_type="antibody", cpu=True, ncpu=2, verbose=False)


def number_chain(sequence: str, scheme: str) -> dict[str, Any]:
    if scheme not in CDR_RANGES:
        raise ValueError("Numbering scheme must be IMGT or Kabat.")
    with _numbering_lock:
        model = _model()
        result = model.number({"chain": sequence})
        if scheme == "kabat":
            result = model.to_scheme("kabat")
    record = result.get("chain")
    if not record or record.get("error") or not record.get("numbering"):
        raise ValueError("ANARCII could not recognize a complete antibody variable domain.")
    chain_type = record.get("chain_type")
    family = "H" if chain_type == "H" else "L" if chain_type in {"K", "L"} else None
    if family is None:
        raise ValueError("ANARCII classified this sequence outside antibody heavy or light chains.")
    start = int(record.get("query_start") or 0)
    residues: list[dict[str, Any]] = []
    sequence_position = start
    cdr_ranges = CDR_RANGES[scheme][family]
    for (position, insertion), aa in record["numbering"]:
        if aa == "-":
            continue
        sequence_position += 1
        label = f"{chain_type}{position}{insertion.strip()}"
        region = next((f"CDR{i}" for i, (lo, hi) in enumerate(cdr_ranges, 1)
                       if lo <= position <= hi), "Framework")
        residues.append({"sequence_position": sequence_position,
                         "number": label, "number_base": position,
                         "aa": aa, "region": region})
    if "".join(r["aa"] for r in residues) not in sequence:
        raise ValueError("Numbering could not be mapped back to the input sequence.")
    cdrs = []
    for i in range(1, 4):
        selected = [r for r in residues if r["region"] == f"CDR{i}"]
        cdrs.append({"name": f"CDR{i}", "sequence": "".join(r["aa"] for r in selected),
                     "start": selected[0]["sequence_position"] if selected else None,
                     "end": selected[-1]["sequence_position"] if selected else None,
                     "number_start": selected[0]["number"] if selected else None,
                     "number_end": selected[-1]["number"] if selected else None})
    return {"chain_type": chain_type, "scheme": scheme,
            "definition": scheme.upper(), "model": "ANARCII",
            "domain_start": residues[0]["sequence_position"],
            "domain_end": residues[-1]["sequence_position"],
            "residues": residues, "cdrs": cdrs}
