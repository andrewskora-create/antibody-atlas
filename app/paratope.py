"""Paratope prediction, used to steer PTM fixes away from antigen contacts.

A PTM liability inside a CDR is the common case, and the substitution that best
removes the chemistry is often at a position that touches antigen. Predicting
which positions those are lets an equivalent fix elsewhere be preferred, and
lets a fix that must sit on a contact be flagged rather than proposed silently.

Prediction uses Paragraph (Chinery et al., *Bioinformatics* 2023), an
equivariant graph neural network over the Fv structure, BSD-3-Clause licensed
with pretrained weights. Paragraph needs an IMGT-numbered Fv, which is what this
project's structure step already produces.

A prediction is not an epitope. Paragraph reports a per-residue probability of
belonging to the paratope; the authors recommend a classifier cut-off of 0.734
to recover the observed number of binding residues. That cut-off is used here
and is reported alongside every call, so a reader can see what "likely contact"
meant.
"""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from typing import Any

from .ptm_precedents import position_order

# Authors' recommended classifier cut-off (Paragraph README): recovers the true
# number of binding residues in the CDR loops plus two residues either side.
CONTACT_CUTOFF = 0.734
METHOD = "Paragraph 1.0.0 (EGNN), pretrained paired-chain weights"
WEIGHTS_FILE = "pretrained_weights.pt"
# 20D one-hot residue type + 2D one-hot chain id, per the pretrained model.
NUM_FEATS = 22
CITATION = ("Chinery, Wahome, Moal, Deane. Paragraph: antibody paratope "
            "prediction using graph neural networks with minimal feature vectors. "
            "Bioinformatics 39(1) 2023. doi:10.1093/bioinformatics/btac732")

LIMITATION = (
    "Predicted paratope, not a measured epitope or a binding assay. Paragraph is "
    "trained on bound antibody-antigen structures and predicts which residues "
    "would contact an antigen, without knowing this antibody's antigen. A "
    "likely-contact position is a reason to prefer an equivalent fix elsewhere "
    "and a reason to test binding, never a prohibition on substitution. A "
    "position predicted not to contact can still be essential to affinity "
    "through conformation."
)


def available() -> dict[str, Any]:
    try:
        import Paragraph  # noqa: F401
    except ImportError as exc:
        return {"available": False, "reason": f"Paragraph is not installed ({exc}). "
                                              f"Install with: pip install "
                                              f"git+https://github.com/oxpig/Paragraph"}
    return {"available": True, "method": METHOD, "cutoff": CONTACT_CUTOFF}


def _band(probability: float) -> str:
    if probability >= CONTACT_CUTOFF:
        return "likely contact"
    if probability >= 0.4:
        return "possible contact"
    if probability >= 0.15:
        return "unlikely contact"
    return "not predicted to contact"


def predict(
    pdb_path: str | Path,
    *,
    heavy_chain_id: str = "H",
    light_chain_id: str | None = "L",
    chain_label_map: dict[str, str] | None = None,
    imgt_to_position: dict[str, dict[str, int]] | None = None,
) -> dict[str, Any]:
    """Predict the paratope of one IMGT-numbered Fv structure.

    ``chain_label_map`` maps PDB chain id to this project's chain label, for the
    usual case where they already agree. ``imgt_to_position`` maps chain label
    and IMGT number to one-based sequence position, so predictions can be keyed
    the way the rest of the pipeline keys residues; without it, results are
    keyed by IMGT number only and the caller must map them.
    """
    status = available()
    if not status["available"]:
        return {"status": "unavailable", "reason": status["reason"],
                "by_chain": {}, "limitation": LIMITATION}

    import torch
    from Paragraph.model import EGNN_Model
    from Paragraph.predict import get_dataloader, evaluate_model

    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        return {"status": "unavailable", "reason": f"No structure at {pdb_path}.",
                "by_chain": {}, "limitation": LIMITATION}

    label_map = chain_label_map or {heavy_chain_id: "H",
                                    **({light_chain_id: "L"} if light_chain_id else {})}

    weights = Path(__import__("Paragraph").__file__).parent / "trained_model" / WEIGHTS_FILE
    if not weights.is_file():
        return {"status": "unavailable", "reason": f"Paragraph weights missing at {weights}.",
                "by_chain": {}, "limitation": LIMITATION}

    try:
        device = torch.device("cpu")
        # Architecture must match the pretrained checkpoint: 20D one-hot residue
        # type plus 2D one-hot chain id, six graph layers, two linear layers.
        net = EGNN_Model(num_feats=NUM_FEATS,
                         graph_hidden_layer_output_dims=[NUM_FEATS] * 6,
                         linear_hidden_layer_output_dims=[10] * 2)
        net.load_state_dict(torch.load(str(weights), map_location=device))
        net = net.to(device)

        with tempfile.TemporaryDirectory() as workdir:
            key_csv = Path(workdir) / "key.csv"
            # Paragraph's key file carries no header row; its own example is a
            # bare "4edw,H,L". A header would be parsed as a structure name.
            with key_csv.open("w", newline="") as handle:
                csv.writer(handle).writerow(
                    [pdb_path.stem, heavy_chain_id, light_chain_id or ""])
            dataloader = get_dataloader(str(key_csv), str(pdb_path.parent))
            frame = evaluate_model(model=net, dataloader=dataloader, device=device)
        rows = frame.to_dict("records")
    except Exception as exc:  # noqa: BLE001 - report, never crash the analysis run
        return {"status": "error", "reason": f"{type(exc).__name__}: {exc}",
                "method": METHOD, "by_chain": {}, "limitation": LIMITATION}

    by_chain: dict[str, dict[int | str, Any]] = {}
    contacts: dict[str, list[str]] = {}
    for row in rows:
        chain_id = str(row.get("chain_id") or row.get("chain_type") or "").strip()
        label = label_map.get(chain_id, chain_id)
        imgt = str(row.get("IMGT") or "").strip()
        try:
            probability = float(row.get("pred") or "nan")
        except ValueError:
            continue
        if probability != probability:  # NaN
            continue
        record = {
            "imgt": imgt,
            "aa": str(row.get("AA") or "").strip(),
            "probability": round(probability, 4),
            "band": _band(probability),
            "likely_contact": probability >= CONTACT_CUTOFF,
        }
        bucket = by_chain.setdefault(label, {})
        position = ((imgt_to_position or {}).get(label) or {}).get(imgt)
        bucket[position if position is not None else f"IMGT{imgt}"] = record
        if record["likely_contact"]:
            contacts.setdefault(label, []).append(imgt)

    return {
        "status": "available" if by_chain else "error",
        "method": METHOD,
        "citation": CITATION,
        "cutoff": CONTACT_CUTOFF,
        "by_chain": by_chain,
        "summary": {
            "scored_residues": sum(len(v) for v in by_chain.values()),
            "likely_contact_count": sum(len(v) for v in contacts.values()),
            # IMGT order, not string order. The previous key was (len, s),
            # which orders equal-width numbers correctly and then puts every
            # insertion last: 112A landed after 114, being five characters
            # against four.
            "likely_contact_imgt": {k: sorted(v, key=position_order)
                                    for k, v in contacts.items()},
            "cutoff": CONTACT_CUTOFF,
            "cutoff_basis": "Author-recommended classifier cut-off (Paragraph README).",
            "keyed_by": "one-based sequence position where mappable, otherwise IMGT number",
        },
        "limitation": LIMITATION,
    }


def imgt_position_map(chains: dict[str, dict]) -> dict[str, dict[str, int]]:
    """Build the IMGT-number to sequence-position map Paragraph output needs.

    ``chains`` follows the analysis record shape. Paragraph reports bare IMGT
    numbers, so the leading chain-type letter this project uses in residue
    labels is stripped.
    """
    out: dict[str, dict[str, int]] = {}
    for label, chain in chains.items():
        numbering = (chain.get("numbering_imgt") or chain.get("numbering")
                     if isinstance(chain, dict) else None)
        if not numbering or numbering.get("scheme") != "imgt":
            continue
        mapping: dict[str, int] = {}
        for residue in numbering.get("residues", []):
            number = str(residue.get("number", ""))
            bare = number.lstrip("HKL") or number
            mapping[bare] = residue["sequence_position"]
        out[label] = mapping
    return out
