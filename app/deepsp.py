"""DeepSP: spatial aggregation and charge descriptors from the Fv sequence alone.

DeepSP (Kalejaye, Chen, Lai et al. 2024; MIT licence) predicts thirty
structure-derived descriptors -- SAP_pos, SCM_neg and SCM_pos, each over the
six CDRs and four larger regions -- directly from a paired Fv sequence, with no
structure prediction step.  The reference values it was trained on come from
molecular-dynamics simulations of 20,530 natural human Fv pairs.

There are two reasons this sits alongside the PROPERMAB descriptors:

  * PROPERMAB is licensed by Regeneron for academic research only.  DeepSP is
    MIT, so an industry user who cannot accept the PROPERMAB terms still gets
    trained-predictor output for the same physical properties.
  * DeepSP needs no structure.  It runs in milliseconds where the
    structure-based route needs roughly fifteen seconds of folding, which
    matters for batch runs.

What it is not: DeepSP predicts what a molecular-dynamics simulation would have
reported, not what an experiment measured.  It inherits whatever the simulation
protocol got wrong, and its own regression error on top of that.  Treat the
numbers as a fast triage signal, not as a measurement.

Implementation note.  DeepSP ships Keras HDF5 weights.  Rather than take a
TensorFlow dependency for three small 1D convolutional networks, the weights
are converted to .npz at setup time (``setup/fetch_deepsp.py``) and evaluated
here in numpy.  The conversion is checked against the sixteen reference
molecules DeepSP publishes alongside the model; ``tests/test_deepsp.py``
re-runs that check and fails if the agreement degrades.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np

MODEL_DIR = Path(__file__).resolve().parent.parent / "models" / "deepsp"
REFERENCE_PATH = Path(__file__).resolve().parent.parent / "data" / "deepsp_reference.json"

HEADS = ("SAP_pos", "SCM_neg", "SCM_pos")
REGIONS = ("CDRH1", "CDRH2", "CDRH3", "CDRL1", "CDRL2", "CDRL3", "CDR", "Hv", "Lv", "Fv")

# The alignment grid is DeepSP's own, reproduced from its seq_preprocessing step.
# Heavy chain: IMGT 1-128 with the CDR3 insertion ladder 111A-111H and
# 112I-112A spliced in between 111 and 112.  Light chain: IMGT 1-127 flat, with
# no insertion columns at all -- DeepSP does not model light-chain insertions.
_H_GRID: tuple[str, ...] = (
    tuple(str(i) for i in range(1, 112))
    + tuple(f"111{c}" for c in "ABCDEFGH")
    + tuple(f"112{c}" for c in "IHGFEDCBA")
    + ("112",)
    + tuple(str(i) for i in range(113, 129))
)
_L_GRID: tuple[str, ...] = tuple(str(i) for i in range(1, 128))

assert len(_H_GRID) == 145 and len(_L_GRID) == 127

_H_INDEX = {label: i for i, label in enumerate(_H_GRID)}
_L_INDEX = {label: i for i, label in enumerate(_L_GRID)}

ALIGNED_LENGTH = len(_H_GRID) + len(_L_GRID)  # 272

_ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
_ALPHABET_INDEX = {c: i for i, c in enumerate(_ALPHABET)}

# Which grid slots belong to which CDR, used only to decide whether a residue
# that fell off the grid should invalidate the prediction.
_H_CDR_SPAN = ((27, 38), (56, 65), (105, 117))
_L_CDR_SPAN = ((27, 38), (56, 65), (105, 117))


class DeepSPUnavailable(RuntimeError):
    """Raised when the converted DeepSP weights are not installed."""


def available() -> bool:
    return all((MODEL_DIR / f"{head}.npz").exists() for head in HEADS)


def _grid_key(residue: dict[str, Any]) -> str:
    """IMGT label without the chain-type prefix, e.g. 'H111A' -> '111A'."""
    number = str(residue["number"])
    return number[1:] if number[:1].isalpha() else number


def _place(residues: Iterable[dict[str, Any]], grid_size: int,
           index: dict[str, int], cdr_spans) -> tuple[list[str], list[dict[str, Any]]]:
    slots = ["-"] * grid_size
    dropped: list[dict[str, Any]] = []
    for residue in residues:
        key = _grid_key(residue)
        slot = index.get(key)
        if slot is None:
            base = residue.get("number_base")
            in_cdr = any(lo <= base <= hi for lo, hi in cdr_spans) if base else False
            dropped.append({"number": residue["number"], "aa": residue["aa"], "in_cdr": in_cdr})
            continue
        slots[slot] = residue["aa"]
    return slots, dropped


def align(heavy: dict[str, Any], light: dict[str, Any]) -> dict[str, Any]:
    """Project a numbered heavy/light pair onto DeepSP's 272-column grid.

    Returns the aligned string plus any residues that had no column.  A residue
    with no column is dropped, which is what DeepSP's own preprocessing does;
    the difference is that this reports it, because a truncated CDR-H3 makes
    the prediction describe a molecule the user did not submit.
    """
    for label, numbering in (("heavy", heavy), ("light", light)):
        if (numbering or {}).get("scheme") != "imgt":
            raise ValueError(f"DeepSP needs IMGT numbering; {label} chain is "
                             f"{(numbering or {}).get('scheme')!r}.")
    h_slots, h_dropped = _place(heavy["residues"], len(_H_GRID), _H_INDEX, _H_CDR_SPAN)
    l_slots, l_dropped = _place(light["residues"], len(_L_GRID), _L_INDEX, _L_CDR_SPAN)
    aligned = "".join(h_slots + l_slots)
    unknown = sorted({c for c in aligned if c not in _ALPHABET_INDEX})
    dropped = ([dict(d, chain="H") for d in h_dropped]
               + [dict(d, chain="L") for d in l_dropped])
    return {
        "aligned": aligned,
        "dropped": dropped,
        "dropped_from_cdr": [d for d in dropped if d["in_cdr"]],
        "unknown_residues": unknown,
    }


def one_hot(aligned: str) -> np.ndarray:
    if len(aligned) != ALIGNED_LENGTH:
        raise ValueError(f"Aligned sequence must be {ALIGNED_LENGTH} columns, got {len(aligned)}.")
    x = np.zeros((ALIGNED_LENGTH, len(_ALPHABET)), dtype=np.float32)
    try:
        columns = [_ALPHABET_INDEX[c] for c in aligned]
    except KeyError as exc:
        raise ValueError(f"DeepSP cannot encode residue {exc.args[0]!r}.") from exc
    x[np.arange(ALIGNED_LENGTH), columns] = 1.0
    return x


# --- inference -------------------------------------------------------------

def _conv1d(x: np.ndarray, kernel: np.ndarray, bias: np.ndarray) -> np.ndarray:
    """Keras Conv1D with 'valid' padding and stride 1."""
    width, _, out_channels = kernel.shape
    n, length, _ = x.shape
    out_len = length - width + 1
    offsets = np.arange(out_len)[:, None] + np.arange(width)[None, :]
    patches = x[:, offsets, :].reshape(n, out_len, -1)
    return patches @ kernel.reshape(-1, out_channels) + bias


def _max_pool1d(x: np.ndarray, pool: int, stride: int) -> np.ndarray:
    n, length, channels = x.shape
    out_len = (length - pool) // stride + 1
    offsets = np.arange(out_len)[:, None] * stride + np.arange(pool)[None, :]
    return x[:, offsets, :].max(axis=2)


@lru_cache(maxsize=len(HEADS))
def _load(head: str) -> tuple[tuple[dict[str, Any], ...], dict[str, np.ndarray]]:
    path = MODEL_DIR / f"{head}.npz"
    if not path.exists():
        raise DeepSPUnavailable(
            f"DeepSP weights for {head} are not installed. Run "
            f"'python setup/fetch_deepsp.py' to download and convert them.")
    blob = np.load(path, allow_pickle=False)
    plan = tuple(json.loads(str(blob["plan"])))
    arrays = {k: blob[k] for k in blob.files if k != "plan"}
    return plan, arrays


def _forward(plan, arrays: dict[str, np.ndarray], x: np.ndarray) -> np.ndarray:
    for step in plan:
        op = step["op"]
        if op == "Conv1D":
            x = _conv1d(x, arrays[step["name"] + ".kernel"], arrays[step["name"] + ".bias"])
        elif op == "Dense":
            x = x @ arrays[step["name"] + ".kernel"] + arrays[step["name"] + ".bias"]
        elif op == "BatchNormalization":
            name = step["name"]
            variance = arrays[name + ".moving_variance"] + step["epsilon"]
            x = ((x - arrays[name + ".moving_mean"]) / np.sqrt(variance)
                 * arrays[name + ".gamma"] + arrays[name + ".beta"])
            continue
        elif op == "MaxPooling1D":
            x = _max_pool1d(x, step["pool_size"], step["strides"])
            continue
        elif op == "Flatten":
            x = x.reshape(x.shape[0], -1)
            continue
        else:  # pragma: no cover - converter rejects anything else
            raise ValueError(f"Unsupported DeepSP layer {op}.")
        activation = step.get("activation", "linear")
        if activation == "relu":
            x = np.maximum(x, 0.0)
        elif activation != "linear":  # pragma: no cover
            raise ValueError(f"Unsupported activation {activation}.")
    return x


def predict_aligned(aligned: Iterable[str]) -> list[dict[str, float]]:
    """Run all three heads over pre-aligned 272-column strings."""
    batch = list(aligned)
    if not batch:
        return []
    x = np.stack([one_hot(s) for s in batch])
    out: list[dict[str, float]] = [{} for _ in batch]
    for head in HEADS:
        plan, arrays = _load(head)
        values = _forward(plan, arrays, x)
        for row, record in zip(values, out):
            for region, value in zip(REGIONS, row):
                record[f"{head}_{region}"] = float(value)
    return out


# --- reference distribution ------------------------------------------------

@lru_cache(maxsize=1)
def _reference() -> dict[str, Any] | None:
    if not REFERENCE_PATH.exists():
        return None
    with REFERENCE_PATH.open() as handle:
        return json.load(handle)


def _percentile(descriptor: str, value: float) -> float | None:
    reference = _reference()
    if not reference:
        return None
    breaks = reference.get("quantiles", {}).get(descriptor)
    if not breaks:
        return None
    # breaks[i] is the value at the i-th percentile, i = 0..100
    position = float(np.searchsorted(np.asarray(breaks, dtype=float), value, side="right"))
    return max(0.0, min(100.0, position))


#: Descriptors worth putting in front of a user, and what a high value means.
#: Restricted to the whole-Fv and whole-CDR rollups plus CDR-H3, because the
#: per-CDR values are noisier and a reader cannot act on thirty numbers.
HIGHLIGHT = (
    ("SAP_pos_Fv", "Spatial aggregation propensity over the whole Fv", "SAP, whole Fv",
     "aggregation", "Hydrophobic surface available for self-association."),
    ("SAP_pos_CDR", "Spatial aggregation propensity across the CDRs", "SAP, CDRs",
     "aggregation", "Hydrophobic patch concentrated in the binding site."),
    ("SAP_pos_CDRH3", "Spatial aggregation propensity of CDR-H3", "SAP, CDR-H3",
     "aggregation", "CDR-H3 is the usual source of hydrophobic self-association."),
    ("SCM_neg_Fv", "Negative spatial charge map over the whole Fv", "SCM negative, whole Fv",
     "non_specific_binding", "Negative surface charge; relevant to viscosity and clearance."),
    ("SCM_pos_Fv", "Positive spatial charge map over the whole Fv", "SCM positive, whole Fv",
     "non_specific_binding", "Positive surface charge drives polyspecificity and fast clearance."),
    ("SCM_pos_CDR", "Positive spatial charge map across the CDRs", "SCM positive, CDRs",
     "non_specific_binding", "Positive charge in the binding site predicts off-target binding."),
)

#: A molecule is flagged when it sits above this percentile of the DeepSP
#: training distribution.  There is no experimental threshold here -- this is a
#: distributional outlier call, so it is deliberately conservative.
FLAG_PERCENTILE = 95.0


def summarize(values: dict[str, float]) -> dict[str, Any]:
    reference = _reference()
    rows = []
    for key, label, short, risk, note in HIGHLIGHT:
        if key not in values:
            continue
        percentile = _percentile(key, values[key])
        rows.append({
            "descriptor": key,
            "label": label,
            "short_label": short,
            "value": round(values[key], 3),
            "percentile": None if percentile is None else round(percentile, 1),
            "flagged": bool(percentile is not None and percentile >= FLAG_PERCENTILE),
            "risk": risk,
            "note": note,
        })
    return {
        "highlights": rows,
        "flagged": [r for r in rows if r["flagged"]],
        "reference": None if not reference else {
            "n": reference.get("n"),
            "description": reference.get("description"),
            "flag_percentile": FLAG_PERCENTILE,
        },
    }


def analyze(heavy: dict[str, Any], light: dict[str, Any]) -> dict[str, Any]:
    """Full DeepSP pass for one paired Fv.

    Raises DeepSPUnavailable if the weights are not installed; the caller is
    expected to treat that as an absent layer rather than an error.
    """
    if not heavy or not light:
        return {"status": "skipped",
                "reason": "DeepSP was trained on paired heavy/light Fv only."}
    alignment = align(heavy, light)
    if alignment["unknown_residues"]:
        return {"status": "skipped",
                "reason": "Sequence contains residues DeepSP cannot encode: "
                          + ", ".join(alignment["unknown_residues"])}
    values = predict_aligned([alignment["aligned"]])[0]
    result: dict[str, Any] = {
        "status": "ok",
        "descriptors": {k: round(v, 4) for k, v in values.items()},
        "alignment": {
            "dropped": alignment["dropped"],
            "dropped_from_cdr": alignment["dropped_from_cdr"],
        },
    }
    result.update(summarize(values))
    if alignment["dropped_from_cdr"]:
        positions = ", ".join(d["number"] for d in alignment["dropped_from_cdr"])
        result["status"] = "unreliable"
        result["warning"] = (
            f"{len(alignment['dropped_from_cdr'])} CDR residue(s) had no column in "
            f"DeepSP's alignment grid and were dropped ({positions}). The grid "
            f"allows seventeen CDR-H3 insertions and none in the light chain, so "
            f"a long CDR-H3 or any light-chain insertion overflows it. These "
            f"numbers describe a truncated molecule, not the one submitted.")
    elif alignment["dropped"]:
        result["note"] = (
            f"{len(alignment['dropped'])} framework residue(s) fell outside DeepSP's "
            f"alignment grid and were dropped, which is what DeepSP's own "
            f"preprocessing does.")
    return result
