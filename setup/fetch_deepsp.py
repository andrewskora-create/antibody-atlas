#!/usr/bin/env python3
"""Fetch the DeepSP model and convert it for numpy inference.

DeepSP is MIT-licensed, so its weights and training data are redistributable.
They are fetched rather than committed because the Keras HDF5 files are 17 MB
and the repository keeps large model weights out of git, exactly as it does for
AntiBERTy.

Three things happen here:

  1. The three Keras architecture JSONs and HDF5 weight files are downloaded
     from the upstream repository and their sha256 recorded.
  2. Each model is converted to a .npz holding a small execution plan plus the
     weight arrays, so inference needs numpy and nothing else.  The conversion
     is verified against the sixteen reference molecules DeepSP publishes with
     the model: if any descriptor disagrees by more than 1e-3 relative, the
     conversion is rejected and nothing is written.
  3. The DeepSP training labels (20,530 Fv pairs, values from molecular-dynamics
     simulation) are reduced to percentile breakpoints so a predicted value can
     be placed in a distribution rather than reported as a bare number.

    python setup/fetch_deepsp.py             # download, convert, verify
    python setup/fetch_deepsp.py --verify    # re-check what is installed

h5py is needed here and nowhere else:  pip install h5py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = ROOT / "models" / "deepsp"
REFERENCE_PATH = ROOT / "data" / "deepsp_reference.json"
PROVENANCE_PATH = MODEL_DIR / "provenance.json"

# Pinned to a commit so a rebuild years from now produces the same model.
COMMIT = "main"
RAW = f"https://raw.githubusercontent.com/Lailabcode/DeepSP/{COMMIT}/"

HEADS = {
    "SAP_pos": ("DeepSP_models/Conv1D_regressionSAPpos.json",
                "DeepSP_models/Conv1D_regression_SAPpos.h5",
                "data/Deep_SAPpos_data.txt"),
    "SCM_neg": ("DeepSP_models/Conv1D_regressionSCMneg.json",
                "DeepSP_models/Conv1D_regression_SCMneg.h5",
                "data/Deep_SCMneg_data.txt"),
    "SCM_pos": ("DeepSP_models/Conv1D_regressionSCMpos.json",
                "DeepSP_models/Conv1D_regression_SCMpos.h5",
                "data/Deep_SCMpos_data.txt"),
}

# Published outputs for sixteen molecules, used to verify the conversion.
REFERENCE_ALIGNED = "seq_aligned_HL2.txt"
REFERENCE_DESCRIPTORS = "DeepSP_descriptors_anarci2.csv"

LICENCE = ("DeepSP is distributed under the MIT License, Copyright (c) 2024 "
           "Lailabcode. See https://github.com/Lailabcode/DeepSP/blob/main/LICENSE")

CITATION = ("Kalejaye L, Wu I-E, Terry T, Lai P-K. DeepSP: Deep learning-based "
            "spatial properties to predict monoclonal antibody stability. "
            "Computational and Structural Biotechnology Journal, 2024. "
            "doi:10.1016/j.csbj.2024.05.029")

SUPPORTED_LAYERS = {"InputLayer", "Dropout", "Conv1D", "Dense",
                    "BatchNormalization", "MaxPooling1D", "Flatten"}


def download(path: str) -> bytes:
    url = RAW + path
    request = urllib.request.Request(url, headers={"User-Agent": "antibody-atlas/deepsp"})
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read()


def build_plan(arch: dict, weights) -> tuple[list[dict], dict[str, np.ndarray]]:
    """Flatten a Keras Sequential config into an ordered list of operations."""
    plan: list[dict] = []
    arrays: dict[str, np.ndarray] = {}
    for layer in arch["config"]["layers"]:
        cls, cfg = layer["class_name"], layer["config"]
        if cls not in SUPPORTED_LAYERS:
            raise ValueError(f"DeepSP architecture changed: unsupported layer {cls}.")
        name = cfg.get("name")
        if cls in ("InputLayer", "Dropout"):
            # Dropout is the identity at inference time.
            continue
        if cls in ("Conv1D", "Dense"):
            if cls == "Conv1D":
                if cfg.get("padding") != "valid" or list(cfg.get("strides", [1])) != [1]:
                    raise ValueError("DeepSP Conv1D is no longer valid/stride-1.")
            arrays[f"{name}.kernel"] = weights[f"{name}/{name}/kernel:0"][()]
            arrays[f"{name}.bias"] = weights[f"{name}/{name}/bias:0"][()]
            plan.append({"op": cls, "name": name,
                         "activation": cfg.get("activation", "linear")})
        elif cls == "BatchNormalization":
            if list(cfg.get("axis", [2])) not in ([2], [-1]):
                raise ValueError("DeepSP BatchNormalization axis changed.")
            for part in ("gamma", "beta", "moving_mean", "moving_variance"):
                arrays[f"{name}.{part}"] = weights[f"{name}/{name}/{part}:0"][()]
            plan.append({"op": cls, "name": name, "epsilon": float(cfg.get("epsilon", 1e-3))})
        elif cls == "MaxPooling1D":
            if cfg.get("padding") != "valid":
                raise ValueError("DeepSP MaxPooling1D padding changed.")
            plan.append({"op": cls, "pool_size": int(cfg["pool_size"][0]),
                         "strides": int(cfg["strides"][0])})
        elif cls == "Flatten":
            plan.append({"op": cls})
    return plan, arrays


def quantiles(values: np.ndarray) -> list[float]:
    return [round(float(v), 6) for v in np.percentile(values, np.arange(0, 101))]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verify", action="store_true",
                        help="re-check installed weights instead of downloading")
    args = parser.parse_args()

    sys.path.insert(0, str(ROOT))
    from app import deepsp  # noqa: E402

    if args.verify:
        if not deepsp.available():
            print("DeepSP weights are not installed.")
            return 1
        return verify(deepsp)

    try:
        import h5py
    except ImportError:
        print("This script needs h5py:  pip install h5py", file=sys.stderr)
        return 2

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    digests: dict[str, str] = {}
    training: dict[str, np.ndarray] = {}

    for head, (arch_path, weight_path, data_path) in HEADS.items():
        print(f"fetching {head} ...", flush=True)
        arch_bytes = download(arch_path)
        weight_bytes = download(weight_path)
        data_bytes = download(data_path)
        for path, blob in ((arch_path, arch_bytes), (weight_path, weight_bytes),
                           (data_path, data_bytes)):
            digests[path] = hashlib.sha256(blob).hexdigest()

        arch = json.loads(arch_bytes)
        with h5py.File(io.BytesIO(weight_bytes), "r") as weights:
            plan, arrays = build_plan(arch, weights)
        np.savez_compressed(MODEL_DIR / f"{head}.npz", plan=json.dumps(plan), **arrays)
        print(f"  wrote {MODEL_DIR / f'{head}.npz'} "
              f"({(MODEL_DIR / f'{head}.npz').stat().st_size / 1e6:.1f} MB)")

        rows = []
        for line in data_bytes.decode().splitlines():
            parts = line.split()
            if len(parts) != 12:
                continue
            rows.append([float(v) for v in parts[2:]])
        training[head] = np.asarray(rows, dtype=float)
        print(f"  training labels: {training[head].shape[0]} molecules")

    deepsp._load.cache_clear()

    counts = {head: int(arr.shape[0]) for head, arr in training.items()}
    if len(set(counts.values())) != 1:
        print(f"DeepSP training files disagree on molecule count: {counts}", file=sys.stderr)
        return 3

    reference = {
        "source": "DeepSP training labels",
        "description": ("Molecular-dynamics reference values for 20,530 natural human "
                        "Fv pairs drawn from the Observed Antibody Space, published as "
                        "DeepSP's training data. This is a natural-repertoire "
                        "distribution, NOT a clinical-stage benchmark: a molecule in "
                        "its upper tail is unusual among natural antibodies, which is "
                        "not the same as being unusual among approved drugs."),
        "caveat": ("Percentiles place a PREDICTED value against a distribution of "
                   "SIMULATED values. DeepSP's own regression error is not included "
                   "in the percentile."),
        "n": next(iter(counts.values())),
        "licence": LICENCE,
        "citation": CITATION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "quantiles": {f"{head}_{region}": quantiles(training[head][:, i])
                      for head in HEADS for i, region in enumerate(deepsp.REGIONS)},
    }
    REFERENCE_PATH.write_text(json.dumps(reference, indent=2) + "\n")
    print(f"wrote {REFERENCE_PATH}")

    PROVENANCE_PATH.write_text(json.dumps({
        "upstream": "https://github.com/Lailabcode/DeepSP",
        "commit_or_ref": COMMIT,
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "licence": LICENCE,
        "citation": CITATION,
        "sha256": digests,
        "converted_by": "setup/fetch_deepsp.py",
        "conversion_note": ("Keras HDF5 weights flattened to an ordered execution plan "
                            "plus float32 arrays so inference runs in numpy. Dropout "
                            "layers are dropped (identity at inference)."),
    }, indent=2) + "\n")
    print(f"wrote {PROVENANCE_PATH}")

    deepsp._reference.cache_clear()
    return verify(deepsp)


def verify(deepsp) -> int:
    """Check the converted models reproduce DeepSP's published descriptors."""
    print("verifying against DeepSP's published reference molecules ...")
    try:
        aligned_blob = download(REFERENCE_ALIGNED).decode()
        descriptor_blob = download(REFERENCE_DESCRIPTORS).decode()
    except Exception as exc:  # pragma: no cover - network
        print(f"  could not fetch reference data ({exc}); skipping verification.")
        return 0

    names, sequences = [], []
    for line in aligned_blob.splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[1]) == deepsp.ALIGNED_LENGTH:
            names.append(parts[0])
            sequences.append(parts[1])
    published = {row["Name"]: row for row in csv.DictReader(io.StringIO(descriptor_blob))}

    predicted = deepsp.predict_aligned(sequences)
    worst, worst_key = 0.0, None
    compared = 0
    for name, values in zip(names, predicted):
        row = published.get(name)
        if not row:
            continue
        for key, got in values.items():
            want = float(row[key])
            relative = abs(got - want) / max(abs(want), 1e-6)
            compared += 1
            if relative > worst:
                worst, worst_key = relative, f"{name}/{key}"

    if not compared:
        print("  no overlapping reference molecules found; cannot verify.")
        return 1
    print(f"  compared {compared} descriptor values across {len(names)} molecules")
    print(f"  worst relative deviation: {worst:.2e} ({worst_key})")
    if worst > 1e-3:
        print("  FAILED: conversion does not reproduce the published output.", file=sys.stderr)
        return 4
    print("  OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
