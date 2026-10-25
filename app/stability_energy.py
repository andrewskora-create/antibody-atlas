"""Local ThermoMPNN folding ΔΔG prediction, separate from local relaxation energies.

The model uses the parent chain backbone. It does not calculate binding energy,
full-antibody thermodynamics, or an energy difference between repacked structures.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
from collections import OrderedDict

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "370f76ec62bd929f7425e311d8df04a0d094990f"
METHOD = "ThermoMPNN"
SIGN = "ΔΔGfold = ΔGfold(mutant) − ΔGfold(parent); negative predicts stabilization, positive destabilization."
LIMITATIONS = [
    "Prediction from the parent backbone and sequence, not a thermodynamic simulation or experimental measurement.",
    "Single-chain folding stability only: excludes antibody chain-interface and antigen-binding effects.",
    "Not calibrated specifically for antibodies; small differences are not evidence of preserved folding.",
    "Local mutant geometry and hydrogen bonds are assessed separately; they do not determine this ΔΔG.",
    "Does not account explicitly for pH, glycosylation, disulfide connectivity, or formulation conditions.",
    "Single substitutions only; values must not be added to predict combinations or epistasis.",
]
AA3 = dict(zip("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(), "ACDEFGHIKLMNPQRSTVWY"))
AA = "ACDEFGHIKLMNPQRSTVWY"
_CACHE: OrderedDict = OrderedDict()
_LOCK = threading.Lock()


def _vendor() -> Path:
    return Path(os.environ.get("ATLAS_THERMOMPNN_DIR", str(ROOT / "vendor" / "thermompnn"))).resolve()


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stability_capability() -> dict:
    try:
        _verified_manifest()
    except (OSError, ValueError) as exc:
        return {"available": False, "method": METHOD, "local": True, "reason": str(exc)}
    return {"available": True, "method": METHOD, "local": True, "reason": None}


def _verified_manifest() -> dict:
    # The setup script is our pinned trust root, not an editable vendor manifest.
    import importlib.util
    spec = importlib.util.spec_from_file_location("atlas_stability_setup", ROOT / "scripts" / "setup_stability.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for filename, expected in module.FILES.items():
        path = _vendor() / filename
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"Missing or modified ThermoMPNN file: {filename}; run scripts/setup_stability.py.")
    return {"commit": COMMIT, "checkpoint_sha256": module.FILES["models/thermoMPNN_default.pt"],
            "backbone_model_sha256": module.FILES["vanilla_model_weights/v_48_020.pt"],
            "repository": "https://github.com/Kuhlman-Lab/ThermoMPNN", "license": "MIT",
            "paper": "https://doi.org/10.1073/pnas.2314853121"}


def _parse_chain(text: str, chain: str) -> dict:
    """Keep exact PDB residue identifiers; fail on incomplete/disconnected chains.

    Numbering gaps are fine when peptide bonds are continuous (e.g. IMGT).
    Insertion codes retain their own positions. Only the first MODEL is used.
    """
    residues: OrderedDict = OrderedDict()
    ended = False
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("TER") and len(line) > 21 and line[21] == chain and residues:
            ended = True
        if line[:6] not in ("ATOM  ", "HETATM") or len(line) < 54 or line[21] != chain:
            continue
        name = line[17:20].strip()
        if line[:6] == "HETATM" and name not in AA3:
            continue  # Solvent/ligands are outside single-chain model scope.
        if ended:
            raise ValueError("Selected chain has multiple TER-separated segments.")
        if name not in AA3:
            raise ValueError("Selected chain contains an unsupported noncanonical residue.")
        rid = str(int(line[22:26])) + line[26].strip()
        residue = residues.setdefault(rid, {"id": rid, "aa": AA3[name], "atoms": {}})
        if residue["aa"] != AA3[name]:
            raise ValueError("Ambiguous residue identity in alternate conformations.")
        atom = line[12:16].strip()
        if atom not in ("N", "CA", "C", "O"):
            continue
        altloc = line[16]
        occupancy = float(line[54:60].strip() or 0.0) if len(line) >= 60 else 0.0
        priority = (altloc == " ", occupancy, altloc == "A")
        xyz = [float(line[start:start+8]) for start in (30, 38, 46)]
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError("Nonfinite backbone coordinate.")
        previous = residue["atoms"].get(atom)
        if previous is None or priority > previous[0]:
            residue["atoms"][atom] = (priority, xyz)
    if not residues:
        raise ValueError("Selected chain was not found in the parent structure.")
    if len(residues) > 1500:
        raise ValueError("Selected chain exceeds the local stability-model size limit.")
    rows = list(residues.values())
    for i, row in enumerate(rows):
        if set(row["atoms"]) != {"N", "CA", "C", "O"}:
            raise ValueError(f"Incomplete backbone at {chain}:{row['id']}; stability not calculated.")
        if i:
            distance = math.dist(rows[i-1]["atoms"]["C"][1], row["atoms"]["N"][1])
            if not 0.8 <= distance <= 2.2:
                raise ValueError("Discontinuous backbone in selected chain; stability not calculated.")
    return {"sequence": "".join(r["aa"] for r in rows), "residue_ids": [r["id"] for r in rows],
            "coordinates": {atom: [r["atoms"][atom][1] for r in rows] for atom in ("N", "CA", "C", "O")}}


def _run_model(parsed: dict) -> dict:
    env = dict(os.environ)
    env.update({"WANDB_MODE": "disabled", "OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2", "PYTHONDONTWRITEBYTECODE": "1"})
    result = subprocess.run([sys.executable, str(ROOT / "scripts" / "thermompnn_infer.py"), str(_vendor())],
                            input=json.dumps(parsed), text=True, capture_output=True, timeout=180, env=env)
    if result.returncode:
        raise RuntimeError("Local ThermoMPNN inference failed; verify the model installation and PyTorch runtime.")
    value = json.loads(result.stdout)
    if value.get("sequence") != parsed["sequence"] or len(value.get("predictions", [])) != len(parsed["sequence"]):
        raise ValueError("Stability model output does not match the input chain.")
    if any(set(row) != set(AA) or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in row.values()) for row in value["predictions"]):
        raise ValueError("Stability model returned incomplete or nonfinite predictions.")
    return value


def predict_stability(pdb_path, chain, pdb_residue, original, replacement) -> dict:
    """Predict one single substitution. Always return a result; never break repack."""
    answer = {"status": "unavailable", "value_kcal_mol": None, "method": METHOD,
              "sign_convention": SIGN, "limitations": list(LIMITATIONS), "provenance": {},
              "basis": "parent_structure_backbone", "units": "kcal/mol"}
    try:
        if not isinstance(chain, str) or len(chain) != 1 or original not in AA or replacement not in AA or len(original) != 1 or len(replacement) != 1:
            raise ValueError("A single chain and canonical parent/mutant residues are required.")
        rid = str(pdb_residue).strip()
        if not re.fullmatch(r"-?\d+[A-Za-z]?", rid):
            raise ValueError("An exact PDB residue identifier, including insertion code, is required.")
        rid = str(int(re.match(r"-?\d+", rid)[0])) + (rid[-1] if rid[-1].isalpha() else "")
        raw = Path(pdb_path).read_bytes()
        parsed = _parse_chain(raw.decode("utf-8"), chain)
        if rid not in parsed["residue_ids"]:
            raise ValueError("Requested PDB residue was not found in the selected chain.")
        index = parsed["residue_ids"].index(rid)
        if parsed["sequence"][index] != original:
            raise ValueError("Parent residue does not match the supplied structure.")
        with _LOCK:
            provenance = _verified_manifest()
            digest = hashlib.sha256(raw).hexdigest()
            key = (digest, chain, provenance["checkpoint_sha256"], str(_vendor()))
            if key not in _CACHE:
                _CACHE[key] = _run_model(parsed)
                while len(_CACHE) > 8:
                    _CACHE.popitem(last=False)
            else:
                _CACHE.move_to_end(key)
            prediction = _CACHE[key]["predictions"][index][replacement]
            runtime = _CACHE[key].get("torch_version", "unknown")
        provenance.update({"input_sha256": digest, "chain": chain, "pdb_residue": rid,
                           "parent_amino_acid": original, "replacement": replacement,
                           "chain_residue_count": len(parsed["sequence"]), "input_structure": "parent",
                           "execution": "local_cpu", "torch_version": runtime, "model_input": "selected_chain_N_CA_C_O_backbone_and_sequence"})
        answer.update(status="available", value_kcal_mol=round(prediction, 3), provenance=provenance)
        if original == "C" or replacement == "C":
            answer["limitations"].append("Cysteine substitutions require special caution: training-assay artifacts and disulfide changes are not reliably modeled.")
    except Exception as exc:
        answer["reason"] = str(exc) if isinstance(exc, (ValueError, RuntimeError)) else "Local stability calculation unavailable: " + type(exc).__name__
    return answer
