"""Detailed local AntiBERTy outputs, retaining the installed Runner's PLL math.

The Runner's aggregate PLL API omits residue values. We use its existing
tokenizer, model and forward method, collecting each masked target only once.
"""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import math
from pathlib import Path

from .numbering import number_chain


INTERPRETATION = (
    "Masked log-probability asks how well AntiBERTy predicts the observed residue when that position is hidden. "
    "Less-negative values indicate greater sequence-model compatibility. The chain mean is the official "
    "AntiBERTy pseudo log-likelihood; the sum depends on sequence length. Pseudo-perplexity is exp(-mean), "
    "so lower means greater model compatibility. Lowest-decile markers identify relatively unusual positions "
    "within this chain (ties included), not validated engineering liabilities. Classification probabilities are "
    "model softmax outputs, not measured ancestry, percent humanness, PTM, immunogenicity or developability risk. "
    "Unusual CDR residues may be required for binding; these outputs do not recommend substitutions."
)


def _checkpoint_provenance(runner) -> dict:
    path = Path(runner.checkpoint_path)
    artifacts = []
    for filename in ("config.json", "model.safetensors"):
        file = path / filename
        if not file.is_file():
            raise ValueError("Detailed AntiBERTy analysis requires a complete local checkpoint.")
        digest = hashlib.sha256()
        with file.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        artifacts.append({"filename": filename, "sha256": digest.hexdigest(), "size_bytes": file.stat().st_size})
    return {"package": "antiberty", "version": importlib.metadata.version("antiberty"),
            "checkpoint_directory": str(path.resolve()), "checkpoint_artifacts": artifacts,
            "execution": "Local inference", "method": "Installed Runner tokenizer and model; one masked pass per residue, special-token logits excluded as in Runner.pseudo_log_likelihood"}


def _summarize_residues(sequence: str, log_probabilities: list[float], numbering: dict) -> tuple[list[dict], list[dict]]:
    if len(sequence) != len(log_probabilities) or not sequence:
        raise ValueError("Residue likelihoods must cover every input residue.")
    if any(not math.isfinite(value) or value > 0 for value in log_probabilities):
        raise ValueError("Invalid residue log-probability returned by AntiBERTy.")
    ordered = sorted(log_probabilities)
    cutoff = ordered[max(0, math.ceil(len(ordered) * 0.10) - 1)]
    numbered = {r["sequence_position"]: r for r in numbering.get("residues", [])}
    residues, groups = [], {}
    for position, (aa, value) in enumerate(zip(sequence, log_probabilities), 1):
        entry = numbered.get(position, {})
        region = entry.get("region", "Unassigned")
        item = {"sequence_position": position, "residue": aa,
                "number": entry.get("number", str(position)), "region": region,
                "masked_log_probability": value, "observed_probability": math.exp(value),
                "rank_within_chain": 1 + sum(other < value for other in log_probabilities),
                "lowest_decile": value <= cutoff}
        residues.append(item)
        groups.setdefault(region, []).append(item)
    regions = []
    for region, members in groups.items():
        values = [r["masked_log_probability"] for r in members]
        mean = sum(values) / len(values)
        regions.append({"region": region, "residue_count": len(members),
                        "mean_log_likelihood": mean, "sum_log_likelihood": sum(values),
                        "pseudo_perplexity": math.exp(-mean),
                        "lowest_decile_count": sum(r["lowest_decile"] for r in members),
                        "minimum_log_likelihood": min(values)})
    return residues, regions


def analyze_antiberty(runner, sequences: dict[str, str], scheme: str = "imgt", batch_size: int = 16) -> dict:
    """Return detailed results per chain, using an already-loaded local runner.

    Caller owns the CPU/thread lock. No model loading, download or network call
    occurs here. This adapter intentionally fails if the installed API changes.
    """
    import torch

    if not sequences or batch_size < 1 or scheme not in {"imgt", "kabat"}:
        raise ValueError("Provide sequences, a positive batch size, and IMGT or Kabat numbering.")
    if any(not sequence or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence)
           for sequence in sequences.values()):
        raise ValueError("Detailed AntiBERTy analysis requires standard amino acids without masks or gaps.")
    module = importlib.import_module(type(runner).__module__)
    mappings = {"species": module.LABEL_TO_SPECIES, "chain": module.LABEL_TO_CHAIN,
                "graft": module.LABEL_TO_GRAFT}
    provenance = _checkpoint_provenance(runner)
    classification, _, _ = runner._forward(list(sequences.values()))
    classes = {}
    for key, mapping in mappings.items():
        logits = getattr(classification, f"{key}_logits")
        if logits.shape[-1] != len(mapping):
            raise ValueError("AntiBERTy classifier labels do not match the installed model output.")
        classes[key] = torch.softmax(logits, dim=-1).cpu().tolist()
    del classification

    result = {}
    for chain_index, (label, sequence) in enumerate(sequences.items()):
        tokens = runner.tokenizer.tokenize(sequence)
        if tokens != list(sequence):
            raise ValueError("AntiBERTy tokenization did not preserve one token per input residue.")
        targets = torch.tensor(runner.tokenizer.encode(tokens)[1:-1], device=runner.device)
        log_probabilities = []
        for start in range(0, len(tokens), batch_size):
            end = min(start + batch_size, len(tokens))
            masked = [tokens[:i] + [module.MASK_CHAR] + tokens[i + 1:] for i in range(start, end)]
            output, _, _ = runner._forward(masked)
            # Select only the masked target, before discarding full batch logits.
            target_logits = output.prediction_logits[
                torch.arange(end - start, device=runner.device),
                torch.arange(start + 1, end + 1, device=runner.device)
            ].clone()
            target_logits[:, runner.tokenizer.all_special_ids] = -float("inf")
            values = torch.log_softmax(target_logits, dim=-1).gather(1, targets[start:end, None]).squeeze(1)
            log_probabilities.extend(values.cpu().tolist())
            del output, target_logits, values
        numbering = number_chain(sequence, scheme)
        residues, regions = _summarize_residues(sequence, log_probabilities, numbering)
        mean = sum(log_probabilities) / len(log_probabilities)
        probabilities = {key: {mapping[i]: float(value) for i, value in enumerate(classes[key][chain_index])}
                         for key, mapping in mappings.items()}
        result[label] = {
            "pseudo_log_likelihood": mean, "mean_log_likelihood": mean,
            "sum_log_likelihood": sum(log_probabilities), "pseudo_perplexity": math.exp(-mean),
            "species": max(probabilities["species"], key=probabilities["species"].get),
            "chain_class": max(probabilities["chain"], key=probabilities["chain"].get),
            "graft_class": max(probabilities["graft"], key=probabilities["graft"].get),
            "class_probabilities": probabilities, "residue_scores": residues,
            "region_summaries": regions, "scheme": scheme, "residue_count": len(sequence),
            "interpretation": INTERPRETATION, "provenance": provenance,
            "residue_flag_definition": "Lowest 10% of observed-residue masked log-probabilities within this chain; ties included. This is a relative model-surprisal marker, not a therapeutic-risk threshold.",
        }
    return result
