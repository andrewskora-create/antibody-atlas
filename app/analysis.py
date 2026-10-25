"""Transparent sequence-level antibody developability descriptors.

These are screening annotations, not measured stability, immunogenicity, or PTMs.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}


def clean_sequence(raw: str, allow_nucleotide: bool = True) -> str:
    """Accept plain sequence or one FASTA record, rejecting ambiguity.

    A nucleotide sequence is translated first. Without that it would pass
    every check here -- A, C, G and T are amino acids -- and fail later with
    a numbering error that says nothing about the real problem. Use
    ``translation.resolve`` directly when the caller needs to report what was
    translated.
    """
    if allow_nucleotide:
        from .translation import looks_like_nucleotide, to_protein
        if looks_like_nucleotide(raw):
            raw = to_protein(raw)["protein"]
    lines = [line.strip() for line in raw.strip().splitlines() if line.strip()]
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    if any(line.startswith(">") for line in lines):
        raise ValueError("Enter one chain per field; multiple FASTA records are not supported.")
    sequence = re.sub(r"\s+", "", "".join(lines)).upper()
    if not sequence:
        raise ValueError("Enter an amino-acid sequence.")
    if len(sequence) > 600:
        raise ValueError("A chain must contain at most 600 amino acids.")
    bad = sorted(set(sequence) - AA)
    if bad:
        raise ValueError(f"Unsupported amino-acid symbols: {', '.join(bad)}. Use the 20 standard residues.")
    return sequence


def _window_scores(seq: str, size: int = 7) -> list[dict[str, Any]]:
    if len(seq) < size:
        return []
    return [
        {
            "start": i + 1,
            "end": i + size,
            "sequence": seq[i:i + size],
            "mean_hydropathy": round(sum(HYDROPATHY[a] for a in seq[i:i + size]) / size, 3),
            "charged_residues": sum(a in "DEKR" for a in seq[i:i + size]),
        }
        for i in range(len(seq) - size + 1)
    ]


def _site(kind: str, start: int, motif: str, interpretation: str) -> dict[str, Any]:
    return {"kind": kind, "start": start, "end": start + len(motif) - 1,
            "motif": motif, "interpretation": interpretation}


def ptm_sites(seq: str, imgt_numbering: dict | None = None) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    if seq and seq[0] in 'QE':
        sites.append(_site('N-terminal pyroglutamate formation', 1, seq[0],
            'Conditional on this being the free mature chain N terminus; common charge heterogeneity, not automatically an engineering liability. Confirm construct and peptide mapping.'))
    # Count or odd/even parity cannot establish whether a cysteine is unpaired.
    if imgt_numbering and imgt_numbering.get('scheme') == 'imgt':
        canonical = {r['sequence_position'] for r in imgt_numbering.get('residues', [])
                     if r.get('number') in {f"{imgt_numbering.get('chain_type')}{n}" for n in (23, 104)} and r.get('aa') == 'C'}
        for position, aa in enumerate(seq, 1):
            if aa == 'C' and position not in canonical:
                sites.append(_site('Noncanonical cysteine', position, aa,
                    'Outside the two canonical IMGT Cys23/Cys104 sites. May form a beneficial additional disulfide; free-thiol status and pairing require experiment. Do not automatically remove.'))
    for i in range(len(seq) - 2):
        tri = seq[i:i + 3]
        if tri[0] == "N" and tri[1] != "P" and tri[2] in "ST":
            sites.append(_site("N-linked glycosylation sequon", i + 1, tri,
                               "Sequence motif only; occupancy requires experiment."))
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        if pair in {"NG", "NS", "NT", "NQ"}:
            sites.append(_site("Asparagine deamidation motif", i + 1, pair,
                               "Context-dependent chemical liability."))
        if pair == 'DP':
            sites.append(_site('Asp-Pro cleavage motif', i + 1, pair,
                'Potential acid-sensitive peptide bond; assess process pH, hold time and fragmentation experimentally. Exposure alone does not predict cleavage.'))
        if pair in {"DG", "DS", "DT"}:
            sites.append(_site("Aspartate isomerization motif", i + 1, pair,
                               "Context-dependent chemical liability."))
    for i, aa in enumerate(seq, 1):
        if aa in "MW":
            sites.append(_site("Oxidation-susceptible residue", i, aa,
                               "Exposure and formulation determine actual oxidation."))
    return sorted(sites, key=lambda x: (x["start"], x["kind"]))


def sequence_metrics(seq: str, ph: float = 7.4, imgt_numbering: dict | None = None) -> dict[str, Any]:
    if not 0 <= ph <= 14:
        raise ValueError("pH must be between 0 and 14.")
    protein = ProteinAnalysis(seq)
    windows = _window_scores(seq)
    patches = [w for w in windows if w["mean_hydropathy"] >= 1.6 and w["charged_residues"] <= 1]
    counts = Counter(seq)
    return {
        "length": len(seq),
        "molecular_weight_da": round(protein.molecular_weight(), 1),
        "isoelectric_point": round(protein.isoelectric_point(), 2),
        "net_charge_at_ph": round(protein.charge_at_pH(ph), 2),
        "ph": ph,
        "gravy": round(protein.gravy(), 3),
        "aromatic_fraction": round(protein.aromaticity(), 3),
        "hydrophobic_windows": windows,
        "hydrophobic_patches": patches,
        "aggregation_proxy": {
            "patch_count": len(patches),
            "definition": "7-residue windows with mean Kyte–Doolittle hydropathy ≥1.6 and ≤1 D/E/K/R residue; sequence-only screening proxy, not an aggregation prediction",
        },
        "cysteines": counts["C"],
        "ptm_sites": ptm_sites(seq, imgt_numbering),
        "composition": dict(sorted(counts.items())),
    }
