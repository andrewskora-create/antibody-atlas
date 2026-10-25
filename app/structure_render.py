"""Render the predicted Fv surface to a PNG without a browser.

The interactive application draws structures with NGL and exports a PNG. A
headless run has no browser, so the report used to ship with an empty box where
the surface view belongs. This renders the same view offline from the predicted
coordinates, so a report produced by `run_panel.py` carries a figure.

The view is oriented paratope-on: the molecule is rotated so the CDR centroid
faces the viewer, because the CDR face is what a developability reader wants to
look at. Colour is relative exposed hydrophobic signal, matching the legend the
report already prints.

This is a depth-sorted sphere rendering, not a solvent-excluded surface. It is
a figure for orientation, never a measurement.
"""
from __future__ import annotations

import math
from io import BytesIO
from pathlib import Path
from typing import Any

from .analysis import HYDROPATHY

# The report's existing legend ramp, cold to hot.
RAMP = ["#28669a", "#6bb4d1", "#e7d779", "#ee9856", "#d64354", "#b91f43"]
CHARGE_RAMP = ["#b91f43", "#e08a9a", "#f0f0ee", "#8fc4dc", "#1f6f9a"]  # acidic -> basic
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.10}
SIDECHAIN_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}


def _parse(pdb_path: str | Path) -> list[dict[str, Any]]:
    atoms = []
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if element == "H":
            continue
        try:
            atoms.append({
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resseq": line[22:27].strip(),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "element": element,
            })
        except ValueError:
            continue
    return atoms


def _residue_signal(atoms, structure_analysis, mode):
    """Per-residue 0-1 signal, from the structural screen where it is available."""
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V",
    }
    exposure: dict[tuple[str, str], float] = {}
    if isinstance(structure_analysis, dict):
        for residue in structure_analysis.get("residues", []):
            value = residue.get("relative_sidechain_sasa")
            if isinstance(value, (int, float)) and math.isfinite(value):
                exposure[(residue.get("chain"), str(residue.get("number", "")).lstrip("HKL"))] = value

    raw: dict[tuple[str, str], float] = {}
    for atom in atoms:
        key = (atom["chain"], atom["resseq"])
        if key in raw:
            continue
        aa = three_to_one.get(atom["resname"])
        if aa is None:
            continue
        rsa = exposure.get(key)
        if rsa is None:
            rsa = 0.5  # exposure unknown: do not pretend it is buried
        if mode == "charge":
            raw[key] = SIDECHAIN_CHARGE.get(aa, 0.0) * max(rsa, 0.05)
        else:
            raw[key] = max(0.0, HYDROPATHY[aa]) / 4.5 * rsa

    if not raw:
        return {}
    if mode == "charge":
        span = max(abs(v) for v in raw.values()) or 1.0
        return {k: 0.5 + v / (2 * span) for k, v in raw.items()}
    # Anchor on a high percentile rather than the single maximum, or one
    # extreme residue flattens everything else to the bottom of the ramp.
    ordered = sorted(raw.values())
    high = ordered[int(0.92 * (len(ordered) - 1))] or (ordered[-1] or 1.0)
    return {k: min(1.0, v / high) for k, v in raw.items()}


def _orient(coords, cdr_keys, keys):
    """Rotate so the CDR face points at the viewer; fall back to PCA."""
    import numpy as np

    centre = coords.mean(axis=0)
    centred = coords - centre
    cdr_mask = np.array([k in cdr_keys for k in keys])
    if cdr_mask.any():
        direction = centred[cdr_mask].mean(axis=0)
    else:
        # No CDR annotation: look down the smallest principal axis.
        direction = np.linalg.svd(centred, full_matrices=False)[2][-1]
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        return centred
    target = direction / norm
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(target, z)
    s, c = np.linalg.norm(v), float(np.dot(target, z))
    if s < 1e-8:
        rotation = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation = np.eye(3) + kmat + kmat @ kmat * ((1 - c) / s**2)
    return centred @ rotation.T


def render_structure_png(
    pdb_path: str | Path,
    structure_analysis: dict | None = None,
    *,
    mode: str = "hydrophobic",
    width_px: int = 1000,
    height_px: int = 620,
) -> bytes | None:
    """PNG bytes of the Fv coloured by exposed hydrophobic signal, or None."""
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError:
        return None

    atoms = _parse(pdb_path)
    if len(atoms) < 50:
        return None

    signal = _residue_signal(atoms, structure_analysis, mode)
    cdr_keys = set()
    if isinstance(structure_analysis, dict):
        for residue in structure_analysis.get("residues", []):
            if str(residue.get("region", "")).startswith("CDR"):
                cdr_keys.add((residue.get("chain"),
                              str(residue.get("number", "")).lstrip("HKL")))

    coords = np.array([a["xyz"] for a in atoms], dtype=float)
    keys = [(a["chain"], a["resseq"]) for a in atoms]
    base = _orient(coords, cdr_keys, keys)

    cmap = LinearSegmentedColormap.from_list(
        "atlas", RAMP if mode != "charge" else CHARGE_RAMP)
    values = np.array([signal.get(k, 0.0) for k in keys], dtype=float)
    radii = np.array([VDW.get(a["element"], 1.6) for a in atoms], dtype=float)

    # Two views: paratope-on, and turned 90 degrees about the vertical axis.
    # One view of a tall Fv leaves most of a landscape figure empty, and the
    # side view shows the face the first one hides.
    angle = math.radians(90.0)
    turn = np.array([[math.cos(angle), 0.0, math.sin(angle)],
                     [0.0, 1.0, 0.0],
                     [-math.sin(angle), 0.0, math.cos(angle)]])
    views = [base, base @ turn.T]

    dpi = 150
    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)
    fig.patch.set_alpha(0.0)

    # One shared scale so the two panels are directly comparable in size.
    span = max(max(np.ptp(v[:, 0]), np.ptp(v[:, 1])) for v in views) / 2 + 3.0
    panel_aspect = (height_px) / (width_px / len(views))

    for index, view in enumerate(views):
        view = view.copy()
        view[:, 0] -= view[:, 0].mean()
        view[:, 1] -= view[:, 1].mean()
        ax = fig.add_axes([index / len(views), 0, 1 / len(views), 1])
        ax.set_xlim(-span, span)
        ax.set_ylim(-span * panel_aspect, span * panel_aspect)
        ax.set_aspect("equal")
        ax.axis("off")

        zs = view[:, 2]
        z_lo, z_span = float(zs.min()), float(np.ptp(zs)) or 1.0
        points_per_angstrom = (width_px / len(views) / dpi * 72.0) / (2 * span)
        order = np.argsort(zs)
        for i in order:
            x, y, z = view[i]
            red, green, blue, _ = cmap(float(np.clip(values[i], 0, 1)))
            depth = 0.45 + 0.55 * ((z - z_lo) / z_span)
            ax.scatter([x], [y], s=(2 * radii[i] * points_per_angstrom) ** 2,
                       c=[(red * depth, green * depth, blue * depth)],
                       linewidths=0.0, edgecolors="none")
        ax.text(0.5, 0.012, "CDR face" if index == 0 else "side view",
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=7.5, color="#5d7180")

    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, transparent=True)
    plt.close(fig)
    return buffer.getvalue()
