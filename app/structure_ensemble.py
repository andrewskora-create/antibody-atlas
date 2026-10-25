"""Fold with several ABodyBuilder2 models and keep the disagreement visible.

A single predicted structure is one hypothesis presented without an error bar.
Where that matters most here is solvent exposure: a PTM site's priority turns
on whether its side chain is exposed, and a loop that one model places facing
out and another places tucked in flips that call. The user sees a confident
"high exposure" either way.

ABodyBuilder2 ships four models and its authors recommend using the ensemble.
TAP2 makes the same recommendation for the same reason, and resolves it the
same way this does: **flag if any model flags.** A site exposed in one model of
three is a site worth looking at, because the cost of checking it is an assay
and the cost of missing it is a liability in the clinic.

What this does NOT do is average coordinates. An averaged structure is not a
structure; it is a smear that no force field would accept and that FreeSASA
would read as something physically impossible. Each model is analysed on its
own and the per-residue results are reconciled afterwards.

The default remains a single model. This triples fold time, and for a triage
pass over a large panel that is the wrong trade; it is the right trade when a
specific site is about to decide an experiment.
"""
from __future__ import annotations

import hashlib
import importlib.util
from functools import lru_cache
from pathlib import Path
from typing import Any

#: ABodyBuilder2's own ensemble. Which of these is actually usable is decided
#: at runtime, not assumed -- see usable_models().
ALL_MODEL_IDS = (1, 2, 3, 4)

#: Exposure bands, matching structure_analysis. Ordered least to most exposed,
#: so "worst across models" is a max over this ordering.
BAND_ORDER = ("unknown", "low", "intermediate", "high")


def _checkpoint_dir() -> Path | None:
    spec = importlib.util.find_spec("ImmuneBuilder")
    if not spec or not spec.origin:
        return None
    return Path(spec.origin).parent / "trained_model"


@lru_cache(maxsize=2)
def usable_models(kind: str = "antibody") -> tuple[tuple[int, ...], dict[int, str]]:
    """Which checkpoints actually load, and why the others do not.

    This check exists because of a real failure: ``antibody_model_3`` in a
    working installation here was 135 MB where models 2 and 4 are 214 MB -- a
    truncated download that had never been noticed, because nothing used
    anything but model 1. Loading it raises deep inside torch with a miniz
    error that would surface to a user as an unexplained fold failure, or, if
    an ensemble silently skipped it, as a quietly smaller ensemble.
    """
    directory = _checkpoint_dir()
    if directory is None:
        return (), {i: "ImmuneBuilder is not installed" for i in ALL_MODEL_IDS}
    import torch
    good: list[int] = []
    problems: dict[int, str] = {}
    for model_id in ALL_MODEL_IDS:
        path = directory / f"{kind}_model_{model_id}"
        if not path.is_file() or path.stat().st_size == 0:
            problems[model_id] = "checkpoint not downloaded"
            continue
        try:
            torch.load(path, map_location="cpu", weights_only=False)
            good.append(model_id)
        except Exception as exc:
            problems[model_id] = f"checkpoint unreadable ({type(exc).__name__})"
    return tuple(good), problems


def plan(requested: int, kind: str = "antibody") -> dict[str, Any]:
    """Decide which models an ensemble of `requested` size will actually use."""
    good, problems = usable_models(kind)
    chosen = list(good[:max(1, requested)])
    return {
        "requested": requested,
        "model_ids": chosen,
        "available": list(good),
        "unusable": {str(k): v for k, v in problems.items()},
        "short": len(chosen) < requested,
        "note": (f"Asked for {requested} models, {len(chosen)} usable. "
                 + "; ".join(f"model {k}: {v}" for k, v in sorted(problems.items()))
                 if len(chosen) < requested else ""),
    }


def fold_models(sequences: dict[str, str], model_ids: list[int], out_dir: Path,
                stem: str) -> list[dict[str, Any]]:
    """Predict one structure per model id. Returns a record per successful model."""
    paired = "L" in sequences
    if paired:
        from ImmuneBuilder import ABodyBuilder2 as Builder
    else:
        from ImmuneBuilder import NanoBodyBuilder2 as Builder
    out_dir.mkdir(parents=True, exist_ok=True)
    produced: list[dict[str, Any]] = []
    for model_id in model_ids:
        path = out_dir / f"{stem}.model{model_id}.pdb"
        try:
            Builder(model_ids=[model_id]).predict(sequences).save(str(path))
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError("the model produced an empty PDB")
            produced.append({"model_id": model_id, "path": path,
                             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        except Exception as exc:
            produced.append({"model_id": model_id, "path": None,
                             "error": f"{type(exc).__name__}: {str(exc)[:200]}"})
    return produced


def _band(exposure: float | None) -> str:
    """The banding in structure_analysis, applied to a reconciled value."""
    if exposure is None:
        return "unknown"
    return "low" if exposure < 0.10 else "intermediate" if exposure < 0.30 else "high"


def consensus(per_model: list[tuple[int, dict[str, Any]]]) -> dict[str, Any]:
    """Reconcile per-residue exposure across models.

    `per_model` is (model_id, analyze_structure result). Returns the worst-case
    exposure per residue, the spread, and the sites where the models disagree
    about the band -- which is the output worth reading, because a site whose
    band is stable across three independent models is a different kind of
    finding from one that is not.
    """
    if not per_model:
        return {"status": "unavailable", "reason": "No model produced an analysis."}

    by_residue: dict[tuple[str, int], dict[str, Any]] = {}
    for model_id, analysis in per_model:
        for residue in analysis.get("residues", []):
            position = residue.get("sequence_position")
            if position is None:
                continue
            key = (residue["chain"], position)
            entry = by_residue.setdefault(key, {
                "chain": residue["chain"], "sequence_position": position,
                "number": residue.get("number"), "aa": residue.get("aa"),
                "region": residue.get("region"), "values": {}})
            entry["values"][model_id] = residue.get("relative_sidechain_sasa")

    rows = []
    unstable = 0
    for entry in by_residue.values():
        measured = [v for v in entry["values"].values() if isinstance(v, (int, float))]
        if not measured:
            rows.append({**{k: entry[k] for k in
                            ("chain", "sequence_position", "number", "aa", "region")},
                         "exposure_max": None, "exposure_min": None, "spread": None,
                         "band": "unknown", "band_stable": None, "models": len(entry["values"])})
            continue
        high, low = max(measured), min(measured)
        bands = {_band(v) for v in measured}
        stable = len(bands) == 1
        unstable += not stable
        rows.append({
            **{k: entry[k] for k in ("chain", "sequence_position", "number", "aa", "region")},
            # Worst case, not mean: this is the flag-if-any rule, applied at the
            # residue level so everything downstream inherits it.
            "exposure_max": round(high, 3), "exposure_min": round(low, 3),
            "spread": round(high - low, 3),
            "band": _band(high), "band_stable": stable,
            "bands_seen": sorted(bands, key=BAND_ORDER.index),
            "models": len(measured),
        })

    rows.sort(key=lambda r: (-(r["spread"] or 0), r["chain"], r["sequence_position"]))
    model_ids = [model_id for model_id, _ in per_model]
    return {
        "status": "ok",
        "model_ids": model_ids,
        "residues": rows,
        "residues_assessed": sum(r["exposure_max"] is not None for r in rows),
        "band_disagreements": unstable,
        "rule": ("Exposure is the maximum across models, not the mean: a site exposed in "
                 "any model is treated as exposed. This follows ABodyBuilder2's own "
                 "ensemble recommendation and TAP2's flag-if-any convention."),
        "caveat": ("Agreement between models is not accuracy. These are four models trained "
                   "on overlapping data; they can agree and all be wrong, and a stable band "
                   "is not experimental evidence of exposure. The disagreement count is also "
                   "not reproducible run to run: ABodyBuilder2's OpenMM refinement and "
                   "hydrogen placement are stochastic, so the same sequence can give a "
                   "different count. Treat it as indicative of how settled the surface is, "
                   "not as a measurement."),
    }


def reconcile_sites(sites: list[dict[str, Any]], consensus_rows: list[dict[str, Any]]
                    ) -> list[dict[str, Any]]:
    """Apply the ensemble's worst-case exposure to PTM sites from model 1.

    The site list comes from the primary model. Each chemistry is reprioritized
    with the same screening rules at the maximum exposure; the original
    single-model exposure and priority remain available for audit.
    """
    index = {(r["chain"], r["sequence_position"]): r for r in consensus_rows}
    out = []
    for site in sites:
        row = index.get((site.get("chain"), site.get("sequence_position")))
        updated = dict(site)
        if row and row.get("exposure_max") is not None:
            updated["single_model_exposure"] = site.get("relative_sidechain_sasa")
            updated["single_model_exposure_band"] = site.get("exposure_band")
            updated["relative_sidechain_sasa"] = row["exposure_max"]
            # A different exposure can cross a screening threshold. Retaining
            # the original priority would silently defeat the flag-if-any rule.
            if all(key in site for key in ("kind", "motif", "start", "residue", "number", "region")):
                from .structure_analysis import _site_priority, _group_physical_sites
                residue = {"chain": site["chain"], "aa": site["residue"],
                           "number": site["number"], "region": site["region"],
                           "relative_sidechain_sasa": row["exposure_max"],
                           "sidechain_sasa_a2": site.get("sidechain_sasa_a2")}
                members = [_site_priority({**site, **flag}, residue, site.get("scheme", "imgt"))
                           for flag in (site.get("flags") or [site])]
                revised = _group_physical_sites(members)[0]
                updated.update(revised)
                updated["single_model_priority"] = site.get("priority")
                updated["single_model_flags"] = site.get("flags")
                updated["relative_sidechain_sasa"] = row["exposure_max"]
            updated["exposure_band"] = row["band"]
            updated["exposure_models"] = row["models"]
            updated["exposure_spread"] = row["spread"]
            updated["exposure_band_stable"] = row["band_stable"]
            if not row["band_stable"]:
                updated["ensemble_note"] = (
                    f"Models disagree on exposure here: {', '.join(row['bands_seen'])} "
                    f"across {row['models']} models (relative side-chain SASA "
                    f"{row['exposure_min']} to {row['exposure_max']}). The highest band is "
                    f"used. Treat the exposure evidence for this site as weak.")
            elif updated["single_model_exposure_band"] != row["band"]:
                updated["ensemble_note"] = (
                    f"Band changed from {updated['single_model_exposure_band']} to "
                    f"{row['band']} once the ensemble was used.")
        out.append(updated)
    return out
