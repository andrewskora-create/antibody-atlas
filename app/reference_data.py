"""Resolution of reference datasets that this project does not redistribute.

Two of the bundled reference corpora are derived from Oxford Thera-SAbDab, which
publishes no licence grant (see DATA_PROVENANCE.md). They are therefore built on
the user's own machine from the upstream CSV rather than shipped in this
repository. This module gives every consumer a single, uniform failure mode so a
fresh clone explains what to run instead of raising a bare FileNotFoundError.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

# name -> (filename, what it is, how to build it)
OPTIONAL_DATASETS: dict[str, tuple[str, str, str]] = {
    "clinical_frameworks": (
        "clinical_frameworks.json",
        "clinical-stage framework and CDR library used for positional precedent",
        "python setup/build_reference_data.py --clinical",
    ),
    "therapeutic_ptm_sequences": (
        "therapeutic_ptm_sequences.json",
        "therapeutic reference chains for documented-PTM site matching",
        "python setup/build_reference_data.py --ptm-sequences",
    ),
    "profiler_comparators": (
        "profiler_comparators.json",
        "protocol-matched therapeutic reference panel",
        "python setup/build_reference_data.py --panel",
    ),
}


class ReferenceDataMissing(RuntimeError):
    """A non-redistributed reference dataset has not been built locally yet."""

    def __init__(self, key: str) -> None:
        filename, description, command = OPTIONAL_DATASETS[key]
        self.key = key
        self.filename = filename
        self.command = command
        super().__init__(
            f"The {description} is not present at data/{filename}.\n"
            f"This project does not redistribute it because its upstream source "
            f"(Oxford Thera-SAbDab) states no licence grant. Build it locally with:\n"
            f"    {command}\n"
            f"See DATA_PROVENANCE.md for the source, terms and verification hashes."
        )


def dataset_path(key: str) -> Path:
    return DATA_DIR / OPTIONAL_DATASETS[key][0]


def available(key: str) -> bool:
    return dataset_path(key).is_file()


def missing_datasets() -> list[str]:
    return [k for k in OPTIONAL_DATASETS if not available(k)]


@lru_cache(maxsize=8)
def load(key: str) -> dict:
    """Load an optional dataset, or raise ReferenceDataMissing with instructions."""
    path = dataset_path(key)
    if not path.is_file():
        raise ReferenceDataMissing(key)
    try:
        return json.loads(path.read_text())
    except ValueError as exc:
        raise RuntimeError(
            f"data/{path.name} is present but not valid JSON ({exc}). "
            f"Rebuild it with: {OPTIONAL_DATASETS[key][2]}"
        ) from exc
