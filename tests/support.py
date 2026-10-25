"""Guards for tests that need material a fresh clone does not have.

Two kinds of thing are deliberately absent from the published repository: the
Thera-SAbDab-derived corpora, which carry no licence grant and are rebuilt
locally, and PROPERMAB, which is Regeneron's and is installed opt-in. Tests
covering those features must SKIP rather than fail, because a clone whose suite
reports failures reads as a broken project when it is behaving exactly as
designed.
"""
from __future__ import annotations

import unittest

from app.reference_data import OPTIONAL_DATASETS, dataset_path


def dataset_available(key: str) -> bool:
    if key not in OPTIONAL_DATASETS:
        raise KeyError(f"Unknown optional dataset '{key}'")
    return dataset_path(key).is_file()


def requires_dataset(key: str):
    """Skip unless a non-redistributed reference dataset has been built."""
    _, description, command = OPTIONAL_DATASETS[key]
    return unittest.skipUnless(
        dataset_available(key),
        f"needs the {description}; build it with: {command}")


def propermab_available() -> bool:
    """Use the app's own check: PROPERMAB is vendored, not pip-installed."""
    from app.propermab import availability
    return bool(availability()["available"])


def requires_propermab():
    """Skip unless PROPERMAB is installed; it is academic-licensed and opt-in."""
    return unittest.skipUnless(
        propermab_available(),
        "needs PROPERMAB (Regeneron, academic research use only); "
        "install it with: python setup/fetch_optional.py")


def antiberty_checkpoint():
    """Where AntiBERTy's weights actually are, or None.

    get_weights() is the package's own resolution order -- an explicit path,
    then $ANTIBERTY_WEIGHTS_DIR, then a directory inside the installed
    package, then a Hugging Face download -- and the same one the Runner
    uses. download=False keeps a test gate off the network.
    """
    try:
        from antiberty import get_weights
        from pathlib import Path
        return Path(get_weights(download=False))
    except Exception:
        return None


def requires_antiberty():
    """Skip unless AntiBERTy's checkpoint resolves.

    This used to test for models/antiberty/model.safetensors, a path this
    project never populates: AntiBERTy carries its own checkpoint. The gate
    was therefore false on every machine, and these tests had never run
    anywhere -- including where AntiBERTy was installed and working.
    """
    return unittest.skipUnless(
        antiberty_checkpoint() is not None,
        "needs AntiBERTy weights; they download on first use, or set "
        "ANTIBERTY_WEIGHTS_DIR to a local checkpoint")


def requires_deepsp():
    """Skip unless the converted DeepSP weights are installed.

    DeepSP is MIT, so unlike PROPERMAB nothing stops a user installing it; the
    weights are simply fetched rather than committed because they are 17 MB.
    """
    from app import deepsp
    return unittest.skipUnless(
        deepsp.available(),
        "needs the DeepSP weights; fetch them with: python setup/fetch_deepsp.py")
