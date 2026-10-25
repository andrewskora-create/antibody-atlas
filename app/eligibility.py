"""Enforce Regeneron's academic-only licence for PROPERMAB inside the app.

Documenting a restriction is not the same as honouring it. PROPERMAB is
licensed to "an academic institution or academically employed full-time
personnel only" and prohibits commercial use. This project ships none of it,
but once a user has installed it, nothing stopped the application from
computing and displaying its descriptors.

So no PROPERMAB descriptor is computed or returned until the user affirms
eligibility in the application. The affirmation is recorded locally with a
timestamp and the exact wording, so there is a record of what was agreed and
when. It is not a legal control and does not verify anything about the user;
it is a deliberate stop that makes the restriction visible at the point of use
rather than buried in a licence file.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# Overridable so a test or a shared install can place it elsewhere.
RECORD = Path(os.environ.get('ATLAS_ELIGIBILITY_FILE', ROOT / '.propermab_eligibility.json'))

REQUIRED_PHRASE = "I am eligible under Regeneron's academic licence"

LICENCE_SUMMARY = (
    "PROPERMAB belongs to Regeneron Pharmaceuticals, Inc. and is not part of this project. "
    "Its licence limits use to “an academic institution or academically employed full-time "
    "personnel only”, for academic research purposes, and prohibits commercial use, "
    "integration into software used for further commercialization, and distribution in any "
    "commercial product or online service."
)
CONSEQUENCE = (
    "Affirming records your statement locally with a timestamp and unlocks the PROPERMAB "
    "surface-patch, spatial-charge-map and Heiden descriptors. It grants you no rights: your "
    "licence is with Regeneron, not with this project, and nonprofit intent alone does not "
    "establish eligibility. If you are at a commercial entity, do not affirm — contact "
    "Regeneron. Everything else in this application is permissively licensed and unaffected."
)


def installed() -> bool:
    return (ROOT / 'vendor' / 'propermab').is_dir() and any((ROOT / 'vendor' / 'propermab').iterdir())


def read() -> dict | None:
    try:
        record = json.loads(RECORD.read_text())
    except (OSError, ValueError):
        return None
    return record if record.get('affirmed') else None


def affirmed() -> bool:
    return read() is not None


def record(statement: str, *, source: str = 'application') -> dict:
    """Store an affirmation. Raises ValueError if the wording does not match."""
    cleaned = ' '.join(str(statement or '').split()).rstrip('.')
    if cleaned.casefold() != REQUIRED_PHRASE.casefold():
        raise ValueError(f'To proceed, type exactly: {REQUIRED_PHRASE}')
    entry = {
        'affirmed': True,
        'statement': REQUIRED_PHRASE,
        'affirmed_at_utc': datetime.now(timezone.utc).isoformat(),
        'source': source,
        'licence': 'Regeneron PROPERMAB academic research use only',
        'note': 'Recorded locally by the user of this installation. Not a verification of '
                'eligibility and not a grant of rights.',
    }
    RECORD.write_text(json.dumps(entry, indent=2))
    return entry


def revoke() -> bool:
    """Withdraw the affirmation; PROPERMAB stops being used immediately."""
    if RECORD.is_file():
        RECORD.unlink()
        return True
    return False


def status() -> dict:
    current = read()
    return {
        'installed': installed(),
        'affirmed': current is not None,
        'affirmed_at_utc': (current or {}).get('affirmed_at_utc'),
        'required_phrase': REQUIRED_PHRASE,
        'licence_summary': LICENCE_SUMMARY,
        'consequence': CONSEQUENCE,
        'blocked_reason': None if current else (
            'PROPERMAB descriptors are withheld until academic eligibility is affirmed.'
            if installed() else
            'PROPERMAB is not installed. Install it with setup/fetch_optional.py --propermab '
            'if you are eligible; everything else works without it.'),
    }


def gate() -> str | None:
    """None when PROPERMAB may be used, otherwise the reason it may not."""
    if not installed():
        return 'PROPERMAB is not installed in this environment.'
    if not affirmed():
        return ('PROPERMAB is installed but academic eligibility has not been affirmed in this '
                'application. Open the PROPERMAB section and affirm, or see USE_TIERS.md.')
    return None
