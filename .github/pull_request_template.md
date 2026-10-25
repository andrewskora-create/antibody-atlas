## What this changes

<!-- One or two sentences. Link the issue if there is one. -->

## Checks

- [ ] `python -m pytest tests/ -q` passes
- [ ] `python tools/make_oss_export.py --out /tmp/check` passes
- [ ] A fresh clone without the optional components still runs — new tests skip
      rather than fail, using the guards in `tests/support.py`

## If this adds or changes a method

- [ ] Attribution added in `app/method_sources.py` and `app/attribution.py`,
      with the limitations written honestly
- [ ] Licence verified against the upstream source, and recorded in
      `THIRD_PARTY_NOTICES.md` / `DATA_PROVENANCE.md` / `COMPLIANCE.md` §4
- [ ] `USE_TIERS.md` updated if this changes what a commercial user can do
- [ ] Any reimplementation is verified against upstream's own published output
      by a test that fails if they disagree

## Anything a reviewer should be sceptical of

<!-- Where you are least confident. This section is more useful than the
     first one. -->
