# Contributing

Thanks for looking. A few things about this project will shape what makes a
good contribution, so they are worth stating before the mechanics.

## What this project is trying to be

A screen that an academic lab can run on its own machine and then defend in a
paper. Two consequences follow, and most review comments trace back to one of
them.

**Every number has to say where it came from.** Each assessment carries its own
attribution line naming the method that produced it, and says explicitly when
the rule is this project's own screening policy rather than a published
threshold. A new feature that emits a number without that is incomplete, even
if the number is correct.

**A missing optional component degrades a section; it does not break a run.**
The clinical corpora are not redistributable, PROPERMAB is academic-only, and
DeepSP is fetched rather than committed. A fresh clone has none of them and
must still screen a molecule end to end, with the absent sections saying what
is absent and printing the command that installs it. Tests for those features
skip, using the guards in `tests/support.py` — they do not fail. A clone whose
suite reports failures reads as a broken project when it is behaving exactly as
designed.

## Before you open a pull request

```bash
python -m pip install -r requirements-dev.txt   # pytest; not a runtime dependency
python -m pytest tests/ -q      # Python regressions
for test in tests/*.cjs; do node "$test"; done  # browser logic regressions (Node.js 22)
python scripts/doctor.py        # what is installed on your machine
python tools/make_oss_export.py --out /tmp/check     # the compliance gate
```

The export step matters even for a small change. It asserts that nothing
restricted leaks into a publishable tree, and that the files the app needs at
runtime are present — a 1.3 MB browser asset once vanished from every export
without a single test noticing, because nothing checked that it arrived.

## Adding a third-party method

This is the most common substantial contribution, and the bar is higher than
"it imports".

1. **Check the licence against the upstream source**, not against what a
   README or a blog post says it is. Record what you actually read. This
   project had FreeSASA recorded as LGPL for a while; it is MIT.
2. **Do not commit weights or corpora you have not cleared.** Fetch them in a
   `setup/` script that records the sha256 of every file it reads, as
   `setup/fetch_deepsp.py` and `setup/fetch_constant_regions.py` do.
3. **Verify any reimplementation against upstream's own published output.**
   If you reimplement, convert or re-host someone's model, there has to be a
   test that fails when your version and theirs disagree. `tests/test_deepsp.py`
   is the pattern: it checks both the numerics and the preprocessing, because
   the preprocessing is where a substitution silently changes the answer.
4. **Add the notice to `THIRD_PARTY_NOTICES.md`, the terms to
   `DATA_PROVENANCE.md`, and the row to `COMPLIANCE.md` §4.** Add it to
   `USE_TIERS.md` if it changes what a commercial user can do.
5. **Write the limitations honestly** in `app/method_sources.py`. The
   `limitations` field is read by users deciding whether to act on a result.
   "Predicts what a simulation would have reported, not what an experiment
   measured" is the kind of sentence that belongs there.

## Style

Comments explain why, not what. If a line looks odd and the reason is a real
constraint — an upstream quirk, a licence, a bug that was once shipped — write
that reason down. Several comments in this codebase exist because something
broke silently once, and they are the main defence against it breaking the same
way again.

No new runtime dependency without a reason that survives the question "can this
be done with what is already here". DeepSP's weights are read in numpy rather
than Keras precisely because three small convolutional networks did not justify
a TensorFlow dependency for every user.

## Things that will not be merged

- Anything that makes a screening heuristic read as a prediction of
  aggregation, clearance, immunogenicity or viscosity.
- An acceptance threshold presented as derived from a reference panel whose
  members were not selected for favourable experimental behaviour. Those panels
  give an observed range, which is not an acceptance interval, and the
  distinction is load-bearing.
- A vendored copy of a component with restricted terms, however convenient.
- Telemetry of any kind.

## Reporting a problem with a result

Please include the input that produced it — or, if you cannot share the
sequence, a public one that reproduces the behaviour. "The tool scored this
wrong" is hard to act on; "ANARCII numbered this position differently from
ANARCI, here is the sequence" is a fix.
