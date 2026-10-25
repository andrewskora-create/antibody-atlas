# Publishing this to GitHub

Four steps. The only one that needs judgement is the first.

## 1. Decide what is public

**Publish from an export, never from a working copy.** The export tool is the
gate: it refuses to write a tree containing restricted sequences or molecule
names, drops the components that are not ours to redistribute, and asserts
that the files the app needs at runtime are present.

```bash
python tools/make_oss_export.py --out ../atlas-export
```

If it exits non-zero, do not push — read what it said.

Always regenerate the export for the exact files you intend to release. A folder name is not proof that its current contents passed the gate.

Then confirm with your own eyes, in whichever export you are publishing:

```bash
python -m pytest tests/ -q      # expect passes and skips, zero failures
git init && git add -A && git status
```

In that `git status`, check that **none** of these appear:

- `data/clinical_frameworks.json`, `data/profiler_comparators.json`,
  `data/therapeutic_ptm_sequences.json` — Thera-SAbDab derived, rebuilt locally
- `.venv/`, `models/`, `vendor/`, `native-build/`, `native-source/`
- anything under `reports/`, `results/`, `out/`, `runs/`
- any `.pdb` other than `tests/fixtures/example_fv.pdb`

`.gitignore` already excludes all of it. The point of looking is that
`.gitignore` only works on files it knows about.

## 2. Set the repository identity

Choose the GitHub owner and repository name. For a fork, update `repository-code` in `CITATION.cff` to that exact identity, verify the citation version/date, and add the resulting repository and Releases links to README. Search the exported tree for unresolved placeholders before publishing. Do not invent a URL before the destination is chosen.

## 3. Create the repository and push

Set shell variables `OWNER` and `NAME` to the selected identity. Create it empty on GitHub — no README, no licence, no `.gitignore`, since this
tree already has all three and GitHub's versions would conflict. Then:

```bash
git add -A
git commit -m "Antibody Atlas: local antibody developability screening"
git branch -M main
git remote add origin https://github.com/$OWNER/$NAME.git
git push -u origin main
```

If you have the `gh` CLI, `gh repo create $NAME --public --source=. --push`
does the create and the push together.

## 4. After the first push

- **Check the Actions tab.** `.github/workflows/tests.yml` runs the suite on
  Python 3.11 and 3.12 and separately re-runs the export gate. The first run
  is the one most likely to fail, usually on a dependency that installs
  differently on a clean runner than on your machine. That is worth knowing.
- **Set the repository description and topics.** Suggested topics:
  `antibody`, `developability`, `protein-engineering`, `bioinformatics`,
  `immunology`, `drug-discovery`.
- **Turn on private vulnerability reporting** (Settings → Code security), which
  is what `SECURITY.md` tells people to use.
- **Add a release.** Use a versioned tag and keep `CITATION.cff` in sync. A maintainer can separately connect the repository to Zenodo for archival releases and DOI assignment.

## Downloads for people who do not use Git

Use a GitHub Release with a versioned ZIP as the prominent download, a tar.gz for Linux, and checksums. Build these assets locally without publishing:

```bash
python tools/build_release.py --version 0.1.0-rc4 --out ../release-downloads
```

This reruns the clean export gate, removes optional PROPERMAB-derived panel descriptors, verifies package contents, preserves launcher executable permissions and produces both archives plus a SHA256SUMS file. It does not include installed environments, reference corpora or model weights. A maintainer must still inspect and test the extracted release. The **Prepare release downloads** Actions workflow produces the same assets as a downloadable workflow artifact, without creating or publishing a GitHub Release.

After the clean-install CI checks pass, create a draft GitHub Release and attach the three files. Its opening instructions should be: download ZIP, unzip, install the listed prerequisites once, double-click **Start Antibody Atlas.command**. Explain that initial dependency/model installation needs internet and several minutes. Include the test platform, known gaps, change notes and source version. Publish only after review. Do not call this a one-click standalone app: it is a source package with guided setup.

A later signed/notarized macOS application could remove the prerequisite and script-approval steps. That requires a separate binary packaging and validation effort; it is not delivered by these source archives.

## What stays off GitHub, permanently

Your own molecules. The export tool checks for them by hash on every run, and
`tools/make_oss_export.py` is the only thing standing between the working copy
and a public tree. Run it every time, including for a one-line change — the
cost is thirty seconds and the failure mode is not recoverable.
