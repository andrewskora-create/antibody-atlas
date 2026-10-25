# Licence and redistribution audit

Audit date: 2026-10-10. Scope: everything in this repository plus every declared
dependency. Every licence below was read from the upstream source on the audit
date, not from recall. Where a licence could not be found, that is recorded as
"none stated" rather than assumed permissive.

This is an engineering audit, not legal advice. Copyright and licence choice are
settled (section 6): the software was developed independently by the author. The
remaining obligations are all third-party ones, plus the pre-push sequence check
in section 6.3.

---

## 1. Summary

| | Count |
|---|---|
| Components cleared for redistribution | 9 |
| Components removed from the published tree | 3 |
| Decisions required before publishing | 2 |
| Proprietary sequence material found in the tree | 0 |

The repository as it stood could not be published as-is. Two components carried
terms incompatible with a public repository, and a third carried GPLv3
obligations that conflicted with a permissive project licence. All three are
needed only by one optional analysis path, so all three now install on the
user's own machine instead of shipping here. Nothing in the scientific core was
lost.

## 2. Proprietary material: clear

The maintainer's own unpublished antibody sequences are covered by a
provisional application that has not published, so the tree is scanned for them
directly before every export. The molecules are not named here: naming them
would leak by description what the scan exists to withhold.

- No unpublished CDR or framework sequence appears anywhere in the tree.
- Antibody names appearing in `run_panel.py` docstrings are placeholders
  (`"QVQL..."`), not sequences.
- The bundled example report embeds no amino-acid sequence; its input is
  identified only by hash.
- Short apparent hits are public germline segments present in
  `human_v_ogrdb.json`, which any near-germline antibody will match. Public
  germline sequence is not proprietary content.

`tools/make_oss_export.py` enforces this mechanically rather than by review. It
holds SHA-256 digests of the restricted fragments, never the fragments
themselves, and slides a window of each recorded length over every amino-acid
run in every file it is about to ship. Absence is therefore proven without the
sequences being present in the checker, which itself ships. The fragment labels
are opaque for the same reason. The export fails rather than warns.

`.gitignore` excludes report and job output directories, so a generated report
containing a real sequence cannot be committed by accident.

## 3. Components removed from the published tree

### 3.1 PROPERMAB — Regeneron, non-commercial, academic-only

`vendor/propermab/` has been removed from the published tree.

The licence restricts the licensee class, not merely the use:

> "You" refers to an academic institution or academically employed full-time personnel only.

and prohibits:

> any integration into other code or software that is used for further commercialization ... distribution of the source code, in whole or in part, or any resulting executables, in any commercial product, or ... use of the source code, in whole or in part, or any resulting executables, in any commercial online service.

Two separate questions follow, and they have different answers.

**Your own use — covered.** The work was performed for Fred Hutchinson Cancer
Center, an academic institution. Academic research use at an academic
institution is what the grant permits. All analysis performed to date stands.
One factual caveat, recorded rather than resolved here: the grant is scoped to
*academic research purposes*, so applying it to a candidate covered by a
provisional application is a question about the purpose of the use rather than
the identity of the institution. That is a licensing judgement for whoever
claims eligibility, not an engineering question, and it does not affect any
result already produced. The in-app gate in app/eligibility.py records an
explicit attestation before the PROPERMAB path runs.

**Redistribution — not permitted here.** Eligibility to use is not eligibility
to publish. A public GitHub repository distributes to everyone, including
commercial entities who are expressly outside the licensee class. The licence
does permit an eligible academic to distribute for academic research purposes,
but a public repository cannot restrict recipients to that class, and the
obligation to pass on the licence terms cannot be discharged by a notice file
that commercial users are free to ignore.

**Resolution.** PROPERMAB is now an opt-in component that the user installs from
Regeneron's own repository. `setup/fetch_optional.py` requires the user to
affirm eligibility, prints the licence, and refuses to proceed without an
explicit acknowledgement. The repository ships no Regeneron code, so your
publication makes no representation about anyone else's eligibility.

**Descriptor values are a separate question.** The reference panel contains
numeric descriptor values that PROPERMAB computed. Those are research outputs,
not the licensed source code, and publishing research outputs is ordinary
practice. A conservative reading is nonetheless available, so the export script
takes `--panel-descriptors permissive-only` to emit a panel with the
PROPERMAB-derived fields dropped. The default keeps them, with provenance
recorded per value; switch to `permissive-only` if you want the conservative
reading.

### 3.2 Thera-SAbDab-derived corpora — no licence grant

Removed from the published tree: `data/clinical_frameworks.json` (12.5 MB, 774
clinical chains), `data/therapeutic_ptm_sequences.json`, and the `sequences`
field of every entry in `data/profiler_comparators.json`.

All three derive from one file:
`TheraSAbDab_SeqStruc_OnlineDownload.csv`, sha256
`cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d`,
retrieved 2026-10-04.

Oxford's Thera-SAbDab pages state **no licence and no terms of use**. They
request citation, say the resource was built "under an open-innovation
agreement" without publishing that agreement's terms, and direct licensing
enquiries to `opig@stats.ox.ac.uk`. Absence of a licence is not permission:
with no grant, there is no default right to redistribute a derived sequence
corpus, and the size of these files makes them a substantial extract rather
than incidental quotation.

**Resolution.** `setup/build_reference_data.py` downloads the CSV from Oxford on
the user's machine and regenerates all three datasets locally, verifying the
download against the recorded sha256 and each rebuilt sequence against its
per-chain hash. The published panel keeps every computed descriptor plus the
per-chain hashes, so a rebuild is provably identical to the one these results
came from. Drug names, isotypes and clinical-stage labels are retained: those
are facts, not expression.

**Recommended, not done without your say-so:** email `opig@stats.ox.ac.uk` and
ask whether redistribution of a derived corpus is permitted. If they say yes in
writing, the datasets can ship directly and setup gets simpler. This is a
one-paragraph email and it is the only route to certainty.

### 3.3 NanoShaper and CGAL — GPLv3

Removed from the published tree: `native-source/nanoshaper/` and
`native-source/cgal42/`.

NanoShaper is GPLv3. CGAL 4.2 is a mix of GPLv3+ and LGPL depending on the
package. Shipping GPLv3 source is permitted, but it would place GPL obligations
on a repository otherwise offered permissively, and would invite the question of
whether the application is a derivative work. Both are needed only by
PROPERMAB's surface-patch path, which is now opt-in, so both leave with it.
`setup/fetch_optional.py` documents the build and leaves it to the user, whose
own GPL obligations then concern only their own machine.

## 4. Components cleared for redistribution

| Component | Licence | Verified at | In repo |
|---|---|---|---|
| ANARCI | BSD-3-Clause | github.com/oxpig/ANARCI | dependency |
| ANARCII | BSD-3-Clause | github.com/oxpig/ANARCII | dependency |
| ImmuneBuilder / ABodyBuilder2 | BSD-3-Clause | github.com/oxpig/ImmuneBuilder | dependency |
| Paragraph (paratope) | BSD-3-Clause, © 2021 L. A. Chinery | github.com/oxpig/Paragraph | dependency |
| AntiBERTy | MIT, © 2023 Jeffrey Ruffolo | bundled `LICENSE.MD` | weights fetched |
| FreeSASA | MIT | github.com/mittinatten/freesasa | dependency |
| DeepSP | MIT, © 2024 Lailabcode | github.com/Lailabcode/DeepSP/blob/main/LICENSE | weights fetched |
| NGL viewer | MIT, © 2014-2017 A. S. Rose | bundled `LICENSE` | bundled |
| OGRDB germline sets | CC0 | ogrdb.airr-community.org | bundled |
| Biopython, FastAPI, ReportLab, PyTorch, OpenMM, PDBFixer | permissive (BSD/MIT/Apache family) | upstream | dependencies |

Note on FreeSASA: it is MIT, not LGPL. An earlier assumption in this project's
notes said LGPL; that was wrong and is corrected here.

Note on AntiBERTy: the 104 MB checkpoint is no longer committed. It exceeds
GitHub's recommended blob size and is fetched at setup. The MIT licence covers
both code and weights, so this is a size decision, not a licence one.

Note on DeepSP: the same reasoning. MIT covers the code, the trained weights
and the published training labels, so all three could be committed. The 17 MB
of Keras checkpoints are fetched at setup instead, and `setup/fetch_deepsp.py`
records the sha256 of every upstream file it reads. What the repository does
carry is derivative and small: the sixteen reference molecules DeepSP publishes
with the model, kept as a test fixture under
`tests/fixtures/deepsp_reference_molecules.json`, and — once built — percentile
breakpoints over the training labels in `data/deepsp_reference.json`, which are
summary statistics rather than the labels. The MIT notice is reproduced in
THIRD_PARTY_NOTICES.md, which is what MIT requires.

The DeepSP weights are converted from Keras HDF5 to numpy at setup. A format
conversion of someone else's trained weights is a derivative work, and MIT
permits it provided the notice travels with it, which it does. No weight value
is altered; the conversion is verified against DeepSP's own published output
before anything is written.

## 5. Literature-derived values — cite, do not mirror

`data/experimental_reference_candidates.json` and
`data/propermab_published_ranges.json` hold small numbers of values read from
Jain et al. 2017 (*PNAS*) Table S1 and from the mAbs supplementary archive for
doi:10.1080/19420862.2025.2474521.

Individual measured values are facts and are not protected by copyright.
Reproducing a substantial portion of a supplementary table is a different
matter. Both files have been checked and hold only the handful of values the
code actually consumes, each carrying its DOI. The upstream archives are not
mirrored. This is the right side of the line and should stay there: if these
files grow, keep them to values in use.

**One item to verify:** `experimental_reference_candidates.json` records
`https://github.com/Graylab/FLAb` as a sequence-correspondence source. FLAb
aggregates datasets under mixed terms, many of them CC BY-NC, CC BY-NC-ND or
CC BY-NC-SA. The field reads as a cross-check citation rather than a copy, and
no FLAb-specific content was found in the file, but confirm that no sequence or
assay value was taken from FLAb rather than from Thera-SAbDab. If any was, it
must be removed: a non-commercial term in a repository offered to industry is
the same problem as PROPERMAB.

## 6. Copyright and licence

**6.1 Copyright holder.**
The author developed this software independently and at no cost. It was not
produced under employment by, or paid engagement with, any institution, so no
work-for-hire or institutional assignment applies. Copyright vests in the
author, and the licence choice is the author's to make. `LICENSE` and `NOTICE`
name the author as holder.

The antibody molecules this software was used to analyse are a separate matter:
that work is a collaboration with Fred Hutchinson Cancer Center and is covered
by a provisional application that has not published. Nothing about the
collaboration affects ownership of this code, and no collaborative or
molecule-specific material is in this repository. If a written collaboration
agreement exists, it should be read as covering the antibody material rather
than this tooling.

**6.2 Licence.**
Apache-2.0 is drafted and recommended: permissive, so academic and industry use
are both unrestricted; it carries an express patent grant and a contribution
term, which matters for a repository adjacent to a patent application; and it is
compatible with every cleared dependency. MIT is the alternative if brevity is
preferred. Permissive in either case keeps the holder's own commercial options
open, since the holder is not bound by the licence it grants others.

**6.3 The one pre-push check that remains.**
Run `tools/make_oss_export.py` before any public push. It refuses to produce an
export containing the provisional antibody sequences, matched by hash over
amino-acid runs rather than by literal, so the verifier does not itself carry
the sequences. This is the only ownership- or disclosure-related gate on
publication.

## 7. What a recipient must do

Nothing, to read the code or run the permissively-licensed core. To use the
optional surface-patch path they must obtain PROPERMAB themselves and satisfy
Regeneron that they are eligible. To use the reference panel they must build it
from Oxford's CSV. `USE_TIERS.md` states this per module for academic and for
industry users, and the application prints the same distinction at startup.
