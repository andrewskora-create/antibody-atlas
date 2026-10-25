# Dataset provenance and terms

One row per dataset: where it came from, what terms it carries, and whether
this repository redistributes it. "Redistributed: no" means the dataset is
built on your machine by `setup/build_reference_data.py`.

---

## Redistributed

### `data/human_v_ogrdb.json`, `data/human_v_ogrdb_numbered.json`
- **Source:** OGRDB AIRR germline sets, *Homo sapiens*, published releases
  (IGH v10, IGK v5, IGL v4) — https://ogrdb.airr-community.org
- **Terms:** CC0. The site states "Information provided on this site is
  published under a CC0 (Creative Commons 0) license." Copyright in individual
  contributions is retained by contributors; CC0 carries no attribution
  condition, though citing the OGRDB paper is good practice.
- **Redistributed:** yes. CC0 permits it.
- **Used for:** germline residue frequency in substitution ranking; humanness.

### `data/therapeutic_ptm_evidence.json`
- **Source:** per-case literature citations with DOIs and FDA label URLs. Holds
  citations and short factual findings, not reproduced text.
- **Terms:** facts and citations; no licence needed.
- **Redistributed:** yes.

### `data/propermab_published_ranges.json`
- **Source:** values read from supplementary CSVs of
  doi:10.1080/19420862.2025.2474521 (*mAbs*), retrieved 2026-10-04.
- **Terms:** individual measured values are facts. Only the values the code
  consumes are held here; the upstream archive is not mirrored.
- **Redistributed:** yes, as a minimal factual extract with DOI attribution.
- **Note:** these are observed published ranges, not author-defined acceptance
  limits, and were computed under a different APBS version and sampling scheme
  than this project's local protocol.

### `data/experimental_reference_candidates.json`
- **Source:** Jain et al. 2017, *PNAS*, Table S1 subset (PMC5293111);
  sequence correspondence cross-checked against Thera-SAbDab and
  https://github.com/Graylab/FLAb.
- **Terms:** measured values are facts; held to those in use.
- **Redistributed:** excluded from the default permissive release while the source audit below remains open. The local working copy may retain this optional extract.
- **Open item:** FLAb aggregates datasets under mixed terms, many CC BY-NC.
  Confirm no sequence or assay value here originated from FLAb rather than
  Thera-SAbDab. See COMPLIANCE.md §5.

### `static/vendor/ngl/`
- **Source:** NGL viewer, © 2014-2017 Alexander S. Rose.
- **Terms:** MIT. Licence retained alongside the code.
- **Redistributed:** yes.

---

## Not redistributed — built locally

All three derive from a single upstream file:

- **URL:** https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv
- **sha256:** `cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d`
- **Retrieved:** 2026-10-04
- **Clinical-stage column snapshot:** Highest_Clin_Trial, February 2025
- **Terms:** **none stated.** Oxford's Thera-SAbDab pages publish no licence and
  no terms of use. They request citation of Raybould, Marks et al., *Nucleic
  Acids Res.* (2019) gkz827, state the resource was built "under an
  open-innovation agreement" without publishing that agreement, and direct
  licensing enquiries to opig@stats.ox.ac.uk. No licence grant means no default
  right to redistribute a derived corpus.

### `data/clinical_frameworks.json` — 774 clinical-stage chains with IMGT numbering
Rebuild: `python setup/build_reference_data.py --clinical`

### `data/therapeutic_ptm_sequences.json` — reference chains for documented-PTM matching
Rebuild: `python setup/build_reference_data.py --ptm-sequences`

### `data/profiler_comparators.json` — protocol-matched reference panel, 36 entries
Rebuild: `python setup/build_reference_data.py --panel`

The published panel is shipped **with sequences removed and everything else
intact**: computed AntiBERTy scores, structure hashes, descriptor values,
drug names, isotypes, clinical-stage labels, and a per-chain sha256 for each
removed sequence. Rebuilding reinserts the sequences and verifies each against
its hash, so a local rebuild is provably identical to the panel these published
results were computed from. Names and isotypes are facts and are retained.

### Not a dataset, but fetched the same way
The AntiBERTy checkpoint (104 MB, MIT) is fetched on first use rather
than committed, because it exceeds GitHub's recommended blob size. Licence
permits redistribution; this is a size decision only.

---

## Verification

```
python setup/build_reference_data.py --verify
```

reports, for each dataset, whether it is present, whether the upstream CSV hash
matches the recorded one, and whether every rebuilt chain matches its recorded
per-chain hash. A hash mismatch means Oxford has updated the CSV since
2026-10-04; the script says so rather than silently producing a different panel,
because panel membership and descriptor ranges would both shift.

## Human immunoglobulin constant regions — UniProtKB

| | |
|---|---|
| Source | UniProtKB (reviewed / Swiss-Prot), release 2026_03 |
| Retrieved | 2026-10-11 |
| Licence | **CC BY 4.0** — <https://creativecommons.org/licenses/by/4.0/> |
| Redistribution | **Permitted.** UniProt applies CC BY 4.0 to all copyrightable parts of its databases and allows copying provided the copyright notice travels with the data. |
| Shipped in-tree | **Yes** — `data/ig_constant_regions.json` |
| File sha256 | `2c2b02b3ea0e105cabab6be3dc2759226ccc7434d19c5b1269845a5d6ad77dd2` |
| Rebuild / verify | `python setup/fetch_constant_regions.py [--verify]` |

Required copyright notice, reproduced with the data:

> UniProtKB data are distributed under the Creative Commons Attribution 4.0 International (CC BY 4.0) License. Copyright 2002-2026 UniProt Consortium. See https://www.uniprot.org/help/license

Entries (9): IGHG1 (P01857-1, 330 aa), IGHG2 (P01859-1, 326 aa), IGHG3 (P01860-1, 377 aa), IGHG4 (P01861-1, 327 aa), IGKC (P01834, 107 aa), IGLC1 (P0CG04, 106 aa), IGLC2 (P0DOY2, 106 aa), IGLC3 (P0DOY3, 106 aa), IGLC7 (A0M8Q6, 106 aa)

This is the one reference dataset in this project that **is** redistributable.
Thera-SAbDab-derived corpora are not, and are rebuilt locally instead; the
contrast is deliberate and each dataset's row states which applies.

### Two traps this dataset contains

**The displayed UniProt isoform for the gamma chains is membrane-bound, not
secreted.** `P01857` as served is 399 aa with a transmembrane helix at 351-371.
The secreted constant region is 330 aa ending `...SLSLSPGK`. Appending the
membrane form to an Fv would add a hydrophobic helix absent from the drug and
corrupt every hydrophobicity and charge figure downstream. The fetcher
therefore requests the `-1` isoform and refuses any heavy chain that does not
end at the CH3 C-terminal lysine or that contains a hydrophobic 20-mer near its
C terminus.

**`P0CG05` and `P0CG06` are demerged.** IGLC2 is now `P0DOY2` and IGLC3 is
`P0DOY3`. The old accessions return an inactive entry with no sequence.

## Therapeutic reference panel — scope change

The panel was 36 named molecules and is now every approved human or humanised
paired whole monoclonal antibody in the Thera-SAbDab February 2025 snapshot:
**123 molecules** under identical criteria and an identical protocol.

The criteria did not change. `scripts/build_profiler_comparators_v2.py`
previously carried a hardcoded list that was validated against those criteria
at build time; it now applies them to every row in the snapshot instead. The
per-molecule protocol — ABodyBuilder2 model 1, the same NanoShaper build
(sha256 recorded in the panel file), the same PROPERMAB adapter — is unchanged,
so the original 36 were reused from checkpoint rather than recomputed, and the
new members are directly comparable to them rather than approximately so.

Why this was worth hours of compute: the old panel's own first stated
limitation was *"deliberately selected therapeutic references; this panel is
not a representative therapeutic population."* A hand-picked subset has no
defined coverage, so the only defensible summary of it was the observed range,
and a percentile over it would have described the selection rather than the
population. The panel is now complete within its definition, which is what
makes a percentile mean something. See `MODULE_GUIDE.md` for how percentiles
are reported and when they are withheld.

What did **not** change: approved antibodies remain a survivorship-filtered
population, panel members were still not chosen for favourable experimental
behaviour, and neither the range nor a percentile is an acceptance interval.

## DeepSP — model weights and reference distribution

**Not redistributed. Fetched and converted locally.**

| | |
|---|---|
| Source | <https://github.com/Lailabcode/DeepSP> |
| Licence | MIT, Copyright (c) 2024 Lailabcode |
| Build with | `python setup/fetch_deepsp.py` |
| Verify with | `python setup/fetch_deepsp.py --verify` |
| Needs | `h5py`, at build time only |
| Cite | Kalejaye et al., *Comput Struct Biotechnol J* (2024), doi:10.1016/j.csbj.2024.05.029 |

MIT would permit redistribution. The weights are fetched rather than committed
for the same reason AntiBERTy's are: 17 MB of binary in a repository that is
otherwise 5 MB of text.

What the build step produces:

- `models/deepsp/{SAP_pos,SCM_neg,SCM_pos}.npz` — the authors' trained weights,
  converted from Keras HDF5 to numpy arrays plus an ordered execution plan.
  The conversion is mechanical: no weight value changes, and dropout layers are
  dropped because they are the identity at inference. Converting rather than
  calling Keras avoids adding a TensorFlow dependency for three small
  convolutional networks.
- `models/deepsp/provenance.json` — the sha256 of every upstream file read,
  the ref it was read at, and the fetch timestamp.
- `data/deepsp_reference.json` — percentile breakpoints over DeepSP's published
  training labels (molecular-dynamics values for 20,530 natural human Fv pairs
  from the Observed Antibody Space). Summary statistics only; the labels
  themselves are not stored.

**Verification.** The build refuses to finish unless the converted models
reproduce the descriptor values DeepSP publishes for its own sixteen reference
molecules to better than 1e-3 relative. `tests/test_deepsp.py` re-runs that
check offline against a committed fixture, and separately checks that ANARCII
places every residue in the same alignment column legacy ANARCI did — the
substitution that would otherwise change the answer silently.

**What the percentiles mean.** They place a *predicted* value against a
distribution of *simulated* values for *natural repertoire* antibodies. That is
three steps from an experimental measurement on a therapeutic, and none of the
three is modelled as error. A molecule in the upper tail is unusual among
natural antibodies; that is not the same as unusual among approved drugs, and
the clinical panel built by `setup/build_reference_data.py` is what answers the
second question.
