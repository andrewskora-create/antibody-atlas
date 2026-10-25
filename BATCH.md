# Screening a panel

```bash
python run_batch.py panel.fasta --out results/ --fold
open results/index.html
```

Built for 20–100 molecules, where a run takes hours and the failure modes that
matter are the ones that waste a day.

## Input

One multi-FASTA. Chains pair by name, and these header conventions all work:

```
>MAB01|H         >MAB01|L
>MAB01_VH        >MAB01_VL
>MAB01_HC        >MAB01_LC
>MAB01 heavy     >MAB01 light
>MAB01:H         >MAB01:K
```

A `name,heavy,light` CSV works too; the extension decides.

**The whole file is validated before any work starts.** Unidentifiable chain
tokens, non-standard residues, duplicate names, light chains with no heavy
chain and out-of-range lengths are all reported at once with line numbers.
A typo on line 40 fails in seconds, not an hour in.

```bash
python run_batch.py panel.fasta --out results/ --validate-only
```

## Running

```
--workers N        parallel molecules (default: cores / 3)
--threads N        threads per worker (default 2)
--fold             predict structure: exposure, patches, paratope, figures
--no-substitutions skip the PTM fix proposals
--no-pdf           data only, no reports
--no-resume        start over instead of continuing
--scheme kabat     CDR convention (default imgt)
```

Three properties matter at this scale:

**Failures are isolated.** A molecule that cannot be numbered is recorded and
the run continues. It never takes the batch down, and the dashboard lists what
failed and why.

**Progress survives.** State is written to `batch_state.json` after every
completion, so a crash at molecule 80 does not restart from 1. Re-running the
same command skips what finished and retries only what did not. Ctrl-C lets
in-flight molecules finish, then stops cleanly.

**Threads are capped per worker.** torch and OpenMM each try to use every core.
Several workers doing that at once is slower than one and will stall a laptop,
so each worker is limited to `--threads` and the default worker count leaves
headroom.

Expect roughly 1–3 minutes per molecule with `--fold`, less without it. The run
prints a rolling estimate.

## Outputs

```
index.html              sortable dashboard, links to everything
batch_summary.csv       one row per molecule, every descriptor, for ranking
<name>_brief.pdf        one page: lights, key numbers, flagged sites
<name>_report.pdf       four pages: evidence, biophysics, precedents, citations
<name>.json             the full record
<name>.pdb              predicted structure (with --fold)
batch_state.json        progress, timings, failures
```

`batch_summary.csv` is the stack-ranking artefact: one row per molecule with
CDR and framework site counts, approved-drug precedent totals, pI, charge,
GRAVY, hydrophobic windows, humanness, patch priorities and predicted contact
counts. Sort it however the decision needs sorting.

## Reading the dashboard

Green means the screen found no qualifying signal in that category. Amber means
it called review, or a sequence descriptor sits outside the range observed
across the approved and Phase-III reference panel. Red means the screen already
rated it high priority — which only CDR sites can be, since framework sites are
capped below high and reported with the number of approved drugs carrying the
same motif at the same position.

"Approved precedent" sums those counts across a molecule's framework sites. A
high number means its framework liabilities are well precedented, which is
usually the end of the question for them.

None of it is a prediction. Reference panels were not selected for favourable
experimental behaviour, so their range is an observed minimum to maximum and a
within-range result is not a pass.

## Two tiers, and what they actually cost

**Sequence-only (default).** Numbering in both schemes, PTM liabilities,
germline comparison, humanness, the per-residue table and the molecule-level
charge and hydrophobicity flags. Add `--isotype IGHG1` (or 2/3/4) to compute pI
and net charge on the assembled full-length H2L2 rather than the Fv.

**Structure (`--fold`).** Everything above plus ABodyBuilder2, relative
side-chain exposure, surface patches and the predicted paratope. Minutes per
molecule rather than seconds.

### Measured, not estimated

On two workers, twelve molecules, sequence-only:

| Mode | Wall time | Per molecule (steady state) |
|---|---|---|
| Default (precedents + substitutions) | 3 min 15 s | 10–20 s |
| `--no-precedents --no-substitutions` | 56 s | 2.4–6.7 s |

Two things to understand before sizing a run:

**Each worker pays a one-off start-up of roughly a minute** — imports, the
numbering model, and the reference corpora, which include a 14 MB clinical
library. In the run above the first molecule on each worker took 69 s and the
tenth took 10 s. Throughput only settles after several molecules per worker, so
a two-molecule batch tells you nothing about the cost of a two-hundred-molecule
one, and spreading a small panel across many workers makes it slower rather
than faster.

**Precedent lookups dominate the steady-state cost.** Warm, a molecule costs
about 6 s, of which roughly 5 s is the documented-PTM and clinical comparator
searches. `--no-precedents` is the right switch for a large triage pass: you
keep the table, the liabilities and the flags, and lose the "how many approved
drugs carry this motif" column. Run the survivors again without it.

Sizing from the numbers above: on four workers a triage pass settles at roughly
2–3 s of wall time per molecule, so a few hundred molecules in well under an
hour. Structure-tier runs are a different order of magnitude — plan those in
molecules per hour, not per minute.
