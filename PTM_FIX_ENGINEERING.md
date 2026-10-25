# Engineering PTM fixes

A developability screen that only reports liabilities leaves you with a list of
problems. This part of the project proposes the substitutions that remove them,
and reports what each one costs you elsewhere in the molecule.

The short version: for every PTM liability it finds, it enumerates all 19
alternative residues, keeps only those that actually remove the chemistry,
checks whether each one creates a new liability, and then ranks what survives by
how conservative it is. Alongside each candidate it reports the effect on pI and
charge, on aggregation proxies, and on predicted antigen contact, so a fix that
looks clean on chemistry but moves the surface is visible as such.

---

## 1. What counts as a liability

From sequence, with IMGT numbering:

| Motif | Chemistry |
|---|---|
| `NG`, `NS`, `NT`, `NQ` | asparagine deamidation |
| `DG`, `DS`, `DT` | aspartate isomerization |
| `N-X-S/T`, X not P | N-linked glycosylation sequon |
| `DP` | acid-sensitive Asp-Pro cleavage |
| `M`, `W` | oxidation-susceptible residue |
| Cys outside IMGT 23 / 104 | noncanonical cysteine |
| N-terminal `Q`/`E` | pyroglutamate formation |

Two of these are reported but never auto-fixed. A noncanonical cysteine may be
forming a deliberate extra disulfide, and N-terminal pyroglutamate is usually
charge heterogeneity rather than an engineering target. Removing either on the
strength of a sequence motif is more likely to break something than fix it.

Where a structure is available, each site also carries relative side-chain SASA
and an exposure band. A buried deamidation motif is a different proposition from
an exposed one, and the exposure is reported rather than folded into a single
score.

## 1a. Framework sites are never high priority

A framework PTM has the same chemistry as a CDR one and a different
consequence. It is not at the binding interface, approved products carry these
motifs routinely, and the ones that matter are a CMC and control-strategy
question rather than a potency question. Escalating a framework site to high
on exposure alone buries the CDR sites that actually decide whether the
molecule works.

So **high priority is reserved for CDR sites.** A framework or unassigned-region
liability is capped at review however exposed it is, and the record carries
`framework_priority_capped` with the reason. The site is still reported in
full — it is not hidden or downgraded in substance, only in rank.

**Every framework site is then searched automatically against approved and
Phase-III drugs.** The clinical library is built from Thera-SAbDab entries at
those two stages only (365 drugs, 774 chains), so a match is a product that
reached at least Phase III carrying the same motif at the same aligned
position, in a chain of at least 90% framework identity. Each site reports how
many approved and how many Phase-III drugs carry it, with names.

This is the most directly relevant evidence available for a framework
liability, and it usually resolves the question. In one screen, all three
framework sites came back well precedented:

| Site | Motif | Exposure | Approved | Phase-III |
|---|---|---|---|---|
| H69 | DS | 70% | 37 | 34 |
| H85 | NS | 51% | 22 | 17 |
| H92 | NS | 47% | 46 | 37 |

Retention is not proof of tolerance — those products have process controls and
release specifications that this screen knows nothing about, and none of this
establishes that the site is unmodified. But a motif carried by 46 approved
antibodies at the same position is not where to spend an engineering cycle, and
the report should not imply otherwise.

## 2. Enumeration, then verification

For each fixable site, all 19 substitutions are generated and the mutated
sequence is **re-scanned with the same detector**. This matters more than it
sounds:

- A candidate that leaves the target chemistry intact is discarded. It is not a
  fix.
- A candidate that removes the target but creates a new motif is kept and
  marked, not hidden. `clean` is true only when nothing new appears.

The textbook example is the one that bit this project during development.
At an `N`-`N`-`G`-style arrangement, the conventional `N→Q` and `N→S` fixes can
each create a *new* deamidation motif with the neighbouring residue, while
`N→A` is clean. You only see that by re-scanning, and the ranking shows it.

## 3. Ranking: conservativeness first, deliberately

Candidates are ordered by:

1. **does it edit a position you asked to preserve** (or a predicted paratope
   contact — see §5) — these go last, never hidden
2. **is it clean** — no new liability
3. **conservativeness** — the mean of a normalised BLOSUM62 exchange score and
   the complement of a physicochemical distance over volume, hydrophobicity,
   charge and polarity
4. **positional precedent** — how often this residue occurs at this IMGT
   position in human germline (OGRDB) and in the clinical-stage library
5. **surface consequence** — the charge and aggregation tie-break from §4
6. residue letter, so ordering is deterministic

Conservativeness leads on purpose. An earlier version of this ranking put
surface descriptors higher and buried the conventional choices — `N55Q`,
`M102L` — beneath exotic substitutions that happened to score well on charge.
That is the wrong answer: a well-precedented conservative fix is the one you
should try first, and a surface heuristic should not displace it. The surface
terms are a tie-break among otherwise-equivalent candidates.

A separate safeguard: a dictionary of conventional fixes per chemistry is
checked explicitly, and if a standard choice does not make the top *n* by score
it is still surfaced, annotated `conventional choice recovered below the top
ranks`. You should never have to go looking for `N→Q`.

## 4. What the fix does to charge, pI and aggregation

Chemistry is not the only axis. `N→D` removes a deamidation motif by
introducing the very residue deamidation would have produced, with a full
negative charge. `N→L` removes it by adding exposed hydrophobic surface. Neither
shows up in a BLOSUM score.

Every candidate therefore carries a `developability_delta`:

**Charge and pI** — per chain and for the paired Fv: isoelectric point, net
charge at the chosen pH, formal charge as (K+R)−(D+E), and the Fv
charge-symmetry parameter (the product of the two chains' formal charges).
Each is reported before, after and as a delta.

**Aggregation proxies** — GRAVY, and the count of hydrophobic sequence windows
under the same rule the parent-molecule screen uses: a 7-residue window with
mean Kyte–Doolittle hydropathy ≥ 1.6 and at most one charged residue. A
candidate that creates a new window is flagged.

**Exposure-weighted hydrophobicity** — the hydropathy change multiplied by the
residue's relative side-chain SASA. Adding a leucine at a buried position is not
the same as adding one on the surface, and without a structure this term is
reported unavailable rather than assumed.

These become plain-language flags shown next to the substitution:

```
H34D   dpI -1.22   fv_csp +3   penalty 0.250
   attention | Shifts this chain's pI by -1.22 (6.46 to 5.24).
   attention | Adds a negative charge.
   note      | Fv charge-symmetry parameter moves +3 (0 to 3).
   good      | Reduces exposed hydrophobic surface (weighted hydropathy -2.27).
```

The bounded 0–1 `penalty` combines pI shift, new hydrophobic windows,
exposure-weighted hydrophobicity and structural risk (proline, free cysteine,
buried charge). It is a tie-break, applied after conservativeness and
precedent. It is not a risk probability.

**Scale caveat, stated plainly.** One residue changes a whole IgG's properties
very little. These deltas rank candidates against each other. They do not
estimate what the molecule will do.

## 5. Staying off the antigen contact surface

Most CDR PTM liabilities sit near the paratope, which is exactly where you least
want to edit. Paratope prediction (Paragraph, an equivariant GNN over the Fv
structure) gives each position a probability of contacting antigen, and the
authors' recommended classifier cut-off of 0.734 defines "likely contact".

Two things follow:

- Every candidate reports the predicted contact probability and band for the
  position it edits.
- Where you have not specified preserved positions yourself, likely-contact
  positions are **preserved automatically**, which pushes candidates that edit
  them to the end of the ranking. They remain visible, with the reason given.
  If you do specify preserved positions for a chain, yours win: a known or
  measured contact is better evidence than a prediction.

So the common case resolves the way you want without configuration. Where a
liability sits off-contact, the ordinary conservative fix comes first. Where it
sits on a predicted contact, the fixes that avoid it are ranked ahead of the
ones that do not, and if the only fix is on the contact you are told so rather
than having it quietly promoted.

A prediction is not an epitope. Paragraph does not know this antibody's antigen;
it predicts which residues would contact one. A position predicted not to
contact can still matter to affinity through conformation. Treat a
likely-contact flag as a reason to prefer an alternative and to test binding,
never as a prohibition.

## 5a. How the paratope prediction actually works

Paragraph is an **E(n)-equivariant graph neural network** over the Fv structure.
The pieces:

- **Nodes.** One per residue, placed at its Cα. Each carries a 22-dimensional
  feature vector: a 20-way one-hot of residue type plus a 2-way one-hot of
  heavy or light chain. That is the entire feature set — no evolutionary
  profile, no physicochemical descriptors, no antigen information.
- **Edges.** Built by spatial proximity between those Cα positions, so the
  graph encodes the folded geometry of the loops rather than their sequence
  order.
- **Equivariance.** The network updates node features using interatomic
  distances rather than raw coordinates, so a rotated or translated structure
  gives identical predictions. No alignment step is needed.
- **Output.** A per-residue probability of belonging to the paratope, from six
  graph layers and two linear layers.

Two consequences matter in practice.

**It only looks at the loops.** The scored region is `IMGT_CDRplus2`: IMGT
25–40, 54–67 and 103–119, the three CDRs plus two residues either side. That
is roughly 74–78 positions. Framework residues are never scored, so a framework
contact cannot be predicted at all, only missed.

**It does not know your antigen.** The model is trained on bound
antibody–antigen structures and learns which positions *tend* to form
paratopes, given the shape and composition of these particular loops. It is
answering "if this antibody bound something, where would it touch", not "where
does this antibody touch *your* antigen". Where an antibody's real paratope is
unusual, the prediction will be confidently conventional and wrong.

The threshold of 0.734 is the authors' recommendation, chosen to recover the
observed number of binding residues rather than to maximise accuracy.

## 5b. What it got right and wrong on a real cocrystal

`validation/paratope_cocrystal.py` runs the predictor on **sacituzumab bound to
Trop2 (PDB 9PI9, 1.56 Å X-ray)** and scores it against the observed contacts.
Run it yourself; it downloads everything it needs.

Why this is a fair test:

- **Held out, and the script proves it.** It fetches Paragraph's own six
  train/validation/test files and checks the entry against all 1,184 PDB
  entries in them. 9PI9 is absent. A hold-out claim you cannot check is just an
  assertion, so the check runs every time rather than being quoted here.
- **The coordinates are good.** At 1.56 Å the side-chain positions are
  reliable, so the contact set is a measurement. Contact sets taken from
  low-resolution cryo-EM are estimates, and a validation is only as good as its
  ground truth.
- **The interface is real.** Trop2 is a protein antigen, not a peptide in a
  groove.
- **Two independent copies.** The asymmetric unit holds two Fab–antigen
  complexes, so the same interface is measured twice.
- **Sacituzumab is an approved therapeutic**, which is the kind of molecule
  this project gets pointed at.

Against the 29 antibody residues within 4.5 Å of Trop2:

| Input | Precision | Recall | MCC | AUC |
|---|---|---|---|---|
| Crystal copy 1 (chains A/B) | 0.88 | 0.81 | 0.76 | 0.938 |
| Crystal copy 2 (chains C/D) | 0.88 | 0.81 | 0.76 | 0.936 |
| ABodyBuilder2 model from sequence | 0.77 | 0.77 | 0.64 | 0.898 |

The two crystal copies agree to within 0.002 AUC and give identical TP/FP/FN,
which bounds how much of this is crystallographic noise: almost none.

**Predicting from a model built from sequence alone costs something real but
not disqualifying** — AUC 0.898 against 0.938, MCC 0.64 against 0.76. That is
the number to believe, because building from sequence is the pipeline's normal
path. Note that the pipeline runs a *single* ABodyBuilder2 model rather than
the four-model ensemble, and the validation does the same so it measures what
users actually get.

Stable across contact definitions: at 4.0 Å precision 0.88 / recall 0.91 /
AUC 0.965, at 5.0 Å precision 0.92 / recall 0.79 / AUC 0.938.

Three observed contacts (H69, L1, L83) lie outside the scored window and are
unreachable by construction, capping recall at 26/29 before the model does
anything.

**By loop:**

| Loop | Scored | Observed contacts | TP | FP | Precision |
|---|---|---|---|---|---|
| H CDR1 | 12 | 2 | 2 | 1 | 0.67 |
| H CDR2 | 12 | 3 | 2 | 2 | 0.50 |
| H CDR3 | 18 | 8 | 7 | 0 | 1.00 |
| L CDR1 | 10 | 6 | 4 | 0 | 1.00 |
| L CDR2 | 7 | 3 | 3 | 0 | 1.00 |
| L CDR3 | 13 | 4 | 3 | 0 | 1.00 |

The heavy-chain CDR3 is called perfectly — 7 of 8 contacts, no false positives
— which is the loop that usually dominates a paratope and the one you most want
protected. Both light-chain CDR1 and CDR3 are clean as well.

**The errors are concentrated in the heavy CDR1 and CDR2**, which between them
produce three of the four false positives at 0.67 and 0.50 precision. That is
the conservative-prediction failure: those loops commonly contact antigen, the
model has only shape and composition to reason from, and here they contribute
less than the average.

The practical consequence is narrow. A false positive makes the tool *more*
conservative — it protects a residue that did not need protecting, so a
substitution gets deprioritised rather than a liability being missed. The cost
is a fix you could have taken; the cost is never a liability you failed to see.

One antibody is a sanity check, not a benchmark, and a 0.88-precision result on
a well-ordered 1.56 Å interface is close to a best case. Expect worse on a
flexible or unusual paratope. When a fix is deprioritised only because of a
predicted contact, that is the case to check: supply `preserve` yourself from
real epitope data when you have it, because a measured contact always beats a
prediction and the code keeps yours when you give it.

## 6. Running it

In the interface, fixes appear inline under each liability — no separate step.
Headless:

```
python run_panel.py --panel panel.json --out results/ --fold --substitutions
```

`--fold` predicts structure, which enables the exposure and paratope terms.
Without it the chemistry, charge and sequence-level aggregation terms still
work and the exposure-dependent ones report as unavailable.

To pin positions yourself, add `preserve` to a panel entry with one-based input
positions:

```json
{"H": "EVQ...", "L": "DIQ...", "preserve": {"H": [55, 57], "L": [91]}}
```

Top-ranked proposals can be emitted as variant specifications for the
comparison runner. Each line is a single substitution: combinations are not
auto-generated, because stacking edits before testing each one confounds the
readout.

## 7. What this does not tell you

Ranking is for conservativeness, precedent and surface consequence. None of it
predicts retained antigen binding, potency, expression, thermostability or
aggregation. Germline and clinical frequencies are observational: a common
residue at a position is not established as tolerated *in this antibody*, and an
absent one is not established as harmful.

Every candidate needs binding and activity testing, and the underlying
liability should be confirmed analytically — peptide mapping for deamidation and
isomerization, forced degradation for oxidation — before you spend a cycle
engineering it away. A motif is a sequence pattern, not an observed
modification.
