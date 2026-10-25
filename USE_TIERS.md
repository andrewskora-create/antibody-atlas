# Who may use what

This project is offered under a permissive licence, so **the core is open to
everyone, academic and commercial alike, with no extra permission needed.** Two
optional components are not ours to license, and those carry conditions set by
their owners, not by us.

Read the two columns below for your own situation. The application prints the
same summary at startup and marks any result that depended on a restricted
component.

## At a glance

| Capability | Academic / non-profit | Commercial / industry |
|---|---|---|
| CDR and framework PTM liability screening | Yes | Yes |
| IMGT / Kabat numbering (ANARCI, ANARCII) | Yes | Yes |
| Fv structure prediction (ABodyBuilder2) | Yes | Yes |
| Solvent exposure and exposure-banded PTM priority | Yes | Yes |
| Sequence descriptors: pI, net charge, GRAVY, charge asymmetry | Yes | Yes |
| Local contact-graph surface patch screen | Yes | Yes |
| Paratope prediction (Paragraph) | Yes | Yes |
| PTM-fix substitution engineering and ranking | Yes | Yes |
| Humanness / germline analysis (OGRDB, CC0) | Yes | Yes |
| AntiBERTy language-model scoring | Yes | Yes |
| Reports and the web interface | Yes | Yes |
| DeepSP spatial aggregation and charge prediction (MIT) | Yes | Yes |
| **Clinical-precedent library and reference panel** | Build locally | Build locally |
| **PROPERMAB surface-patch descriptors** | Yes, if you qualify | **No — needs a licence from Regeneron** |

Nothing in the first block requires permission from anyone. That block is the
whole scientific workflow: liabilities in, ranked fixes out, with structure,
exposure, charge, aggregation proxies and paratope context.

## The two conditional components

### Clinical-precedent library and reference panel — both tiers build locally

Used for: positional precedent in substitution ranking, documented-PTM site
matching, and the protocol-matched therapeutic comparison panel.

These derive from Oxford Thera-SAbDab, which publishes no licence, so we do not
redistribute them. Everyone builds them on their own machine:

```
python setup/build_reference_data.py --all
```

This downloads Oxford's CSV, verifies it against a recorded hash, and
regenerates the datasets, checking each rebuilt chain against its per-chain
hash so your copy is provably the same input these published results used. It
takes a few minutes and needs no credentials. Oxford asks that you cite
Raybould, Marks et al., *Nucleic Acids Res.* (2019) gkz827.

Without this step the application still runs. Precedent columns read
"not available — reference data not built" rather than failing.

### PROPERMAB descriptors — academic only, and not from us

Used for: the published-protocol charge and hydrophobic surface-patch
descriptors, the spatial charge map, and the Heiden hydrophobic score.

PROPERMAB is Regeneron's, under a licence that limits use to "an academic
institution or academically employed full-time personnel only" and prohibits
commercial use, integration into software used for commercialization, and
distribution in any commercial product or online service.

We ship none of it. If you qualify:

```
python setup/fetch_optional.py --propermab
```

which prints Regeneron's licence, requires you to affirm eligibility, then
fetches it from Regeneron's own repository. It also needs a GPLv3 NanoShaper
build and APBS/PDB2PQR; the script documents both and leaves the build to you,
so any GPL obligation stays on your machine.

**If you are at a company, do not run this.** Contact Regeneron for terms.
Everything else here is yours to use freely, and the alternatives below cover
most of what the descriptors are used for.

## If you cannot use PROPERMAB

Since the DeepSP integration this is a much smaller gap than it was.

**DeepSP (MIT) predicts the same spatial properties from sequence alone.** It
is a convolutional network trained on molecular-dynamics simulations of 20,530
natural human Fv pairs, and it reports spatial aggregation propensity
(SAP_pos) and positive and negative spatial charge maps (SCM_pos, SCM_neg)
over each CDR and over the whole Fv. Install it with:

```
python setup/fetch_deepsp.py
```

There is no eligibility question and no commercial restriction: DeepSP is MIT,
and the weights are the authors' own, converted from Keras to numpy so no
TensorFlow dependency is added. The conversion is verified against the sixteen
reference molecules the authors publish with the model, and `tests/test_deepsp.py`
re-runs that check, including the alignment step, where this project uses
ANARCII and DeepSP used legacy ANARCI.

Two honest caveats. DeepSP predicts what a simulation would have reported, not
what an experiment measured, so it carries the simulation's error as well as
its own. And its percentiles place you against a *natural repertoire*, not
against approved drugs — being unusual among natural antibodies is not the same
as being unusual among therapeutics. Where PROPERMAB gives you a value computed
under a published protocol and compared with a clinical panel, DeepSP gives you
a prediction compared with a repertoire. Those answer different questions.

Beyond DeepSP, you still get, from permissively-licensed code only:

- the local contact-graph patch screen, which finds charged and hydrophobic
  patches from structure and relative side-chain exposure and reports count,
  summed exposed area and CDR adjacency;
- full sequence-level charge analysis, including pI, net charge at a chosen pH,
  per-chain charge, charge asymmetry and the Fv charge-symmetry parameter;
- aggregation-relevant hydrophobicity, including GRAVY and exposed hydrophobic
  area;
- paratope prediction, so charge and PTM findings can be read against predicted
  antigen contacts.

What you lose is the comparison of patch *areas* against a therapeutic panel
computed under the identical published protocol. That is a comparability claim,
not a measurement. The Therapeutic Antibody Profiler at
`opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap` is a free alternative
for charge-patch context, and is what we would suggest reaching for first.

## What none of this establishes

Every priority this software emits is a screening heuristic for ordering
experiments. None is a probability of aggregation, clearance, immunogenicity or
high viscosity, and a within-range result on any panel is not a pass. Reference
panels here were not selected for favourable experimental behaviour, so their
range is an observed minimum to maximum, never an acceptance interval.
