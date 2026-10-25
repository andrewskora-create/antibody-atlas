# Antibody Atlas — module guide

This guide describes what each module calculates and what its output can support. The application is local and intended for research triage, not clinical decisions or assay replacement. After analysis starts, enabled structure, AntiBERTy, and PROPERMAB steps run through the workflow; the overall assessment waits until selected steps reach terminal status. Failed or interrupted steps can be retried manually.

## Analysis modules

| Module | What it does | How to read it |
|---|---|---|
| Sequence validation and numbering | Accepts one VH and/or VL domain; ANARCII numbers residues using IMGT or Kabat and assigns CDRs. ImmuneBuilder also calls legacy ANARCI for structure preparation; HMMER (`hmmscan`) is a system prerequisite for that path, not for ANARCII itself. | Residue positions and CDR boundaries depend on the selected scheme. Insertions retain lettered positions. |
| Framework identity (“humanness”) | Compares numbered IMGT FR1–FR3 positions with the closest functional human V allele in the bundled OGRDB set. Shows coverage, numerator/denominator, mismatches, alternatives, and sequence fingerprint. | A transparent framework identity percentage. It excludes CDRs and FR4; it is not OASis, a calibrated humanness measure, or clinical immunogenicity prediction. |
| Sequence descriptors | Reports charge, pI, molecular mass, mean hydropathy, hydrophobic windows, CDR lengths, and selected PTM-susceptible motifs/residues. | Descriptors and motif flags are screening signals. They do not demonstrate PTM occupancy/rate, aggregation, expression, stability, or viscosity. Hydrophobic windows use a disclosed 7-residue Kyte–Doolittle rule, not a trained assay predictor. |
| Structure prediction and viewer | ImmuneBuilder predicts a local Fv model; paired sequences use ABodyBuilder2, heavy-only input uses NanoBodyBuilder2. One model by default; three is selectable. NGL displays cartoon, atoms, and surfaces. | Predicted coordinates are hypotheses; flexible loops, especially CDR-H3, can be uncertain. This is not a full-length antibody or Fc model and does not predict antigen binding. |
| Structure ensemble (optional) | Folds the Fv with several of ABodyBuilder2's four models and reconciles relative side-chain exposure across them by **maximum**, not mean, so a site exposed in any model is treated as exposed. Reports every residue where the models disagree on the exposure band. Checkpoints are verified before use; a corrupt one is named, not silently skipped. | Costs roughly one fold per model. Agreement is not accuracy — the models share training data and can agree while all being wrong, and the disagreement count is not reproducible run to run because ABodyBuilder2's refinement is stochastic. Flag-if-any deliberately trades specificity for sensitivity, so expect more sites to check, not more real liabilities. Coordinates are never averaged. |
| Hydrogen-bond context | Donor/acceptor geometry from the structure by standard distance and angle cut-offs, annotated with burial and local packing. Reports how many side-chain hydrogen bonds a position holds, which is what a substitution there would have to replace. Runs on every fold (sub-second). | Counts and geometry, never an energy: two substitutions that each break one buried bond cannot be ranked here. Hydrogen bonds live on side chains, and side chains are the least reliable part of a predicted backbone, so every finding carries the structure it came from. Losing an exposed bond is usually tolerable; losing a buried one is the case worth stopping for. |
| Modelled substitution (optional, per substitution) | Builds the replacement side chain with PDBFixer, relaxes it and its neighbours against a rigid environment with OpenMM, and reports hydrogen bonds lost and gained plus any clash introduced. Answers what the molecule looks like *after* the change, rather than inferring it by deleting the wild-type side chain — which systematically over-penalises conservative swaps. About three seconds, so it is asked for per substitution rather than run over every candidate. | **No free energy is reported and none should be inferred.** Force-field potential energies of two different sequences are not comparable, so a kcal/mol here would be dominated by chemistry unrelated to stability. No solvent model is used; that is affordable only because the environment is frozen, so the side chain relaxes into a fixed pocket rather than the structure re-equilibrating. One relaxed conformation is not a rotamer search, so a substitution tolerated in a different rotamer can read as a loss. The wild type is relaxed by the identical route and an identity control establishes the noise floor; off-site bond changes are reported so the signal-to-noise of each answer is visible. A position carrying an IMGT insertion code cannot be addressed unambiguously and is refused rather than guessed. |
| DeepSP spatial descriptors (optional) | Predicts spatial aggregation propensity and positive/negative spatial charge maps over each CDR and the whole Fv from the paired sequence alone, in milliseconds and with no structure step. MIT, so unlike PROPERMAB it carries no commercial restriction. | Predicts what a molecular-dynamics simulation would have reported, not what an experiment measured. Percentiles place it against 20,530 natural repertoire Fvs, which is not a clinical-stage benchmark. A CDR-H3 that overflows DeepSP's alignment grid is reported unreliable rather than scored. |
| Surface screen | FreeSASA and a local side-chain contact graph display relative hydrophobic/positive/negative exposure and clusters. PTM sites can be mapped onto the model and assessed for relative exposure. | Colors and areas are method-specific screening values, not calibrated aggregation or electrostatic-potential scores. PTM motifs are not proof of modification or clinical risk. Unmapped/unknown evidence must remain unknown. |
| AntiBERTy | Runs the bundled AntiBERTy model locally as a workflow step when enabled for chain/species/graft classification and pseudo-log-likelihood. | A model-compatibility score, not an assay, calibrated developability score, humanness percentage, or VH/VL pairing score. |
| TAP handoff | Offers a manual route to Oxford's Therapeutic Antibody Profiler. The app can copy/open the handoff but does not submit sequences. | Local cluster counts are not official TAP PSH/PPC/PNC/SFvCSP scores. Pasting/submitting on the external site is a separate external action. |

## Variant comparison

After an original analysis, compare up to eight named candidates using substitutions at 1-based input positions or full VH/VL sequences. Source residues and chain/length compatibility are checked. The app reports gained, lost, and retained PTM flags, sequence descriptor changes, framework identity, and fresh structure/AntiBERTy/PROPERMAB workflow results when those steps are enabled and available. Structural measures require candidate-specific coordinates and compatible mapping; they are not copied from the parent sequence.

These are calculated deltas for prioritizing experiments. They do not predict a candidate's measured PK, concentration behavior, stability, aggregation, or CMC readiness, and they do not automatically propose mutations.

## PROPERMAB (optional academic research module)

The vendored source is the Regeneron PROPERMAB implementation, with attribution and its upstream academic-use license retained in [vendor/propermab/LICENSE.md](vendor/propermab/LICENSE.md) and source notice in [vendor/propermab/UPSTREAM.txt](vendor/propermab/UPSTREAM.txt). Before using or redistributing it, the receiving group must read and satisfy that included license. This documentation does not grant rights or interpret eligibility for a particular user. The optional Python requirements are in `requirements-propermab.txt`. Install them through the launcher with `python3.11 run.py --setup-only --with-propermab` (substitute `python3.12` as needed). Run `python3.12 run.py --doctor` after installing system prerequisites; doctor only diagnoses and does not install requirements. `--warmup` fetches the ABodyBuilder2/NanoBodyBuilder2 model weights without an antibody sequence.

The base app works without this optional environment. PROPERMAB calculates upstream Fv descriptors and, when native dependencies are available, surface-patch features. The full patch path needs compatible APBS, PDB2PQR, multivalue, and NanoShaper. The Apple Silicon NanoShaper executable from the development machine is excluded from the portable ZIP. The ZIP includes optional `native-source/` (NanoShaper 0.7.8 and patched CGAL 4.2 source with licenses) plus `scripts/build_native.sh`. Before building on macOS, install Xcode Command Line Tools and run `brew install cmake boost gmp mpfr`. On Ubuntu 24.04, install `cmake libboost-all-dev libgmp-dev libmpfr-dev` in addition to the Python/build requirements above. On Unix-like systems, build the binary with `scripts/build_native.sh` and point `PROPERMAB_NANOSHAPER` to the resulting absolute path (for example, `native-build/nanoshaper/NanoShaper`). The build helper was tested from an extracted release on the development Apple Silicon Mac; other recipient platforms are unvalidated; install the required native dependencies first. Patch features may be unavailable while other modules continue to work. Windows users should build/run under WSL2.

The shipped install has been verified on an Apple Silicon Mac only. macOS Intel, Linux, and Windows/WSL2 have not been verified. Use WSL2 for Windows; native Windows is untested. APBS version differences may change numeric results; cross-version equivalence is not claimed. PROPERMAB outputs here are descriptors, not validated PK, viscosity, or solubility predictions. No fitted endpoint checkpoints are bundled, and unavailable calculations are not treated as zero.

## Using an experimental structure

Everything downstream of folding takes a structure path, but until recently
only the folding step wrote one, so a crystal or cryo-EM model could not be
used. `POST /api/structure` accepts a PDB or mmCIF file and produces a job in
the same shape a fold does, so substitutions, PROPERMAB, the paratope
predictor and the viewer all work against it unchanged. The difference is
`structure_source`, which is `experimental` rather than `predicted` and
travels with every finding.

Chains are classified by numbering them: ANARCII recognises heavy and light
variable domains and does not recognise anything else, so the Fv is separated
from the antigen without being told which is which. Heavy and light are paired
by proximity, which also separates copies in an asymmetric unit. Constant
domains are dropped from the Fv by keeping only residues the variable-domain
numbering covers.

Three things are reported rather than smoothed over. Unresolved density
becomes a gap, never a modelled-in loop and never a silently shortened chain.
Multiple copies are listed and one is chosen, never averaged — an averaged
structure is not a structure. And `POST /api/structure/inspect` describes a
file without committing to it, so a copy can be chosen before anything runs.

**With an antigen present, epitope distance is measured.** Paragraph predicts
the *paratope* — positions on the antibody likely to contact something — which
is an inference. A complex gives the minimum heavy-atom distance from each Fv
residue to the antigen, which is a measurement of that structure. The two are
reported separately because they are different claims. What a measured
epitope still is not: necessarily the epitope in solution, since crystal
contacts and the construct used can both differ.

## Hydrogen bonds and side-chain repacking

The substitution ranking scores conservativeness, precedent and predicted
contact — none of which is geometry, so a buried asparagine holding a
framework hydrogen bond scored the same as a solvent-facing one.

`app/hbonds.py` finds donors and acceptors by the conventional distance and
angle criteria and annotates each bond with burial and local packing. Burial
is what makes it actionable: losing a solvent-exposed bond is usually
tolerable, losing a buried one often is not. It reports counts, not an energy.

`app/repack.py` answers the question the inventory cannot. Inferring
consequence by subtraction assumes any replacement lacking the donor or
acceptor loses the bond, which over-penalises conservative swaps. Repacking
builds the mutant with PDBFixer, relaxes the neighbourhood with OpenMM and
re-measures. On the test Fv at an aspartate with six side-chain bonds, Glu,
Asn and Gln form replacement bonds while Ser and Thr do not — the latter
carry a hydroxyl that could donate but are too short to reach. Subtraction
rates all five identically.

Three things make the number trustworthy, and they are worth checking if you
modify this module. The wild type goes through the identical pipeline, so the
difference is the mutation rather than the protocol. Only a shell around the
site moves; the rest is frozen by zeroing masses, which is both faster and a
cleaner comparison than restraints. And `identity_control` — mutate a residue
to itself — reports the noise floor, which is currently zero bonds changed.

No free energy is reported, deliberately. Force-field potential energies of
two different sequences are not comparable, so a plausible-looking kcal/mol
would be dominated by chemistry unrelated to stability.

There is no implicit solvent, which is a real compromise and the reason is
measured: GBn2 took 421 s for twenty minimisation steps on a 3,483-atom Fv,
the cheaper OBC and HCT variants 41–50 s, vacuum with a cutoff 1.4 s. In
vacuum charged side chains normally collapse into spurious salt bridges; what
makes it defensible here is that the environment is rigid, so the only motion
available is a side chain relaxing a few ångström into a fixed pocket.

On a predicted backbone this models one thing on top of another. The `source`
label exists so that stays visible, and an experimental structure is where
this stops compounding two errors.

## The therapeutic reference panel

The panel is every approved human or humanised paired whole monoclonal
antibody in the Thera-SAbDab February 2025 snapshot, scored locally under one
protocol. Selection is programmatic, not a curated list: the panel is the
complete set under its own definition, which is what makes a position within
it interpretable.

That matters for how a comparison is reported. Previously the only honest
summary was the observed minimum to maximum, because a percentile over a
hand-picked subset would have described the curation. Now each descriptor also
carries a percentile, withheld below 60 references on the reasoning that with
n references the finest distinction expressible is 100/n points and the tails
turn on one or two molecules.

A percentile is strictly more informative than within/above/below range, which
cannot distinguish a molecule sitting with the bulk from one pressed against
the wall. It remains a **position**, not a verdict: approved antibodies are a
survivorship-filtered population, the panel was not selected for favourable
experimental behaviour, and nothing here establishes that any position is
acceptable.

## Herceptin ±10% comparison: what it means

The patch assessment can be set to Off, 5%, 10%, or 20%; ±10% is the default. The app compares three categories separately (hydrophobic, positive-charge, negative-charge), using both total and CDR-adjacent PROPERMAB patch areas against compatible trastuzumab (Herceptin) reference values. Complete residue mapping/exposure and matching structure provenance are required; a signal above the reference maximum prevents the favorable comparative call. At ±10%, passing this computational rule can produce **“Low comparative risk — Herceptin-like.”**

This is a chosen computational screening policy. It is not a laboratory assay comparison, an assay-derived threshold, or evidence of equivalent HIC, aggregation, viscosity, PK, safety, or clinical behavior. The tolerance is not validated as predictive. Published reference ranges are observed data ranges, not acceptance limits. The rule does not change PTM priorities. JSON and the concise PDF record the selected tolerance and comparison outcome.

## Report

The PDF is generated locally when the workflow reaches terminal status. It is limited to three pages and includes an abbreviated input identifier, a hydrophobic-surface image and a short two-paragraph summary on page one, followed by prioritized PTM findings and selected precedent/next-step content. It makes missing assay or lineage evidence explicit. Full input/reference/PDB hashes and detailed provenance are available in the app and JSON; the PDF does not contain the full hashes. In the normal workflow the report becomes available after folding and the selected profiler steps finish.

Reports and downloaded PDB/JSON files are saved by the browser to its configured download location, outside automatic app cleanup. Review the report's methods and missing-data statements before sharing or using its conclusions.

## Local privacy and limits

The app binds to localhost (`127.0.0.1`), has no login, and does not send input sequences to remote analysis APIs. Data are kept in server memory for the session; generated PDB intermediates are removed on normal shutdown. Model weights for structure prediction may download on first use. The bundled AntiBERTy weights are read from disk. The app can be reached by local processes on the computer, and downloaded artifacts may enter browser backups or cloud sync. Hosted assistant chats and external tools have separate privacy behavior.

This version does not predict antigen/epitope binding, expression yield, measured aggregation, viscosity, stability, or clinical ADA risk. Single-model structural output has no ensemble confidence estimate. Sequence-window and surface/PTM thresholds are not calibrated to assays. Extra cysteine counts do not establish disulfide pairing or free thiols. Experimental measurement and expert review remain necessary for decisions that depend on these properties.
