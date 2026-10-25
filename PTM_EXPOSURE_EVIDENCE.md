# PTM exposure evidence and screening rules

Evidence review: 4 October 2026. Scope: primary antibody studies of methionine/tryptophan oxidation, asparagine deamidation, aspartate isomerization, and N-linked glycosylation, plus solvent-accessibility method documentation. This is a targeted literature review, not a preregistered systematic review or meta-analysis. Different endpoints, stressors, structures, and SASA definitions prevent pooling the reported cutoffs.

## What “uncalibrated” means

The app's screening categories have not been fitted and independently validated against measured site-specific modification rates or frequencies under specified conditions. High priority means prioritize experimental follow-up. It does not mean a known probability of modification, a percentage modified, an expected shelf life, or clinical risk. CDR location describes possible functional impact separately from chemical susceptibility.

## What exposure measures

FreeSASA calculates solvent-accessible surface area (SASA) on the supplied structure. Relative side-chain SASA is side-chain accessible area divided by the corresponding reference area. Multiplication by 100 produces the displayed exposure percentage. It is not the percentage of the antibody's total surface. FreeSASA relative values may exceed 100% in unusual conformations or at chain ends. Absolute Å², whole-residue relative SASA, side-chain relative SASA, atom-specific accessibility, and molecular-dynamics averages are different descriptors; their cutoffs cannot be exchanged. [FreeSASA documentation](https://freesasa.github.io/python/classes.html)

## Primary evidence

| Chemistry | Study and result | Consequence for this app |
|---|---|---|
| Met oxidation | Yang et al. (2017), 121 clinical-stage antibody variable regions expressed in IgG1/kappa constructs, 0.1% H₂O₂ for 24 h: fractional side-chain SASA around 9% (8–11% range) separated segment-level oxidation counts. Buried H3 Met included false negatives. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5419077/) | Exposure supports follow-up prioritization. Neither burial nor a threshold transferred between normalization methods establishes stability. |
| Met oxidation | Sankar et al.: incorporation of local contacts and fluctuations improved prediction over SASA alone. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6284603/) | Static accessibility omits local interactions and dynamics. |
| Met oxidation | Tavella et al. (2022) compared side-chain SASA, water coordination, and nearby hydroxyl groups. Their >15% relative-SASA boundary used a different normalization and was selected for agreement in the studied datasets; partially buried reactive Met remained challenging. [Study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0279689) | The published boundary is not a universal product-quality limit. This app does not implement the WCN-OH dynamics method. |
| Trp oxidation | Sharma et al.: 38 Trp across 17 mAbs; >30% **MD-averaged side-chain SASA**, approximately 80 Å², identified 13/14 labile and 20/24 nonlabile sites under AAPH stress. Liability was >35% oxidation under that experiment. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4284567/) | A promising context-specific classifier; a static model's 30% is not equivalent to this published MD descriptor. |
| Trp oxidation | Lenka et al. (2026), 187 mAbs/264 CDR Trp: SASA >50 Å² alone had sensitivity 91% but specificity 46%. Only 50/163 exposed sites were reactive. Electrostatic context improved specificity. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13475319/) | Exposure alone cannot support a high-probability call. A local charge patch is not the study's calculated electrostatic potential. |
| Photooxidation | Delmar et al. (2021) combined solvent exposure with environmental features. Cross-validated AUCs were 0.926 for Met and 0.860 for Trp, with separate holdout evaluation. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8060516/) | Peroxide, radical, and light-induced oxidation must be distinguished. |
| Asn deamidation | Yan et al. (2018), five IgG1/five IgG4: a 10 Å² SASA boundary appeared within a decision tree also using structure, flexibility, and reaction geometry. Exposed stable and buried reactive sites occurred. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6152450/) | 10 Å² is not 10% exposure. Preserve motif concern even for buried sites. |
| Asn/Asp | Hoffmann et al. (2024), 57 mAbs and 452 motifs: static models using side-chain ASA and the next residue's ASA achieved balanced accuracies of 62% for Asn and 55% for Asp; the broader method incorporates dynamic and chemical descriptors. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10984128/) | Exposure alone is weak evidence, especially for Asp. Adjacent backbone accessibility and conformational flexibility matter. |
| N-linked glycosylation | Experimental sequon insertion showed that sites buried in the final fold can be glycosylated during biosynthesis. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5562741/) Native accessibility can affect subsequent glycan processing, a different endpoint. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4142333/) | Never downrank a sequon from final folded-state burial. Occupancy remains unresolved. |

## App rules, explicitly heuristic

Exposure bins are **low <10%, intermediate 10–<30%, high ≥30%**. They organize screening and display; they are not universal literature-validated risk cutoffs.

- NG/NS deamidation and DG/DS isomerization flags remain **review** when exposure is low or intermediate, and become **high screening priority** at ≥30% side-chain exposure.
- Other selected Asn/Asp motifs and Met/Trp oxidation sites are **review** at ≥10%, otherwise **lower structural concern**. Exposure alone never gives oxidation a high call.
- N-linked glycosylation sequons are **review**, independent of final-fold exposure.
- Missing exposure or unreliable mapping is unknown/unassessed, never evidence of low concern.
- Overlapping chemistry flags at one residue count as one physical site; all mechanisms remain visible.
- Candidate-level follow-up is recommended for any high site or at least three review sites. Review is recommended for any review site or at least five lower-concern sites. Otherwise no priority site is identified. These count thresholds are also app heuristics, not a combined probability.

Colors: **red = high screening priority; amber = review; green = lower structural concern; gray = unassessed**. Green does not establish stability. A site without a structural measurement can retain a sequence-driven review call while its exposure remains explicitly unknown.

## Experimental interpretation

Use site-resolved peptide mapping under intended storage and relevant forced-degradation conditions to establish modification rate and occupancy. Link modified species to binding/potency and stability where appropriate. Formulation, pH, temperature, duration, oxidant/light, and expression system must accompany any measured result. Calibration would require matching descriptors and endpoints, independent held-out validation, and uncertainty estimates.
