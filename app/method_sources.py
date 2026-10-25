"""Human-readable method provenance; references verified 2026-10-04.

Runtime versions/checkpoints are recorded separately. Citation is not evidence
that a method ran, and papers do not validate the app's combined heuristics.
"""

METHOD_SOURCES = [
    {
        "id": "antiberty", "name": "AntiBERTy",
        "role": "Local antibody language-model species/chain classification and mean masked-residue pseudo log-likelihood when explicitly run.",
        "publication_title": "Deciphering antibody affinity maturation with language models and weakly supervised learning (Ruffolo, Gray and Sulam, 2021; preprint)",
        "publication_url": "https://arxiv.org/abs/2112.07782",
        "software_url": "https://github.com/jeffreyruffolo/AntiBERTy",
        "limitations": "Sequence-model likelihood and species classification are not percent human framework identity, immunogenicity, affinity, or clinical developability probabilities. AntiBERTy is distinct from AntiBERTa.",
    },
    {
        "id": "immunebuilder", "name": "ImmuneBuilder / ABodyBuilder2",
        "role": "Local paired antibody variable-domain structure prediction; ImmuneBuilder also provides NanoBodyBuilder2 for heavy-only domains.",
        "publication_title": "ImmuneBuilder: Deep-Learning models for predicting the structures of immune proteins (Abanades et al., 2023)",
        "publication_url": "https://www.nature.com/articles/s42003-023-04927-7",
        "software_url": "https://github.com/oxpig/ImmuneBuilder",
        "limitations": "Predicted coordinates and ensemble uncertainty are not experimental structures. Loop/side-chain errors and conformational changes affect surface and PTM screens; Fv predictions do not establish full-IgG behavior.",
    },
    {
        "id": "anarcii", "name": "ANARCII",
        "role": "Local antigen-receptor numbering used to assign the selected IMGT or Kabat positions and regions.",
        "publication_title": "ANARCII: A Generalised Language Model for Antigen Receptor Numbering (Greenshields-Watson et al., 2025; preprint)",
        "publication_url": "https://doi.org/10.1101/2025.04.16.648720",
        "software_url": "https://github.com/oxpig/ANARCII",
        "limitations": "Numbering conventions define different CDR boundaries. Numbering does not identify an experimentally established antigen-contacting paratope; ambiguous/unassigned residues require review.",
    },
    {
        "id": "freesasa", "name": "FreeSASA",
        "role": "Local solvent-accessible surface area, including relative side-chain exposure for PTM screening and local patch descriptors.",
        "publication_title": "FreeSASA: An open source C library for solvent accessible surface area calculation (Mitternacht, 2016)",
        "publication_url": "https://doi.org/10.12688/f1000research.7931.1",
        "software_url": "https://freesasa.github.io/",
        "limitations": "SASA is geometric exposure in the supplied structure. It is not percent modified or a PTM probability; relative side-chain, whole-residue and atom-specific SASA are different quantities.",
    },
    {
        "id": "ogrdb", "name": "OGRDB human germline reference",
        "role": "Reference human V alleles for the app's exact framework identity calculation; the result records the bundled locus release and reference hash.",
        "publication_title": "OGRDB: a reference database of inferred immune receptor genes (Lees et al., 2020)",
        "publication_url": "https://doi.org/10.1093/nar/gkz822",
        "software_url": "https://ogrdb.airr-community.org/germline_sets/Human",
        "limitations": "The app's FR1–FR3 identity metric is its own calculation against this reference, not a validated clinical ADA or immunogenicity prediction. Allele choice and denominator affect the percentage; CDRs and FR4 are excluded.",
    },
    {
        "id": "ngl", "name": "NGL Viewer",
        "role": "Local browser rendering of coordinates, surfaces, side chains, heat maps and PTM markers.",
        "publication_title": "NGL viewer: web-based molecular graphics for large complexes (Rose et al., 2018)",
        "publication_url": "https://doi.org/10.1093/bioinformatics/bty419",
        "software_url": "https://github.com/nglviewer/ngl",
        "limitations": "Rendering does not validate the model. Patch heat colors and PTM marker colors are assigned by this app, not an NGL scientific risk predictor.",
    },
    {
        "id": "tap", "name": "Therapeutic Antibody Profiler (official TAP)",
        "role": "Optional external Oxford service; the local app provides a handoff rather than executing official TAP.",
        "publication_title": "Five computational developability guidelines for therapeutic antibody profiling (Raybould et al., 2019)",
        "publication_url": "https://doi.org/10.1073/pnas.1810576116",
        "software_url": "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap",
        "limitations": "Local exposed-cluster counts are not TAP PSH, PPC, PNC or SFvCSP scores. Official submission sends the pasted sequence to Oxford. Reference-distribution flags do not guarantee drug success.",
    },
    {
        "id": "propermab", "name": "PROPERMAB",
        "role": "Official local sequence and structural descriptor calculations when explicitly run; per-run results identify completed and unavailable features.",
        "publication_title": "PROPERMAB: an integrative framework for in silico prediction of antibody developability using machine learning (Li et al., 2025)",
        "publication_url": "https://doi.org/10.1080/19420862.2025.2474521",
        "software_url": "https://github.com/regeneron-mpds/propermab",
        "limitations": "The separate Atlas contact-graph display is not PROPERMAB. Local descriptor output does not constitute a fitted druggability model; no trained outcome checkpoint is supplied. Academic research license applies. Per-run availability identifies the APBS/NanoShaper scope.",
    },
    {
        "id": "ptm_heuristic", "name": "App PTM screening rules",
        "role": "App-authored sequence motifs plus relative side-chain exposure, with separate potential CDR functional impact and physical-residue deduplication.",
        "publication_title": "Local heuristic; no independent validation publication",
        "publication_url": "",
        "software_url": "",
        "limitations": "The 10%/30% exposure bins and high/review/lower categories are disclosed screening heuristics, not calibrated probabilities, measured PTM levels or literature-validated universal cutoffs. Supporting literature does not validate this combined rule set.",
        "evidence_sources": [
            {"title": "Yang et al. 2017: Met oxidation and side-chain SASA", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5419077/"},
            {"title": "Sharma et al. 2014: MD-based Trp oxidation and Asp isomerization", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4284567/"},
            {"title": "Yan et al. 2018: structure-based Asn deamidation", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6152450/"},
            {"title": "Hoffmann et al. 2024: deamidation/isomerization model comparison", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10984128/"},
            {"title": "N-glycosylation and protein folding (2017)", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5562741/"},
        ],
    },
]

METHOD_SOURCES.append({
    'id':'biopython', 'name':'Biopython ProtParam and Kyte-Doolittle hydropathy',
    'role':'Local sequence molecular mass, isoelectric point, aromatic fraction and GRAVY calculations; the app uses the Kyte-Doolittle scale for its sliding-window display.',
    'publication_title':'A simple method for displaying the hydropathic character of a protein (Kyte and Doolittle, 1982)',
    'publication_url':'https://doi.org/10.1016/0022-2836(82)90515-0',
    'software_url':'https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html',
    'limitations':'These are sequence descriptors. The app-authored 7-residue hydrophobic-window threshold and simple side-chain pKa charge model do not predict measured aggregation, viscosity, solubility, or PK.'
})

METHOD_SOURCES.append({
    'id':'therapeutic_precedents', 'name':'Therapeutic PTM evidence and clinical sequence retention',
    'role':'Local sequence-specific comparison with literature-documented PTMs and a separate Phase III/approved therapeutic sequence snapshot.',
    'publication_title':'Thera-SAbDab: the Therapeutic Structural Antibody Database (Raybould et al., 2020)',
    'publication_url':'https://doi.org/10.1093/nar/gkz827',
    'software_url':'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/',
    'limitations':'Clinical sequence retention does not establish modification, tolerance, current regulatory status, or engineering history. Source stage snapshot is February 2025. Human V-gene CDR identity is a protein-sequence match, not proven germline ancestry. Each experimental observation has its own citation and conditions.'
})


METHOD_SOURCES.append({
    'id': 'additional_liabilities', 'name': 'Additional chemical and cysteine flags',
    'role': 'Asp-Pro cleavage, conditional N-terminal Q/E cyclization, and cysteines outside canonical IMGT Cys23/Cys104.',
    'publication_title': 'Experimental and numbering references for additional liability flags',
    'limitations': 'Motif presence is not a reaction rate. Asp-Pro and noncanonical cysteines receive review priority independent of SASA. Additional cysteines may be stabilizing disulfides, not free thiols. Pyroglutamate requires a free mature terminus and is lower-priority characterization context, not generic adverse function.',
    'evidence_sources': [
        {'title': 'Detecting aspartate isomerization and backbone cleavage by NMR (2021)', 'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7897204/'},
        {'title': 'N-terminal Glu to pyroglutamate conversion in vivo (2011)', 'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3064176/'},
        {'title': 'IMGT conserved numbering: Cys23 and Cys104', 'url': 'https://imgt.org/IMGTindex/numbering.php'},
        {'title': 'Noncanonical intradomain disulfide linkages can stabilize variable domains (2023)', 'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10641266/'}
    ]
})


METHOD_SOURCES.append({
    'id': 'paragraph', 'name': 'Paragraph',
    'role': 'Local paratope prediction over the predicted Fv, used to report predicted antigen contact per position and to prefer PTM fixes away from the binding surface.',
    'publication_title': 'Paragraph: antibody paratope prediction using graph neural networks with minimal feature vectors (Chinery, Wahome, Moal and Deane, 2023)',
    'publication_url': 'https://doi.org/10.1093/bioinformatics/btac732',
    'software_url': 'https://github.com/oxpig/Paragraph',
    'limitations': 'Predicted paratope, not a measured epitope. The model does not know this antibody’s antigen and scores only IMGT CDR+/-2 positions, so framework contacts cannot be predicted. The 0.734 classifier cut-off is the authors’ recommendation for recovering the observed number of binding residues, not an accuracy optimum. A local validation against a held-out therapeutic antibody cocrystal is in PTM_FIX_ENGINEERING.md.',
})

METHOD_SOURCES.append({
    'id': 'substitution_ranking', 'name': 'Conservative substitution ranking',
    'role': 'Ranking of candidate PTM fixes by BLOSUM62 exchange score, a physicochemical distance, human germline residue frequency and clinical-stage residue frequency at the aligned position.',
    'publication_title': 'Amino acid substitution matrices from protein blocks (Henikoff and Henikoff, 1992)',
    'publication_url': 'https://doi.org/10.1073/pnas.89.22.10915',
    'limitations': 'Conservativeness and positional precedent are not predictions of retained binding, potency, expression or stability. Germline and clinical frequencies are observational: a common residue is not established as tolerated in this antibody, and an absent one is not established as harmful.',
    'evidence_sources': [
        {'title': 'OGRDB human germline reference (Lees et al., 2020)', 'url': 'https://doi.org/10.1093/nar/gkz822'},
        {'title': 'Thera-SAbDab clinical-stage sequences (Raybould et al., 2020)', 'url': 'https://doi.org/10.1093/nar/gkz827'},
    ],
})

METHOD_SOURCES.append({
    'id': 'charge_clearance', 'name': 'Charge, surface and pharmacokinetic context',
    'role': 'Interpretation of net charge, pI and charge-patch descriptors in relation to nonspecific binding, clearance and high-concentration viscosity.',
    'publication_title': 'In silico selection of therapeutic antibodies for development: viscosity, clearance, and chemical stability (Sharma et al., 2014)',
    'publication_url': 'https://doi.org/10.1073/pnas.1421779112',
    'limitations': 'The charge–clearance relationship is non-monotonic: both strongly positive and strongly negative variants have shown faster clearance, so a low pI is not automatically favourable. These descriptors inform clearance and viscosity risk; none of them predicts it, and viscosity is concentration-dependent and must be measured.',
    'evidence_sources': [
        {'title': 'Datta-Mannan et al. 2016 — surface charge engineering, self-association and clearance', 'url': 'https://www.nature.com/articles/srep38644'},
        {'title': 'Rai et al. 2023 — electrostatic surface and antibody viscosity', 'url': 'https://www.nature.com/articles/s41598-023-28841-4'},
        {'title': 'Jain et al. 2017 — biophysical properties of clinical-stage antibodies', 'url': 'https://doi.org/10.1073/pnas.1616408114'},
    ],
})

METHOD_SOURCES.append({
    'id': 'clinical_pk_charge',
    'name': 'Human clearance against estimated pI and poly-specificity',
    'role': 'Direction and interpretation of the pharmacokinetic charge flag, using the published '
            'association between estimated isoelectric point and measured linear clearance in '
            'clinical-stage antibodies.',
    'publication_title': 'Identifying biophysical assays and properties that enrich for slow clearance '
                         'in clinical-stage therapeutic antibodies (Grinshpun et al., 2021)',
    'publication_url': 'https://doi.org/10.1080/19420862.2021.1932230',
    'limitations': 'The published analysis used 64 full-length IgG molecules and measured human linear '
                   'clearance. This project computes pI from variable domains only, so the value is a '
                   'rank-order proxy and not the molecule pI the authors thresholded on. The authors '
                   'also identify a poly-specificity assay as the stronger discriminator, and that is an '
                   'experiment this project cannot substitute for. Clearance is additionally governed by '
                   'FcRn interaction, target-mediated disposition and anti-drug antibodies, none of which '
                   'is visible from a variable-domain sequence. No clearance value is predicted.',
})

METHOD_SOURCES.append({
    'id': 'clinical_biophysical_landscape',
    'name': 'Clinical-stage biophysical property boundaries',
    'role': 'Context for whether a descriptor sits inside the range spanned by antibodies that reached '
            'advanced clinical stages.',
    'publication_title': 'Biophysical properties of the clinical-stage antibody landscape (Jain et al., 2017)',
    'publication_url': 'https://doi.org/10.1073/pnas.1616408114',
    'limitations': 'The published boundaries come from twelve measured assays on 137 isotype-matched '
                   'full-length IgG1 molecules. This project has not measured any of those assays and does '
                   'not reproduce their thresholds; the reference is cited for the concept of empirical '
                   'clinical-stage limits, not as a validation of any value computed here.',
})

METHOD_SOURCES.append({
    'id': 'colloidal_ph_margin',
    'name': 'Formulation pH margin from the isoelectric point',
    'role': 'Basis for the colloidal-stability flag: net charge approaches zero near the isoelectric '
            'point, reducing electrostatic repulsion between molecules.',
    'publication_title': 'Physical stability of proteins in aqueous solution: mechanism and driving forces '
                         'in nonnative protein aggregation (Chi et al., 2003)',
    'publication_url': 'https://doi.org/10.1023/a:1025771421906',
    'limitations': 'The pH-versus-pI margin is a colloidal-stability heuristic, not an aggregation '
                   'prediction. Real behaviour depends on ionic strength, buffer species, excipients, '
                   'concentration and conformational stability, none of which is modelled. The margin is '
                   'computed from a variable-domain pI, which is not the full molecule pI. This project '
                   'sets the margin bands as its own screening policy.',
})

METHOD_SOURCES.append({
    'id': 'ig_constant_regions',
    'name': 'Human immunoglobulin constant regions (UniProtKB)',
    'role': 'Secreted heavy (IGHG1-4) and light (IGKC, IGLC1/2/3/7) constant-region sequences, '
            'appended to the variable domains so molecule-level charge descriptors are computed '
            'on a full-length H2L2 rather than on the Fv alone.',
    'publication_title': 'UniProt: the Universal Protein Knowledgebase',
    'publication_url': 'https://doi.org/10.1093/nar/gkae1010',
    'software_url': 'https://www.uniprot.org/',
    'limitations': 'Sequences are the UniProt canonical secreted isoform for one allotype per '
                   'subclass. Allotypes differ at charged positions and shift the computed pI. No '
                   'glycan is modelled, so sialylated glycoforms that lower the real pI are absent, '
                   'and no therapeutic Fc engineering (IgG4 S228P, silencing, half-life or charge '
                   'variants) is applied. Assembly is sequence concatenation, not a structural '
                   'model. UniProtKB data are CC BY 4.0.',
})

METHOD_SOURCES.append({
    'id': 'deepsp',
    'name': 'DeepSP',
    'role': 'Sequence-only prediction of spatial aggregation propensity (SAP_pos) and spatial '
            'charge maps (SCM_neg, SCM_pos) over the CDRs and the whole Fv. Used as a fast, '
            'MIT-licensed alternative to the structure-based descriptor route, which needs '
            'folding and, for the PROPERMAB descriptors, a Regeneron academic licence.',
    'publication_title': 'DeepSP: Deep learning-based spatial properties to predict monoclonal '
                         'antibody stability (Kalejaye, Wu, Terry and Lai, 2024)',
    'publication_url': 'https://doi.org/10.1016/j.csbj.2024.05.029',
    'software_url': 'https://github.com/Lailabcode/DeepSP',
    'limitations': 'DeepSP predicts what a molecular-dynamics simulation would have reported, not '
                   'what an experiment measured, so it carries the simulation protocol\'s error '
                   'as well as its own regression error. It was trained on natural human Fv pairs '
                   'from the Observed Antibody Space, so the training distribution is a natural '
                   'repertoire rather than a clinical-stage panel. Its alignment grid holds at '
                   'most seventeen CDR-H3 insertions and no light-chain insertions; a molecule '
                   'that overflows it is reported as unreliable rather than scored. Weights are '
                   'the authors\' own, converted from Keras to numpy and verified against the '
                   'sixteen reference molecules they publish with the model.',
})

METHOD_SOURCES.append({
    'id': 'structure_ensemble',
    'name': 'ABodyBuilder2 model ensemble',
    'role': 'Folding the same Fv with several of ABodyBuilder2\'s four trained models and '
            'taking the worst-case relative side-chain exposure across them, so a PTM site '
            'that is exposed in any model is treated as exposed.',
    'publication_title': 'ImmuneBuilder: Deep-Learning models for predicting the structures of '
                         'immune proteins (Abanades et al., 2023)',
    'publication_url': 'https://www.nature.com/articles/s42003-023-04927-7',
    'software_url': 'https://github.com/oxpig/ImmuneBuilder',
    'limitations': 'Agreement between models is not accuracy. The four models share training '
                   'data and can agree while all being wrong, so a stable exposure band is not '
                   'experimental evidence. Coordinates are never averaged: each model is '
                   'analysed separately and only the per-residue exposure is reconciled. The '
                   'flag-if-any rule raises sensitivity at the cost of specificity by design, '
                   'which means more sites to check, not more real liabilities. Checkpoints are '
                   'verified before use and a corrupt one is reported rather than skipped.',
})

METHOD_SOURCES.append({
    'id': 'repack',
    'name': 'Side-chain repacking (PDBFixer + OpenMM)',
    'role': 'Builds the substituted side chain, relaxes it and its immediate neighbours '
            'against a rigid environment, and reports the hydrogen bonds lost and gained '
            'and any steric clash introduced. Answers what the molecule looks like after '
            'a substitution, rather than inferring the consequence by removing the '
            'wild-type side chain.',
    'publication_title': 'OpenMM 7: Rapid development of high performance algorithms for '
                         'molecular dynamics (Eastman et al., 2017)',
    'publication_url': 'https://doi.org/10.1371/journal.pcbi.1005659',
    'software_url': 'https://github.com/openmm/pdbfixer',
    'limitations': 'No free energy is reported and none should be inferred: force-field '
                   'potential energies of two different sequences are not comparable, so a '
                   'kcal/mol here would be dominated by chemistry unrelated to stability. '
                   'Hydrogen-bond counts and clashes are what this measures. No solvent '
                   'model is used, which is affordable only because the environment is held '
                   'rigid; the mutated side chain relaxes into a fixed pocket rather than '
                   'the structure re-equilibrating. A single relaxed conformation is not a '
                   'rotamer search, so a substitution that would be tolerated in a different '
                   'rotamer can read as a loss. On a predicted backbone this models one '
                   'thing on top of another, and the result carries its structure source for '
                   'that reason. Off-site bond changes are reported as protocol noise so the '
                   'signal-to-noise of each answer is visible.',
})
