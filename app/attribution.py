"""Which published method produced each assessment, and where to read it.

Every conclusion this project prints is produced by a specific method, most of
them other people's. A reader deciding whether to act on a priority call needs
to see what produced it without digging through a methods appendix, so each
assessment carries its own attribution line and link.

Citation is not evidence that a method ran, and no cited paper validates this
project's combined heuristics. Where a call is this project's own screening
policy rather than a published method, the entry says so.
"""
from __future__ import annotations

from .method_sources import METHOD_SOURCES

BY_ID = {source['id']: source for source in METHOD_SOURCES}

# assessment -> (what produced it, source ids, whether the rule is ours)
ASSESSMENTS: dict[str, dict] = {
    'ptm_priority': {
        'label': 'PTM site detection and priority',
        'basis': 'Sequence motifs with IMGT numbering and relative side-chain exposure',
        'ids': ['ptm_heuristic', 'anarcii', 'immunebuilder', 'freesasa', 'additional_liabilities'],
        'local_rule': 'Exposure banding and priority tiers are this project’s screening policy, '
                      'not a published threshold. Framework sites are capped below high.',
    },
    'framework_comparator': {
        'label': 'Approved and Phase-III comparator search',
        'basis': 'Exact motif at the aligned position in clinical-stage therapeutic sequences',
        'ids': ['therapeutic_precedents', 'anarcii'],
        'local_rule': 'Library is restricted to Approved and Phase-III entries. Retention is not '
                      'proof of tolerance.',
    },
    'documented_ptm': {
        'label': 'Documented PTM precedent',
        'basis': 'Literature-reported modification at a homologous site',
        'ids': ['therapeutic_precedents'],
        'local_rule': 'Each case carries its own citation and experimental conditions.',
    },
    'folding_stability': {
        'label': 'Predicted folding stability',
        'basis': 'ThermoMPNN single-substitution estimate from parent-chain coordinates',
        'ids': ['thermompnn'],
        'local_rule': 'Not a measured energy or binding prediction; it does not reorder the initial shortlist.',
    },
    'surface_patches': {
        'label': 'Surface patch findings',
        'basis': 'Atlas contact-graph clustering of exposed side chains; separate PROPERMAB results only when installed and run',
        'ids': ['freesasa', 'propermab', 'immunebuilder', 'charge_clearance'],
        'local_rule': 'Follow-up tiers are this project’s policy. The reference panel range is an '
                      'observed minimum to maximum, not an acceptance interval.',
    },
    'structure_ensemble': {
        'label': 'Ensemble exposure',
        'basis': 'Worst-case relative side-chain exposure across several ABodyBuilder2 models',
        'ids': ['structure_ensemble', 'immunebuilder', 'freesasa'],
        'local_rule': 'Flag-if-any follows TAP2\u2019s convention and ABodyBuilder2\u2019s own '
                      'ensemble recommendation. Which models are usable is checked at runtime; '
                      'the band thresholds are this project\u2019s screening policy.',
    },
    'repack': {
        'label': 'Modelled substitution geometry',
        'basis': 'Mutant side chain built and relaxed against a rigid environment',
        'ids': ['repack', 'pdbfixer', 'openmm', 'immunebuilder', 'freesasa'],
        'local_rule': 'The wild type is relaxed by the identical route, so the difference '
                      'reflects the substitution rather than the protocol; an identity '
                      'control establishes the noise floor. No stability change is '
                      'reported, by design.',
    },
    'hydrogen_bonds': {
        'label': 'Hydrogen-bond context',
        'basis': 'Donor-acceptor geometry with burial and local packing',
        'ids': ['repack', 'freesasa', 'immunebuilder'],
        'local_rule': 'Standard distance and angle cut-offs, not this project\u2019s. Burial '
                      'bands match the exposure bands used elsewhere. Counts and geometry '
                      'only \u2014 no energy, so two substitutions that each break one '
                      'buried bond cannot be ranked here.',
    },
    'deepsp': {
        'label': 'Sequence-only spatial aggregation and charge descriptors',
        'basis': 'Convolutional network over the IMGT-aligned paired Fv sequence',
        'ids': ['deepsp', 'anarcii'],
        'local_rule': 'Percentiles place a predicted value against DeepSP\u2019s own training '
                      'labels, which are simulated values for 20,530 natural human Fv pairs. '
                      'That is a natural-repertoire backdrop, not a clinical-stage benchmark, '
                      'and the 95th-percentile flag is this project\u2019s screening policy.',
    },
    'full_length_assembly': {
        'label': 'Full-length constant-region assembly',
        'basis': 'Secreted human constant regions appended to the variable domains as H2L2',
        'ids': ['ig_constant_regions', 'biopython'],
        'local_rule': 'One allotype per subclass, no glycan, no Fc engineering. The assembly makes '
                      'charge descriptors comparable with published full-length analyses; it does '
                      'not make them predictions.',
    },
    'biophysical': {
        'label': 'Biophysical and formulation descriptors',
        'basis': 'Sequence descriptors plus structure-derived charge, hydrophobicity and viscosity terms',
        'ids': ['ptm_heuristic', 'propermab', 'charge_clearance'],
        'local_rule': 'Descriptors inform clearance, viscosity and self-association risk; none predicts it.',
    },
    'paratope': {
        'label': 'Predicted antigen contact',
        'basis': 'Equivariant graph neural network over the predicted Fv',
        'ids': ['paragraph', 'immunebuilder'],
        'local_rule': 'Cut-off 0.734 is the authors’ recommendation. Prediction is not a measured epitope.',
    },
    'substitutions': {
        'label': 'PTM fix proposals and ranking',
        'basis': 'Exhaustive substitution with rescan, ranked by conservativeness and positional precedent',
        'ids': ['substitution_ranking', 'ogrdb', 'therapeutic_precedents'],
        'local_rule': 'Ordering, the surface tie-break and the conventional-fix recovery rule are this '
                      'project’s policy.',
    },
    'humanness': {
        'label': 'Human framework identity',
        'basis': 'Exact FR1–FR3 identity against human germline V alleles',
        'ids': ['ogrdb', 'anarcii'],
        'local_rule': 'This identity metric is this project’s calculation, not a validated '
                      'immunogenicity prediction.',
    },
    'language_model': {
        'label': 'Antibody language-model scoring',
        'basis': 'Masked-residue pseudo log-likelihood',
        'ids': ['antiberty'],
        'local_rule': 'Likelihood is not affinity, immunogenicity or developability.',
    },
    'structure': {
        'label': 'Predicted structure and exposure',
        'basis': 'Fv structure prediction with solvent-accessible surface area',
        'ids': ['immunebuilder', 'freesasa', 'anarcii', 'anarci', 'hmmer'],
        'local_rule': 'Predicted coordinates are a hypothesis, especially for CDR-H3.',
    },
}


def sources_for(assessment: str) -> list[dict]:
    """Full source records behind one assessment."""
    entry = ASSESSMENTS.get(assessment)
    return [BY_ID[i] for i in entry['ids'] if i in BY_ID] if entry else []


def citation(assessment: str) -> dict | None:
    """Compact attribution for one assessment: names, primary link, caveat."""
    entry = ASSESSMENTS.get(assessment)
    if not entry:
        return None
    sources = sources_for(assessment)
    primary = next((s for s in sources if s.get('publication_url')), None)
    return {
        'assessment': assessment,
        'label': entry['label'],
        'basis': entry['basis'],
        'methods': [s['name'] for s in sources],
        'local_rule': entry.get('local_rule'),
        'primary_url': (primary or {}).get('publication_url'),
        'primary_title': (primary or {}).get('publication_title'),
        'references': [
            {'name': s['name'], 'title': s.get('publication_title'),
             'url': s.get('publication_url'), 'software_url': s.get('software_url')}
            for s in sources
        ],
    }


def line(assessment: str, *, max_methods: int = 4) -> str:
    """Compact attribution for a fixed-layout report block."""
    entry = citation(assessment)
    if not entry:
        return ''
    methods = entry['methods'][:max_methods]
    if len(entry['methods']) > max_methods:
        methods.append(f"and {len(entry['methods']) - max_methods} more")
    text = f"Method: {entry['basis']}. {', '.join(methods)}."
    if entry.get('local_rule'):
        text += f" {entry['local_rule']}"
    return text


def all_citations() -> dict[str, dict]:
    return {name: citation(name) for name in ASSESSMENTS}
