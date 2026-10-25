"""Local deterministic synthesis; clinical precedent does not recalibrate PTM chemistry."""
import math
from .patch_assessment import assess_patches


def assessed(site):
    value = site.get('relative_sidechain_sasa')
    return isinstance(value, (int, float)) and math.isfinite(value) and value >= 0


def focused(site):
    oxidation_only = all('oxidation' in str(f.get('kind') or f.get('chemistry', '')).lower()
                         for f in (site.get('flags') or [site]))
    return assessed(site) and (site.get('priority') == 'high' or
        (site.get('priority') == 'review' and not (site.get('region') == 'Framework' and oxidation_only)))


def imgt_position(number):
    """Strip the IMGT chain-type prefix, leaving the canonical position."""
    value = str(number or '')
    return value[1:] if value[:1] in {'H', 'K', 'L'} and value[1:2].isdigit() else value


def site_label(chain, number, residue=None, position=None):
    """An unambiguous site label.

    IMGT prefixes a position with the chain TYPE -- H heavy, K kappa, L lambda
    -- and those letters collide with the one-letter amino-acid codes. 'K56'
    reads as lysine-56 when it means kappa chain, IMGT position 56, and the
    residue there may be anything. Worse, prepending the user's own chain label
    to an already-prefixed number produced strings like 'LK56'.

    So the label leads with the chain as the user supplied it and the residue
    actually present, and reports the IMGT position separately and explicitly.
    """
    imgt = imgt_position(number)
    if residue and position is not None:
        return f"{chain} {residue}{position} (IMGT {imgt})"
    return f"{chain} IMGT {imgt}"


def comparator_summary(site, data):
    """Approved and Phase-III drugs carrying this exact motif at this position.

    The clinical library is built from Thera-SAbDab entries at Approved or
    Phase-III only, so every match is a drug that reached at least Phase III
    carrying this liability. Retention is not proof of tolerance, but for a
    framework site it is the most directly relevant evidence available.
    """
    index = {(s_.get('query_chain'), s_.get('position'), s_.get('kind')): s_
             for s_ in (data.get('clinical_framework_precedents') or {}).get('sites', [])}
    position = site.get('sequence_position', site.get('start'))
    entries = [index[key] for key in
               ((site.get('chain'), position, f.get('kind')) for f in (site.get('flags') or [site]))
               if key in index]
    if not entries:
        return {'searched': False, 'approved': 0, 'phase_iii': 0, 'drugs': [],
                'note': 'Not searched: the site is outside the numbered region the library covers.'}
    matches = [m for entry in entries for m in entry.get('matches', [])]
    approved = [m for m in matches if m.get('stage') == 'Approved']
    phase_iii = [m for m in matches if m.get('stage') == 'Phase-III']
    drugs, seen = [], set()
    for match in approved + phase_iii:
        if match['drug'] not in seen:
            seen.add(match['drug']); drugs.append(f"{match['drug']} ({match['stage']})")
    return {'searched': True, 'approved': len({m['drug'] for m in approved}),
            'phase_iii': len({m['drug'] for m in phase_iii}), 'drugs': drugs[:6],
            'note': ('No approved or Phase-III drug in the library carries this motif at this '
                     'position with a comparable framework.' if not matches else
                     'Retention in clinical-stage products does not prove PTM tolerance.')}


def build_assessment(data):
    structure = data.get('structure_analysis') or {}
    patch_risk = assess_patches(data)
    all_sites = structure.get('ptm_sites', [])
    sites = sorted(filter(focused, all_sites), key=lambda s: (s.get('priority') != 'high',
                   not str(s.get('region', '')).startswith('CDR'), s.get('chain', ''), -(s.get('relative_sidechain_sasa') or 0), s.get('sequence_position', s.get('start', 0))))
    # Every site carries its approved/Phase-III comparator result, so a framework
    # liability is always reported next to the drugs that already carry it.
    for site in sites:
        site['comparators'] = comparator_summary(site, data)
    unknown = len(structure.get('unassessed_ptm_sites', [])) + sum(not assessed(s) for s in all_sites)
    if not structure:
        unknown = len({(label, s['start']) for label, chain in data.get('chains', {}).items()
                       for s in chain.get('metrics', {}).get('ptm_sites', [])})
    counts = dict(high=sum(s.get('priority') == 'high' for s in sites),
                  medium=sum(s.get('priority') == 'review' for s in sites), unknown=unknown,
                  deferred=sum(assessed(s) and not focused(s) for s in all_sites))
    state = 'Partial' if structure and (unknown or structure.get('mapping_warnings')) else 'Structure assessed' if structure else 'Preliminary'
    examples = '; '.join(
        f"{site_label(s['chain'], s.get('number'), s.get('residue'), s.get('sequence_position', s.get('start')))}"
        f" ({s.get('region')}, {s.get('motif', s.get('kind', 'PTM'))})" for s in sites[:3])
    first = (f"The structure-informed screen identifies {counts['high']} high- and {counts['medium']} medium-priority PTM sites"
             + (f", led by {examples}. " if examples else '. ')) if structure else 'Structure and exposure results are unavailable; PTM and surface conclusions remain preliminary. '
    if unknown: first += f"{unknown} PTM sites remain unassessed. "
    if structure:
        names={'hydrophobic':'Hydrophobic','positive':'Positive-charge','negative':'Negative-charge'}
        first += ' '.join(f"{names[p['kind']]} patches: {('low comparative risk (Herceptin-like)' if p.get('herceptin_like') else dict(high='high follow-up priority', review='medium follow-up priority', lower='lower concern', unknown='unresolved')[p['priority']])}{' on the incomplete model' if p.get('provisional') else ''}." for p in patch_risk['categories'])+' '
    anti = [label for label, chain in data.get('chains', {}).items() if chain.get('antiberty')]
    pm = data.get('propermab')
    if not all(c.get('antiberty') for c in data.get('chains',{}).values()) or (set(data.get('chains',{}))=={'H','L'} and not pm):
        state = 'Preliminary · profilers incomplete' if state!='Partial' else 'Partial · profilers incomplete'
    if pm and pm.get('structure_partial'):
        first += 'PROPERMAB structural comparisons are withheld because the model is incomplete. '
    if pm and pm.get('terminal_repair'):
        first += ('Rebuilt terminal coordinates are included and remain uncertain. ' if pm.get('structure_selection')=='terminal_repaired' else 'Terminal repair was not selected; the original model is retained. ')
    first += 'These calls prioritize testing; they do not predict modification rates or functional loss.'
    keys = {(s.get('chain'), s.get('sequence_position', s.get('start')), f.get('kind'))
            for s in sites for f in (s.get('flags') or [s])}
    clinical = [s for s in (data.get('clinical_framework_precedents') or {}).get('sites', [])
                if (s.get('query_chain'), s.get('position'), s.get('kind')) in keys and s.get('matches')]
    documented = []
    for case in (data.get('ptm_precedents') or {}).get('cases', []):
        for match in case.get('matches', []):
            if (match.get('query_chain'), match.get('query_sequence_position'), case.get('chemistry')) in keys:
                documented.append((case, match))
    second = ''
    if clinical:
        s = clinical[0]; m = s['matches'][0]
        second += f"{site_label(s['query_chain'], s.get('number'), s['motif'][:1], s['position'])} {s['motif']} has a retained-motif match in {m['drug']} ({m['stage']}; catalog snapshot). Confirm the liability before mutating solely for motif presence. "
    if documented:
        c, m = documented[0]
        second += f"A matched site was experimentally modified in {c['antibody']} ({c['site']}); experimental conditions and effects are summarized in the cited precedent. "
    if not second:
        second = 'No matched therapeutic evidence is available for the focused sites in the local catalog. ' if structure else 'Therapeutic context will be summarized for priority sites after exposure analysis. '
    second += 'Clinical retention does not prove PTM tolerance. Confirm priority sites by peptide mapping and binding/activity testing before engineering. Use HIC and self-interaction assays to test surface concerns, then assess viscosity at the intended concentration. PK and CMC readiness remain unmeasured.'
    return dict(status=state, patch_risk=patch_risk, paragraphs=[first, second], counts=counts, sites=sites,
                clinical_sites=clinical, documented_matches=[dict(case=c, match=m) for c, m in documented])
