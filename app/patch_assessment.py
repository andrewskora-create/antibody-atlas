"""Evidence-bounded patch triage. No clinical probability or borrowed TAP cutoff."""
import math

SOURCES = [
    {'title':'Raybould et al. 2019 — Therapeutic Antibody Profiler', 'url':'https://doi.org/10.1073/pnas.1810576116'},
    {'title':'Datta-Mannan et al. 2016 — surface engineering, self-association and clearance', 'url':'https://www.nature.com/articles/srep38644'},
    {'title':'Rai et al. 2023 — electrostatic surface and antibody viscosity', 'url':'https://www.nature.com/articles/s41598-023-28841-4'},
]
RULES = ('Local screening policy, not validated risk thresholds: lower = no qualifying local patches and no compatible PROPERMAB upper-range signal; medium = qualifying local patches, or a nonzero compatible PROPERMAB patch area; high = at least one compatible PROPERMAB patch-area descriptor above the maximum of its small selected reference panel. High means prioritize follow-up, not predicted failure. Incomplete mapping or exposure with observed patches receives provisional medium follow-up; absence of a signal on an incomplete model remains unknown. Whole-input risk and reference comparisons remain unresolved. CDR proximity informs mutation caution, not automatic escalation. No TAP cutoffs are applied.')
CONFIG = {
 'hydrophobic': ('Hydrophobic patches', 'Potential self-association, aggregation and formulation liabilities.', 'Measure HIC, self-interaction and aggregation under relevant stress; assess viscosity at the intended concentration.', {'hyd_patch_area','hyd_patch_area_cdr'}),
 'positive': ('Positive-charge patches', 'Potential nonspecific interactions and associated clearance liabilities; PK cannot be inferred from patch size.', 'Test nonspecific binding and self-interaction; characterize PK experimentally and assess pH/salt dependence.', {'pos_patch_area','pos_patch_area_cdr'}),
 'negative': ('Negative-charge patches', 'Potential charge-mediated self-association and viscosity; negative charge alone is not adverse.', 'Measure self-interaction and concentration-dependent viscosity across relevant pH/salt conditions. Full-IgG/Fc interactions are not modeled here.', {'neg_patch_area','neg_patch_area_cdr'}),
}


def assess_patches(data):
    s=data.get('structure_analysis') or {}
    coverage=s.get('mapping_coverage_percent') or {}
    expected=set(data.get('chains') or coverage)
    complete=bool(s and expected and expected <= set(coverage)) and all(
        isinstance(coverage[k],(int,float)) and math.isfinite(coverage[k]) and coverage[k]>=100 for k in expected)
    complete=complete and not s.get('mapping_warnings')
    if any(r.get('aa')!='G' and (not isinstance(r.get('relative_sidechain_sasa'),(int,float)) or not math.isfinite(r['relative_sidechain_sasa'])) for r in s.get('residues', [])):
        complete=False
    pm=data.get('propermab') or {}
    # Only use comparator evidence attached to this same, fully mapped model.
    same_model=bool(pm.get('input_pdb_sha256') and pm['input_pdb_sha256']==data.get('provenance',{}).get('structure_pdb_sha256'))
    compatible=complete and same_model and not pm.get('structure_partial')
    tolerance=data.get('assessment_policy',{}).get('herceptin_tolerance_percent',10)
    if isinstance(tolerance,bool) or not isinstance(tolerance,(int,float)) or not math.isfinite(tolerance) or not 0 <= tolerance <= 20:
        tolerance=10
    exposure_complete=bool(expected)
    for chain in expected:
        sequence=(data.get('chains',{}).get(chain) or {}).get('sequence','')
        rows=[r for r in s.get('residues',[]) if r.get('chain')==chain]
        positions=[r.get('sequence_position') for r in rows]
        exposure_complete=exposure_complete and bool(sequence) and len(rows)==len(sequence) and set(positions)==set(range(1,len(sequence)+1))
        exposure_complete=exposure_complete and all(isinstance(r.get('sequence_position'),int) and 1<=r['sequence_position']<=len(sequence) and r.get('aa')==sequence[r['sequence_position']-1] and (r.get('aa')=='G' or isinstance(r.get('relative_sidechain_sasa'),(int,float)) and math.isfinite(r['relative_sidechain_sasa'])) for r in rows)
    results=[]
    for kind,(label,concern,action,keys) in CONFIG.items():
        patches=[p for p in s.get('patches',[]) if p.get('kind')==kind]
        refs=[b for b in pm.get('summary',{}).get('benchmarks',[]) if compatible and b.get('key') in keys and b.get('reference_n',0)>0 and b.get('query_relation') in ('above_range','within_range','below_range')]
        above=[b for b in refs if b['query_relation']=='above_range']
        measured_area=any(same_model and f.get('key') in keys and isinstance(f.get('value'),(int,float)) and math.isfinite(f['value']) and f['value']>0 for f in pm.get('features',[]))
        positive_area=measured_area or any(isinstance(b.get('query_value'),(int,float)) and b['query_value']>0 for b in refs)
        observed='high' if above else 'review' if patches or positive_area else 'lower'
        priority=observed if complete else 'review' if patches or positive_area else 'unknown'
        anchors={b['key']:b for b in refs if b.get('reference_anchor')}
        def near_anchor(b):
            baseline=b['reference_anchor'].get('value'); value=b.get('query_value')
            return all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) for v in (baseline,value)) and baseline>0 and abs(value-baseline)<=baseline*tolerance/100+1e-9
        herceptin_like=bool(tolerance and compatible and exposure_complete and not above and keys <= anchors.keys() and all(near_anchor(anchors[k]) for k in keys))
        if herceptin_like: priority='lower'
        reasons=[]
        largest=max((p.get('exposed_sidechain_sasa_sum_a2',0) for p in patches),default=0)
        if s: reasons.append(f'{len(patches)} local patches; largest summed side-chain accessible area {largest:g} Å²; {sum(bool(p.get("near_cdr")) for p in patches)} near a CDR.')
        if above:
            reasons.append('Above compatible reference maximum: '+', '.join(b.get('label',b['key']) for b in above)+'. The reference minimum–maximum is not an acceptance range; panel members were not selected for favorable experimental behavior.')
        elif herceptin_like:
            reasons.append(f'Both total and CDR-adjacent patch areas are within ±{tolerance:g}% of Herceptin under the compatible local calculation protocol: '+', '.join(f"{anchors[k].get('label',k)} {100*(anchors[k]['query_value']-anchors[k]['reference_anchor']['value'])/anchors[k]['reference_anchor']['value']:+.1f}%" for k in sorted(keys))+'. Low comparative risk is a user-selected screening policy, not an experimentally validated cutoff or a prediction of HIC, PK or viscosity.')
        elif refs: reasons.append('These PROPERMAB patch areas are at or below the largest values in the matched therapeutic reference panel. This is descriptive similarity, not a favorable experimental result or a low-risk call.')
        else: reasons.append('No compatible PROPERMAB patch comparison is available.')
        if positive_area and not patches: reasons.append('PROPERMAB detects nonzero patch area even though the local contact-graph screen detects no qualifying clusters; these are different surface definitions.')
        if not complete:
            missing=[]
            for chain_label,info in s.get('chain_mapping',{}).items():
                mapped={r.get('sequence_position') for r in s.get('residues',[]) if r.get('chain')==chain_label}
                positions=sorted(set(range(1,info.get('input_residues',0)+1))-mapped)
                if positions:missing.append(chain_label+' input positions '+', '.join(map(str,positions)))
            reasons.append(('Unmodeled residues: '+'; '.join(missing)+'. ' if missing else '')+'Observed patch findings remain usable provisionally. Whole-input risk is unresolved; missing residues and altered termini can change charge and surface patches.')
        elif not patches and not positive_area: reasons.append('No patch signal meets this screen’s criteria; other liabilities remain possible.')
        results.append(dict(kind=kind,label=label,priority=priority,provisional=not complete,observed_priority=observed if s else 'unknown',
                            herceptin_like=herceptin_like,herceptin_tolerance_percent=tolerance,
                            call=('Low comparative risk — Herceptin-like' if herceptin_like else 'Provisional · medium follow-up priority' if not complete and priority=='review' else {'high':'High follow-up priority','review':'Medium follow-up priority','lower':'Lower concern in this screen','unknown':'Partial assessment — no assignable priority'}[priority]),
                            count=len(patches) if s else None,largest_area_a2=largest if s else None,
                            reasons=reasons,concern=concern,action=action,comparisons=refs))
    return dict(categories=results,definition=RULES+f' Herceptin-like low comparative risk requires both total and CDR-adjacent areas within ±{tolerance:g}% of the compatible Herceptin reference, complete structure coverage, and no above-reference maximum signal. This configurable policy does not alter PTM priorities; zero disables it.',herceptin_tolerance_percent=tolerance,sources=SOURCES,complete=complete,
                limitation='Heuristic follow-up priorities, not probabilities of aggregation, clearance or high viscosity. Single Fv conformation; formulation, concentration, full-IgG context and experimental outcomes are not modeled.')
