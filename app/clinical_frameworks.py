"""Local retrieval of retained framework and CDR PTM motifs in historical clinical sequences.

Sequence retention is not evidence that the PTM occurred or was clinically acceptable.
"""
from functools import lru_cache
from pathlib import Path
import json
from .numbering import number_chain
from .analysis import ptm_sites
from .germline_cdr import annotate_germline_cdr
from .ptm_precedents import _canonical, _position_order

REFERENCE_PATH = Path(__file__).resolve().parent.parent / 'data' / 'clinical_frameworks.json'
THRESHOLD = 90.0
EVIDENCE = 'Retained residue and motif in deposited clinical sequence; engineering history unknown. PTM occurrence is unproven.'
NOTE = ('Historical highest-stage snapshot from TheraSAbDab, February 2025; Approved is not restricted to FDA approval, and Phase III does not mean currently active. '
        'Only same H/K/L subtype, identical complete motif at aligned canonical IMGT positions, and the same selected region (Framework or CDR1/2/3) in both sequences are required. '
        'Framework identity covers FR1–FR4 using the union of numbered framework positions; gaps are mismatches and CDRs and unnumbered flanks are excluded. '
        'Framework-site retrieval requires at least 90% framework identity, a retrieval rule rather than a biological risk cutoff. CDR-site retrieval requires an exact whole selected CDR sequence at the same canonical positions (union denominator; gaps mismatch), without a framework identity threshold; frameworks may differ. A clinical CDR match does not by itself prove germline origin. All qualifying matches are shown. '
        'Retained sequence does not establish an observed PTM, acceptable modification level, engineering history, or clinical suitability. No match is not evidence of absence.')

@lru_cache(maxsize=1)
def _reference_data():
    from .reference_data import load
    return load('clinical_frameworks')

def _maps(canonical, selected):
    by_seq = {r['sequence_position']:r for r in selected['residues']}
    all_positions = {_canonical(r):r for r in canonical['residues']}
    framework = {key:r for key,r in all_positions.items() if by_seq.get(r['sequence_position'],{}).get('region')=='Framework'}
    return all_positions, framework, by_seq

def clinical_framework_precedents(chains, scheme='imgt'):
    if scheme not in {'imgt','kabat'}:raise ValueError('Numbering scheme must be IMGT or Kabat.')
    dataset=_reference_data()
    result={k:dataset.get(k) for k in ('version','source_url','source_snapshot','source_sha256')}
    result.update(note=NOTE,retrieval_threshold_percent=THRESHOLD,coverage=dict(dataset['coverage']),sites=[],skipped_chains=[])
    references=[]
    for reference in dataset['references']:
        canonical=reference['numbering']['imgt']; selected=reference['numbering'][scheme]
        references.append((reference,canonical,*_maps(canonical,selected)))
    for chain,value in chains.items():
        sequence=value['sequence'] if isinstance(value,dict) else value
        try:
            canonical=number_chain(sequence,'imgt')
            selected=value.get('numbering') if isinstance(value,dict) else None
            if not selected or selected.get('scheme')!=scheme:selected=canonical if scheme=='imgt' else number_chain(sequence,scheme)
        except ValueError:
            result['skipped_chains'].append({'chain':chain,'reason':'Antibody numbering unavailable'})
            continue
        if chain not in {'H','L'} or ('H' if canonical['chain_type']=='H' else 'L') != chain:
            result['skipped_chains'].append({'chain':chain,'reason':'Input label does not match numbered chain family'})
            continue
        query,qfr,qselected=_maps(canonical,selected)
        q_byseq={r['sequence_position']:r for r in canonical['residues']}
        for flag in ptm_sites(sequence, canonical):
            position=flag['start']; motif=flag['motif']
            region=qselected.get(position,{}).get('region')
            if region not in {'Framework','CDR1','CDR2','CDR3'}:continue
            qregion={k:r for k,r in query.items() if qselected.get(r['sequence_position'],{}).get('region')==region}
            motif_rows=[q_byseq.get(p) for p in range(position,position+len(motif))]
            if any(r is None for r in motif_rows):continue
            keys=[_canonical(r) for r in motif_rows]
            site={'query_chain':chain,'position':position,'number':qselected[position]['number'],'query_imgt_number':q_byseq[position]['number'],'query_chain_type':canonical['chain_type'],'query_region':region,'motif':motif,'kind':flag['kind'],'matches':[]}
            if region!='Framework':
                site['germline_cdr']=annotate_germline_cdr(site,canonical,selected)
            for ref,rcanonical,rmap,rfr,rselected in references:
                if rcanonical['chain_type']!=canonical['chain_type']:continue
                aligned=[rmap.get(k) for k in keys]
                if any(r is None for r in aligned) or ''.join(r['aa'] for r in aligned)!=motif:continue
                rp=aligned[0]['sequence_position']
                if flag['kind']=='N-terminal pyroglutamate formation' and (position!=1 or rp!=1):continue
                if [r['sequence_position'] for r in aligned]!=list(range(rp,rp+len(motif))):continue
                if rselected.get(rp,{}).get('region')!=region:continue
                fr_union=set(qfr)|set(rfr)
                matches=sum(k in qfr and k in rfr and qfr[k]['aa']==rfr[k]['aa'] for k in fr_union)
                identity=100*matches/len(fr_union) if fr_union else 0
                rregion={k:r for k,r in rmap.items() if rselected.get(r['sequence_position'],{}).get('region')==region}
                region_union=set(qregion)|set(rregion)
                region_matches=sum(k in qregion and k in rregion and qregion[k]['aa']==rregion[k]['aa'] for k in region_union)
                if region=='Framework' and identity<THRESHOLD:continue
                if region!='Framework' and region_matches!=len(region_union):continue
                match_type=('Exact framework' if matches==len(fr_union) else 'Similar framework and aligned site') if region=='Framework' else f'Exact {region} sequence'
                all_keys=sorted(set(query)|set(rmap),key=_position_order)
                anchor=all_keys.index(keys[0]); window=all_keys[max(0,anchor-5):anchor+6]
                qc=''.join(query[k]['aa'] if k in query else '-' for k in window)
                rc=''.join(rmap[k]['aa'] if k in rmap else '-' for k in window)
                site['matches'].append({'drug':ref['drug'],'stage':ref['stage'],'development_status':ref['development_status'],'reference_chain':ref['reference_chain'],'reference_number':rselected[rp]['number'],'reference_imgt_number':aligned[0]['number'],'reference_position':rp,'reference_region':region,
                    'region_identity_percent':round(100*region_matches/len(region_union),2),'region_matches':region_matches,'region_denominator':len(region_union),
                    'framework_identity_percent':round(identity,2),'framework_matches':matches,'framework_denominator':len(fr_union),'match_type':match_type,
                    'query_context':qc,'reference_context':rc,'local_matches':sum(a==b and a!='-' for a,b in zip(qc,rc)),'local_positions':len(window),'local_numbered_positions':window,'context_anchor_index':window.index(keys[0]),'evidence_note':EVIDENCE})
            site['matches'].sort(key=lambda m:(m['match_type']!='Exact framework',m['stage']!='Approved',-m['framework_identity_percent'],m['drug'],m['reference_chain']))
            result['sites'].append(site)
    result['coverage'].update(query_sites=len(result['sites']),query_framework_sites=sum(s['query_region']=='Framework' for s in result['sites']),query_cdr_sites=sum(s['query_region']!='Framework' for s in result['sites']),matched_sites=sum(bool(s['matches']) for s in result['sites']),total_matches=sum(len(s['matches']) for s in result['sites']))
    return result
