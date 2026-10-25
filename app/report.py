"""Four-page decision summary. Full evidence remains in the local application."""
from datetime import datetime, timezone
from .assessment import build_assessment, imgt_position, site_label
from .attribution import line as method_line, citation as method_citation
from html import escape
from io import BytesIO
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph

NAVY='#163449'; TEAL='#2a7187'; MUTED='#5d7180'; PALE='#edf5f6'
RISK={'high':'#a8505b','review':'#94601f','lower':'#397d6e','unassessed':MUTED}
ORDER={'high':0,'review':1,'unassessed':2,'lower':3}


LIABILITY_SHORT = {
    'Asparagine deamidation motif': 'Deamidation',
    'Aspartate isomerization motif': 'Isomerization',
    'Oxidation-susceptible residue': 'Oxidation',
    'N-linked glycosylation sequon': 'Glyc. sequon',
    'N-terminal pyroglutamate formation': 'Pyro-Glu',
    'Asp-Pro cleavage motif': 'Asp-Pro',
    'Noncanonical cysteine': 'Noncanon. Cys',
}


def priority_sites(data):
    structure=data.get('structure_analysis') or {}
    sites=list(structure.get('ptm_sites',[]))+list(structure.get('unassessed_ptm_sites',[]))
    if not sites:
        for label,chain in data.get('chains',{}).items():
            numbered={r['sequence_position']:r for r in chain['numbering']['residues']}
            for flag in chain['metrics'].get('ptm_sites',[]):
                r=numbered.get(flag['start'],{})
                sites.append({**flag,'chain':label,'sequence_position':flag['start'],
                              'number':r.get('number',str(flag['start'])),'region':r.get('region','Unassigned'),
                              'priority':'unassessed'})
    # Count physical residues once even when multiple chemistry flags overlap.
    unique={}
    for site in sites:
        key=(site.get('chain'),site.get('sequence_position',site.get('start')))
        if key not in unique or ORDER.get(site.get('priority','unassessed'),2)<ORDER.get(unique[key].get('priority','unassessed'),2):
            unique[key]=site
    return sorted(unique.values(),key=lambda s:(ORDER.get(s.get('priority','unassessed'),2),
                  not str(s.get('region','')).startswith('CDR'),-(s.get('relative_sidechain_sasa') or 0),
                  s.get('chain',''),s.get('sequence_position',s.get('start',0))))


def make_report(data, title='Antibody developability report', structure_image=None):
    """Exactly five fixed-layout pages, bounded selections, explicit omissions."""
    output=BytesIO(); c=canvas.Canvas(output,pagesize=(612,792))
    c.setTitle(str(title)[:100]); c.setAuthor('Antibody Atlas')
    chains=data.get('chains',{}); structure=data.get('structure_analysis') or {}
    assessment=build_assessment(data)
    sites=assessment['sites']; counts={k:sum(s.get('priority','unassessed')==k for s in sites) for k in ORDER}
    pm=data.get('propermab') or {}; provenance=data.get('provenance') or {}
    scheme=str(data.get('scheme','imgt')).upper()

    def text(value,x,y,width=516,height=40,size=10,color=NAVY,bold=False,truncate=True):
        # Bound each block without compressing the type or spilling onto another page.
        value=str(value).replace('\u2011','-'); value=value.replace('\n',' ')
        style=ParagraphStyle('block',fontName='Helvetica-Bold' if bold else 'Helvetica',fontSize=size,
                             leading=size*1.35,textColor=colors.HexColor(color))
        def paragraph(s):return Paragraph(escape(s),style)
        p=paragraph(value); _,h=p.wrap(width,height)
        if h>height:
            if not truncate: raise ValueError('Assessment exceeds its reserved report space')
            lo,hi=0,len(value)
            while lo<hi:
                mid=(lo+hi+1)//2
                if paragraph(value[:mid]+'...').wrap(width,height)[1]<=height:lo=mid
                else:hi=mid-1
            p=paragraph(value[:lo].rstrip()+'...');_,h=p.wrap(width,height)
        p.drawOn(c,x,y-h)
        return h

    def heading(label,y):
        text(label,48,y,height=24,size=14,bold=True)
        c.setStrokeColor(colors.HexColor('#dce7e9'));c.line(48,y-28,564,y-28)

    def page(number,label,total=5):
        c.setFillColor(colors.HexColor(TEAL));c.rect(0,780,612,12,fill=1,stroke=0)
        text('ANTIBODY ATLAS  /  DECISION SUMMARY',48,758,height=15,size=9,color=TEAL,bold=True)
        text(label,48,729,height=35,size=23)
        c.setStrokeColor(colors.HexColor('#dce7e9'));c.line(48,43,564,43)
        text(f'Local research screening | {scheme} | {number} / {total}',48,31,height=15,size=8,color=MUTED)
        text('Full evidence and detailed results remain in the app.',280,31,width=284,height=15,size=8,color=MUTED)

    def row(values,y,widths,height=34,header=False):
        c.setFillColor(colors.HexColor(NAVY if header else '#f2f7f8'));c.rect(48,y-height,516,height,fill=1,stroke=0)
        x=48
        for val,w in zip(values,widths):
            text(val,x+7,y-7,width=w-14,height=height-10,size=9,color='#ffffff' if header else NAVY,bold=header);x+=w

    def compact_row(values,y,widths,height=22,header=False,tint=None):
        fill=NAVY if header else (tint or '#f2f7f8')
        c.setFillColor(colors.HexColor(fill));c.rect(48,y-height,516,height,fill=1,stroke=0)
        x=48
        for val,w in zip(values,widths):
            text(val,x+5,y-5,width=w-10,height=height-10,size=8,
                 color='#ffffff' if header else NAVY,bold=header);x+=w

    page(1,'Developability overview')
    text(str(title)[:100],48,683,height=31,size=12,bold=True)
    text(f"{datetime.now(timezone.utc).strftime('%Y-%m-%d')} | Input ID: {provenance.get('input_sha256','not recorded')[:16]} | Charge pH: {data.get('ph','not recorded')}",48,646,height=20,size=9,color=MUTED)
    heading('Overall assessment · ' + assessment['status'].lower(),620)
    used=text(assessment['paragraphs'][0],48,579,height=140,size=9.5,truncate=False)
    y=579-used-12
    text(assessment['paragraphs'][1],48,y,height=y-354,size=9.5,truncate=False)
    text('Predicted hydrophobic surface',48,333,width=300,height=22,size=13,bold=True)
    if structure_image:
        c.drawImage(ImageReader(BytesIO(structure_image)),48,107,width=302,height=198,preserveAspectRatio=True,anchor='c',mask='auto')
    else:
        c.setFillColor(colors.HexColor(PALE));c.roundRect(48,115,302,192,8,fill=1,stroke=0)
        text('Structure image unavailable. Export from the app after local structure prediction to include the surface view.',68,230,width=262,height=70,size=11,color=MUTED)
    stops=['#28669a','#6bb4d1','#e7d779','#ee9856','#d64354','#b91f43']
    for i,color in enumerate(stops):
        c.setFillColor(colors.HexColor(color));c.rect(48+i*30,99,30,7,fill=1,stroke=0)
    text('Lower',48,95,width=90,height=13,size=8,color=MUTED)
    text('Higher signal',162,95,width=100,height=13,size=8,color=MUTED)
    text('Relative exposed hydrophobic signal; colors are not measured aggregation risk or the PROPERMAB patch surface.',48,79,width=302,height=30,size=8,color=MUTED)
    text('Variable-domain properties',370,333,width=194,height=22,size=11,bold=True)
    for i,(label,chain) in enumerate(list(chains.items())[:2]):
        h=chain.get('humanness',{});m=chain['metrics'];top=300-i*98
        text(f"{label} chain · human framework {h.get('percent','-')}%",370,top,width=194,height=28,size=10,bold=True)
        text(f"{h.get('nearest_human_v_gene','Unavailable')} | {h.get('matches','-')}/{h.get('framework_positions','-')} matches",370,top-28,width=194,height=28,size=9)
        text(f"pI {m.get('isoelectric_point','-')} | Charge {m.get('net_charge_at_ph','-')}\nGRAVY {m.get('gravy','-')}",370,top-56,width=194,height=29,size=9)
    text('Identity uses IMGT FR1-FR3 against OGRDB human germline alleles (Lees et al., 2020). This local screen guides follow-up; clinical acceptance limits are not established.',370,112,width=194,height=58,size=8,color=MUTED)
    c.showPage()

    page(2,'Priority sites and surface findings')
    text(f"Focused PTMs: {assessment['counts']['high']} high | {assessment['counts']['medium']} medium | {assessment['counts']['deferred']} deferred | {assessment['counts']['unknown']} unassessed",48,682,height=30,size=11,bold=True)
    shown=sites[:6]
    text(f'Showing {len(shown)} of {len(sites)} unique flagged residues, ranked by screening priority, CDR location and exposure. Site is chain, residue and input position; the second column gives the region and the canonical {scheme} position. Framework sites are capped below high and searched against approved and Phase-III drugs.',48,650,height=42,size=9,color=MUTED)
    widths=[66,84,122,42,92,110]
    row(['Site',f'Region / {scheme}','Chemistry / motif','Exp.','Approved / Ph3','Priority'],598,widths,31,True)
    for i,s in enumerate(shown):
        value=s.get('relative_sidechain_sasa');exp=f'{value:.0%}' if isinstance(value,(int,float)) else '?'
        kind=str(s.get('kind','PTM flag')).replace('Asparagine deamidation motif','Deamidation').replace('Aspartate isomerization motif','Isomerization').replace('Oxidation-susceptible residue','Oxidation').replace('N-linked glycosylation sequon','N-glycosylation')
        priority=s.get('priority','unassessed')
        comp=s.get('comparators') or {}
        if not comp.get('searched'):cell='not searched'
        elif comp.get('approved') or comp.get('phase_iii'):cell=f"{comp.get('approved',0)} appr / {comp.get('phase_iii',0)} P3"
        else:cell='none found'
        # Chain and residue first: an IMGT label like 'K56' means kappa chain,
        # position 56, not lysine. Residue and input position remove the clash.
        site_cell=f"{s.get('chain','?')} {s.get('residue','')}{s.get('sequence_position',s.get('start','?'))}"
        region_cell=f"{s.get('region','Unknown')} · {imgt_position(s.get('number'))}"
        row([site_cell,region_cell,f"{kind}: {s.get('motif','')}",exp,cell,priority.upper()],567-i*34,widths,34)
        c.setFillColor(colors.HexColor(RISK.get(priority,MUTED)));c.rect(48,533-i*34,3,34,fill=1,stroke=0)
    if not shown:text('No assessed high or medium sites in the focused screen. Deferred and unassessed sites remain in the app; this does not exclude liabilities.',48,556,height=38)
    table_bottom=567-len(shown)*34 if shown else 510
    text(method_line('ptm_priority'),48,table_bottom-6,width=516,height=34,size=7.5,color=MUTED)
    url=(method_citation('ptm_priority') or {}).get('primary_url')
    cmp_cite=method_citation('framework_comparator') or {}
    text('Comparator search: '+ (cmp_cite.get('basis') or '') + '. ' + (cmp_cite.get('primary_title') or ''),
         48,table_bottom-38,width=516,height=22,size=7.5,color=TEAL)
    if cmp_cite.get('primary_url'):
        c.linkURL(cmp_cite['primary_url'],(48,table_bottom-60,516,table_bottom-38),relative=0)

    heading('Surface findings',table_bottom-72)
    base=table_bottom-112
    for i,risk in enumerate(assessment['patch_risk']['categories']):
        priority={'high':'HIGH','review':'MEDIUM','lower':'LOWER','unknown':'UNKNOWN'}[risk['priority']]
        if risk.get('herceptin_like'):priority='LOW COMPARATIVE (Herceptin-like)'
        if risk.get('provisional') and risk['priority']=='review':priority+=' (provisional)'
        label={'hydrophobic':'Hydrophobic','positive':'Positive charge','negative':'Negative charge'}[risk['kind']]
        area=risk['largest_area_a2']
        detail=f"{label}: {priority}. {risk['count']} patches; largest {area:g} A^2." if area is not None else f"{risk['label']}: UNKNOWN - structure unavailable."
        if risk['priority']=='review' and risk['count']==0:
            detail=f"{label}: {priority}. PROPERMAB area detected; 0 local clusters."
        text(detail,48,base-i*21,height=20,size=10)
    patch_cite=method_citation('surface_patches') or {}
    text(method_line('surface_patches'),48,base-70,width=516,height=34,size=7.5,color=MUTED)
    if patch_cite.get('primary_url'):
        text(patch_cite.get('primary_title') or '',48,base-106,width=516,height=20,size=7.5,color=TEAL)
        c.linkURL(patch_cite['primary_url'],(48,base-126,516,base-106),relative=0)
    text('Risk key: HIGH = prioritized predicted chemistry/exposure concern (CDR only); REVIEW = follow-up warranted; LOWER = lower concern in this screen; UNASSESSED = insufficient structural evidence. Exposure informs priority but does not predict PTM rate. Within reference min-max alone is not a pass.',48,base-130,width=516,height=44,size=8,color=MUTED)
    c.showPage()

    # ---------------- biophysical and formulation developability -------------
    # Vertical budget: table top 690, footer rule 43. Rows are 25pt so the 9pt
    # text (12.15pt leading) fits inside the 15pt inner box without truncating.
    page(3,'Biophysical and formulation properties')
    benchmarks={b.get('key'):b for b in (pm.get('summary') or {}).get('benchmarks',[]) if b.get('key')}
    features={f.get('key'):f for f in pm.get('features',[]) if f.get('key')}

    def descriptor(key):
        bench=benchmarks.get(key); feature=features.get(key)
        value=(bench or {}).get('query_value',(feature or {}).get('value'))
        if not isinstance(value,(int,float)):return None
        out={'value':f"{value:,.4g}",'range':'-','position':'no panel'}
        # Helvetica drops the superscript glyphs, so units go in as ASCII,
        # matching the 'A^2' already used on the surface-findings page.
        unit=(bench or feature or {}).get('unit') or ''
        unit=unit.replace('\u00c5\u00b2','A^2').replace('\u00b2','^2').replace('\u00c5','A')
        if unit:out['value']+=f' {unit}'
        if bench and isinstance(bench.get('min'),(int,float)) and bench.get('reference_n'):
            out['range']=f"{bench['min']:g} to {bench['max']:g} (n={bench['reference_n']})"
            out['position']={'within_range':'within','above_range':'ABOVE',
                             'below_range':'below'}.get(bench.get('query_relation'),'-')
            # A percentile says more than within/above/below: 'within' covers both
            # a molecule sitting with the bulk and one pressed against the wall.
            # Only present once the panel is large enough for it to mean anything.
            percentile=bench.get('percentile')
            if isinstance(percentile,(int,float)):
                out['position']=('ABOVE' if out['position']=='ABOVE'
                                 else f"P{percentile:g}")
                if out['position'].startswith('P') and (percentile>=90 or percentile<=10):
                    out['position']+=' tail'
        return out

    def per_chain(metric):
        parts=[f"{name} {chain['metrics'][metric]:g}" for name,chain in list(chains.items())[:2]
               if isinstance(chain.get('metrics',{}).get(metric),(int,float))]
        return {'value':' | '.join(parts),'range':'-','position':'no panel'} if parts else None

    flags=data.get('molecule_flags') or {}
    fv=flags.get('fv') or {}
    full=flags.get('basis')=='full_length'
    basis=('Assembled '+str((flags.get('assembly') or {}).get('format','IgG'))
           if full else 'Variable-domain molecule')
    def molecular_value(key, fallback):
        value=fv.get(key)
        return ({'value':f'{value:g}', 'range':'-', 'position':'no panel'}
                if isinstance(value,(int,float)) else per_chain(fallback))
    GROUPS=[('Charge',[((basis if isinstance(fv.get('pi'),(int,float)) else 'Variable-domain chain')+' pI',molecular_value('pi','isoelectric_point')),
                       ((basis if isinstance(fv.get('net_charge_assay_ph'),(int,float)) else 'Variable-domain chain')+' charge (e)',molecular_value('net_charge_assay_ph','net_charge_at_ph')),
                       ('Fv formal charge',descriptor('fv_charge')),
                       ('Charge symmetry (VH x VL)',descriptor('fv_csp')),
                       ('VH minus VL charge',descriptor('fv_chml'))]),
            ('Hydrophobicity',[('GRAVY',per_chain('gravy')),
                               ('Apolar accessible area',descriptor('hyd_asa')),
                               ('Heiden hydrophobic score',descriptor('heiden_score'))]),
            ('Aggregation and self-association',[
                               ('Hydrophobic patch area',descriptor('hyd_patch_area')),
                               ('Hydrophobic patch area, CDR',descriptor('hyd_patch_area_cdr')),
                               ('Positive patch area',descriptor('pos_patch_area')),
                               ('Negative patch area',descriptor('neg_patch_area'))]),
            ('High-concentration behaviour',[
                               ('Spatial charge map (SCM)',descriptor('scm')),
                               ('Negative patch area, CDR',descriptor('neg_patch_area_cdr'))])]

    # DeepSP predicts the same spatial properties from sequence alone. It is
    # grouped separately rather than merged into the rows above because its
    # reference distribution is a different thing: natural repertoire Fvs, not
    # the clinical panel, and a predicted value rather than a computed one.
    deepsp=data.get('deepsp') or {}
    if deepsp.get('status') in {'ok','unreliable'}:
        deepsp_rows=[]
        for highlight in deepsp.get('highlights',[]):
            percentile=highlight.get('percentile')
            deepsp_rows.append((highlight.get('short_label') or highlight['label'][:44],{
                'value':f"{highlight['value']:,.4g}",
                'range':'-' if percentile is None
                        else f"{(deepsp.get('reference') or {}).get('n','?'):,} repertoire Fvs",
                'position':'-' if percentile is None
                           else ('ABOVE' if highlight.get('flagged') else f"P{percentile:g}")}))
        if deepsp_rows:
            GROUPS.append(('Predicted from sequence (DeepSP)',deepsp_rows))

    widths=[214,112,106,84]
    y=690
    row(['Descriptor','Value','Panel range','Position'],y,widths,25,True); y-=25
    shown=0
    for group,items in GROUPS:
        rows=[(label,info) for label,info in items if info]
        if not rows:continue
        c.setFillColor(colors.HexColor(PALE));c.rect(48,y-18,516,18,fill=1,stroke=0)
        text(group,55,y-3,width=500,height=14,size=8.5,color=TEAL,bold=True); y-=18
        for label,info in rows:
            shown+=1
            row([label,info['value'],info['range'],info['position']],y,widths,25)
            if info['position']=='ABOVE':
                c.setFillColor(colors.HexColor(RISK['high']));c.rect(48,y-25,3,25,fill=1,stroke=0)
            y-=25
        y-=5

    if not shown:
        text('No biophysical descriptors were computed for this run.',48,y,height=28,size=10);y-=34

    has_percentile=any(isinstance(b.get('percentile'),(int,float)) for b in benchmarks.values())
    note=('Panel range is the observed minimum to maximum across the matched therapeutic reference set under an '
          'identical protocol. It is not an acceptance interval; panel members were not selected for favourable '
          'experimental behaviour. ABOVE means prioritise follow-up, not predicted failure.')
    if has_percentile:
        note+=(' Pnn is the position within that panel: P90 means higher than 90% of the reference molecules. '
               'A percentile is a position among approved antibodies, not a verdict on acceptability, and '
               '"tail" marks the outer tenth at either end.')
    if deepsp.get('status') in {'ok','unreliable'}:
        note+=(' DeepSP rows are predicted from sequence, not computed from a structure, and their position is a '
               'percentile of natural repertoire Fv pairs rather than of the clinical panel. ABOVE there means the '
               '95th percentile of that repertoire distribution.')
    if not benchmarks:
        note=('Structure-derived surface, charge-distribution and viscosity descriptors require the optional '
              'PROPERMAB component, which is not installed here. ')
        note+=('The DeepSP rows above predict spatial aggregation and charge from sequence alone and need no '
               'structure and no academic licence, but they are predictions of a simulation result, not '
               'measurements.' if deepsp.get('status') in {'ok','unreliable'} else
               'Only sequence-level charge and hydrophobicity are shown. See USE_TIERS.md for alternatives, '
               'including the free Therapeutic Antibody Profiler and the MIT-licensed DeepSP, which this project '
               'can install with: python setup/fetch_deepsp.py')
    if deepsp.get('status')=='unreliable':
        text('DeepSP prediction is NOT usable for this molecule: '+str(deepsp.get('warning','')),
             48,282,width=516,height=34,size=8,color=colors.HexColor(RISK['high']))
    text(note,48,240,width=516,height=38,size=8,color=MUTED)
    bio_cite=method_citation('biophysical') or {}
    bio_line=method_line('biophysical')
    if deepsp.get('status') in {'ok','unreliable'}:
        bio_line+=(' DeepSP rows: convolutional network over the IMGT-aligned paired Fv '
                   'sequence (DeepSP, ANARCII).')
    text(bio_line,48,200,width=516,height=34,size=7.5,color=MUTED)
    link_cite=bio_cite
    if not benchmarks and deepsp.get('status') in {'ok','unreliable'}:
        link_cite=method_citation('deepsp') or bio_cite
    if link_cite.get('primary_url'):
        text(link_cite.get('primary_title') or '',48,170,width=516,height=20,size=7.5,color=TEAL)
        c.linkURL(link_cite['primary_url'],(48,150,516,170),relative=0)

    heading('What these do and do not establish',146)
    text('Charge and surface descriptors inform clearance, viscosity and self-association risk; none of them predicts '
         'it. Both high positive and high negative net charge have been associated with faster clearance, so a low pI '
         'is not automatically favourable. Viscosity is concentration-dependent and must be measured at the intended '
         'formulation concentration. Confirm with HIC, self-interaction (AC-SINS or equivalent), concentration-dependent '
         'viscosity, accelerated and freeze-thaw stability, and a PK study. No thermal stability, solubility, '
         'immunogenicity or expression prediction is made here.',48,112,width=516,height=64,size=8.5)
    c.showPage()

    page(4,'Relevant precedents and next steps')
    clinical=data.get('clinical_framework_precedents') or {}
    cases=[item['case'] for item in assessment['documented_matches']]
    clinical_sites=assessment['clinical_sites']
    heading('Evidence relevant to this sequence',678)
    entries=[]
    for case in cases[:2]:
        m=next(item['match'] for item in assessment['documented_matches'] if item['case'] is case)
        detail=(f"Documented PTM: {case['antibody']} {case['site']}; matches input "
                f"{site_label(m.get('query_chain',''), m.get('query_number'))} {m['query_region']}. {case['observation']}")
        if case.get('terminal_requirement'):detail+=' Only applies to a free mature-chain N-terminus.'
        entries.append((detail,case.get('url'),case.get('citation') or 'Published study'))
    for site in clinical_sites:
        if len(entries)>=3:break
        m=site['matches'][0];g=site.get('germline_cdr') or {}
        detail=(f"Retained clinical motif: input "
                f"{site_label(site.get('query_chain',''), site.get('number'), site['motif'][:1], site['position'])} "
                f"{site.get('query_region','')} {site['motif']} matches {m['drug']} ({m['stage']}, source snapshot). "
                f"{m['match_type']}; framework identity {m['framework_identity_percent']}%.")
        if g.get('status')=='exact_match':detail+=' CDR also matches human V-gene germline reference.'
        entries.append((detail,clinical.get('source_url'),
                        'Thera-SAbDab clinical-stage sequences (Raybould et al., 2020)'))

    # Flow each entry by its measured height instead of a fixed slot, so short
    # entries do not leave a gap before their own source link.
    y=640
    if not entries:
        y-=text('No matching clinical/experimental evidence is available in the local catalog for the focused sites. Absence of a match is not proof of safety.',48,y,height=46)+14
    for i,(detail,url,label) in enumerate(entries,1):
        used=text(f'{i}. {detail}',48,y,width=516,height=58,size=9.5)
        y-=used+3
        if url:
            text(label,62,y,width=502,height=14,size=8,color=TEAL)
            c.linkURL(url,(62,y-14,564,y),relative=0)
            y-=16
        y-=10

    y-=4
    y-=text(f"Selected {len(entries)} examples; {len(cases)} documented PTM matches and "
            f"{clinical.get('coverage',{}).get('total_matches',0)} clinical site-reference matches in the app. "
            f"Clinical-stage library covers Approved and Phase-III entries only; snapshot February 2025. "
            f"Retention or germline identity does not prove tolerance or an absence of historical mutation.",
            48,y,width=516,height=42,size=8.5,color=MUTED)

    heading('Recommended next steps',y-22)
    y-=56
    actions=[]
    if sites:actions.append('Confirm priority PTMs by site-resolved peptide mapping under formulation and relevant stress conditions; assess binding and activity before choosing mutations.')
    if not structure:actions.append('Run local structure prediction and exposure analysis before using structural PTM or surface conclusions.')
    else:actions.append('Check predicted exposed patches with self-interaction and aggregation assays; measure viscosity and solubility at the intended formulation concentration.')
    actions.append('Establish stability, potency and formulation data. PK and CMC readiness require experimental evidence beyond this sequence and structure screen.')
    for i,action in enumerate(actions[:3]):
        y-=text(f'{i+1}. {action}',48,y,width=516,height=40,size=9.5)+9

    names=[]
    for key in ('structure','ptm_priority','framework_comparator','surface_patches','biophysical','paratope','substitutions','humanness','language_model'):
        for source in (method_citation(key) or {}).get('methods',[]):
            if source not in names:names.append(source)
    text('Methods used in this report: '+'; '.join(names)+'.',48,132,width=516,height=46,size=7.5,color=MUTED)
    text('Each assessment above carries its own attribution. Full citations, limitations and source links are in the local app and in MODULE_GUIDE.md. Citation is not evidence that a method ran, and no cited paper validates this report\u2019s combined heuristics.',48,80,width=516,height=32,size=7.5,color=MUTED)
    c.showPage()

    # Page 5: the consolidated per-residue table. This is the report's point of
    # contact with the thing the app is built around, so it is a page of its
    # own rather than a column somewhere. A full Fv is ~230 rows and will not
    # fit, so only rows worth attention are printed and the omission is stated.
    page(5,'Per-residue assessment')
    table=data.get('residue_table') or {}
    rows=table.get('rows') or []
    if not rows:
        try:
            from .residue_table import build_residue_table
            table=build_residue_table(chains,structure=structure,paratope=data.get('paratope'))
            rows=table['rows']
        except Exception:
            rows=[]

    if not rows:
        text('The per-residue table could not be built for this record. It requires canonical IMGT '
             'numbering for every chain.',48,660,width=516,height=40,size=10,color=MUTED)
    else:
        counts=table.get('verdict_counts',{})
        missing=table.get('missing_layers') or []
        text(f"{table.get('notable_rows',0)} of {len(rows)} residues are worth attention  |  "
             f"{counts.get('preserve',0)} preserve · {counts.get('evaluate',0)} evaluate · "
             f"{counts.get('watch',0)} watch · {counts.get('clear',0)} clear",
             48,682,height=26,size=11,bold=True)
        alt=next((r.get('alt_scheme') for r in rows if r.get('alt_scheme')),None)
        text('Every row is keyed to the canonical IMGT position'
             + (f', with {alt} alongside.' if alt else '.')
             + (f' Not yet joined: {", ".join(missing)}; without a structure every liability stays '
                'at watch, so this deliberately over-reports.' if missing else ''),
             48,656,width=516,height=26,size=8.5,color=MUTED)

        widths=[40,32,44,104,62,54,50,56,74]
        header=['Chain','AA','IMGT',('Region / '+alt) if alt else 'Region',
                'Germline','Exposure','Contact','Call','Liability']
        y=622
        compact_row(header,y,widths,header=True); y-=22
        # Severity order, not sequence order. The page holds ~21 rows and a full
        # Fv has more notable ones than that, so printing in sequence order fills
        # the page with the heavy chain and can drop an 'evaluate' on the light chain
        # off a decision summary entirely.
        rank={'evaluate':0,'preserve':1,'watch':2,'clear':3}
        notable=sorted((r for r in rows if r.get('notable')),
                       key=lambda r:(rank.get(r.get('verdict'),4), r.get('chain',''),
                                     r.get('sequence_position',0)))
        limit=21
        tints={'evaluate':'#fdf1ef','preserve':'#eef3f5','watch':'#fdf7ec','clear':'#f2f7f8'}
        for r in notable[:limit]:
            region=r.get('imgt_region','')
            if r.get('alt_region') and r['alt_region']!=region:
                region=f"{region} / {r['alt_region']}"
            if r.get('germline_aa') is None: germline='-'
            elif r.get('germline_match'): germline='match'
            else: germline=f"{r['germline_aa']} differs"
            exposure='-' if r.get('exposure') is None else f"{round(r['exposure']*100)}%"
            contact='-' if r.get('paratope_score') is None else f"{r['paratope_score']:.2f}"
            liability='; '.join(LIABILITY_SHORT.get(l.get('kind'),l.get('kind',''))
                                + (f" {l['motif']}" if l.get('motif') else '')
                                for l in (r.get('liabilities') or [])) or '-'
            compact_row([r.get('chain',''),r.get('aa',''),r.get('imgt',''),region,germline,
                         exposure,contact,r.get('verdict',''),liability],
                        y,widths,tint=tints.get(r.get('verdict')))
            y-=22

        omitted=len(notable)-limit
        note=(f"{omitted} further notable residue{'s' if omitted!=1 else ''} and "
              f"{len(rows)-len(notable)} unremarkable ones are not shown; the full table is in the app. "
              if omitted>0 else
              f"{len(rows)-len(notable)} unremarkable residues are not shown; the full table is in the app. ")
        legend=table.get('verdict_legend') or {}
        text(note+'Germline compares against the assigned human V allele; a difference inside CDR1/2 is '
             'expected for a matured antibody and CDR3 has no V-gene reference at all.',
             48,y-8,width=516,height=34,size=8,color=MUTED)
        text('  |  '.join(f"{k}: {v}" for k,v in legend.items() if k in ('preserve','evaluate','watch')),
             48,y-46,width=516,height=44,size=7.5,color=MUTED)

    c.showPage();c.save()
    return output.getvalue()
