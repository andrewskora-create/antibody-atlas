/* Conclusions use existing computed evidence; language-model scores never change PTM risk. */
window.outcomeEscape=t=>String(t??'').replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
// Display triage only: preserve all computed chemistry and exposure calls.
window.ptmFocusSite=function(site){
 const flags=site.flags?.length?site.flags:[site];
 const oxidationOnly=flags.every(f=>/oxidation/i.test(f.kind||f.chemistry||''));
 const assessed=Number.isFinite(site.relative_sidechain_sasa)&&site.relative_sidechain_sasa>=0;
 return assessed&&(site.priority==='high'||(site.priority==='review'&&!(site.region==='Framework'&&oxidationOnly)));
};
window.ptmFocusedSites=structure=>(structure?.ptm_sites||[]).filter(window.ptmFocusSite);
window.ptmMatchFocused=function(structure,chain,position,kind){
 return window.ptmFocusedSites(structure).some(s=>s.chain===chain&&Number(s.sequence_position??s.start)===Number(position)&&
   (s.flags?.length?s.flags:[s]).some(f=>f.kind===kind));
};
window.ptmFocusNote='Default: high and medium (review) screening priorities. Framework oxidation is deferred unless high priority; this is an engineering display preference, not evidence of safety. Lower and unknown sites remain available under Show all. Risk categories are not modification probabilities.';
window.ptmFocusCounts=function(structure){
 const sites=window.ptmFocusedSites(structure);
 return {high:sites.filter(s=>s.priority==='high').length,medium:sites.filter(s=>s.priority==='review').length,
 unknown:(structure?.unassessed_ptm_sites||[]).length+(structure?.ptm_sites||[]).filter(s=>!Number.isFinite(s.relative_sidechain_sasa)||s.relative_sidechain_sasa<0).length,
 deferred:(structure?.ptm_sites||[]).filter(s=>Number.isFinite(s.relative_sidechain_sasa)&&s.relative_sidechain_sasa>=0&&!window.ptmFocusSite(s)).length};
};
window.ptmOutcome=function(chain,label,structure){
 const flags=chain.metrics?.ptm_sites||[];
 if(!structure)return flags.length?'PTM motifs are present. Run structure analysis to assess exposure and follow-up priority.':'No screened PTM motifs detected; structure-based assessment has not run.';
 const allSites=(structure.ptm_sites||[]).filter(s=>s.chain===label);
 const sites=allSites.filter(window.ptmFocusSite);
 const high=sites.filter(s=>s.priority==='high'),review=sites.filter(s=>s.priority==='review');
 const unknown=[...(structure.unassessed_ptm_sites||[]),...allSites.filter(s=>!Number.isFinite(s.relative_sidechain_sasa)||s.relative_sidechain_sasa<0)].filter(s=>s.chain===label);
 const deferred=allSites.filter(s=>Number.isFinite(s.relative_sidechain_sasa)&&s.relative_sidechain_sasa>=0&&!window.ptmFocusSite(s)).length;
 const targets=[...high,...review].slice(0,3).map(s=>`${s.number} (${s.region}${Number.isFinite(s.relative_sidechain_sasa)?`, ${Math.round(s.relative_sidechain_sasa*100)}% exposed`:''})`);
 let outcome=high.length?`${high.length} high-priority PTM site(s)${review.length?` and ${review.length} at medium priority`:""}; confirm modification and functional impact before engineering.`:review.length?`${review.length} PTM site(s) have medium screening priority; confirm the predicted liabilities experimentally.`:'No high or medium PTM sites in the focused view.';
 if(targets.length)outcome+=` Prioritize ${targets.join('; ')}${high.length+review.length>3?'; remaining sites in the structure details':''}.`;
 if(deferred)outcome+=` ${deferred} lower-priority or framework-oxidation site(s) deferred; Show all to inspect.`;
 if(unknown.length)outcome+=` ${unknown.length} additional site(s) lack reliable structural assessment.`;
 if(structure.mapping_warnings?.length)outcome+=' Structural mapping has limitations; see the structure details.';
 return outcome;
};
window.antiConclusion=function(anti){
 if(!anti)return 'Not run.';
 const b=(anti.summary?.benchmarks||[]).find(b=>b.key==='pseudo_log_likelihood');
 const relation={below_range:'Less typical',within_range:'Within the range',above_range:'More typical'}[b?.query_relation];
 return relation?`${relation} ${b.query_relation==='within_range'?'of':'than'} the ${b.reference_n} matching reference chains in the model comparison.`:'Sequence scored; no compatible reference comparison is available.';
};
window.propermabAvailability=function(result){
 const sequenceKeys=new Set(['n_charged_res_fv','vh_charge','vl_charge','fv_charge','fv_csp']);
 const structural=(result.features||[]).filter(f=>!sequenceKeys.has(f.key));
 const missing=(result.unavailable||[]).filter(x=>!/hydrophobic moment/i.test(x.feature||''));
 const hasPatches=Boolean(result.patch_pipeline);
 return {structuralCount:structural.length,hasPatches,missing,
   status:result.structure_selection==='terminal_repaired'?'Full-input repaired structure selected':result.structure_partial?'Partial structure — terminal residues omitted':!structural.length?'Sequence descriptors only — structure unavailable':!hasPatches?'Partial result — surface patches unavailable':'Structure and surface patches complete'};
};
window.propermabSurfaceReadout=function(result){
 const e=window.outcomeEscape;
 const descriptors=[
  ['hyd_patch_area','Hydrophobic patch area','patch'],
  ['hyd_patch_area_cdr','Hydrophobic patch area near CDR','patch'],
  ['hyd_asa','Apolar accessible area','structure'],
  ['pos_patch_area','Positive patch area','patch'],
  ['neg_patch_area','Negative patch area','patch']
 ];
 const features=new Map((result.features||[]).map(f=>[f.key,f]));
 const unavailable=result.unavailable||[];
 const rows=descriptors.map(([key,label,group])=>{
  const feature=features.get(key);
  const value=Number.isFinite(feature?.value)?`${feature.value.toFixed(1)} ${feature.unit||'Å²'}`:null;
  const anchor=feature?.benchmark?.reference_anchor;
  const relative=Number.isFinite(anchor?.relative_difference_percent)?` (${anchor.relative_difference_percent>=0?'+':''}${anchor.relative_difference_percent.toFixed(1)}%)`:'';
  const reason=unavailable.find(u=>u.feature===label)?.reason||
   unavailable.find(u=>u.feature===(group==='patch'?'APBS/NanoShaper surface patches':'Structure descriptors'))?.reason||
   'Descriptor was not returned by this run.';
  return `<div class="pm-surface-item"><span>${e(feature?.label||label)}</span><strong>${value===null?'Unavailable':e(value)}</strong>${value===null?`<small>${e(reason)}</small>`:anchor?`<small class="herceptin-reference">Herceptin: ${anchor.value.toFixed(1)} ${e(feature.unit||'Å²')} · input Δ ${anchor.delta>=0?'+':''}${anchor.delta.toFixed(1)}${relative}</small>`:''}</div>`;
 });
 return `<section class="pm-surface-readout" aria-label="PROPERMAB surface descriptors"><h5>Hydrophobic and charged surface</h5><div class="pm-surface-grid">${rows.join('')}</div><p class="footnote">Areas are computed descriptors in Å². Missing calculations are unavailable, not zero.</p></section>`;
};
window.propermabConclusion=function(result){
 const e=window.outcomeEscape;
 const availability=window.propermabAvailability(result);
 const hic=result.summary?.published_context?.herceptin;
 const herceptin=hic?`<p class="herceptin-reference"><b>Herceptin benchmark:</b> published study HIC ${hic.hic_rt_min.toFixed(1)} min. This is recombinant IgG1 in the cited assay; your HIC retention time is not predicted. Amber markers identify its locally computed hydrophobicity and charge descriptors. Proximity supports comparison, not a developability pass. <a href="${e(hic.source_url)}" target="_blank" rel="noopener noreferrer">Source</a></p>`:'';
 const diagnostics=availability.missing.length?`<div class="mapping-warning"><b>${e(availability.status)}</b>${availability.missing.map(x=>`<p><b>${e(x.feature)}:</b> ${e(x.reason)}</p>`).join('')}${!availability.structuralCount&&!result.input_pdb_sha256?'<p>Run “Predict and analyze structure” above, wait for it to finish, then rerun PROPERMAB. The automatic workflow runs PROPERMAB after folding.</p>':''}</div>`:'';
 const scored=(result.summary?.benchmarks||[]).filter(b=>b.query_relation&&b.query_relation!=='unavailable');
 const outside=scored.filter(b=>b.query_relation!=='within_range');
 const standout=outside.filter(b=>b.query_relation==='above_range'&&(/patch_area/.test(b.key)||b.key==='hyd_asa')).slice(0,2);
 let conclusion=!result.features?.length?'No PROPERMAB descriptors are available.':!scored.length?'Descriptors were calculated, but no compatible reference comparison is available.':outside.length?`${outside.length} of ${scored.length} compared descriptors differ from the therapeutic reference panel. Review these differences before deciding on engineering.`:'No descriptor falls outside the available reference spans. This is not a developability pass.';
 const next=standout.length?`Main follow-up: ${standout.map(b=>b.label.toLowerCase()).join(' and ')} exceed the reference span. Check self-interaction, aggregation and formulation behavior experimentally.`:outside.length?'No single direction establishes risk. Review the outlying descriptors with the PTM and surface findings.':'PK, maximum concentration and CMC readiness still require experimental data.';
 const missing=result.unavailable?.some(x=>/structure|patch/i.test(x.feature))?' Structural descriptors are incomplete.':'';
 return `<section class="profiler-summary"><h4>PROPERMAB conclusion</h4>${diagnostics}${window.propermabSurfaceReadout(result)}${herceptin}${window.propermabRepairReadout(result)}${result.structure_note?`<p class="mapping-warning">${e(result.structure_note)}</p>`:""}<p>${e(conclusion+missing)}</p><p>${e(next)}</p><p class="footnote">Reference range = observed minimum to maximum among matched therapeutic antibodies. Within range is not a developability pass; this panel is not selected for favorable experimental behavior.</p></section>`;
};

// Sequence cards must show the same assessment as their matching structural site.
window.sequencePtmRisk=function(flag,chain,structure){
 if(!structure)return {priority:'unassessed',label:'Unknown · structure not run',note:'Run structure analysis to assess exposure. A sequence motif alone has no assigned structural priority.'};
 const site=(structure.ptm_sites||[]).find(s=>s.chain===chain&&Number(s.sequence_position??s.start)===Number(flag.start)&&
   (s.flags?.length?s.flags:[s]).some(f=>f.kind===flag.kind));
 if(!site||!Number.isFinite(site.relative_sidechain_sasa)||site.relative_sidechain_sasa<0)
   return {priority:'unassessed',label:'Unknown · exposure unavailable',note:'This site is not classified as low risk.'};
 const labels={high:'High priority',review:'Medium · review',lower:'Lower structural concern'};
 const grouped=site.flags?.length>1;
 return {priority:site.priority,label:labels[site.priority]||'Unknown',
   note:`${(site.relative_sidechain_sasa*100).toFixed(1)}% side-chain exposed. ${grouped?'Priority is for the combined flags at this residue. ':''}${site.evidence_note||'Heuristic screening priority; not a measured modification rate.'}`};
};
window.renderSequencePtmRisk=function(flag,chain,structure){
 const risk=window.sequencePtmRisk(flag,chain,structure),escape=window.outcomeEscape;
 return `<small class="sequence-ptm-risk"><span class="priority ${escape(risk.priority)}">${escape(risk.label)}</span><span>${escape(risk.note)}</span></small>`;
};

window.propermabRepairReadout=function(result){
 const r=result.terminal_repair;if(!r)return '';
 const e=window.outcomeEscape,format=v=>Number.isFinite(v)?v.toFixed(1):'Unavailable';
 return `<section class="propermab-repair-readout"><h4>Terminal repair · sensitivity comparison</h4><p>${e(r.status==='complete'?'Local reconstruction and comparison complete.':'Repair comparison unavailable.')}</p><p>${e(r.note||'')}</p>${r.status==='complete'?`${r.pdb_url?`<a href="${e(r.pdb_url)}" class="download-link">Download repaired PDB hypothesis ↓</a>`:''}<p class="footnote">Rebuilt input positions: ${e(Object.entries(r.omitted_positions||{}).filter(([,v])=>v.length).map(([c,v])=>`${c}: ${v.join(', ')}`).join('; '))}</p><details><summary>Original versus repaired measurements</summary><table class="audit-table"><thead><tr><th>Descriptor</th><th>Original</th><th>Repaired</th><th>Change</th></tr></thead><tbody>${(r.metrics||[]).map(m=>`<tr><td>${e(m.label)} ${e(m.unit)}</td><td>${format(m.original)}</td><td>${format(m.repaired)}</td><td>${format(m.delta)}</td></tr>`).join('')}</tbody></table><p class="footnote">${e(r.method||'')} Missing calculations remain unavailable.</p></details>`:`<p class="mapping-warning">${e(r.error||'Local repair did not finish. Original-model results remain available.')}</p>`}</section>`;
};
