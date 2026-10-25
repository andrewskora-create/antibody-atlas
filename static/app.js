const $ = (id) => document.getElementById(id);
const exampleH = 'EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS';
const exampleL = 'VIWMTQSPSSLSASVGDRVTITCQASQDIRFYLNWYQQKPGKAPKLLISDASNMETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLPFTFGPGTKVDFK';
let report = null, selected = 'H', activeJob = null, activeAntiBertyJob = null, activePropermabJob = null, structureView = null;
let nglStage = null, nglComponent = null, nglMarkerComponent = null, detailMode = 'surface', ptmFilter = 'priority';
let herceptinTolerance=10;
let profilerStates={anti:'not_run',pm:'not_run'};
let analysisGeneration=0, foldGeneration=0, assessmentGeneration=0, foldState='idle';
const HEAT_STOPS = [[0,'#28669a'],[0.22,'#6bb4d1'],[0.42,'#e7d779'],[0.62,'#ee9856'],[0.80,'#d64354'],[1,'#b91f43']];
const PTM_COLORS = {high:'#f16778',review:'#f2bc64',lower:'#72c49b'};

function heatColor(value){
  const v=Math.max(0,Math.min(1,Number(value)||0));
  for(let i=1;i<HEAT_STOPS.length;i++){
    if(v>HEAT_STOPS[i][0])continue;
    const [lo,a]=HEAT_STOPS[i-1],[hi,b]=HEAT_STOPS[i],t=(v-lo)/(hi-lo);
    const channels=[1,3,5].map(j=>Math.round(parseInt(a.slice(j,j+2),16)*(1-t)+parseInt(b.slice(j,j+2),16)*t));
    return `#${channels.map(x=>x.toString(16).padStart(2,'0')).join('')}`;
  }
  return HEAT_STOPS.at(-1)[1];
}
function colorNumber(hex){return parseInt(hex.slice(1),16);}
function colorRgb(hex){return [1,3,5].map(j=>parseInt(hex.slice(j,j+2),16)/255);}
function siteKey(site){return `${site.chain}:${site.start}`;}
function markerNumbers(sites){return new Map([...new Set(sites.map(siteKey))].map((key,i)=>[key,i+1]));}

function safe(text) { return String(text).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;'); }
function metric(name, value, unit='') { return `<div class="metric"><span>${safe(name)}</span><strong>${safe(value)} <em>${safe(unit)}</em></strong></div>`; }
function showError(message) { $('error').textContent = message; $('error').hidden = false; }
function clearError() { $('error').hidden = true; $('error').textContent = ''; }

function renderChain(label) {
  selected = label;
  const chain = report.chains[label], m = chain.metrics, n = chain.numbering, h = chain.humanness;
  $('chain-tabs').innerHTML = Object.keys(report.chains).map(key => `<button type="button" role="tab" aria-selected="${key===label}" class="${key===label?'active':''}" data-chain="${key}">${key==='H'?'Heavy':'Light'} · ${report.chains[key].numbering.chain_type}</button>`).join('');
  $('chain-tabs').querySelectorAll('button').forEach(button => button.onclick = () => renderChain(button.dataset.chain));
  $('metric-grid').innerHTML = metric('Percent humanness · FR identity', `${h.percent}%`) + metric(`Net charge at pH ${m.ph}`, m.net_charge_at_ph) + metric('Isoelectric point', m.isoelectric_point) + metric('Mean hydropathy', m.gravy) + metric('Molecular weight', (m.molecular_weight_da/1000).toFixed(1), 'kDa') + metric('Hydrophobic windows', m.aggregation_proxy.patch_count) + metric('PTM flags', m.ptm_sites.length) + metric('Domain length', m.length, 'aa');
  $('numbering-label').textContent = `${n.model} · ${n.scheme.toUpperCase()} CDR boundaries`;
  const map = new Map(n.residues.map(r => [r.sequence_position, r]));
  $('sequence-map').innerHTML = [...chain.sequence].map((aa,i) => {const r=map.get(i+1); const cls=r ? r.region.toLowerCase() : 'unassigned'; return `<span class="residue ${cls}" title="${safe(r ? `${r.number} · ${r.region} · sequence ${i+1}` : `Unassigned · sequence ${i+1}`)}">${aa}</span>`}).join('');
  $('cdr-grid').innerHTML = n.cdrs.map(c => `<div class="cdr"><span>${c.name} · ${n.scheme.toUpperCase()}</span><strong>${safe(c.sequence || 'No residues assigned')}</strong><small>${c.number_start ? `${safe(c.number_start)}–${safe(c.number_end)} · ` : ''}${c.sequence.length} aa</small></div>`).join('');
  renderHydropathy(m.hydrophobic_windows);
  $('aggregation-note').textContent = `${m.aggregation_proxy.patch_count} flagged windows. ${m.aggregation_proxy.definition}.`;
  renderPtmFocusControls();
  const regionAt=new Map(n.residues.map(r=>[r.sequence_position,r.region]));
  for(const [elementId,group] of [['cdr-ptm-list','CDR'],['framework-ptm-list','Framework'],['unassigned-ptm-list','Unassigned']]){
    const sites=m.ptm_sites.filter(site=>ptmFilter==='all'||window.ptmMatchFocused(report.structure_analysis,label,site.start,site.kind)).filter(site=>group==='CDR'?String(regionAt.get(site.start)||'Unassigned').startsWith('CDR'):(regionAt.get(site.start)||'Unassigned')===group);
    $(elementId).innerHTML=sites.length?sites.map(site=>`<div class="ptm-row"><span><b>${safe(site.kind)}</b><br>${safe(regionAt.get(site.start)||'Unassigned')} · ${safe(site.interpretation)}${window.renderSequencePtmRisk(site,label,report.structure_analysis)}${site.precedent_ids?.length?`<small class="precedent-note">${safe(precedentNames(site.precedent_ids))} · aligned site + motif</small>`:""}</span><code>${safe(site.motif)} · ${site.start}</code></div>`).join(''):'<p class="footnote">No selected flags in this region.</p>';
  }
  const auditRows=(h.candidate_alleles||[]).map(c=>`<tr><td>${safe(c.name)}</td><td>${c.percent}%</td><td>${c.matches}/${c.framework_positions}</td><td>${c.coverage_percent}%</td></tr>`).join('');
  const selectedAllele=(h.candidate_alleles||[]).find(c=>c.name===h.nearest_human_v_gene);
  $('human-detail').innerHTML = `<div class="human-main"><strong>${h.percent}%</strong><div>Closest allele: <b>${safe(h.nearest_human_v_gene)}</b><br>${h.matches}/${h.framework_positions} exact IMGT FR1–FR3 matches · ${h.coverage_percent}% position coverage</div></div><div class="bar"><span style="width:${h.percent}%"></span></div><p class="footnote">Input: ${h.input_length} aa · local SHA-256 prefix <code>${safe(h.input_sha256_prefix)}</code>. Use this fingerprint to confirm two reports analyzed the same sequence.</p><table class="audit-table"><thead><tr><th>Human V allele</th><th>Identity</th><th>Matches</th><th>Coverage</th></tr></thead><tbody>${auditRows}</tbody></table><details><summary>Show differences from the closest allele</summary><p class="footnote">${safe((selectedAllele?.mismatches||[]).map(x=>`${x.number}: ${x.query_aa} vs ${x.reference_aa||'gap'}`).join(' · ')||'No differences')}</p></details><div class="footnote">${safe(h.definition)} Reference: OGRDB ${safe(h.locus)} release ${h.reference_release}. ${safe(h.caveat)}</div>`;
  const anti=chain.antiberty;
  window.renderAntiBertyPanel(report.chains,label,$('antiberty-result'),renderChain,report.structure_analysis);
}

function renderHydropathy(windows) {
  if (!windows.length) { $('hydro-chart').textContent = 'Sequence too short for window analysis.'; return; }
  const width=500, height=140, pad=14, vals=windows.map(w=>w.mean_hydropathy), min=Math.min(-3,...vals), max=Math.max(3,...vals);
  const x=i=>pad+i*(width-2*pad)/Math.max(1, vals.length-1), y=v=>height-pad-(v-min)*(height-2*pad)/(max-min);
  const points=vals.map((v,i)=>`${x(i).toFixed(1)},${y(v).toFixed(1)}`).join(' ');
  const threshold=y(1.6), zero=y(0);
  $('hydro-chart').innerHTML=`<svg viewBox="0 0 ${width} ${height}" preserveAspectRatio="none" role="img" aria-label="Sliding-window Kyte-Doolittle hydropathy"><line x1="0" y1="${zero}" x2="${width}" y2="${zero}" stroke="#d4e3e7"/><line x1="0" y1="${threshold}" x2="${width}" y2="${threshold}" stroke="#e2aa58" stroke-dasharray="4 4"/><polyline points="${points}" fill="none" stroke="#2a7187" stroke-width="2.4" vector-effect="non-scaling-stroke"/></svg>`;
}

async function analyze() {
  const generation=++analysisGeneration; $('fold').disabled=true;$('propermab').disabled=true;$('antiberty').disabled=true; activeJob=null; activeAntiBertyJob=null; activePropermabJob=null; ++foldGeneration; ++assessmentGeneration;
  report=null;foldState='idle';profilerStates={anti:'not_run',pm:'not_run'};window.setVariantBaseline?.(null);$('results').hidden=true;$('empty').hidden=false;
  clearError(); $('analyze').disabled=true; $('analyze').innerHTML='Analyzing…';
  try {
    const response=await fetch('/api/analyze',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({heavy:$('heavy').value,light:$('light').value,scheme:$('scheme').value,ph:Number($('ph').value)})});
    const data=await response.json(); if(generation!==analysisGeneration)return; if(!response.ok) throw Error(data.detail || 'Analysis failed.');
    report=data; window.setVariantBaseline?.(data); ptmFilter='priority'; activePropermabJob=null; $('propermab').disabled=false; $('propermab-status').textContent=''; $('propermab-result').innerHTML=''; window.renderExternalLookup?.(report,$('external-lookup')); renderPrecedents(); renderMethodSources(); activeJob=null; activeAntiBertyJob=null; structureView=null; $('fold').disabled=false; $('antiberty').disabled=false; $('empty').hidden=true; $('results').hidden=false; $('pdb-link').hidden=true; $('structure-canvas').hidden=true; $('structure-3d').hidden=true; $('detail-controls').hidden=true; $('model-legend').hidden=true; $('model-selection').hidden=true; $('viewer-note').hidden=true; $('surface-controls').hidden=true; $('structure-report').hidden=true;$('structure-hypothesis').hidden=true; $('fold-status').textContent='';$('antiberty-status').textContent='';
    if(nglStage){nglStage.removeAllComponents();nglComponent=null;nglMarkerComponent=null;}
    const cdrLength=Object.values(data.chains).reduce((total,c)=>total+c.numbering.cdrs.reduce((sum,loop)=>sum+loop.sequence.length,0),0);
    $('tap-context').textContent=`Local descriptors: total ${data.scheme.toUpperCase()} CDR length ${cdrLength} aa, plus this app's hydrophobic and charge cluster counts after folding. These counts are not TAP PSH/PPC/PNC/SFvCSP values or TAP threshold calls.`;
    $('copy-tap-heavy').disabled=!data.chains.H;
    $('copy-tap-light').disabled=!data.chains.L;
    $('tap-copy-status').textContent=data.chains.H&&data.chains.L?'Opening TAP alone sends no sequence. Submission on TAP is your choice.':'Official TAP requires paired heavy and light variable domains.';
    renderChain(data.chains.H?'H':'L');
    foldState=data.chains.H?'pending':'unsupported';refreshAssessment();
    if(data.chains.H)void fold();else {$('fold-status').textContent='Heavy chain required; light-chain-only structure is unavailable.';startProfilers();}
  } catch(err) {showError(err.message);} finally {$('analyze').disabled=false;$('analyze').innerHTML='Analyze antibody <span>→</span>';}
}

function downloadJson() {
  const blob=new Blob([JSON.stringify(report,null,2)],{type:'application/json'});
  const url=URL.createObjectURL(blob), link=document.createElement('a');link.href=url;link.download='antibody-analysis.json';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}

async function downloadPdf(){
  clearError(); $('export-pdf').disabled=true;
  try{
    const exportId=report.analysis_id, request=assessmentRequest();
    const structure_image=await captureReportStructure();
    if(report.analysis_id!==exportId||JSON.stringify(assessmentRequest())!==JSON.stringify(request))throw Error('Analysis changed during export; export again.');
    const response=await fetch('/api/report',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({title:'Antibody developability report',...request,structure_image})});
    if(!response.ok){const error=await response.json();throw Error(error.detail||'PDF export failed.');}
    const blob=await response.blob(),url=URL.createObjectURL(blob),link=document.createElement('a');
    link.href=url;link.download='antibody-developability-report.pdf';link.click();setTimeout(()=>URL.revokeObjectURL(url),1500);
  }catch(err){showError(err.message);}finally{$('export-pdf').disabled=false;}
}

// Render a dedicated local surface view: the user's orientation and mode stay untouched.
async function captureReportStructure(){
  if(!report?.structure_analysis)return null;
  if(!window.NGL)throw Error('Structure image renderer is unavailable. Refresh and export again.');
  const data=report.structure_analysis, hash=report.provenance.structure_pdb_sha256;
  const url=$('pdb-link').getAttribute('href');
  if(!url||!hash)throw Error('Selected structure is not ready for export.');
  const host=document.createElement('div');
  Object.assign(host.style,{position:'fixed',left:'-2000px',top:'0',width:'900px',height:'680px'});
  document.body.appendChild(host);
  let stage,scheme;
  try{
    const response=await fetch(url);if(!response.ok)throw Error('Selected PDB unavailable.');
    const bytes=await response.arrayBuffer();
    const digest=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))).map(b=>b.toString(16).padStart(2,'0')).join('');
    if(digest!==hash)throw Error('Selected PDB changed; refresh and export again.');
    stage=new NGL.Stage(host,{backgroundColor:'white',sampleLevel:1});
    const component=await stage.loadFile(new Blob([bytes]),{ext:'pdb'});
    const heat=new Map(data.residues.map(r=>[`${r.chain}:${r.pdb_number}`,r.surface_heat?.hydrophobic||0]));
    scheme=NGL.ColormakerRegistry.addScheme(function(){this.atomColor=atom=>colorNumber(heatColor(heat.get(`${atom.chainname}:${atom.resno}${String(atom.inscode||'').trim()}`)||0));},'Report hydrophobic surface');
    component.addRepresentation('surface',{sele:'protein',colorScheme:scheme,opacity:1,surfaceType:'av',useWorker:false});
    await new Promise(resolve=>stage.tasks.onZeroOnce(resolve));
    component.autoView(0);
    const blob=await stage.makeImage({factor:1,antialias:true,trim:false,transparent:false});
    const png=await new Promise((resolve,reject)=>{const reader=new FileReader();reader.onload=()=>resolve(reader.result);reader.onerror=()=>reject(Error('Could not encode structure image.'));reader.readAsDataURL(blob);});
    return {png,structure_sha256:hash};
  }finally{stage?.dispose();if(scheme)NGL.ColormakerRegistry.removeScheme(scheme);host.remove();}
}

async function runAntiBerty(){
  clearError();if(!report)return;
  const analysisId=report.analysis_id, generation=analysisGeneration;
  activeAntiBertyJob=null;for(const c of Object.values(report.chains))delete c.antiberty;
  profilerStates.anti='running';refreshAssessment();
  $('antiberty').disabled=true;$('antiberty-status').textContent='Running model on this computer…';
  try{
    const response=await fetch('/api/antiberty',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({heavy:report.chains.H?.sequence||'',light:report.chains.L?.sequence||'',analysis_id:report.analysis_id})});
    const job=await response.json();if(report?.analysis_id!==analysisId||generation!==analysisGeneration)return;if(!response.ok)throw Error(job.detail||'AntiBERTy could not start.');
    activeAntiBertyJob=job.id;await pollAntiBerty(job.id);
  }catch(err){if(report?.analysis_id!==analysisId||generation!==analysisGeneration)return;showError(err.message);profilerStates.anti='failed';$('antiberty-status').textContent='Failed · retry available';$('antiberty').disabled=false;refreshAssessment();}
}

async function pollAntiBerty(id){
  if(id!==activeAntiBertyJob)return;
  try{
    const response=await fetch(`/api/antiberty/${id}`),job=await response.json();if(id!==activeAntiBertyJob)return;if(!response.ok)throw Error(job.detail||'AntiBERTy status unavailable.');
    if(job.status==='complete'){
      for(const [label,result] of Object.entries(job.results)){if(report.chains[label])report.chains[label].antiberty=result;}
      report.provenance.antiberty_model_version=job.model_version;
      profilerStates.anti='complete';$('antiberty-status').textContent='Complete · local AntiBERTy';$('antiberty').disabled=false;renderChain(selected);refreshAssessment();
    }else if(job.status==='failed'){throw Error(`AntiBERTy failed: ${job.error}`);}
    else{$('antiberty-status').textContent='Scoring local sequences…';setTimeout(()=>pollAntiBerty(id),3000);}
  }catch(err){if(id!==activeAntiBertyJob)return;showError(err.message);profilerStates.anti='failed';$('antiberty-status').textContent='Failed · retry available';$('antiberty').disabled=false;refreshAssessment();}
}

async function fold() {
  clearError(); if(!report?.chains.H) {showError('A heavy variable domain is required for structure prediction.');return;}
  const analysisId=report.analysis_id, generation=++foldGeneration;
  activeJob=null;activePropermabJob=null;profilerStates.pm='not_run';delete report.propermab;delete report.structure_analysis;
  for(const key of ['structure_model','structure_model_version','structure_pdb_sha256','original_structure_pdb_sha256','structure_selection','structure_preparation'])delete report.provenance[key];
  delete report.original_structure_analysis;structureView=null;if(nglStage)nglStage.removeAllComponents();nglComponent=null;nglMarkerComponent=null;
  for(const id of ['structure-hypothesis','pdb-link','structure-report','structure-3d','structure-canvas','surface-controls','detail-controls','model-legend','model-selection','viewer-note'])$(id).hidden=true;
  $('propermab-result').innerHTML='';$('propermab-status').textContent='Waiting for local structure analysis…';$('propermab').disabled=true;
  foldState='pending';renderChain(selected);renderPrecedents();refreshAssessment();
  $('fold').disabled=true;$('fold-status').textContent='Submitting local prediction…';
  try {
    const response=await fetch('/api/fold',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({heavy:report.chains.H.sequence,light:report.chains.L?.sequence||'',ph:report.ph,scheme:report.scheme,analysis_id:analysisId})});
    const job=await response.json();if(report?.analysis_id!==analysisId||generation!==foldGeneration)return;
    if(!response.ok)throw Error(job.detail||'Could not start folding.');
    activeJob=job.id; await pollFold(job.id);
  }catch(err){if(report?.analysis_id!==analysisId||generation!==foldGeneration)return;foldState='failed';showError(err.message);$('fold-status').textContent='Prediction failed · retry available';$('fold').disabled=false;$('propermab').disabled=false;$('propermab-status').textContent='Structure unavailable; sequence descriptors only.';startProfilers();}
}

async function pollFold(id) {
  if(id!==activeJob)return;
  try {
    const response=await fetch(`/api/fold/${id}`), job=await response.json();if(id!==activeJob)return;if(!response.ok)throw Error(job.detail||'Status unavailable.');
    if(job.status==='complete'){
      foldState=job.structure_analysis?'complete':'partial';
      $('fold-status').textContent=`${job.structure_analysis?'Complete':'PDB available; exposure analysis incomplete'} · ${job.method}`;$('fold').disabled=false;$('propermab').disabled=false;$('propermab-status').textContent=job.structure_analysis?'Ready for structure descriptors.':'Exposure unavailable; inspect prediction error.';
      Object.assign(report.provenance,{structure_model:job.method,structure_model_version:job.model_version,structure_pdb_sha256:job.model_sha256});
      $('pdb-link').href=job.pdb_url;$('pdb-link').hidden=false;$('structure-hypothesis').hidden=false;
      if(job.structure_analysis){report.structure_analysis=job.structure_analysis;renderStructuralReport(job.structure_analysis);renderChain(selected);renderPrecedents();drawStructure(job.structure_analysis);}
      else if(job.analysis_error){if(job.structure_trace)drawStructure(job.structure_trace);showError(job.analysis_error);}
      startProfilers();
      if(job.structure_analysis||job.structure_trace)await loadAtomicStructure(job.pdb_url,job.structure_analysis||job.structure_trace);
    }else if(job.status==='failed'){
      throw Error(`Structure prediction failed: ${job.error}`);
    }else{
      $('fold-status').textContent=job.status==='queued'?'Queued…':'Predicting and screening on this computer…';
      setTimeout(()=>pollFold(id),3000);
    }
  }catch(err){if(id!==activeJob)return;foldState='failed';showError(err.message);$('fold-status').textContent='Prediction failed · retry available';$('fold').disabled=false;$('propermab').disabled=false;$('propermab-status').textContent='Structure unavailable; sequence descriptors only.';startProfilers();}
}

function renderStructuralReport(data) {
  $('surface-controls').hidden=false; $('structure-report').hidden=false;
  const counts=data.patch_counts, risk=data.overall_ptm_risk;
  const focused=window.ptmFocusCounts(data),focusCall=focused.high?'High-priority sites':focused.medium?'Medium-priority sites':focused.unknown?'Assessment incomplete':'No sites in focused view';
  const markerIndex=markerNumbers(data.ptm_sites);
  $('structure-summary').innerHTML=`<div><strong>${counts.hydrophobic}</strong><span>Hydrophobic patches</span></div><div><strong>${counts.positive}</strong><span>Positive patches</span></div><div><strong>${counts.negative}</strong><span>Negative patches</span></div><div class="risk-summary"><strong>${safe(focusCall)}</strong><span>Focused PTM follow-up screen</span></div>`;
  $('patch-list').innerHTML=data.patches.length?data.patches.map(p=>`<button type="button" class="patch-card ${safe(p.kind)}" data-patch="${safe(p.id)}"><span class="patch-top"><b>${safe(p.kind)} · ${p.member_count} residues</b><strong>${p.exposed_sidechain_sasa_sum_a2} Å²</strong></span><small>${safe(p.residues.map(r=>`${r.aa}${r.number}`).join(', '))}</small>${p.near_cdr?'<em>Near CDR</em>':''}</button>`).join(''):'<p class="footnote">No exposed clusters met the screening thresholds.</p>';
  $('patch-list').querySelectorAll('button').forEach(button=>button.onclick=()=>{const p=data.patches.find(x=>x.id===button.dataset.patch);setStructureMode(p.kind,p.id);});
  const coverage=Object.entries(data.mapping_coverage_percent||{}).map(([chain,value])=>`${chain} ${value}%`).join(' · ');
  $('ptm-risk').innerHTML=`<div class="risk-call ${focused.high?'high':focused.medium?'review':focused.unknown?'incomplete':'lower'}"><span>${safe(focusCall)}</span><small>${focused.high} high · ${focused.medium} medium · ${focused.deferred} deferred · ${focused.unknown} unknown</small></div>${Object.entries(report.chains).map(([label,chain])=>`<p><b>${label==='H'?'Heavy':'Light'}:</b> ${safe(window.ptmOutcome(chain,label,data))}</p>`).join('')}<p class="footnote"><b>Evidence:</b> ${safe(risk.evidence_status)} · ${safe(data.scheme.toUpperCase())} mapping coverage: ${safe(coverage)}</p>${data.mapping_warnings?.length?`<p class="mapping-warning">${safe(data.mapping_warnings.join(' '))}</p>`:''}<details class="risk-method"><summary>How PTM risk is defined</summary><p>${safe(window.ptmFocusNote)}</p><p>Full unfiltered screen: ${safe(risk.call)}.</p><p>${safe(risk.definition)}</p><p>${safe(risk.calibration_note||"Screening priorities are not validated modification probabilities or rates.")}</p><p>Red: high screening priority. Amber: medium (review). Green: lower structural concern, not absence of risk. Gray: structure unassessed. CDR location describes potential functional impact separately.</p><p>Exposure is side-chain solvent-accessible area divided by its FreeSASA reference area; values can exceed 100% in unusual conformations. It is not the fraction of the antibody surface or the percentage of molecules modified.</p></details>`;
  for(const [elementId,group] of [['structural-ptm-cdr','CDR'],['structural-ptm-framework','Framework'],['structural-ptm-unassigned','Unassigned']]){
    const belongs=s=>group==='CDR'?s.region.startsWith('CDR'):s.region===group;
    const assessed=data.ptm_sites.filter(s=>ptmFilter==='all'||window.ptmFocusSite(s)).filter(belongs).map(s=>`<button type="button" class="struct-site ${safe(s.priority)}" data-site="${siteKey(s)}"><span class="site-index">${markerIndex.get(siteKey(s))}</span><span class="site-body"><b>${safe(s.kind)}</b><small>${safe(s.chain)} · ${safe(s.number)} · sequence ${s.sequence_position} · ${safe(s.scheme.toUpperCase())} ${safe(s.region)} · ${safe(s.exposure_band?`${s.exposure_band} exposure`:"Exposure assessed")}</small><span class="exposure-readout"><span class="exposure-track" aria-hidden="true"><i style="width:${Math.min(100,Math.max(0,s.exposure_percent??s.relative_sidechain_sasa*100))}%"></i></span><strong>${s.exposure_percent==null&&s.relative_sidechain_sasa==null?"Unknown":`${(s.exposure_percent??s.relative_sidechain_sasa*100).toFixed(1)}%`} side-chain exposed</strong></span><small class="site-reasons">${safe(s.evidence_note||s.modification_likelihood)} · Impact: ${safe(s.functional_impact)} · ${safe(s.evidence_status)}</small>${s.precedent_ids?.length?`<small class="precedent-note">${safe(precedentNames(s.precedent_ids))} · aligned site + motif</small>`:""}</span><em class="priority ${safe(s.priority)}">${safe(s.priority==='review'&&Number.isFinite(s.relative_sidechain_sasa)?'Medium · review':s.risk_label||s.priority)}</em></button>`).join('');
    const unassessed=(ptmFilter==='all'?(data.unassessed_ptm_sites||[]):[]).filter(belongs).map(s=>`<div class="struct-site unassessed"><span><b>${safe(s.kind)}</b><small>${safe(s.chain)} · ${safe(s.number)} · sequence ${s.sequence_position} · ${safe(s.scheme.toUpperCase())} ${safe(s.region)} · no reliable structural mapping</small>${s.precedent_ids?.length?`<small class="precedent-note">${safe(precedentNames(s.precedent_ids))} · aligned site + motif; structure unassessed</small>`:""}</span><em class="priority unassessed">unassessed</em></div>`).join('');
    $(elementId).innerHTML=assessed+unassessed||'<p class="footnote">No selected PTM sites in this region.</p>';
  }
  document.querySelectorAll('.structural-ptm-list button').forEach(button=>button.onclick=()=>setStructureMode('ptm',button.dataset.site));
  $('structure-method').textContent=`${data.method}. ${data.surface_heat_definition||''} Patch area is summed exposed side-chain SASA, not a mesh area. Charge uses approximate residue pKa at pH ${data.ph}. Single predicted structure; PTM priorities are screening annotations, not measured modification rates.`;
  $('surface-controls').querySelectorAll('button').forEach(button=>button.onclick=()=>setStructureMode(button.dataset.mode,null));
}

function drawStructure(data) {
  if(data.residues.length<10)return;
  const canvas=$('structure-canvas');canvas.hidden=false;
  const initialMode=data.surface_heat_definition?'hydrophobic':data.ptm_sites.length?'ptm':'trace';
  structureView={data,mode:initialMode,focus:null,ax:.3,ay:.5,pointer:null,heatSchemes:{}};
  canvas.onpointerdown=e=>{structureView.pointer={x:e.clientX,y:e.clientY};canvas.setPointerCapture(e.pointerId)};
  canvas.onpointermove=e=>{if(!structureView?.pointer)return;structureView.ay+=(e.clientX-structureView.pointer.x)*.012;structureView.ax+=(e.clientY-structureView.pointer.y)*.012;structureView.pointer={x:e.clientX,y:e.clientY};paintStructure()};
  canvas.onpointerup=()=>{if(structureView)structureView.pointer=null;};
  setStructureMode(initialMode,null);
}

function setStructureMode(mode,focus) {
  if(!structureView)return;
  structureView.mode=mode;structureView.focus=focus;
  if(nglComponent){
    detailMode=mode==='ptm'?'cartoon':mode==='trace'?'sidechains':'surface';
    $('detail-controls').querySelectorAll('button').forEach(button=>button.classList.toggle('active',button.dataset.detail===detailMode));
  }
  $('surface-controls').querySelectorAll('button').forEach(button=>button.classList.toggle('active',button.dataset.mode===mode));
  $('patch-list').querySelectorAll('button').forEach(button=>button.classList.toggle('active',button.dataset.patch===focus));
  document.querySelectorAll('.structural-ptm-list button').forEach(button=>button.classList.toggle('active',button.dataset.site===focus));
  renderModelLegend();
  paintStructure();
  updateNglHighlights();
}

function renderModelLegend(){
  if(!structureView)return;
  const {data,mode,focus}=structureView;
  $('model-legend').hidden=false;
  const heat=mode==='hydrophobic'||mode==='positive'||mode==='negative';
  $('heat-legend').hidden=!heat;$('ptm-legend').hidden=mode!=='ptm';
  const names={hydrophobic:'Hydrophobic surface',positive:'Positive charge',negative:'Negative charge'};
  $('legend-title').textContent=heat?`${names[mode]} · relative display intensity`:mode==='ptm'?'PTM site markers':'Predicted structure';
  $('legend-explainer').textContent=heat?'Blue = little exposed patch signal. Red = the strongest exposed cluster in this visual scale. Compare within this model and mode only; color is not a measured aggregation or PK risk.':mode==='ptm'?'Numbered pins match the site cards. Red = high screening priority; amber = medium (review); green = lower structural concern. These are follow-up categories, not modification probabilities. Focused view defers framework oxidation. Unknown sites are counted below; Show all reveals their details.':'Rotate and zoom the predicted structure.';
  $('ptm-filter').hidden=mode!=='ptm';
  $('ptm-filter').querySelectorAll('button').forEach(button=>button.classList.toggle('active',button.dataset.ptmFilter===ptmFilter));
  let label='Select a patch or PTM site below to focus the model.';
  if(focus&&mode==='ptm'){
    const sites=data.ptm_sites.filter(s=>siteKey(s)===focus),site=sites[0];
    if(site)label=`Site ${markerNumbers(data.ptm_sites).get(focus)} · ${site.chain} ${site.number} · ${sites.map(s=>s.kind).join(' + ')} · ${site.risk_label||site.priority.toUpperCase()} · ${site.relative_sidechain_sasa==null?"unknown":(site.relative_sidechain_sasa*100).toFixed(1)+"%"} side-chain exposed`;
  }else if(focus){
    const patch=data.patches.find(p=>p.id===focus);
    if(patch)label=`${patch.kind} patch · ${patch.member_count} residues · ${patch.exposed_sidechain_sasa_sum_a2} Å² summed exposed side-chain area${patch.near_cdr?' · near CDR':''}`;
  }
  $('model-selection').textContent=label;
  $('model-selection').hidden=false;
}

function paintStructure() {
  if(!structureView)return;
  const {data,mode,focus,ax,ay}=structureView, atoms=data.residues, canvas=$('structure-canvas'), ctx=canvas.getContext('2d');
  const ratio=devicePixelRatio||1,w=canvas.clientWidth,h=canvas.clientHeight;canvas.width=w*ratio;canvas.height=h*ratio;ctx.scale(ratio,ratio);ctx.fillStyle='#0d2535';ctx.fillRect(0,0,w,h);
  const center=[0,1,2].map(i=>atoms.reduce((s,a)=>s+a.ca_xyz[i],0)/atoms.length);
  const scale=Math.max(...atoms.map(a=>Math.hypot(...a.ca_xyz.map((x,i)=>x-center[i]))));
  const project=a=>{const [x,y,z]=a.ca_xyz.map((v,i)=>v-center[i]),rx=x*Math.cos(ay)+z*Math.sin(ay),rz=-x*Math.sin(ay)+z*Math.cos(ay),ry=y*Math.cos(ax)-rz*Math.sin(ax);return {x:w/2+rx/scale*Math.min(w,h)*.39,y:h/2+ry/scale*Math.min(w,h)*.39,z:rz,atom:a}};
  const points=atoms.map(project);ctx.lineWidth=2.2;
  for(let i=1;i<points.length;i++){if(points[i].atom.chain!==points[i-1].atom.chain)continue;ctx.strokeStyle=points[i].atom.chain==='H'?'#497a89':'#9b8064';ctx.beginPath();ctx.moveTo(points[i-1].x,points[i-1].y);ctx.lineTo(points[i].x,points[i].y);ctx.stroke();}
  const key=a=>`${a.chain}:${a.sequence_position}`;let marked=new Map();
  if(mode==='ptm'){
    for(const site of data.ptm_sites){if(focus&&focus!==siteKey(site))continue;if(!focus&&ptmFilter==='priority'&&!window.ptmFocusSite(site))continue;marked.set(siteKey(site),PTM_COLORS[site.priority]);}
  }else if(mode!=='trace'){
    for(const atom of atoms)if(atom.surface_heat)marked.set(key(atom),heatColor(atom.surface_heat[mode]));
  }
  const numbers=mode==='ptm'?markerNumbers(data.ptm_sites):null;
  for(const p of points.filter(p=>marked.has(key(p.atom))).sort((a,b)=>a.z-b.z)){ctx.beginPath();ctx.arc(p.x,p.y,mode==='ptm'?10:4.5,0,Math.PI*2);ctx.fillStyle=marked.get(key(p.atom));ctx.fill();if(mode==='ptm'){ctx.strokeStyle='#f8fcfd';ctx.lineWidth=1.5;ctx.stroke();ctx.fillStyle='#fff';ctx.font='bold 10px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillText(String(numbers.get(key(p.atom))||''),p.x,p.y);}}
  ctx.fillStyle='#bad7df';ctx.font='12px Arial';ctx.fillText(`${mode==='trace'?'predicted backbone':mode==='ptm'?'PTM sites':mode+' patches'} · drag to rotate`,14,h-14);
}
function pdbSelection(residue){
  if(!residue?.pdb_number)return null;
  const match=String(residue.pdb_number).match(/^(-?\d+)([A-Za-z]?)$/);
  return match?`${match[1]}${match[2]?`^${match[2]}`:''}:${residue.chain}`:null;
}
function focusContextSelection(focused,data){
  if(!focused.length)return null;
  const center=[0,1,2].map(i=>focused.reduce((sum,r)=>sum+r.sidechain_xyz[i],0)/focused.length);
  const nearby=data.residues.filter(r=>Math.hypot(...r.ca_xyz.map((value,i)=>value-center[i]))<=15);
  const tokens=[...new Set(nearby.map(pdbSelection).filter(Boolean))];
  return tokens.length?`(${tokens.join(' or ')})`:null;
}

function heatScheme(mode,data){
  if(structureView.heatSchemes[mode])return structureView.heatSchemes[mode];
  const heatByPdb=new Map(data.residues.map(r=>[`${r.chain}:${r.pdb_number}`,r.surface_heat?.[mode]||0]));
  const scheme=NGL.ColormakerRegistry.addScheme(function(){
    this.atomColor=function(atom){
      const position=`${atom.resno}${String(atom.inscode||'').trim()}`;
      return colorNumber(heatColor(heatByPdb.get(`${atom.chainname}:${position}`)||0));
    };
  },`${mode} relative surface heat`);
  structureView.heatSchemes[mode]=scheme;
  return scheme;
}

function addPtmPins(data,focus,lookup){
  const grouped=new Map();
  for(const site of data.ptm_sites){
    const key=siteKey(site);
    if(focus&&key!==focus)continue;
    if(!focus&&ptmFilter==='priority'&&!window.ptmFocusSite(site))continue;
    const previous=grouped.get(key);
    if(!previous||({high:3,review:2,lower:1}[site.priority]>{high:3,review:2,lower:1}[previous.priority]))grouped.set(key,site);
  }
  if(!grouped.size)return;
  const center=[0,1,2].map(i=>data.residues.reduce((sum,r)=>sum+r.ca_xyz[i],0)/data.residues.length);
  const numbers=markerNumbers(data.ptm_sites),shape=new NGL.Shape('Numbered PTM sites'),labels=[];
  for(const [key,site] of grouped){
    const residue=lookup.get(key);if(!residue)continue;
    const anchor=residue.sidechain_xyz;
    const outward=anchor.map((value,i)=>value-center[i]);
    const length=Math.hypot(...outward)||1;
    const tip=anchor.map((value,i)=>value+outward[i]/length*(focus?3.5:2.5));
    const color=colorRgb(PTM_COLORS[site.priority]);
    shape.addCylinder(anchor,tip,color,0.17,`Site ${numbers.get(key)} stem`);
    shape.addSphere(tip,color,focus?1.45:1.05,`Site ${numbers.get(key)} ${site.number} ${site.kind}`);
    labels.push({tip,number:numbers.get(key),priority:site.priority});
  }
  nglMarkerComponent=nglStage.addComponentFromObject(shape);
  nglMarkerComponent.addRepresentation('buffer');
  for(const label of labels){
    const badge=document.createElement('span');
    badge.className=`model-pin-label ${label.priority}`;
    badge.textContent=String(label.number);
    nglMarkerComponent.addAnnotation(new NGL.Vector3(...label.tip),badge,{offsetY:-4});
  }
}

function updateNglHighlights(){
  if(!nglComponent||!structureView)return;
  if(nglMarkerComponent){nglStage.removeComponent(nglMarkerComponent);nglMarkerComponent=null;}
  nglComponent.removeAllRepresentations();
  const {data,mode,focus}=structureView,lookup=new Map(data.residues.map(r=>[`${r.chain}:${r.sequence_position}`,r]));
  const heat=mode==='hydrophobic'||mode==='positive'||mode==='negative';
  nglComponent.addRepresentation('cartoon',{colorScheme:mode==='trace'?'chainname':'uniform',color:mode==='trace'?undefined:'#a9bdc7',opacity:detailMode==='surface'?0.20:0.74});
  if(detailMode==='sidechains')nglComponent.addRepresentation('ball+stick',{sele:'sidechainAttached',colorScheme:heat?heatScheme(mode,data):'element',scale:0.34,opacity:heat?0.94:0.45});
  if(detailMode==='surface')nglComponent.addRepresentation('surface',{sele:'protein',colorScheme:heat?heatScheme(mode,data):'uniform',color:heat?undefined:'#aebfc8',opacity:heat?0.89:0.30,surfaceType:'av',useWorker:false});
  let focusSelection=null,zoomSelection=null;
  if(focus&&mode==='ptm'){
    const residue=lookup.get(focus),selection=pdbSelection(residue);
    if(selection){nglComponent.addRepresentation('ball+stick',{sele:selection,colorScheme:'uniform',color:'#ffffff',scale:0.65});focusSelection=selection;zoomSelection=focusContextSelection([residue],data);}
  }else if(focus&&heat){
    const residues=data.patches.filter(p=>p.id===focus).flatMap(p=>p.residues.map(r=>lookup.get(`${r.chain}:${r.sequence_position}`)));
    const tokens=[...new Set(residues.map(pdbSelection).filter(Boolean))];
    if(tokens.length){focusSelection=`(${tokens.join(' or ')})`;zoomSelection=focusContextSelection(residues,data);nglComponent.addRepresentation('ball+stick',{sele:focusSelection,colorScheme:heatScheme(mode,data),scale:0.58});}
  }
  if(mode==='ptm')addPtmPins(data,focus,lookup);
  if(zoomSelection||focusSelection)nglComponent.autoView(zoomSelection||focusSelection,450);
}

function setDetailMode(mode){
  detailMode=mode;
  $('detail-controls').querySelectorAll('button').forEach(button=>button.classList.toggle('active',button.dataset.detail===mode));
  updateNglHighlights();
}

async function loadAtomicStructure(url,data){
  const generation=foldGeneration;
  if(!window.NGL){$('viewer-note').textContent='Atomic viewer unavailable; the backbone view remains below.';$('viewer-note').hidden=false;return;}
  try{
    const response=await fetch(url);if(!response.ok)throw Error('PDB unavailable');
    const pdb=await response.text();if(generation!==foldGeneration)return;
    $('structure-3d').hidden=false;
    if(!nglStage)nglStage=new NGL.Stage('structure-3d',{backgroundColor:'#0d2535',sampleLevel:1});
    nglStage.removeAllComponents();
    nglMarkerComponent=null;
    const component=await nglStage.loadFile(new Blob([pdb],{type:'text/plain'}),{ext:'pdb'});
    if(generation!==foldGeneration){nglStage.removeComponent(component);return;}nglComponent=component;
    setDetailMode(detailMode);
    nglComponent.autoView();
    $('detail-controls').hidden=false;$('structure-canvas').hidden=true;
    $('viewer-note').textContent='Drag to rotate; scroll to zoom. Heat colors show relative exposed patch signal. PTM pins are numbered to match the site list; select a site to zoom in.';
    $('viewer-note').hidden=false;
  }catch(err){if(generation!==foldGeneration)return;nglComponent=null;$('structure-3d').hidden=true;$('structure-canvas').hidden=false;$('viewer-note').textContent=`Atomic viewer unavailable (${err.message}); backbone view remains available.`;$('viewer-note').hidden=false;}
}
window.addEventListener('resize',()=>{paintStructure();if(nglStage)nglStage.handleResize();});

$('example').onclick=()=>{$('heavy').value=exampleH;$('light').value=exampleL;};
$('copy-tap-heavy').onclick=()=>copyTapChain('H');
$('copy-tap-light').onclick=()=>copyTapChain('L');
$('analyze').onclick=analyze;$('export').onclick=downloadJson;$('export-pdf').onclick=downloadPdf;$('fold').onclick=fold;$('antiberty').onclick=runAntiBerty;$('propermab').onclick=runPropermab;
$('detail-controls').querySelectorAll('button').forEach(button=>button.onclick=()=>setDetailMode(button.dataset.detail));
$('ptm-filter').querySelectorAll('button').forEach(button=>button.onclick=()=>{setPtmFilter(button.dataset.ptmFilter);});
fetch('/api/health').then(r=>r.json()).then(h=>{$('health').textContent=`Local · numbering ${h.numbering_model?'ready':'unavailable'} · folding ${h.structure_model?'ready':'unavailable'}`}).catch(()=>{$('health').textContent='Local service unavailable'});

async function copyTapChain(label){
  const sequence=report?.chains[label]?.sequence;
  if(!sequence)return;
  try{
    await navigator.clipboard.writeText(sequence);
    $('tap-copy-status').textContent=`${label==='H'?'VH':'VL'} copied locally. Paste it into TAP only if you want to submit it to Oxford.`;
  }catch{
    $('tap-copy-status').textContent='Clipboard access was blocked. Copy the chain from the input field instead.';
  }
}

async function runPropermab(){
  if(!report)return;
  if(foldState==='pending')return;
  const analysisId=report.analysis_id, generation=foldGeneration;
  activePropermabJob=null;delete report.propermab;$('propermab-result').innerHTML='';
  profilerStates.pm='running';refreshAssessment();
  clearError();$('propermab').disabled=true;$('propermab-status').textContent='Running locally…';
  try{
    const response=await fetch('/api/propermab',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({analysis_id:report.analysis_id,fold_job_id:(foldState==='complete'||foldState==='partial')?activeJob:null})});
    const job=await response.json();if(report?.analysis_id!==analysisId||generation!==foldGeneration)return;if(!response.ok)throw Error(job.detail||'PROPERMAB could not start.');
    activePropermabJob=job.id;await pollPropermab(job.id);
  }catch(err){if(report?.analysis_id!==analysisId||generation!==foldGeneration)return;showError(err.message);profilerStates.pm='failed';$('propermab-status').textContent='Could not run';$('propermab').disabled=false;refreshAssessment();}
}
async function pollPropermab(id){
  if(id!==activePropermabJob)return;
  try{
    const response=await fetch(`/api/propermab/${id}`),job=await response.json();
    if(id!==activePropermabJob)return;
    if(!response.ok)throw Error(job.detail||'PROPERMAB status unavailable.');
    if(job.status==='complete'){
      report.propermab=job.results;const result=job.results;
      if(result.structure_selection==='terminal_repaired'&&result.selected_structure_analysis){
        const originalHash=report.provenance.original_structure_pdb_sha256||report.provenance.structure_pdb_sha256;
        if(result.original_input_pdb_sha256!==originalHash||result.input_pdb_sha256!==result.terminal_repair?.repaired_pdb_sha256)throw Error('Repaired model provenance does not match this analysis.');
        report.original_structure_analysis=report.original_structure_analysis||report.structure_analysis;
        report.structure_analysis=result.selected_structure_analysis;
        Object.assign(report.provenance,{original_structure_pdb_sha256:originalHash,structure_pdb_sha256:result.input_pdb_sha256,structure_selection:'terminal_repaired'});
        $('fold-status').textContent='Complete · repaired full-input model selected';
        $('pdb-link').href=result.terminal_repair.pdb_url;$('pdb-link').hidden=false;
        renderStructuralReport(report.structure_analysis);renderChain(selected);renderPrecedents();drawStructure(report.structure_analysis);
        void loadAtomicStructure(result.terminal_repair.pdb_url,report.structure_analysis);
      }
      profilerStates.pm='complete';refreshAssessment();
      $('propermab-status').textContent=`${window.propermabAvailability(result).status} · ${result.features.length} descriptors`; $('propermab').disabled=false;
      $('propermab-result').innerHTML=`${window.propermabConclusion(result)}<details><summary>Supporting descriptors, comparisons and methods</summary>${window.renderProfilerSummary(result.summary)}<p class="footnote">${safe(result.scope)} · ${safe(result.numbering_scheme)} · sequence-charge pH ${safe(result.ph)}</p>${result.patch_pipeline?`<p class="footnote"><b>${safe(result.patch_pipeline)}</b></p><p class="mapping-warning">${safe(result.patch_conditions)}</p>`:""}<div class="descriptor-scroll"><table class="audit-table"><thead><tr><th>Descriptor</th><th>Value</th></tr></thead><tbody>${result.features.map(f=>`<tr><td>${safe(f.label)}</td><td>${safe(typeof f.value==='number'?Number(f.value.toPrecision(6)):f.value)} ${safe(f.unit||'')}</td></tr>`).join('')}</tbody></table></div>${result.unavailable.length?`<details><summary>${result.unavailable.length} unavailable descriptor groups</summary>${result.unavailable.map(u=>`<p class="footnote"><b>${safe(u.feature)}</b>: ${safe(u.reason)}</p>`).join('')}</details>`:''}<p class="footnote">${safe(result.interpretation)} · Source ${safe(result.source_commit)}</p></details>`;
    }else if(job.status==='failed'){throw Error(job.error||'PROPERMAB failed.');}
    else{setTimeout(()=>pollPropermab(id),2000);}
  }catch(err){if(id!==activePropermabJob)return;showError(err.message);profilerStates.pm='failed';$('propermab-status').textContent='Failed';$('propermab').disabled=false;refreshAssessment();}
}
function precedentNames(ids){
  const cases=report?.ptm_precedents?.cases||[];
  return ids.map(id=>cases.find(p=>p.id===id)).filter(Boolean).map(p=>`${p.antibody} ${p.site}`).join('; ');
}
function renderPrecedents(){
  const clinical=report?.clinical_framework_precedents;
  const sites=(clinical?.sites||[]).filter(s=>ptmFilter==='all'||window.ptmMatchFocused(report.structure_analysis,s.query_chain,s.position,s.kind));
  const filteredClinical=clinical?{...clinical,sites,coverage:{...clinical.coverage,query_sites:sites.length,matched_sites:sites.filter(s=>s.matches.length).length,total_matches:sites.reduce((n,s)=>n+s.matches.length,0)}}:null;
  window.renderClinicalFrameworks(filteredClinical,$('clinical-frameworks'));
  const original=report?.ptm_precedents;
  const cases=(original?.cases||[]).map(p=>({...p,matches:(p.matches||[]).filter(m=>ptmFilter==='all'||window.ptmMatchFocused(report.structure_analysis,m.query_chain,m.query_sequence_position,p.chemistry))})).filter(p=>p.matches.length);
  const catalog=original?{...original,cases}:null;
  renderPtmFocusControls();
  if(!catalog){$('ptm-precedents').textContent='No sequence comparison attached. Analyze again.';return;}
  const context=(text,m)=>{const i=m.context_anchor_index,n=m.matched_motif?.length||1;return `${text.slice(0,i)}[${text.slice(i,i+n)}]${text.slice(i+n)}`;};
  const matches=p=>(p.matches||[]).map(m=>`<div class="precedent-match"><strong>${safe(m.match_type)}</strong><p>Your ${safe(m.query_number)} · sequence ${m.query_sequence_position} · ${safe(m.query_scheme.toUpperCase())} ${safe(m.query_region)} → ${safe(p.antibody)} ${safe(p.site)} · reference ${safe(m.reference_number)} ${safe(m.reference_region)} (${safe(m.reference_scheme.toUpperCase())}); canonical IMGT ${safe(m.reference_imgt_number)}</p><p>Numbered-chain identity: <b>${m.identity_percent}%</b> (${m.identity_matches}/${m.identity_denominator} positions). Local context: ${m.local_matches}/${m.local_positions} identical aligned positions.</p><pre class="precedent-alignment">Query  ${safe(context(m.query_context,m))}
Drug   ${safe(context(m.reference_context,m))}</pre><p class="footnote">Residue/motif and numbered position match. Differences in surrounding sequence, exposure, formulation and expression can change the outcome.</p></div>`).join('');
  $('ptm-precedents').innerHTML=`<p class="footnote">${safe(catalog.note)}</p><p><b>${catalog.cases.length} displayed documented PTM precedent${catalog.cases.length===1?"":"s"}</b> · ${catalog.searched_cases??0} curated cases searched locally · ${original.cases.length-catalog.cases.length} matched cases hidden by the current filter</p>${catalog.unavailable_cases?`<p class="mapping-warning">${catalog.unavailable_cases} catalog entries lack a usable verified reference and were not searched.</p>`:""}${!catalog.cases.length?'<p class="mapping-warning">No documented matches for the current PTM selection. Run structure analysis to establish priorities, or choose Show all PTMs and precedents. This does not establish that the sequence has no PTM liability or no precedent elsewhere.</p>':''}${catalog.cases.map(p=>`<article class="precedent-card"><h4>${safe(p.antibody)} · ${safe(p.site)}</h4>${matches(p)}<details><summary>Observed PTM, conditions and sources</summary><p>${safe(p.observation)}</p>${p.terminal_requirement?`<p class="mapping-warning">${safe(p.terminal_requirement)}</p>`:""}${p.retention_evidence?`<p><b>Retention evidence:</b> ${safe(p.retention_evidence)}</p>`:""}<p><b>Conditions:</b> ${safe(p.conditions)}</p><p><b>Functional evidence:</b> ${safe(p.functional_effect)}</p><p><b>Region context:</b> ${safe(p.region_detail)}. ${safe(p.numbering)} · ${safe(p.region_basis)}</p><p>${safe(p.evidence_status)}</p><p><a href="${safe(p.url)}" target="_blank" rel="noopener noreferrer">${safe(p.citation)}</a>${(p.additional_sources||[]).map(x=>` · <a href="${safe(x.url)}" target="_blank" rel="noopener noreferrer">${safe(x.citation)}</a>`).join('')}</p><p><a href="${safe(p.status_url)}" target="_blank" rel="noopener noreferrer">${safe(p.drug_status)}</a></p>${p.sequence_source_url?`<p><a href="${safe(p.sequence_source_url)}" target="_blank" rel="noopener noreferrer">Therapeutic reference sequence source</a></p>`:''}</details></article>`).join('')}`;
}

function renderMethodSources(){
  $('method-sources').innerHTML='<p class="footnote">Installed or cited does not mean executed. Optional model results are shown only after you run them. External TAP and sequence lookup require your own handoff and submission.</p>'+(report?.sources||[]).map(s=>`<details class="precedent-card"><summary><b>${safe(s.name)}</b> · ${safe(s.version)}</summary><p>${safe(s.role)}</p><p>${safe(s.implementation||'')}</p><p>${safe(s.limitations)}</p>${s.publication_url?`<p><a href="${safe(s.publication_url)}" target="_blank" rel="noopener noreferrer">${safe(s.publication_title)}</a></p>`:`<p>${safe(s.publication_title)}</p>`}${s.software_url?`<a href="${safe(s.software_url)}" target="_blank" rel="noopener noreferrer">Official software / source</a>`:''}${(s.evidence_sources||[]).map(e=>`<p><a href="${safe(e.url)}" target="_blank" rel="noopener noreferrer">${safe(e.title)}</a></p>`).join('')}${(s.artifacts||[]).map(a=>`<p class="source-hash">${safe(a.file)} · ${safe(a.status)}<br>SHA-256: ${safe(a.sha256||'Unavailable')}</p>`).join('')}</details>`).join('');
}

function renderPtmFocusControls(){
  for(const [id,anchor] of [['ptm-focus-controls','cdr-ptm-list'],['precedent-focus-controls','clinical-frameworks']]){
    let box=$(id);
    if(!box){box=document.createElement('div');box.id=id;(id==='ptm-focus-controls'?$(anchor).previousElementSibling:$(anchor)).before(box);}
    const c=window.ptmFocusCounts(report?.structure_analysis);
    box.innerHTML=`<button type="button" class="secondary" aria-pressed="${ptmFilter==='all'}">${ptmFilter==='all'?'Focus on medium + high PTMs':'Show all PTMs and precedents'}</button><p class="footnote">${safe(window.ptmFocusNote)}</p><p class="footnote">${report?.structure_analysis?`${c.high} high · ${c.medium} medium in focus · ${c.deferred} deferred · ${c.unknown} unknown`:'Structure assessment has not run: PTM priorities are unknown. Predict a structure above the local profilers, or show all sequence flags.'}</p>`;
    box.querySelector('button').onclick=()=>setPtmFilter(ptmFilter==='all'?'priority':'all');
  }
}
function setPtmFilter(value){
  ptmFilter=value;
  if(report){renderChain(selected);renderPrecedents();}
  if(report?.structure_analysis){renderStructuralReport(report.structure_analysis);if(structureView){structureView.focus=null;renderModelLegend();paintStructure();updateNglHighlights();}}
}

function assessmentRequest(){return {herceptin_tolerance_percent:herceptinTolerance,analysis_id:report.analysis_id,fold_job_id:foldState==='complete'||foldState==='partial'?activeJob:null,antiberty_job_id:Object.values(report.chains).some(c=>c.antiberty)?activeAntiBertyJob:null,propermab_job_id:report.propermab?activePropermabJob:null};}
async function refreshAssessment(){
 delete report.assessment;
 $('export-pdf').disabled=true;
 $('patch-risk-assessment').innerHTML='<p class="footnote">Patch-risk assessment will be available after the workflow finishes.</p>';
 const id=++assessmentGeneration,analysisId=report.analysis_id;
 const terminal=new Set(['complete','failed','unsupported']);
 const ready=(terminal.has(foldState)||foldState==='partial')&&terminal.has(profilerStates.anti)&&terminal.has(profilerStates.pm);
 if(!ready){
  $('assessment-status').textContent='Waiting for workflow completion';
  const labels={not_run:'waiting',running:'running locally',complete:'complete',failed:'failed',unsupported:'not supported',pending:'running locally',partial:'PDB available; exposure incomplete',idle:'waiting'};
  $('assessment-text').textContent=`Structure: ${labels[foldState]||foldState}. AntiBERTy: ${labels[profilerStates.anti]}. PROPERMAB: ${labels[profilerStates.pm]}. The overall assessment will be generated when all steps finish.`;
  return;
 }
 const statuses={pending:'Preliminary · structure running locally',failed:'Preliminary · prediction failed',partial:'Partial · exposure unavailable',unsupported:'Preliminary · heavy chain required'};
 $('assessment-status').textContent=statuses[foldState]||'Updating assessment…';
 $('assessment-text').textContent='Compiling local evidence…';
 try{
  const response=await fetch('/api/assessment',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(assessmentRequest())});
  const data=await response.json();if(id!==assessmentGeneration||report?.analysis_id!==analysisId)return;
  if(!response.ok)throw Error(data.detail||'Assessment unavailable');
  report.assessment=data;report.assessment_policy={herceptin_tolerance_percent:data.patch_risk?.herceptin_tolerance_percent??herceptinTolerance};$('export-pdf').disabled=false;$('assessment-status').textContent=statuses[foldState]||(Object.values(profilerStates).includes('running')?'Preliminary · local profilers running':Object.values(profilerStates).includes('failed')?'Workflow complete · a profiler failed':data.status==='Structure assessed'?'Workflow complete · assessment ready':'Workflow complete · '+data.status);
  renderPatchRisk(data.patch_risk);
  $('assessment-text').innerHTML=data.paragraphs.map(p=>`<p>${safe(p)}</p>`).join('');
 }catch(err){if(id!==assessmentGeneration)return;$('assessment-status').textContent='Assessment unavailable';$('assessment-text').textContent=err.message;}
}

function renderPatchRisk(risk){
 const root=$('patch-risk-assessment');if(!root||!risk)return;
 root.innerHTML=`<p class="footnote">Risk here means follow-up priority, not a probability of failure.</p><label>Herceptin comparison tolerance <select id="herceptin-tolerance">${[0,5,10,20].map(t=>`<option value="${t}" ${t===herceptinTolerance?'selected':''}>${t?`±${t}%`:'Off'}</option>`).join('')}</select></label>${risk.categories.map(p=>`<section class="patch-risk-card ${safe(p.priority)}"><button type="button" class="text-button" data-risk-mode="${safe(p.kind)}">${safe(p.label)} ↗</button><p><span class="priority ${safe(p.priority==='unknown'?'unassessed':p.priority)}">${safe(p.call)}</span></p><p>${safe(p.concern)}</p><details><summary>Why this call and what to test</summary>${p.reasons.map(r=>`<p>${safe(r)}</p>`).join('')}${p.comparisons.map(b=>`<p>${safe(b.label)}: ${safe(Number(b.query_value).toFixed(1))} ${safe(b.unit)}; therapeutic reference min–max ${safe(Number(b.min).toFixed(1))}–${safe(Number(b.max).toFixed(1))} (n=${safe(b.reference_n)}).</p>${b.reference_anchor?`<p class="herceptin-reference"><b>Herceptin: ${b.reference_anchor.value.toFixed(1)} ${safe(b.unit)}</b> · input difference ${b.reference_anchor.delta>=0?'+':''}${b.reference_anchor.delta.toFixed(1)} ${safe(b.unit)}. Comparative policy: both total and CDR-adjacent areas within ±${safe(risk.herceptin_tolerance_percent)}% qualify as Herceptin-like; this is not an experimentally validated cutoff.</p>`:''}`).join('')}${p.comparisons.length?'<p class="footnote">Range = lowest to highest calculated value among compatible therapeutic references. Within range means similar on this descriptor, not proven well behaved. Members were not selected for low aggregation or viscosity.</p>':''}<p><b>Follow-up:</b> ${safe(p.action)}</p></details></section>`).join('')}<details><summary>Patch-risk definitions and sources</summary><p>${safe(risk.definition)}</p><p>${safe(risk.limitation)}</p>${risk.sources.map(s=>`<p><a href="${safe(s.url)}" target="_blank" rel="noopener noreferrer">${safe(s.title)}</a></p>`).join('')}</details>`;
 root.querySelector('#herceptin-tolerance').onchange=e=>{herceptinTolerance=Number(e.target.value);void refreshAssessment();};
 root.querySelectorAll('[data-risk-mode]').forEach(b=>b.onclick=()=>setStructureMode(b.dataset.riskMode,null));
}

function startProfilers(){
 if(!report)return;
 if(profilerStates.anti!=='running'&&!Object.values(report.chains).every(c=>c.antiberty))void runAntiBerty();
 else if(Object.values(report.chains).every(c=>c.antiberty))profilerStates.anti='complete';
 if(report.chains.H&&report.chains.L){if(profilerStates.pm!=='running')void runPropermab();}
 else {profilerStates.pm='unsupported';$('propermab').disabled=true;$('propermab-status').textContent='Requires paired heavy and light chains.';}
 refreshAssessment();
}
