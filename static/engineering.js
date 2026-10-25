/* Engineering workspace: coordinates and evidence stay attached to their own variant. */
(function(root){
'use strict';
const finite = v => typeof v === 'number' && Number.isFinite(v);
const esc = v => String(v == null ? '' : v).replace(/[&<>"']/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const num = v => finite(v) ? (Math.abs(v)<.005?0:v).toFixed(2) : 'Not assessed';
const signed = v => finite(v) ? `${v>=.005?'+':''}${(Math.abs(v)<.005?0:v).toFixed(2)}` : 'Not assessed';
const AA = 'ACDEFGHIKLMNPQRSTVWY';
const HYDRO = {I:4.5,V:4.2,L:3.8,F:2.8,C:2.5,M:1.9,A:1.8,G:-.4,T:-.7,S:-.8,W:-.9,Y:-1.3,P:-1.6,H:-3.2,E:-3.5,Q:-3.5,D:-3.5,N:-3.5,K:-3.9,R:-4.5};
let context = {}, generation = 0, panels = new Set();
const localModels = new Map();
function residueSelection(chain,residue,atom){
  if(!/^[A-Za-z0-9]+$/.test(String(chain)) || !/^-?\d+[A-Za-z]?$/.test(String(residue))) return null;
  const m = /^(-?\d+)([A-Za-z]?)$/.exec(String(residue));
  if(atom && !/^[A-Za-z0-9]+$/.test(String(atom))) return null;
  return `${m[1]}${m[2]?`^${m[2]}`:''}:${chain}${atom?`.${atom}`:''}`;
}
function safeURL(url){ return typeof url==='string' && /^\/(?!\/)[A-Za-z0-9_/.?=&%-]+$/.test(url) ? url : null; }
function metricRow(label,m,unit='',evidence='Calculated'){
  const available=m && m.status==='available' && finite(m.baseline) && finite(m.candidate) && finite(m.delta);
  return `<tr><th scope="row">${esc(label)}${unit?` <small>${esc(unit)}</small>`:''}</th><td>${available?num(m.baseline):'Not assessed'}</td><td>${available?num(m.candidate):'Not assessed'}</td><td>${available?signed(m.delta):'—'}</td><td>${esc(available?evidence:m?.reason||'Both values are required.')}</td></tr>`;
}
function metricsHTML(comp){
  const molecule=comp.molecule||{}, rows=[];
  rows.push(metricRow('Overall charge',molecule.metrics?.net_charge_at_ph,'e'));
  rows.push(metricRow('Overall pI',molecule.metrics?.isoelectric_point));
  for(const [chain,data] of Object.entries(comp.chains||{})){
    rows.push(metricRow(`${chain} input-chain charge`,data.metrics?.net_charge_at_ph,'e'));
    rows.push(metricRow(`${chain} input-chain pI`,data.metrics?.isoelectric_point));
    rows.push(metricRow(`${chain} chain GRAVY`,data.metrics?.gravy));
  }
  for(const [chain,data] of Object.entries(molecule.per_chain||{})){
    if(molecule.isotype && molecule.isotype!=='fv') rows.push(metricRow(`${chain} assembled-chain charge × ${data.copies??'?'}`,data.metrics?.net_charge_at_ph,'e per chain'));
  }
  for(const [kind,value] of Object.entries(comp.structure?.patch_areas||{})) rows.push(metricRow(`${kind.replace(/_/g,' ')} patch area`,value,'Å²','Predicted coordinates'));
  for(const [kind,value] of Object.entries(comp.structure?.patch_counts||{})) rows.push(metricRow(`${kind.replace(/_/g,' ')} patches`,value,'','Predicted coordinates'));
  return `<h4>What changed?</h4><p class="eng-note">pH ${esc(molecule.ph??comp.ph??'not recorded')} · ${esc(String(molecule.basis||molecule.format||'Molecular basis not recorded').replace(/\.$/,''))}. Charge and pI are sequence calculations; surface descriptors require each variant’s own coordinates.</p><div class="eng-scroll"><table class="eng-table"><thead><tr><th>Property</th><th>Parent</th><th>Variant</th><th>Change</th><th>Evidence</th></tr></thead><tbody>${rows.join('')}</tbody></table></div>${comp.structure?.status!=='available'?`<p class="eng-note">Surface comparison: not assessed. ${esc(comp.structure?.reason||'Both variants need their own compatible structure analysis.')}</p>`:''}${comp.structure?.aggregate_reason?`<p class="eng-note">${esc(comp.structure.aggregate_reason)}</p>`:''}${(molecule.assumptions||[]).map(n=>`<p class="eng-note">${esc(n)}</p>`).join('')}`;
}
function propertyIndex(analysis){
  const result=new Map();
  for(const r of analysis?.residues||[]) if(r.pdb_number!=null && r.chain) result.set(`${r.chain}:${r.pdb_number}`,r);
  return result;
}
function atomProperty(atom,index,mode){
  const r=index.get(`${atom.chainname}:${atom.resno}${atom.inscode||''}`);
  if(!r || !finite(r.relative_sidechain_sasa)) return 0x9aa0a6;
  if(mode==='charge'){
    if(!finite(r.sidechain_charge_at_ph)) return 0x9aa0a6;
    const q=r.sidechain_charge_at_ph*Math.max(0,Math.min(1,r.relative_sidechain_sasa));
    return q>.12?0x1f6f9a:q<-.12?0xb91f43:0xdfe6e8;
  }
  if(!finite(HYDRO[r.aa])) return 0x9aa0a6;
  const h=Math.max(0,HYDRO[r.aa])/4.5*Math.max(0,Math.min(1,r.relative_sidechain_sasa));
  return h>.45?0xb91f43:h>.28?0xee9856:h>.14?0xe7d779:0x3f7fa6;
}
function resolveMutation(m,side,analysis){
  const mapped=m[`${side}_structure`];
  if(mapped) return residueSelection(mapped.chain,mapped.residue);
  const r=(analysis?.residues||[]).find(r=>r.chain===m.chain && r.sequence_position===m.position);
  return r?residueSelection(r.chain,r.pdb_number):null;
}
function removePanel(panel){
  if(!panel||panel.dead) return;
  panel.dead=true;
  panel.observer?.disconnect();
  panel.cleanups?.forEach(f=>f());
  panel.stages?.forEach(s=>{try{s.dispose();}catch(_){}});
  panel.schemes?.forEach(id=>{try{root.NGL.ColormakerRegistry.removeScheme(id);}catch(_){}});
  panels.delete(panel);
}
function reset(){ generation++; localModels.clear(); for(const panel of [...panels]) removePanel(panel); context={}; }
function configure(value){ context={...value}; }
async function createViewer(host,pair,mutations,kind){
  const panel={host,stages:[],components:[],schemes:[],cleanups:[],dead:false}; panels.add(panel);
  const stamp=generation;
  host.innerHTML=`<div class="eng-toolbar"><div role="group" aria-label="Structure arrangement"><button type="button" data-arrangement="paired" aria-pressed="true">Side by side</button><button type="button" data-arrangement="overlay" aria-pressed="false">Overlay</button></div><label>Color <select data-color><option value="identity">Parent / variant</option value="charge">Exposed charge</option><option value="hydrophobic">Exposed hydrophobicity</option></select></label><label>Style <select data-style><option value="cartoon">Cartoon + changed residues</option><option value="surface">Molecular surface</option></select></label><label><input type="checkbox" data-sync checked> Link views</label><button type="button" data-fit>Fit structures</button></div><div class="eng-view-grid"><section><h5><span class="eng-dot parent"></span>Parent</h5><div class="eng-canvas" data-view="parent" aria-label="Parent molecular structure"></div></section><section><h5><span class="eng-dot variant"></span>Variant</h5><div class="eng-canvas" data-view="variant" aria-label="Variant molecular structure"></div></section></div><div class="eng-overlay" hidden><h5>Parent / variant overlay</h5><div class="eng-canvas" data-view="overlay" aria-label="Aligned parent and variant molecular structures"></div></div><p class="eng-legend"></p><p class="eng-note" data-view-status role="status">Loading ${esc(kind)} coordinates…</p><p class="eng-note">Changed residues use ball-and-stick. Both coordinate sets are fitted to the parent by Cα alignment. Structural differences can include model variability.</p>`;
  const status=host.querySelector('[data-view-status]');
  if(!root.NGL && context.loadNGL){try{await context.loadNGL();}catch(error){status.textContent=`Could not load structure viewer: ${error.message}`;return panel;}}
  if(panel.dead||stamp!==generation)return panel;
  const urls=[safeURL(pair.baseline?.pdb_url),safeURL(pair.candidate?.pdb_url)];
  if(!urls.every(Boolean) || !root.NGL){status.textContent='Coordinates are not available for both structures. Run structure prediction for this comparison.';return panel;}
  const downloads=root.document?.createElement('p');
  if(downloads){downloads.className='eng-note';downloads.innerHTML=`<a href="${esc(urls[0])}" download>Download parent coordinates</a> · <a href="${esc(urls[1])}" download>Download variant coordinates</a><br>3D scope follows the coordinate files; an assembled-IgG charge calculation does not add constant-domain coordinates.`;host.appendChild(downloads);}
  const analyses=[pair.baseline.structure_analysis,pair.candidate.structure_analysis];
  panel.analyses=analyses;
  const indexes=analyses.map(propertyIndex);
  const propertyReady=indexes.every(i=>i.size>0);
  host.querySelectorAll('[data-color] option').forEach(o=>{if(o.value!=='identity')o.disabled=!propertyReady;});
  try{
    for(const name of ['parent','variant','overlay']) panel.stages.push(new root.NGL.Stage(host.querySelector(`[data-view="${name}"]`),{backgroundColor:'#f8fafb',quality:'medium'}));
    const [a,b,c,d]=await Promise.all([panel.stages[0].loadFile(urls[0],{ext:'pdb'}),panel.stages[1].loadFile(urls[1],{ext:'pdb'}),panel.stages[2].loadFile(urls[0],{ext:'pdb'}),panel.stages[2].loadFile(urls[1],{ext:'pdb'})]);
    if(panel.dead || generation!==stamp){removePanel(panel);return panel;}
    panel.components=[a,b,c,d];
    b.superpose(a,true); d.superpose(c,true);
    panel.selections=analyses.map((analysis,i)=>mutations.map(m=>resolveMutation(m,i?'candidate':'baseline',analysis)).filter(Boolean));
    panel.draw=async()=>{
      if(panel.dead||panel.busy)return;panel.busy=true;
      const controls=[...host.querySelectorAll('[data-color], [data-style]')];controls.forEach(c=>c.disabled=true);
      const mode=host.querySelector('[data-color]').value,style=host.querySelector('[data-style]').value;
      panel.components.forEach(component=>component.removeAllRepresentations());
      panel.schemes.forEach(id=>root.NGL.ColormakerRegistry.removeScheme(id));panel.schemes=[];
      panel.components.forEach((component,i)=>{
        const side=i%2,colour=side?'#087f83':'#65758c';
        let scheme=colour;
        if(mode!=='identity'){
          scheme=root.NGL.ColormakerRegistry.addScheme(function(){this.atomColor=atom=>atomProperty(atom,indexes[side],mode);});
          panel.schemes.push(scheme);
        }
        component.addRepresentation(style,{color:scheme,opacity:i>=2?.7:1,surfaceType:'av',useWorker:true});
        const selection=panel.selections[side].join(' or ');
        if(selection) component.addRepresentation('ball+stick',{sele:selection,color:side?'#b9690c':'#8b3866',radiusScale:.8});
      });
      host.querySelector('.eng-legend').innerHTML=mode==='identity'?'<span class="eng-dot parent"></span>Parent <span class="eng-dot variant"></span>Variant · orange / plum: changed residues':mode==='charge'?'<span class="eng-swatch" style="background:#b91f43"></span>Negative &lt; −0.12 <span class="eng-swatch" style="background:#dfe6e8"></span>Near neutral <span class="eng-swatch" style="background:#1f6f9a"></span>Positive &gt; +0.12 · charge × exposure, not electrostatic potential. Grey: not assessed.':'<span class="eng-swatch" style="background:#3f7fa6"></span>≤ 0.14 <span class="eng-swatch" style="background:#e7d779"></span>&gt; 0.14 <span class="eng-swatch" style="background:#ee9856"></span>&gt; 0.28 <span class="eng-swatch" style="background:#b91f43"></span>&gt; 0.45 · positive Kyte–Doolittle hydropathy / 4.5 × exposure. Grey: not assessed.';
      status.textContent='Drawing structures…';
      try{await Promise.all(panel.stages.map(stage=>new Promise((resolve,reject)=>{
        if(!stage.tasks?.onZeroOnce)return resolve();
        const timer=setTimeout(()=>reject(new Error('Surface drawing timed out; choose cartoon or retry.')),60000);
        const finish=()=>{clearTimeout(timer);resolve();};
        panel.cleanups.push(finish);stage.tasks.onZeroOnce(finish);
      })));
      if(!panel.dead)status.textContent=`${kind}. ${propertyReady?'Each surface map uses its own structure analysis.':'Surface property maps are not assessed for these coordinates.'}`;
      }catch(error){if(!panel.dead)status.textContent=error.message;}
      finally{panel.busy=false;if(!panel.dead)controls.forEach(c=>c.disabled=false);}
    };
    await panel.draw();
    if(panel.dead||generation!==stamp)return panel; panel.stages.forEach(s=>s.autoView(0));
    let syncing=false;
    panel.stages.slice(0,2).forEach((stage,i)=>{
      const sync=()=>{if(syncing||panel.dead||!host.querySelector('[data-sync]').checked)return;syncing=true;try{panel.stages[1-i].viewerControls.orient(stage.viewerControls.getOrientation(),true);}finally{syncing=false;}};
      stage.viewerControls.signals.changed.add(sync);
      panel.cleanups.push(()=>stage.viewerControls.signals.changed.remove(sync));
    });
    // Share a complete initial orientation so linked rotations have the same center and scale.
    panel.stages[1].viewerControls.orient(panel.stages[0].viewerControls.getOrientation(),true);
    host.querySelector('[data-color]').onchange=panel.draw;host.querySelector('[data-style]').onchange=panel.draw;
    host.querySelector('[data-fit]').onclick=()=>{panel.stages.forEach(s=>s.autoView(0));panel.stages[1].viewerControls.orient(panel.stages[0].viewerControls.getOrientation(),true);};
    host.querySelectorAll('[data-arrangement]').forEach(button=>button.onclick=()=>{
      const overlay=button.dataset.arrangement==='overlay';host.querySelector('.eng-view-grid').hidden=overlay;host.querySelector('.eng-overlay').hidden=!overlay;
      host.querySelectorAll('[data-arrangement]').forEach(b=>b.setAttribute('aria-pressed',String(b===button)));
      panel.stages.forEach(s=>s.handleResize());
      if(overlay) panel.stages[2].autoView(0);
    });
    if(root.ResizeObserver){panel.observer=new root.ResizeObserver(()=>panel.stages.forEach(s=>s.handleResize()));panel.observer.observe(host);}
  }catch(error){status.textContent=`Could not display coordinates: ${error.message}. Other comparison evidence remains available.`;}
  return panel;
}
function siteList(items){if(!Array.isArray(items))return 'Not assessed';return items.length?items.map(s=>`${esc(s.chain)}:${esc(s.position)} ${esc(s.kind)}${s.motif?` (${esc(s.motif)})`:''}`).join('; '):'None detected';}
async function comparisonEnergyHTML(comp){
  const mutations=comp.mutations||[];
  if(mutations.length!==1)return '<p class="eng-note">Combined folding ΔΔG is not calculated. Single-substitution estimates must not be added to predict combinations.</p>';
  const m=mutations[0],analysis=context.analysis;
  if(!analysis || comp.baseline_analysis_id!==analysis.analysis_id || comp.ph!==analysis.ph ||
      analysis.chains?.[m.chain]?.sequence?.[m.position-1]!==m.from || !context.foldId)
    return '<p class="eng-note">Folding ΔΔG is not available for this comparison’s matched parent structure.</p>';
  const key=JSON.stringify({chain:m.chain,position:Number(m.position),to:m.to,foldId:context.foldId,ph:Number(analysis.ph)});
  const pending=localModels.get(key);
  if(!pending)return '<p class="eng-note">No single-substitution ΔΔG result is available yet. Model this substitution in Explore &amp; engineer to obtain a local stability estimate, then rerun this comparison.</p>';
  try{
    const result=await pending;
    return '<h4>Single-substitution stability estimate</h4><p class="eng-note">From the original parent structure used in Explore &amp; engineer; the same estimate shown in the substitution shortlist. The independently predicted variant structures above provide separate surface comparisons.</p>'+freeEnergyHTML(result.free_energy);
  }catch(error){return '<p class="eng-note">Single-substitution ΔΔG is unavailable: '+esc(error.message)+'</p>';}
}
async function renderComparison(state,target){
  const stamp=generation;
  const energies=await Promise.all((state.results||[]).map(r=>r.status==='complete'?comparisonEnergyHTML(r.comparison||{}):''));
  if(stamp!==generation)return;
  for(const p of [...panels]) if(target.contains(p.host))removePanel(p);
  target.innerHTML=`<div class="eng-workspace"><p class="eng-note">Canonical amino-acid substitutions · compare the complete variant, including combinations · no automatic advancement decision</p>${(state.results||[]).map((r,i)=>{
    if(r.status!=='complete')return `<section class="eng-candidate"><h3>${esc(r.name)}</h3><p role="status">${esc(r.error||'Comparison not complete.')}</p></section>`;
    const c=r.comparison||{},mutations=(c.mutations||[]).map(m=>`${m.chain}:${m.from}${m.position}${m.to}`).join(', ');
    return `<section class="eng-candidate" data-candidate="${i}"><header><div><h3>${esc(r.name)}</h3><p class="eng-mutations">${esc(mutations)}</p></div><span class="eng-tag">${r.model_status?.structure?.status==='complete'?'Variant structure predicted':'Sequence comparison'}</span></header><div data-coordinate-host></div>${metricsHTML(c)}${energies[i]}<details class="eng-detail"><summary>Chemical liability changes</summary><dl><dt>Removed motifs</dt><dd>${siteList(c.ptm?.removed)}</dd><dt>Introduced motifs</dt><dd>${siteList(c.ptm?.added)}</dd><dt>Retained motifs</dt><dd>${siteList(c.ptm?.retained)}</dd></dl><p class="eng-note">Motifs indicate potential chemistry; they do not measure modification rates.</p></details>${(c.engineering_notes||[]).map(n=>`<p class="eng-caution">${esc(n)}</p>`).join('')}</section>`;
  }).join('')}<p class="eng-note">${esc(state.note||'')}</p></div>`;
  (state.results||[]).forEach((r,i)=>{
    if(r.status!=='complete')return;
    const host=target.querySelector(`[data-candidate="${i}"] [data-coordinate-host]`);
    const base=state.baseline_models?.structure||{},mut=r.model_status?.structure||{};
    const pair={baseline:{...base,structure_analysis:base.structure_analysis||state.baseline_structure_analysis},candidate:{...mut,structure_analysis:mut.structure_analysis||r.structure_analysis}};
    if(safeURL(base.pdb_url)&&safeURL(mut.pdb_url)){
      host.innerHTML='<button type="button" class="eng-open-viewer">View parent and variant structures</button>';
      host.querySelector('button').onclick=()=>createViewer(host,pair,r.comparison?.mutations||[],'Independently predicted structures');
    }else host.innerHTML=`<p class="eng-empty">${esc(mut.reason||'Select “Predict structure for each variant” when running the comparison to visualize new coordinates.')}</p>`;
  });
}
function bondRows(bonds,label){return bonds.map((item,i)=>{const b=item.candidate||item;return `<tr><td><button type="button" data-bond="${label}:${i}">${esc(b.donor_chain)}:${esc(b.donor_resname)}${esc(b.donor_resid)}.${esc(b.donor_atom)} → ${esc(b.acceptor_chain)}:${esc(b.acceptor_resname)}${esc(b.acceptor_resid)}.${esc(b.acceptor_atom)}</button></td><td>${item.baseline?`${num(item.baseline.distance_a)} → `:''}${num(b.distance_a)}</td><td>${item.baseline?`${num(item.baseline.angle_deg)} → `:''}${num(b.angle_deg)}</td><td>${esc(b.kind)}</td></tr>`;}).join('');}
function freeEnergyHTML(value){
  const f=value||{};
  if(f.status!=='available'||!finite(f.value_kcal_mol))return `<p class="eng-caution">ΔΔG not calculated${f.status==='failed'?' — calculation failed':''}. ${esc(f.reason||'A dedicated stability calculation is unavailable for this substitution.')}</p>`;
  const limits=Array.isArray(f.limitations)?f.limitations:[f.limitations].filter(Boolean);
  return `<div class="eng-evidence-strip"><span><b>${signed(f.value_kcal_mol)} kcal/mol</b> predicted folding ΔΔG</span></div><p class="eng-note">${esc(f.method||'Method not recorded')}${f.version?' · '+esc(f.version):''}. ${esc(f.sign_convention||'Positive values predict destabilization; negative values predict stabilization.')}</p>${limits.map(x=>`<p class="eng-caution">${esc(x)}</p>`).join('')}<p class="eng-note">A model estimate, not measured stability or antigen-binding free energy. Small differences should not determine a selection on their own.</p>`;
}
function modelEvent(detail){if(root.dispatchEvent&&root.CustomEvent)root.dispatchEvent(new root.CustomEvent('atlas:local-model',{detail}));}
async function modelLocal({chain,position,to,foldId=context.foldId,ph=context.analysis?.ph??7.4}){
  if(!foldId)throw new Error('Predict or upload a matched parent structure first.');
  if(!/^[HL]$/.test(chain)||!Number.isInteger(Number(position))||Number(position)<1||!AA.includes(to)||String(to).length!==1)throw new Error('Choose a valid chain, sequence position and canonical replacement.');
  const identity={chain,position:Number(position),to,foldId,ph:Number(ph)},key=JSON.stringify(identity),stamp=generation;
  let pending=localModels.get(key);
  if(!pending){
    pending=(async()=>{
      const response=await root.fetch('/api/repack',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({fold_job_id:foldId,chain,position:Number(position),replacement:to,ph:Number(ph)})});
      const result=await response.json();
      if(!response.ok)throw new Error(typeof result.detail==='string'?result.detail:'Local modeling failed.');
      return result;
    })();
    localModels.set(key,pending);
  }
  modelEvent({...identity,state:'modeling'});
  try{
    const result=await pending;
    if(stamp===generation)modelEvent({...identity,state:'complete',result});
    return result;
  }catch(error){
    if(localModels.get(key)===pending)localModels.delete(key);
    if(stamp===generation)modelEvent({...identity,state:'failed',error:error.message});
    throw error;
  }
}
function localEvidenceHTML(result){
  const h=result.hydrogen_bonds||{},cl=result.clashes||{};
  return `<h4>Local interactions</h4><div class="eng-evidence-strip"><span><b>${finite(h.lost_at_site)?h.lost_at_site:'Not assessed'}</b> bonds lost</span><span><b>${finite(h.gained_at_site)?h.gained_at_site:'Not assessed'}</b> bonds gained</span><span><b>${finite(cl.introduced)?cl.introduced:'Not assessed'}</b> clash count increase</span></div><p class="eng-note">Parent clashes: ${finite(cl.wild_type)?cl.wild_type:'not assessed'} → variant: ${finite(cl.mutant)?cl.mutant:'not assessed'}. This is the net increase in close-contact count, not a list of newly introduced clash identities. These are geometric contacts, not binding or stability measurements.</p><p class="eng-note">${esc(result.chemical_role_assumptions||'Hydrogen-bond atom-role and protonation assumptions were not recorded for this result.')}</p>${['lost','gained','retained'].map(label=>`<details class="eng-detail"${label!=='retained'?' open':''}><summary>${label[0].toUpperCase()+label.slice(1)} hydrogen bonds</summary>${Array.isArray(h[label])?(h[label].length?`<div class="eng-scroll"><table class="eng-table"><thead><tr><th>Donor → acceptor · click to highlight</th><th>Distance Å</th><th>Angle °</th><th>Type</th></tr></thead><tbody>${bondRows(h[label],label)}</tbody></table></div>`:'<p>None detected.</p>'):'<p>Not assessed.</p>'}</details>`).join('')}<p class="eng-note">Background rearrangement: ${finite(result.protocol_noise?.bonds_changed_away_from_the_site)?result.protocol_noise.bonds_changed_away_from_the_site:'not assessed'} bonds changed away from the site. ${esc(result.protocol_noise?.note||'')}</p><p class="eng-note">${esc(result.method||'')}</p>${freeEnergyHTML(result.free_energy)}`;
}
async function previewLocal({chain,position,to},target){
  target=target||root.document.getElementById('engineering-local-preview');
  if(!target)throw new Error('Local preview container is missing.');
  if(!/^[HL]$/.test(chain)||!Number.isInteger(Number(position))||Number(position)<1||!AA.includes(to)||String(to).length!==1){target.textContent='Choose a chain, input sequence position, and canonical replacement residue.';return;}
  if(!context.foldId){target.textContent='Predict or upload a matched parent structure before previewing a substitution.';return;}
  const stamp=generation,requestId=Symbol();target._atlasRequest=requestId;
  for(const p of [...panels]) if(target.contains(p.host))removePanel(p);
  target.innerHTML='<p role="status" class="eng-empty">Building the replacement side chain and relaxing parent and variant under the same conditions…</p>';
  target.scrollIntoView?.({block:'start',behavior:'auto'});
  try{
    const result=await modelLocal({chain,position,to});
    if(generation!==stamp||target._atlasRequest!==requestId)return;
    const pair=result.structures||{baseline:{pdb_url:result.parent_pdb_url},candidate:{pdb_url:result.mutant_pdb_url}};
    const selection=pair.selection||{chain:result.chain,residue:result.position};
    target.innerHTML=`<section class="eng-workspace eng-candidate"><header><div><h3>${esc(chain)}:${esc(result.from)}${esc(position)}${esc(to)}</h3><p>Local substitution model</p></div><span class="eng-tag">Locally relaxed coordinates</span></header><p class="eng-note">Starting coordinates: ${esc(result.source||'not recorded')}. Both parent and variant were locally relaxed; this is not an independent fold prediction or an experimentally solved mutant.</p><div data-local-coordinates></div><div data-local-evidence>${localEvidenceHTML(result)}</div></section>`;
    const panel=await createViewer(target.querySelector('[data-local-coordinates]'),pair,[{chain,position:Number(position),baseline_structure:selection,candidate_structure:selection}],'Locally relaxed parent and variant');
    if(generation!==stamp||target._atlasRequest!==requestId||panel.dead)return;
    target.querySelectorAll('[data-bond]').forEach(button=>button.onclick=()=>{
      const [kind,index]=button.dataset.bond.split(':'),entry=result.hydrogen_bonds[kind][Number(index)];
      if(!panel.components.length)return;
      panel.components.forEach((component,i)=>{
        const bond=(i%2===0?entry.baseline:entry.candidate)||entry;
        const donor=residueSelection(bond.donor_chain,bond.donor_resid,bond.donor_atom),acceptor=residueSelection(bond.acceptor_chain,bond.acceptor_resid,bond.acceptor_atom);
        if(!donor||!acceptor)return;
        (component.reprList||[]).filter(r=>r.name==='eng-bond').forEach(r=>component.removeRepresentation(r));
        if((kind==='lost'&&i%2===1)||(kind==='gained'&&i%2===0))return;
        const rep=component.addRepresentation('distance',{atomPair:[[donor,acceptor]],color:'#a24c13',labelColor:'#55300d',labelSize:1.5,name:'eng-bond'});
        const focus=component.addRepresentation('ball+stick',{sele:`${donor} or ${acceptor}`,color:'#a24c13',name:'eng-bond'});component.autoView(`${donor} or ${acceptor}`,0);
      });
    });
  }catch(error){if(generation===stamp&&target._atlasRequest===requestId)target.innerHTML=`<p class="eng-error" role="alert">${esc(error.message)}</p>`;}
}
function focusResidue({chain,position}){
  for(const panel of panels) panel.components.forEach((c,i)=>{const r=(panel.analyses?.[i%2]?.residues||[]).find(r=>r.chain===chain&&r.sequence_position===Number(position));const sele=r?residueSelection(r.chain,r.pdb_number):null;if(sele)c.autoView(sele,0);});
}
const api={configure,reset,renderComparison,previewLocal,modelLocal,focusResidue,_test:{finite,num,signed,esc,residueSelection,safeURL,metricRow,metricsHTML,comparisonEnergyHTML,propertyIndex,atomProperty,resolveMutation,localEvidenceHTML,freeEnergyHTML,createViewer}};
root.AtlasEngineering=api;
if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
