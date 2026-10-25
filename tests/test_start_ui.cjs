// Exercise the actual served inline UI, not the legacy static/app.js.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync('static/start.html','utf8');
let source=html.match(/<script>([\s\S]*?)<\/script>/)[1];
source=source.replace(/\(async function boot\(\)\{[\s\S]*?\}\)\(\);/, '').replace('wireBatch();','');
const elements=new Map();
const el=id=>{ if(!elements.has(id)) elements.set(id,{id,value:'',files:[],innerHTML:'',textContent:'',disabled:false,
  dataset:{},listeners:{},classList:{add(){},remove(){}},focus(){this.focused=true;},scrollIntoView(){},
  querySelectorAll(){return [];},querySelector(){return null;},closest(){return null;},
  addEventListener(type,fn){this.listeners[type]=fn;},setAttribute(k,v){this[k]=v;}});return elements.get(id);};
const ctx=vm.createContext({console,Map,Set,Blob,URL,setTimeout,clearTimeout,
  document:{querySelector:s=>el(s.slice(1)),getElementById:el,querySelectorAll:()=>[],documentElement:{dataset:{}}},
  window:{addEventListener(){},removeEventListener(){},matchMedia:()=>({matches:false,addEventListener(){}})},
  fetch:async()=>({ok:true,blob:async()=>new Blob(['MODEL'])})});
vm.runInContext(source,ctx);
const run=code=>vm.runInContext(code,ctx);
// Render report functions for absent, complete, positive and partial structures.
ctx.mf={flags:[{id:'hydrophobic_patches',label:'Sequence hydrophobic windows',verdict:'amber',action:'run the structure step'},
  {id:'cdr_hydrophobicity',label:'CDR hydrophobicity',verdict:'amber'},
  {id:'charge_asymmetry',label:'Sequence charge imbalance',verdict:'amber'},
  {id:'net_charge',label:'Net charge',verdict:'green'}],fv:{},tradeoffs:[]};
let summary=run('renderMoleculeFlags(mf,{})');
assert(!summary.includes('Sequence hydrophobic windows'));
assert(!summary.includes('CDR hydrophobicity'));
assert(!summary.includes('Sequence charge imbalance'));
assert(summary.includes('composition-based estimates'));
for(const priority of ['lower','review','unknown']){
  ctx.risk={complete:priority!=='unknown',categories:[{kind:'hydrophobic',label:'Hydrophobic patches',priority,
    count:priority==='unknown'?null:priority==='lower'?0:2,largest_area_a2:120,action:'Measure HIC and aggregation under stress.'}]};
  let rendered=run('renderShortlist({},mf,risk,{})');
  assert(!rendered.includes('Sequence hydrophobic windows'));
  assert(!rendered.includes('run the structure step'));
  if(priority==='lower') {assert(rendered.includes('Baseline panel'));assert(rendered.includes('No additional targeted follow-up'));}
  if(priority==='review') {assert(rendered.includes('2 qualifying spatial patches'));assert(rendered.includes('HIC against reference antibodies'));}
  if(priority==='unknown') {assert(rendered.includes('coverage is incomplete'));assert(rendered.includes('Computational follow-up'));}
}
assert(run('renderShortlist({},mf,null,{})').includes('Surface properties not assessed'));
ctx.risk={complete:false,categories:[{kind:'hydrophobic',label:'Hydrophobic patches',priority:'review',count:1,provisional:true,action:'Measure HIC.'}]};
assert(run('renderShortlist({},mf,risk,{})').includes('Observed patches are provisional'));
// Guidance must match chemistry and preserve every selected site.
ctx.sites={rows:Array.from({length:8},(_,i)=>({chain:'H',aa:'N',sequence_position:i+1,imgt:i+1,imgt_region:'CDR1',exposure:0.4,verdict:'evaluate',paratope_contact:i===0,liabilities:[{kind:'Asparagine deamidation motif'},{kind:'Oxidation-susceptible residue'}]}))};
let plan=run('renderShortlist(sites,{}, {complete:true,categories:[{kind:"hydrophobic",priority:"lower",count:0}]},{})');
assert(plan.includes('N8'));assert(plan.includes('thermal/pH stress'));assert(plan.includes('Targeted oxidative stress'));
assert(plan.includes('If engineering a predicted contact'));assert(plan.includes('HDX-MS'));assert(plan.includes('cryo-EM'));
assert(!plan.includes('Thermal stress followed by mass spectrometry for every liability'));
ctx.sites.rows[0].exposure=null;
plan=run('renderShortlist({rows:[sites.rows[0]]},{},{categories:[]},{})');
assert(plan.includes('resolve structural coverage'));assert(!plan.includes('Targeted oxidative stress'));
ctx.sites.rows[0].exposure=0.4;ctx.sites.rows[0].paratope_contact=false;
ctx.sites.rows[0].liabilities=[{kind:'N-linked glycosylation sequon'}];
assert(run('renderShortlist({rows:[sites.rows[0]]},{},{categories:[]},{})').includes('no generic forced-stress assay'));
for(const kind of ['positive','negative']){
 ctx.patch={complete:true,categories:[{kind,priority:'review',count:1,largest_area_a2:100}]};
 const p=run('renderShortlist({}, {}, patch,{})');
 assert(p.includes(kind==='positive'?'Heparin-affinity chromatography':'DLS interaction parameter'));
 if(kind==='negative'){assert(!p.includes('First assay:</b> Heparin'));assert(p.includes('Direct viscosity'));}
}
assert(run('liabilityAssays([{kind:"Aspartate isomerization motif"}])[0].first').includes('isomer-resolving'));
(async()=>{
  el('ph').value='0'; assert.equal(run('inputPH()'),0);
  for(const bad of ['', 'NaN', '-1','15']) { el('ph').value=bad; assert.throws(()=>run('inputPH()'),/0 to 14/); }
  el('ph').value='6.5';assert.equal(run('inputPH()'),6.5);
  ctx.job={id:'j',status:'complete',rows:[{name:'high',pI:10},{name:'low',pI:2},{name:'missing',pI:null}],uiSort:{key:'pI',direction:1}};
  let table=run('renderBatchRows(job)');
  assert(table.indexOf('data-name="low"') < table.indexOf('data-name="high"'));
  assert(table.indexOf('data-name="high"') < table.indexOf('data-name="missing"'));
  assert(table.includes('aria-sort="ascending"'));
  assert(run("structureViewerNote('experimental')").includes('supplied structure'));
  assert(!run("structureViewerNote('experimental')").includes('Predicted coordinates'));
  // New result must attach a fresh viewer, and a stale PDB response must not load into it.
  const stages=[]; let pendingFetch;
  ctx.NGL={Stage:function(host){this.host=host;this.disposed=false;this.loads=0;
    this.dispose=()=>this.disposed=true;this.removeAllComponents=()=>{};
    this.loadFile=async()=>{this.loads++;return {addRepresentation(){},autoView(){}};};
    this.signals={clicked:{removeAll(){},add(){}}}; stages.push(this);},
    ColormakerRegistry:{addScheme:()=> 'colors',removeScheme(){}}};
  ctx.window.NGL=ctx.NGL;
  run('loadNGL=async()=>{}');
  const firstHost=el('viewer');
  await run("showStructure('first',null,'plain','cartoon',false)");
  assert.equal(stages.length,1); assert.equal(stages[0].host,firstHost);
  run('resetStructureViewer()'); elements.delete('viewer');
  await run("showStructure('second',null,'plain','cartoon',false)");
  assert.equal(stages.length,2); assert(stages[0].disposed);assert.equal(stages[1].host,el('viewer'));
  let ready, resolved=false;
  stages[1].tasks={onZeroOnce(callback){ready=callback;}};
  const surface=run("showStructure('surface',null,'plain','surface',false)").then(()=>resolved=true);
  await new Promise(r=>setImmediate(r));
  assert.equal(resolved,false,'Surface is not ready while NGL tasks remain');
  ready();await surface;assert.equal(resolved,true);
  delete stages[1].tasks;
  ctx.fetch=()=>new Promise(r=>pendingFetch=r);
  const stale=run("showStructure('stale',null,'plain','cartoon',false)");
  await Promise.resolve();await Promise.resolve();
  const loads=stages[1].loads;
  run('resetStructureViewer()');
  pendingFetch({ok:true,blob:async()=>new Blob(['STALE'])});await stale;
  assert.equal(stages[1].loads,loads);
  // Rapid control clicks during an unfinished draw must not change its legend.
  const panel=el('panel');panel.open=true;
  el('viewer-host').closest=selector=>selector==='[hidden]'?null:panel;
  const modeButton=el('mode');modeButton.dataset.mode='charge';
  const controls=[modeButton]; panel.querySelectorAll=()=>controls;
  const modes=el('viewer-modes');modes.querySelectorAll=()=>controls;
  let finish;const requests=[];
  ctx.showStructure=(...args)=>{requests.push(args);return new Promise(r=>finish=r);};
  run("wireTools({}, {}, 'x', null, 'experimental')");
  assert(modeButton.disabled);
  modes.listeners.click({target:{closest:()=>modeButton}});
  assert.equal(requests.length,1);assert.equal(requests[0][2],'hydrophobic');
  assert(el('viewer-legend').innerHTML.includes('Exposed hydrophobicity'));
  finish();await new Promise(r=>setImmediate(r));
  assert.equal(modeButton.disabled,false);
  modes.listeners.click({target:{closest:()=>modeButton}});
  assert.equal(requests.length,2);assert.equal(requests[1][2],'charge');finish();
  await new Promise(r=>setImmediate(r));
  // Opening a batch molecule restores its exact settings and discards an old uploaded structure.
  const batchButton=el('batch-button');batchButton.dataset.name='B';
  el('batch-results').querySelectorAll=selector=>selector==='.batch-open'?[batchButton]:[];
  el('structure-file').value='old-model.pdb';el('ph').value='8';el('isotype').value='IGHG4';
  let screens=0;el('run').click=()=>screens++;
  ctx.api=async()=>({chains:{H:'PUBLIC_H',L:'PUBLIC_L'}});
  ctx.batchJob={id:'j',options:{scheme:'kabat',ph:0,isotype:'IGHG1'}};
  run('wireBatchRowActions(batchJob)');await batchButton.onclick();
  assert.equal(el('structure-file').value,'');assert.equal(el('ph').value,0);
  assert.equal(el('scheme').value,'kabat');assert.equal(el('isotype').value,'IGHG1');assert.equal(screens,1);
  el('run').disabled=true;el('heavy').value='CURRENT';await batchButton.onclick();
  assert.equal(el('heavy').value,'CURRENT');assert.equal(screens,1);
  console.log('Served UI: pH, numeric batch sorting, source labels, stage replacement, stale load and rapid control regression checks passed.');
})().catch(e=>{console.error(e);process.exitCode=1});

// Missing-atom counts depend on structure coverage, not repeated render callbacks.
vm.runInContext(`{
 const atoms = [{chainname:'H',resno:1},{chainname:'H',resno:2},{chainname:'H',resno:3}];
 const structure = {eachAtom: fn => atoms.forEach(fn)};
 const index = new Map([['H:1',{aa:'A',rsa:.4}],['H:2',{aa:'G'}]]);
 if(countUnassessedAtoms(structure,index,'hydrophobic') !== 1) throw Error('Incorrect missing atom count');
 if(countUnassessedAtoms(structure,index,'plain') !== 0) throw Error('Plain mode cannot imply exposure coverage');
}`, ctx);
// Coordinates may use IMGT numbering and insertion codes; input indexes are not PDB keys.
vm.runInContext(`{
 const index=residueIndex({residues:[{chain:'H',sequence_position:52,pdb_number:'57',aa:'S',number:'H57',relative_sidechain_sasa:.5},{chain:'H',sequence_position:99,pdb_number:'112A',aa:'Y',number:'H112A',relative_sidechain_sasa:.2}],ptm_sites:[{chain:'H',sequence_position:52,priority:'high'}]});
 if(index.has('H:52') || index.get('H:57').sequence_position!==52 || index.get('H:57').flagged!=='high')throw Error('Sequence/PDB mapping mixed');
 const structure={eachAtom:fn=>fn({chainname:'H',resno:112,inscode:'A'})};
 if(countUnassessedAtoms(structure,index,'hydrophobic')!==0)throw Error('Insertion code lost');
}`,ctx);
