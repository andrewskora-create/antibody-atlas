'use strict';
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const eng=require('../static/engineering.js');
const html=fs.readFileSync('static/start.html','utf8');
let source=html.match(/<script>([\s\S]*?)<\/script>/)[1].replace(/\(async function boot\(\)\{[\s\S]*?\}\)\(\);/, '').replace('wireBatch();','');
const elements=new Map(),listeners={};
const el=id=>{if(!elements.has(id))elements.set(id,{id,innerHTML:'',textContent:'',value:'',files:[],dataset:{},classList:{add(){},remove(){}},querySelector(){return null;},querySelectorAll(){return [];},addEventListener(){},setAttribute(){}});return elements.get(id);};
const window={addEventListener(name,fn){listeners[name]=fn;},removeEventListener(){},matchMedia:()=>({matches:false,addEventListener(){}})};
const ctx=vm.createContext({console,Map,Set,Blob,URL,setTimeout,clearTimeout,window,document:{querySelector:s=>el(s.slice(1)),getElementById:el,querySelectorAll:()=>[],documentElement:{dataset:{}}}});
vm.runInContext(source,ctx);const run=s=>vm.runInContext(s,ctx);
const candidates='ACDEFGHIKLMNPQRSTVWY'.split('').filter(x=>x!=='S').map((aa,i)=>({mutation:`S52${aa}`,replacement_aa:aa,variant_spec:`H:S52${aa}`,conservativeness:{score:.8},selection_score:{score:.4,residue_similarity:.8,clinical_frequency:0,basis:'similarity_clinical_blend',reason:'Equal-weight'},clinical_precedent:{percent:0,observed:0,total:20,snapshot:'fixture'},germline_precedent:{},clean:true}));
ctx.d={chain:'H',position:52,original_aa:'S',candidates,candidate_count:19,top_model_candidates:candidates.slice(0,5).map(c=>c.variant_spec)};
const rendered=run('renderPositionResult(d)');
assert.equal((rendered.match(/class="preview-local"/g)||[]).length,19,'All candidates keep individual modeling');
assert.match(rendered,/Residue similarity/);assert.match(rendered,/Clinical reference frequency/);assert.match(rendered,/0 \/ 20/);assert.match(rendered,/Not modeled/);
assert.match(rendered,/50% residue similarity and 50% clinical reference frequency/);
ctx.c={model_state:'complete',local_model:{hydrogen_bonds:{lost_at_site:0,gained_at_site:1,retained_at_site:2},clashes:{introduced:0},free_energy:{status:'available',value_kcal_mol:0,method:'Test prediction'}}};
let cells=run('positionModelCells(c)');assert.match(cells.bonds,/−0 \/ \+1/);assert.match(cells.energy,/0.00/);assert.doesNotMatch(cells.energy,/Not calculated/);
ctx.c.local_model.free_energy={status:'available',value_kcal_mol:NaN};assert.match(run('positionModelCells(c)').energy,/Not calculated/);
ctx.c.local_model.free_energy={status:'failed',reason:'Engine missing'};assert.match(run('positionModelCells(c)').energy,/Calculation failed/);
assert.match(eng._test.freeEnergyHTML({status:'available',value_kcal_mol:-.75,method:'Test method',limitations:['Predicted only']}),/-0.75 kcal\/mol/);
assert.match(eng._test.freeEnergyHTML({status:'available',value_kcal_mol:NaN}),/not calculated/);
assert.match(eng._test.localEvidenceHTML({free_energy:{status:'available',value_kcal_mol:1.2,method:'Test method'}}),/\+1.20 kcal\/mol/);
(async()=>{
 let calls=[],inFlight=0,maxInFlight=0;
 window.AtlasEngineering={async modelLocal(request){calls.push(request);inFlight++;maxInFlight=Math.max(maxInFlight,inFlight);await Promise.resolve();inFlight--;listeners['atlas:local-model']({detail:{...request,state:'complete',result:{hydrogen_bonds:{lost_at_site:0,gained_at_site:1},free_energy:{status:'unavailable'}}}});return {};}};
 run("positionModelLookup={generation:viewerGeneration,chain:'H',position:52,foldId:'parent-a',ph:7.4,data:d}");
 await run('modelPositionTopFive(positionModelLookup)');
 assert.equal(calls.length,5);assert.equal(maxInFlight,1);assert.deepEqual(calls.map(c=>c.to),candidates.slice(0,5).map(c=>c.replacement_aa));
 assert(candidates.slice(0,5).every(c=>c.model_state==='complete'));assert.equal(candidates[5].model_state,undefined);
 assert.match(el('pos-model-status').textContent,/5 completed/);
 // Wrong source or pH must never overwrite a visible row.
 candidates[0].model_state='complete';
 listeners['atlas:local-model']({detail:{chain:'H',position:52,to:candidates[0].replacement_aa,foldId:'parent-a',ph:6,state:'failed'}});
 assert.equal(candidates[0].model_state,'complete');
 // Cancel the queue while a request is in flight: do not start another candidate.
 let release;calls=[];
 window.AtlasEngineering.modelLocal=request=>{calls.push(request);return new Promise(resolve=>release=resolve);};
 const pending=run('modelPositionTopFive(positionModelLookup)');run('positionModelLookup=null');release({});await pending;assert.equal(calls.length,1);
 // Same generation cancellation also applies when the entire report is replaced.
 run("positionModelLookup={generation:viewerGeneration,chain:'H',position:52,foldId:'parent-a',ph:7.4,data:d}");calls=[];
 const replacing=run('modelPositionTopFive(positionModelLookup)');run('viewerGeneration++');release({});await replacing;assert.equal(calls.length,1);
 // Automatic and manual requests for the same model share a single request.
 eng.reset();eng.configure({foldId:'parent-a',analysis:{ph:7.4}});
 let fetchCount=0,resolveFetch;global.fetch=()=>{fetchCount++;return new Promise(resolve=>resolveFetch=resolve);};
 const identity={chain:'H',position:52,to:'T'};
 const a=eng.modelLocal(identity),b=eng.modelLocal(identity);
 assert.equal(fetchCount,1);resolveFetch({ok:true,json:async()=>({free_energy:{status:'available',value_kcal_mol:0}})});
 assert.deepEqual(await a,await b);await eng.modelLocal(identity);assert.equal(fetchCount,1);
 // Changed pH/parent form separate cache keys; failure is retryable.
 global.fetch=async()=>{fetchCount++;return {ok:false,json:async()=>({detail:'Unavailable model'})};};
 await assert.rejects(eng.modelLocal({...identity,ph:6}),/Unavailable model/);
 await assert.rejects(eng.modelLocal({...identity,ph:6}),/Unavailable model/);assert.equal(fetchCount,3);
 console.log('Top-five model UI passed: all 19 controls, explicit missing values, sequential queue, cancellation, source/pH identity, cache reuse and retry.');
})().catch(e=>{console.error(e);process.exitCode=1;});
