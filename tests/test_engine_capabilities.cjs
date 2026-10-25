'use strict';
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync('static/start.html','utf8');
const boot=html.match(/\(async function boot\(\)\{[\s\S]*?\}\)\(\);/)[0];
async function render(caps,features){
 const elements=new Map(),get=s=>{if(!elements.has(s))elements.set(s,{});return elements.get(s);};
 const ctx=vm.createContext({api:async url=>url==='/api/health'?{structure_model:true,features}:caps,$:get,esc:s=>String(s),caps:null});
 await vm.runInContext(boot,ctx);return elements;
}
(async()=>{
 const base={endpoints:['/api/structure','/api/repack'],features:{top_five_models:true,stability_prediction:true},clinical_precedent_available:true};
 let el=await render({...base,stability_energy:{available:true,method:'Dedicated method'}});
 assert.match(el.get('#engines').textContent,/stability ready/);assert.match(el.get('#engines').title,/Dedicated method/);assert(!el.has('#error-slot'));
 el=await render({...base,stability_energy:{available:false,reason:'No licensed engine configured'}});
 assert.match(el.get('#engines').textContent,/stability unavailable/);assert.match(el.get('#engines').title,/No licensed engine/);assert(!el.has('#error-slot'),'An optional missing energy engine alone is not a stale server');
 el=await render({...base,features:{}});assert.match(el.get('#error-slot').innerHTML,/automatic top-five models/);assert.match(el.get('#error-slot').innerHTML,/newer than the server/);
 el=await render({...base,features:undefined},{top_five_models:true});assert(!el.has('#error-slot'),'Health features supported as fallback');
 console.log('Engine capabilities passed: stability availability and stale top-five backend detection.');
})().catch(e=>{console.error(e);process.exitCode=1;});
