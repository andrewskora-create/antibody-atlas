const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
const elements=new Map();const el=id=>{if(!elements.has(id))elements.set(id,{querySelectorAll:()=>[],hidden:false,textContent:'',innerHTML:''});return elements.get(id)};
const ctx=vm.createContext({window:{addEventListener:()=>{}},document:{getElementById:el},fetch:async()=>({json:async()=>({})}),setTimeout:()=>{},console,Blob});
vm.runInContext(fs.readFileSync('static/app.js','utf8'),ctx);
let requests=0;ctx.fetch=async()=>{requests++;return {ok:true,json:async()=>({status:'Structure assessed',paragraphs:['Final assessment'],patch_risk:null})}};
(async()=>{
 vm.runInContext("report={analysis_id:'public',chains:{H:{},L:{}},provenance:{}};foldState='pending';profilerStates={anti:'not_run',pm:'not_run'};",ctx);
 await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,0);assert(el('export-pdf').disabled);
 vm.runInContext("foldState='complete';profilerStates={anti:'running',pm:'running'}",ctx);
 await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,0);
 vm.runInContext("profilerStates.anti='complete'",ctx);await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,0);
 vm.runInContext("profilerStates.pm='complete'",ctx);await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,1);assert.equal(el('export-pdf').disabled,false);
 assert(el('assessment-text').innerHTML.includes('Final assessment'));
 vm.runInContext("profilerStates.pm='running'",ctx);await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,1);assert.equal(vm.runInContext('report.assessment',ctx),undefined);
 vm.runInContext("profilerStates.pm='failed'",ctx);await vm.runInContext('refreshAssessment()',ctx);assert.equal(requests,2);assert(el('assessment-status').textContent.includes('failed'));
 console.log('Assessment waits for all steps; retries invalidate prior output; failed terminal state is explicit.');
})().catch(e=>{console.error(e);process.exitCode=1});
