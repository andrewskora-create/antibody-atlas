'use strict';
const assert=require('node:assert/strict');
const eng=require('../static/engineering.js');
const h=eng._test;
assert.equal(h.num(0),'0.00');
assert.equal(h.num(null),'Not assessed');
assert.equal(h.num(NaN),'Not assessed');
assert.equal(h.signed(-1.25),'-1.25');
assert.equal(h.residueSelection('H','112A','NZ'),'112^A:H.NZ');
assert.equal(h.residueSelection('H','112B'),'112^B:H');
assert.equal(h.residueSelection('H','112 or all'),null);
assert.equal(h.safeURL('//external.example/x.pdb'),null);
assert.equal(h.safeURL('https://external.example/x.pdb'),null);
assert.equal(h.safeURL('/api/repack/abc/parent.pdb'),'/api/repack/abc/parent.pdb');
const parent={residues:[{chain:'H',pdb_number:'112A',sequence_position:110,aa:'K',relative_sidechain_sasa:1,sidechain_charge_at_ph:1}]};
const mutant={residues:[{chain:'H',pdb_number:'112A',sequence_position:110,aa:'E',relative_sidechain_sasa:1,sidechain_charge_at_ph:-1}]};
const atom={chainname:'H',resno:112,inscode:'A'};
assert.equal(h.atomProperty(atom,h.propertyIndex(parent),'charge'),0x1f6f9a);
assert.equal(h.atomProperty(atom,h.propertyIndex(mutant),'charge'),0xb91f43);
assert.equal(h.atomProperty(atom,new Map(),'charge'),0x9aa0a6);
assert.equal(h.atomProperty({...atom,inscode:'B'},h.propertyIndex(parent),'charge'),0x9aa0a6);
assert.equal(h.resolveMutation({chain:'H',position:110},'baseline',parent),'112^A:H');
assert.equal(h.resolveMutation({chain:'H',position:110},'baseline',null),null);
assert.equal(h.resolveMutation({baseline_structure:{chain:'L',residue:'51A'}},'baseline',null),'51^A:L');
const zero=h.metricRow('Charge',{status:'available',baseline:0,candidate:0,delta:0});
assert.match(zero,/0\.00/);assert.doesNotMatch(zero,/Not assessed/);
assert.match(h.metricRow('Charge',{status:'unavailable',baseline:0,candidate:0,delta:0,reason:'Missing coordinates'}),/Missing coordinates/);
assert.match(h.metricsHTML({ph:7.4,chains:{L:{metrics:{net_charge_at_ph:{status:'available',baseline:1,candidate:2,delta:1}}}}}),/L input-chain charge/);
const evidence=h.localEvidenceHTML({hydrogen_bonds:{lost:[],gained:[],lost_at_site:0,gained_at_site:0},clashes:{wild_type:0,mutant:0,introduced:0}});
assert.match(evidence,/<b>0<\/b> bonds lost/);assert.match(evidence,/Retained hydrogen bonds<\/summary><p>Not assessed/);
assert.match(evidence,/ΔΔG not calculated/);
const retained=h.localEvidenceHTML({hydrogen_bonds:{retained:[{baseline:{distance_a:3,angle_deg:150},candidate:{distance_a:2.8,angle_deg:155}}]}});assert.match(retained,/3.00 → 2.80/);
assert.equal(h.esc('<img src=x>'),'&lt;img src=x&gt;');
console.log('Engineering workspace helper contracts passed: identity, insertion codes, missing evidence, matched property maps, and charge coverage.');
// Exercise the actual async viewer lifecycle against the bundled NGL API contract.
(async()=>{
 const stages=[],disposed=[],removed=[],schemes=[];
 function component(url){return {url,reprList:[],structure:{},superpose(other,aligned){this.alignment={other,aligned};},removeAllRepresentations(){this.reprList=[];},addRepresentation(kind,options){const r={kind,options};this.reprList.push(r);return r;},autoView(){}};}
 global.NGL={ColormakerRegistry:{addScheme(builder){const item={};builder.call(item);schemes.push(item);return `scheme${schemes.length}`;},removeScheme(id){removed.push(id);}},Stage:class{
  constructor(el){this.element=el;stages.push(this);this.handlers=[];this.viewerControls={getOrientation:()=>({owner:this}),orient:(value,silent)=>{this.orientation=value;if(!silent)this.handlers.slice().forEach(f=>f());},signals:{changed:{add:f=>this.handlers.push(f),remove:f=>this.handlers.splice(this.handlers.indexOf(f),1)}}};}
  async loadFile(url){const c=component(url);(this.components??=[]).push(c);return c;}
  autoView(){}handleResize(){}dispose(){disposed.push(this);}
 }};
 const elements={};
 const arrangement=['paired','overlay'].map(value=>({dataset:{arrangement:value},setAttribute(key,value){this[key]=value;}}));
 const host={innerHTML:'',querySelector(selector){return elements[selector]??=(selector==='[data-color]'?{value:'identity'}:selector==='[data-style]'?{value:'cartoon'}:selector==='[data-sync]'?{checked:true}:{});},querySelectorAll(selector){return selector==='[data-arrangement]'?arrangement:[];}};
 const pair={baseline:{pdb_url:'/api/base.pdb',structure_analysis:parent},candidate:{pdb_url:'/api/mutant.pdb',structure_analysis:mutant}};
 const storedNGL=global.NGL;delete global.NGL;let loadCount=0;eng.configure({loadNGL:async()=>{loadCount++;global.NGL=storedNGL;}});
 const panel=await h.createViewer(host,pair,[{chain:'H',position:110}],'Predicted');
 assert.equal(loadCount,1);assert.equal(stages.length,3);assert.equal(panel.components.length,4);
 assert.equal(panel.components[1].alignment.other,panel.components[0]);
 assert.equal(panel.components[3].alignment.other,panel.components[2]);
 assert.equal(panel.components[0].reprList[1].options.sele,'112^A:H');
 stages[0].handlers[0]();assert.equal(stages[1].orientation.owner,stages[0]);
 elements['[data-sync]'].checked=false;stages[1].orientation=null;stages[0].handlers[0]();assert.equal(stages[1].orientation,null);
 elements['[data-color]'].value='charge';await panel.draw();
 assert.equal(schemes[0].atomColor(atom),0x1f6f9a);assert.equal(schemes[1].atomColor(atom),0xb91f43);
 await panel.draw();assert.equal(removed.length,4,'Changing maps releases previous color schemes');
 arrangement[1].onclick();assert.equal(elements['.eng-view-grid'].hidden,true);assert.equal(elements['.eng-overlay'].hidden,false);
 eng.reset();assert.equal(disposed.length,3);assert.equal(stages[0].handlers.length,0);
 // Delayed loading after reset must dispose, never attach stale coordinates.
 let resolveLoad;global.NGL.Stage=class {constructor(){this.viewerControls={};} loadFile(){return new Promise(resolve=>{(this.resolvers??=[]).push(resolve);resolveLoad??=[];resolveLoad.push(resolve);});}dispose(){}};
 const pending=h.createViewer(host,pair,[],'Predicted');eng.reset();resolveLoad.forEach(resolve=>resolve(component('/late.pdb')));const cancelled=await pending;assert.equal(cancelled.dead,true);
 console.log('Engineering viewer lifecycle passed: owned coordinate URLs, alignment, synchronized camera, map separation, disposal, stale-load cancellation.');
})().catch(error=>{console.error(error);process.exitCode=1;});
