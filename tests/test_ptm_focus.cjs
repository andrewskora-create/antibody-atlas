const vm=require('node:vm'),fs=require('node:fs'),assert=require('node:assert/strict');
const context={window:{}};vm.createContext(context);vm.runInContext(fs.readFileSync('static/outcome_summary.js','utf8'),context);
const focus=context.window.ptmFocusSite;
for(const region of ['CDR1','CDR2','CDR3']){
 for(const aa of ['W','M'])assert.equal(focus({kind:'Oxidation-susceptible residue',motif:aa,region,priority:'review',relative_sidechain_sasa:0.25}),true,`${region} ${aa} oxidation retained`);
}
assert.equal(focus({kind:'Oxidation-susceptible residue',motif:'W',region:'Framework',priority:'review',relative_sidechain_sasa:0.25}),false);
assert.equal(focus({kind:'Oxidation-susceptible residue',motif:'W',region:'CDR1',priority:'review',relative_sidechain_sasa:null}),false,'Unknown is not assumed assessed');
console.log('CDR Met/Trp focus checks passed.');
