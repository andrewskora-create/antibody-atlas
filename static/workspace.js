/* Result navigation and residue selection. Analysis stays local and unchanged. */
(function(root){
'use strict';
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const names={overview:'Overview',engineer:'Explore & engineer',compare:'Compare variants',evidence:'Evidence & methods'};
let context=null, host=null, active='overview';
function positionFor(analysis,chain,value,basis){
 const c=analysis?.chains?.[chain]; if(!c) throw Error('This chain is not part of the current analysis.');
 const raw=String(value).trim().toUpperCase();
 if(basis==='imgt'){
  const rows=(analysis.residue_table?.rows||[]).filter(r=>r.chain===chain&&String(r.imgt).replace(/^[HKL]/,'').toUpperCase()===raw);
  if(rows.length!==1) throw Error('IMGT position is absent or ambiguous in this chain. Choose sequence numbering or select a residue from the table.');
  return Number(rows[0].position ?? rows[0].sequence_position);
 }
 if(!/^\d+$/.test(raw)||Number(raw)<1||Number(raw)>c.sequence.length) throw Error(`Enter a sequence position from 1 to ${c.sequence.length}.`);
 return Number(raw);
}
function paneFor(id){
 if(['sec-explore-substitutions-at-a-position','sec-liability-substitution-proposals','sec-chemical-liabilities','sec-predicted-binding-surface','sec-3d-structure','sec-all-residues','sec-predicted-antigen-contact-surface'].includes(id))return 'engineer';
 if(['sec-molecule-level-properties','sec-surface-patches','sec-what-to-measure-next'].includes(id))return 'overview';
 return 'evidence';
}
function activate(view,focus=false){
 if(!host||!names[view])return; active=view;
 host.querySelectorAll('[data-work-pane]').forEach(p=>{p.hidden=p.dataset.workPane!==view;});
 host.querySelectorAll('[data-work-tab]').forEach(b=>{const on=b.dataset.workTab===view;b.setAttribute('aria-selected',String(on));b.tabIndex=on?0:-1;if(on&&focus)b.focus();});
 const jump=host.querySelector('#sec-jump');if(jump){for(const o of jump.options)o.hidden=!document.getElementById(o.value)?.closest(`[data-work-pane="${view}"]`);}
 window.dispatchEvent(new CustomEvent('atlas:viewchange',{detail:{view}}));
}
function reveal(id){const target=document.getElementById(id);if(!target)return;const p=target.closest('[data-work-pane]');if(p)activate(p.dataset.workPane);if(target.tagName==='DETAILS')target.open=true;target.scrollIntoView({block:'start',behavior:'auto'});}
function selectPosition(chain,position,navigate=true){
 if(!context)return;
 const c=context.analysis.chains?.[chain];if(!c||!Number.isInteger(Number(position))||position<1||position>c.sequence.length)return;
 const rows=context.table?.rows||[];const r=rows.find(r=>r.chain===chain&&Number(r.position??r.sequence_position)===Number(position));
 document.getElementById('pos-chain').value=chain;document.getElementById('pos-number').value=position;
 document.getElementById('pos-basis').value='sequence';
 const box=document.getElementById('selected-residue');
 if(box)box.innerHTML=`<strong>${esc(chain)} ${esc(c.sequence[position-1])}${esc(position)}</strong><span>IMGT ${esc(r?.imgt??'not assigned')} · ${esc(r?.imgt_region??'region not assigned')}</span><span>${esc((r?.liabilities||[]).map(x=>x.kind||x.motif).join('; ')||'No sequence liability flagged at this position')}</span>`;
 host.querySelectorAll('[data-seq-position]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.seqChain===chain&&Number(b.dataset.seqPosition)===Number(position))));
 if(navigate){activate('engineer');const p=document.getElementById('sec-explore-substitutions-at-a-position');if(p)p.open=true;}
}
function mount(results,ctx){
 context=ctx;host=results;active='overview';
 if(ctx.table)context.analysis={...ctx.analysis,residue_table:ctx.table};
 const original=Array.from(results.children);
 results.classList.add('atlas-workspace'); document.body.classList.add('has-workspace');
 const shell=document.createElement('div');shell.className='workspace-heading';
 const chainLengths=Object.entries(ctx.analysis.chains||{}).map(([k,c])=>`${k} ${c.sequence.length} aa`).join(' · ');
 shell.innerHTML=`<div><p class="workspace-eyebrow">CANDIDATE WORKSPACE · LOCAL ANALYSIS</p><h2>From sequence to a testable variant</h2><p>${esc(chainLengths)} · ${esc(String(ctx.analysis.scheme||'imgt').toUpperCase())} · pH ${esc(ctx.analysis.ph)} · ${esc(ctx.isotype||'Fv only')}</p></div><button type="button" class="ghost" id="edit-candidate">Edit inputs</button>`;
 results.replaceChildren(shell);
 const nav=document.createElement('div');nav.className='workspace-tabs';nav.setAttribute('role','tablist');nav.setAttribute('aria-label','Candidate workspace');
 nav.innerHTML=Object.entries(names).map(([k,n])=>`<button type="button" role="tab" data-work-tab="${k}" id="work-tab-${k}" aria-controls="work-pane-${k}" aria-selected="${k==='overview'}" tabindex="${k==='overview'?0:-1}">${n}</button>`).join('');results.append(nav);
 const panes={};for(const k of Object.keys(names)){const p=document.createElement('section');p.id=`work-pane-${k}`;p.dataset.workPane=k;p.setAttribute('role','tabpanel');p.setAttribute('aria-labelledby',`work-tab-${k}`);p.hidden=k!=='overview';panes[k]=p;results.append(p);}
 const intro=document.createElement('section');intro.className='workspace-intro';intro.innerHTML=`<div><h3>Choose the next experiment</h3><p>This screen identifies computational concerns and candidate changes. Manufacturing readiness still needs expression, stability, formulation and retained-function evidence.</p></div><div class="workspace-actions"><button type="button" class="primary" data-open-work="engineer">Explore substitutions</button><button type="button" class="ghost" data-open-work="compare">Compare variants</button></div>`;panes.overview.append(intro);
 for(const el of original){
  if(el.classList.contains('secnav')){panes.evidence.append(el);continue;}
  panes[el.matches('details.panel')?paneFor(el.id):'overview'].append(el);
 }
 const compare=document.getElementById('tool-variants');if(compare){compare.open=true;panes.compare.append(compare);}
 const local=document.createElement('div');local.id='engineering-local-preview';
 const lookup=document.getElementById('sec-explore-substitutions-at-a-position');if(lookup){lookup.open=true;panes.engineer.prepend(lookup);lookup.append(local);}
 const viewer=document.getElementById('sec-3d-structure');if(viewer&&lookup)lookup.after(viewer);
 const selected=document.createElement('div');selected.className='selected-residue';selected.id='selected-residue';selected.setAttribute('role','status');selected.innerHTML='<strong>Select a residue</strong><span>Enter a position, choose a sequence letter, or select a site in the structure.</span>';
 if(lookup)lookup.querySelector('.p-body').prepend(selected);
 const seq=document.createElement('details');seq.className='sequence-picker';seq.innerHTML='<summary>Select from the sequence</summary>'+Object.entries(ctx.analysis.chains||{}).map(([chain,c])=>`<p>${esc(chain)} chain · sequence positions</p><div class="sequence-strip">${Array.from(c.sequence).map((aa,i)=>`<button type="button" data-seq-chain="${esc(chain)}" data-seq-position="${i+1}" aria-pressed="false" title="${esc(chain)} ${aa}${i+1}" aria-label="Select ${esc(chain)} ${aa}${i+1}">${aa}<small>${i+1}</small></button>`).join('')}</div>`).join('');if(lookup)lookup.querySelector('.p-body').append(seq);
 const evidence=document.createElement('section');evidence.className='card evidence-gap';evidence.innerHTML='<h3>Experimental evidence</h3><p>No experimental assay results were entered in this screen. These remain unassessed:</p><ul><li>Expression yield, recovery and product purity</li><li>Monomer content, aggregation and stability under relevant conditions</li><li>Self-association, solubility and viscosity at the intended concentration</li><li>Binding and functional activity of the parent and each proposed variant</li></ul><p class="note">A supplied structure provides coordinates; it does not supply these measurements. Missing results are not a pass.</p>';panes.overview.append(evidence);
 const input=document.getElementById('input-card');if(input){input.hidden=true;document.getElementById('edit-candidate').onclick=()=>{input.hidden=!input.hidden;if(!input.hidden)input.scrollIntoView({block:'start'});};}
 nav.onclick=e=>{const b=e.target.closest('[data-work-tab]');if(b)activate(b.dataset.workTab);};
 nav.onkeydown=e=>{const keys=Object.keys(names),i=keys.indexOf(active);let n;if(e.key==='ArrowRight')n=(i+1)%keys.length;else if(e.key==='ArrowLeft')n=(i+keys.length-1)%keys.length;else if(e.key==='Home')n=0;else if(e.key==='End')n=keys.length-1;if(n!==undefined){e.preventDefault();activate(keys[n],true);}};
 if(results._workspaceClick)results.removeEventListener('click',results._workspaceClick);
 results._workspaceClick=e=>{const open=e.target.closest('[data-open-work]');if(open)activate(open.dataset.openWork);const s=e.target.closest('[data-seq-position]');if(s)selectPosition(s.dataset.seqChain,Number(s.dataset.seqPosition));};results.addEventListener('click',results._workspaceClick);
 activate('overview');
}
root.AtlasWorkspace={mount,activate,reveal,selectPosition,positionFor,paneFor};
if(typeof module!=='undefined')module.exports={positionFor,paneFor};
})(typeof window==='undefined'?globalThis:window);
