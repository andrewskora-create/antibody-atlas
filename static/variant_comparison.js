(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const MAX_CANDIDATES = 8;
  const AA = /^[ACDEFGHIKLMNPQRSTVWY]+$/;
  let baseline = null;
  let activeJob = null;
  let revision = 0;
  let mode = 'mutations';
  const escapeHtml = value => String(value ?? '').replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;');
  const fmt = value => value == null ? 'Unavailable' : typeof value === 'number' ? Number(value.toFixed(3)).toString() : String(value);
  const detail = value => typeof value === 'string' ? value : JSON.stringify(value);

  function setMode(next) {
    mode = next;
    $('variant-mutations-input').hidden = next !== 'mutations';
    $('variant-sequences-input').hidden = next !== 'sequences';
    for (const [id, selected] of [['variant-mode-mutations', next === 'mutations'], ['variant-mode-sequences', next === 'sequences']]) {
      $(id).classList.toggle('active', selected);
      $(id).setAttribute('aria-pressed', String(selected));
    }
    clearError();
  }

  function clearError() { $('variant-error').hidden = true; $('variant-error').textContent = ''; }
  function fail(message) { $('variant-error').textContent = message; $('variant-error').hidden = false; }
  function sequence(chain) { return baseline?.chains?.[chain]?.sequence || ''; }
  function parseMutationRow(line, index) {
    const parts = line.split('|').map(x => x.trim());
    if (parts.length !== 2 || !parts[0] || !parts[1]) throw Error(`Line ${index}: use name | H:N35Q,L:M20L.`);
    const name = parts[0];
    const mutations = parts[1].split(',').map(token => {
      const match = /^([HL]):([ACDEFGHIKLMNPQRSTVWY])([1-9]\d*)([ACDEFGHIKLMNPQRSTVWY])$/i.exec(token.trim());
      if (!match) throw Error(`Line ${index}: invalid substitution “${token.trim()}”. Use H:N35Q.`);
      const [, chainRaw, fromRaw, positionRaw, toRaw] = match;
      const chain = chainRaw.toUpperCase(), from_aa = fromRaw.toUpperCase(), to_aa = toRaw.toUpperCase(), position = Number(positionRaw);
      const reference = sequence(chain);
      if (!reference) throw Error(`Line ${index}: the original has no ${chain} chain.`);
      if (position > reference.length) throw Error(`Line ${index}: ${chain} has ${reference.length} residues; position ${position} is outside the input sequence.`);
      if (reference[position - 1] !== from_aa) throw Error(`Line ${index}: ${chain}:${from_aa}${position}${to_aa} does not match the original residue ${reference[position - 1]} at input position ${position}.`);
      if (from_aa === to_aa) throw Error(`Line ${index}: ${chain}:${from_aa}${position}${to_aa} makes no change.`);
      return {chain, position, from_aa, to_aa};
    });
    if (new Set(mutations.map(m => `${m.chain}:${m.position}`)).size !== mutations.length) throw Error(`Line ${index}: the same position appears more than once.`);
    return {name, mutations};
  }

  function parseSequenceRow(line, index) {
    const parts = line.split('|').map(x => x.trim());
    if (parts.length < 2 || parts.length > 3 || !parts[0]) throw Error(`Line ${index}: use name | VH | VL.`);
    const [name, heavyRaw, lightRaw = ''] = parts;
    const heavy = heavyRaw.toUpperCase().replace(/\s/g, ''), light = lightRaw.toUpperCase().replace(/\s/g, '');
    for (const chain of ['H', 'L']) {
      const candidate = chain === 'H' ? heavy : light, reference = sequence(chain);
      if (!!candidate !== !!reference) throw Error(`Line ${index}: the candidate must have the same ${chain} chain set as the original.`);
      if (candidate && !AA.test(candidate)) throw Error(`Line ${index}: ${chain} contains a nonstandard amino acid or FASTA header.`);
      if (candidate.length !== reference.length) throw Error(`Line ${index}: ${chain} length must match the original (${reference.length} aa); substitutions only are supported.`);
    }
    if (heavy === sequence('H') && light === sequence('L')) throw Error(`Line ${index}: the candidate is identical to the original.`);
    return {name, heavy, light};
  }

  function candidates() {
    const input = $(mode === 'mutations' ? 'variant-mutations' : 'variant-sequences').value;
    const rows = input.split(/\r?\n/).map(x => x.trim()).filter(Boolean);
    if (!rows.length) throw Error('Enter at least one candidate.');
    if (rows.length > MAX_CANDIDATES) throw Error(`Enter at most ${MAX_CANDIDATES} candidates per batch.`);
    const parsed = rows.map((row, i) => mode === 'mutations' ? parseMutationRow(row, i + 1) : parseSequenceRow(row, i + 1));
    if (new Set(parsed.map(x => x.name.toLowerCase())).size !== parsed.length) throw Error('Candidate names must be unique within a batch.');
    return parsed;
  }

  function metricCell(label, item) {
    if (!item || item.status === 'unavailable' || item.candidate == null || item.baseline == null) {
      return `<tr><th>${escapeHtml(label)}</th><td colspan="3" class="variant-unavailable">Unavailable${item?.reason ? ` · ${escapeHtml(item.reason)}` : ''}</td></tr>`;
    }
    const delta = item.delta == null ? '—' : `${item.delta > 0 ? '+' : ''}${fmt(item.delta)}`;
    return `<tr><th>${escapeHtml(label)}</th><td>${escapeHtml(fmt(item.baseline))}</td><td>${escapeHtml(fmt(item.candidate))}</td><td>${escapeHtml(delta)}</td></tr>`;
  }

  function siteText(site) {
    const location = `${site.chain || '?'}:${site.position ?? site.start ?? '?'}`;
    return `${location} ${site.kind || 'PTM'}${site.motif ? ` (${site.motif})` : ''}`;
  }

  function transitionText(label, transition) {
    if (!transition || transition.status === 'unavailable') return `${label}: unassessed`;
    if (transition.status === 'unchanged') return '';
    return `${label}: ${fmt(transition.baseline)} → ${fmt(transition.candidate)}`;
  }

  function retainedText(site) {
    const changes = [transitionText('priority', site.priority_transition), transitionText('exposure band', site.exposure_transition), transitionText('exposure (%)', site.exposure_percent_transition)].filter(Boolean);
    if (site.baseline_motif !== site.candidate_motif) changes.unshift(`motif: ${site.baseline_motif || '?'} → ${site.candidate_motif || '?'}`);
    return `${siteText(site)}${changes.length ? ` — ${changes.join('; ')}` : ' — no modeled transition'}`;
  }

  function renderPtm(ptm) {
    if (!ptm) return '<p class="variant-unavailable">PTM comparison unavailable.</p>';
    const lines = [];
    for (const [key, label] of [['added', 'Added'], ['removed', 'Removed']]) {
      const sites = ptm[key] || [];
      lines.push(`<div class="variant-ptm-group"><strong>${label} (${sites.length})</strong><span>${sites.length ? sites.map(site => escapeHtml(siteText(site))).join(' · ') : 'None'}</span></div>`);
    }
    const retained = ptm.retained || [];
    const changed = retained.filter(site => site.priority_transition?.status === 'changed' || site.exposure_transition?.status === 'changed' || site.exposure_percent_transition?.status === 'changed' || site.baseline_motif !== site.candidate_motif);
    lines.push(`<div class="variant-ptm-group"><strong>Retained (${retained.length})</strong><span>${changed.length} with a motif or modeled context change; ${retained.length - changed.length} unchanged or unassessed; expand for details.</span></div>`);
    if (retained.length) lines.push(`<details class="variant-retained"><summary>Show retained PTM sites, including unchanged oxidation flags</summary><ul>${retained.map(site => `<li>${escapeHtml(retainedText(site))}</li>`).join('')}</ul></details>`);
    return lines.join('');
  }

  function renderStructure(structure) {
    if (!structure || structure.status !== 'available') return `<p class="variant-unavailable">Structural comparison unavailable${structure?.reason ? `: ${escapeHtml(structure.reason)}` : '. Run separate folds for the original and each candidate.'}</p>`;
    const patches = Object.entries(structure.patch_counts || {}).map(([kind, count]) => `${kind}: ${fmt(count.baseline)} → ${fmt(count.candidate)} (${count.delta > 0 ? '+' : ''}${fmt(count.delta)})`);
    const areas = Object.entries(structure.patch_areas || {}).map(([kind, area]) => `${kind} exposed patch area: ${fmt(area.baseline)} → ${fmt(area.candidate)} Å²`);
    const coverage = Object.entries(structure.mapping_coverage_percent || {}).map(([chain, value]) => `${chain} mapping coverage: ${fmt(value.baseline)}% → ${fmt(value.candidate)}%`);
    return `<p>${escapeHtml([...patches, ...areas, ...coverage].join(' · ') || 'Separate folds available; no comparable structural descriptors were returned.')}</p>${structure.aggregate_reason ? `<p class="variant-unavailable">${escapeHtml(structure.aggregate_reason)}</p>` : ''}`;
  }

  function renderModelStatuses(models, title) {
    if (!models || !Object.keys(models).length) return '';
    return `<div class="variant-model-status"><strong>${escapeHtml(title)}</strong>${Object.entries(models).map(([model, entry]) => {
      const status = typeof entry === 'object' && entry !== null ? (entry.status || entry.state || 'Unavailable') : entry;
      const reason = typeof entry === 'object' && entry !== null ? (entry.reason || entry.error || '') : '';
      return `<span>${escapeHtml(model)}: ${escapeHtml(status)}${reason ? ` · ${escapeHtml(reason)}` : ''}</span>`;
    }).join('')}</div>`;
  }


  function renderPropermab(result) {
    if(!result || result.status === 'unavailable') return `<p class="variant-unavailable">PROPERMAB comparison unavailable: ${escapeHtml(result?.reason || 'Not run for both inputs.')}</p>`;
    const features=result.features||[], missing=result.unavailable||[];
    return `<p>${features.length} compatible descriptors compared · ${missing.length} unavailable. Changes are not calibrated improvements.</p><details><summary>PROPERMAB descriptor changes</summary><div class="variant-table-wrap"><table><thead><tr><th>Descriptor (unit)</th><th>Original</th><th>Candidate</th><th>Change</th></tr></thead><tbody>${features.map(f=>metricCell(`${f.name||f.key} (${f.unit||'unitless'})`,f)).join('')}</tbody></table></div>${missing.length?`<ul>${missing.map(f=>`<li>${escapeHtml(f.name||f.key)}: ${escapeHtml(f.reason)}</li>`).join('')}</ul>`:''}</details>`;
  }

  function renderComparison(comparison) {
    if (!comparison) return '<p class="variant-unavailable">Comparison unavailable.</p>';
    const mutationLabels = (comparison.mutations || []).map(m => `${m.chain}:${m.from || m.from_aa || '?'}${m.position}${m.to || m.to_aa || '?'}${m.baseline_region ? ` (${m.baseline_region})` : ''}`);
    const metricNames = [['isoelectric_point', 'Isoelectric point'], ['net_charge_at_ph', 'Net charge at pH'], ['gravy', 'Mean hydropathy'], ['humanness_percent', 'Framework identity (%)'], ['antiberty_mean_log_likelihood', 'AntiBERTy mean log likelihood']];
    const chainSections = Object.entries(comparison.chains || {}).map(([chain, value]) => `<div class="variant-chain"><h5>${chain === 'H' ? 'Heavy' : 'Light'} chain</h5><div class="variant-table-wrap"><table><thead><tr><th>Descriptor</th><th>Original</th><th>Candidate</th><th>Change</th></tr></thead><tbody>${metricNames.map(([key, label]) => metricCell(label, value.metrics?.[key])).join('')}</tbody></table></div></div>`).join('');
    const limitations = (comparison.limitations || []).map(x => `<li>${escapeHtml(x)}</li>`).join('');
    const validation = (comparison.next_validation || []).map(x => `<li>${escapeHtml(x)}</li>`).join('');
    return `<p class="variant-mutation-list"><strong>${comparison.mutation_count ?? mutationLabels.length} substitution(s):</strong> ${escapeHtml(mutationLabels.join(', ') || 'None reported')}</p>${(comparison.engineering_notes||[]).map(note=>`<p class="variant-unavailable">${escapeHtml(note)}</p>`).join('')}${chainSections}<h5>PTM motif changes</h5>${renderPtm(comparison.ptm)}<h5>Structural context</h5>${renderStructure(comparison.structure)}<h5>PROPERMAB</h5>${renderPropermab(comparison.propermab)}${comparison.conclusion ? `<p class="variant-conclusion">${escapeHtml(comparison.conclusion)}</p>` : ''}${limitations ? `<details><summary>Limits of this comparison</summary><ul>${limitations}</ul></details>` : ''}${validation ? `<details><summary>Suggested experimental checks</summary><ul>${validation}</ul></details>` : ''}`;
  }

  function renderResults(results, baselineModels) {
    const sorted = [...results].sort((a, b) => (a.comparison?.mutation_count ?? Infinity) - (b.comparison?.mutation_count ?? Infinity));
    $('variant-results').innerHTML = `<p class="variant-results-note">Original versus each candidate · single substitutions shown first. Changes describe local calculations, not a composite developability score.</p>${renderModelStatuses(baselineModels, 'Original model availability')}${sorted.map((item, i) => `<article class="variant-result"><h4>${escapeHtml(item.name || `Candidate ${i + 1}`)} <span>${item.status === 'complete' ? 'Compared' : escapeHtml(item.status || 'Unavailable')}</span></h4>${renderModelStatuses(item.model_status, 'Candidate model availability')}${item.error ? `<p class="variant-unavailable">${escapeHtml(item.error)}</p>` : renderComparison(item.comparison)}</article>`).join('')}`;
  }

  async function poll(jobId, currentRevision) {
    while (jobId === activeJob && currentRevision === revision) {
      await new Promise(resolve => setTimeout(resolve, 1800));
      if (jobId !== activeJob || currentRevision !== revision) return;
      const response = await fetch(`/api/variants/${encodeURIComponent(jobId)}`);
      const job = await response.json();
      if (jobId !== activeJob || currentRevision !== revision) return;
      if (!response.ok) throw Error(job.detail || 'Variant job status unavailable.');
      $('variant-status').textContent = `${job.status === 'running' ? 'Analyzing' : job.status === 'queued' ? 'Queued' : job.status}${job.progress ? ` · ${job.progress}` : ''}`;
      if (job.status === 'failed') throw Error(job.error || 'Variant comparison failed.');
      if (job.status === 'complete') {
        renderResults(job.results || [], job.baseline_models);
        $('variant-run').disabled = false;
        $('variant-status').textContent = `${(job.results || []).length} candidate(s) processed`;
        return;
      }
    }
  }

  async function run() {
    clearError();
    const currentRevision = ++revision;
    try {
      if (!baseline?.analysis_id) throw Error('Analyze the original antibody before comparing candidates.');
      if ($('variant-run-propermab').checked && (!sequence('H') || !sequence('L') || baseline.ph < 4 || baseline.ph > 10)) throw Error('PROPERMAB requires paired VH/VL and the original analysis pH must be 4–10.');
      const variants = candidates();
      $('variant-run').disabled = true;
      $('variant-results').innerHTML = '';
      $('variant-status').textContent = 'Submitting local comparison…';
      const response = await fetch('/api/variants', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({analysis_id: baseline.analysis_id, variants, run_structure: $('variant-run-structure').checked, run_antiberty: $('variant-run-antiberty').checked, run_propermab: $('variant-run-propermab').checked})});
      const job = await response.json();
      if (currentRevision !== revision) return;
      if (!response.ok) throw Error(job.detail || 'Could not start variant comparison.');
      activeJob = job.id || job.job_id;
      if (!activeJob) throw Error('Variant service returned no job ID.');
      await poll(activeJob, currentRevision);
    } catch (error) {
      if (currentRevision !== revision) return;
      fail(error.message || 'Variant comparison failed.');
      $('variant-status').textContent = '';
      $('variant-run').disabled = false;
    }
  }

  window.setVariantBaseline = report => {
    revision++;
    baseline = report;
    activeJob = null;
    $('variant-run').disabled = !report;
    $('variant-results').innerHTML = '';
    $('variant-status').textContent = report ? 'Original ready for comparison' : 'Analyze an original sequence first';
    clearError();
  };
  $('variant-mode-mutations').addEventListener('click', () => setMode('mutations'));
  $('variant-mode-sequences').addEventListener('click', () => setMode('sequences'));
  $('variant-run').addEventListener('click', run);
  $('variant-run-propermab').addEventListener('change', event => {
    if (event.target.checked) $('variant-run-structure').checked = true;
    $('variant-run-structure').disabled = event.target.checked;
  });
})();
