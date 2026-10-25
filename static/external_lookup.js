/* Optional, per-use handoff. No query sequence is placed in a URL or sent by this app. */
(() => {
  'use strict';
  const destination = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/';
  const cleanups = new WeakMap();
  function element(tag, text, className) {
    const node = document.createElement(tag);
    if (text) node.textContent = text;
    if (className) node.className = className;
    return node;
  }

  window.renderExternalLookup = function renderExternalLookup(report, container) {
    if (!container) return;
    cleanups.get(container)?.();
    container.replaceChildren();
    // Capture the analyzed result, never the editable input or a mutable global.
    const chains = Object.entries(report?.chains || {})
      .filter(([key, value]) => ['H', 'L'].includes(key) && typeof value.sequence === 'string' && value.sequence.length)
      .map(([key, value]) => ({ key, sequence: value.sequence, label: key === 'H' ? 'VH' : 'VL' }));
    let currentDialog = null;
    let disposed = false;
    cleanups.set(container, () => { disposed = true; currentDialog?.remove(); });
    const launch = element('button', 'Find therapeutic sequence precedents (external)');
    launch.type = 'button';
    launch.disabled = !chains.length;
    const status = element('p', 'Optional manual handoff to Oxford Thera-SAbDab. No external sequence search has been run.', 'note');
    status.setAttribute('role', 'status');
    container.append(launch, status);

    launch.addEventListener('click', () => {
      if (disposed) return;
      currentDialog?.remove();
      const dialog = element('dialog');
      currentDialog = dialog;
      dialog.style.cssText = 'max-width:580px;width:calc(100% - 40px);border:1px solid #abbac7;border-radius:16px;padding:28px;color:#172c3a;background:#fff;box-shadow:0 20px 80px #10283a55';
      const title = element('h2', 'Review external search handoff');
      title.id = `external-lookup-title-${Date.now()}`;
      dialog.setAttribute('aria-labelledby', title.id);
      const description = element('p', 'Destination: University of Oxford, Thera-SAbDab (opig.stats.ox.ac.uk). Approval copies one selected chain to your clipboard and opens its blank official search page. The app sends no sequence to Oxford.');
      const instruction = element('p', 'To transmit the sequence, you must paste it into the matching Heavy or Light chain field and click Search Thera-SAbDab on Oxford’s website. That search leaves your computer and is handled by Oxford. You can cancel here or leave the external page without submitting.');
      const label = element('label', 'Chain to copy from this analysis: ');
      const select = element('select');
      select.setAttribute('aria-label', 'Chain to copy');
      chains.forEach((chain, index) => {
        const option = element('option', `${chain.label} · ${chain.sequence.length} residues`);
        option.value = String(index);
        select.append(option);
      });
      label.append(select);
      const note = element('p', 'This prepares a search; it does not produce matches. Search similarity and clinical-stage records are context, not evidence that a shared PTM is safe. Copying uses your device’s clipboard settings.', 'note');
      const controls = element('div');
      controls.style.cssText = 'display:flex;gap:12px;flex-wrap:wrap;margin-top:20px';
      const approve = element('button');
      approve.type = 'button';
      const setLabel = () => { approve.textContent = `Approve: copy ${chains[Number(select.value)].label} and open Oxford`; };
      select.addEventListener('change', setLabel);
      setLabel();
      const cancel = element('button', 'Cancel');
      cancel.type = 'button';
      cancel.addEventListener('click', () => dialog.close());
      dialog.addEventListener('close', () => { dialog.remove(); if (currentDialog === dialog) currentDialog = null; });
      approve.addEventListener('click', async () => {
        if (disposed) return;
        const chain = chains[Number(select.value)];
        approve.disabled = true;
        // Clipboard call occurs only in this explicit approval handler.
        let copyTask;
        try { copyTask = navigator.clipboard.writeText(chain.sequence); }
        catch (error) { copyTask = Promise.reject(error); }
        window.open(destination, '_blank', 'noopener,noreferrer');
        dialog.close();
        try {
          await copyTask;
          if (!disposed) status.textContent = `${chain.label} copied locally; requested Oxford’s blank search page. Paste into its ${chain.key === 'H' ? 'Heavy' : 'Light'} chain field and submit only if you approve sending it. If the browser blocked the new tab, allow popups and repeat this handoff. No external search result has been retrieved.`;
        } catch (error) {
          if (!disposed) status.textContent = 'Clipboard access was blocked. The blank Oxford page was requested, but no sequence was copied or submitted. Repeat this handoff after enabling clipboard access, or copy your input manually.';
        }
      });
      controls.append(approve, cancel);
      dialog.append(title, description, instruction, label, note, controls);
      container.append(dialog);
      dialog.showModal();
    });
  };
})();
