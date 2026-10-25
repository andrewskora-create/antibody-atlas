"""One-page developability brief, with a link to the full report.

A four-page report is right when you are deciding about one molecule. It is the
wrong unit when you are ranking twenty. This renders a single page that answers
"does this one need attention, and where" in a few seconds, and links to the
detailed report for anyone who wants the evidence.

Every signal here is taken from the assessments the full report already makes.
No new threshold is invented for the brief: a light is amber here because the
screen already called that category review, not because the brief applied a
rule of its own.
"""
from __future__ import annotations

from datetime import datetime, timezone
from html import escape
from io import BytesIO

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

from .assessment import build_assessment, imgt_position
from .attribution import line as method_line

NAVY = '#163449'; TEAL = '#2a7187'; MUTED = '#5d7180'; PALE = '#edf5f6'
GREEN = '#397d6e'; AMBER = '#94601f'; RED = '#a8505b'
LIGHT = {'good': GREEN, 'review': AMBER, 'high': RED, 'unknown': MUTED}
LABEL = {'good': 'OK', 'review': 'REVIEW', 'high': 'ATTENTION', 'unknown': 'NO DATA'}


def _descriptor(pm, key):
    for feature in (pm.get('features') or []):
        if feature.get('key') == key and isinstance(feature.get('value'), (int, float)):
            return feature['value']
    return None


def _panel_position(pm):
    """How many compared descriptors sit outside the reference panel range."""
    marks = [b for b in (pm.get('summary') or {}).get('benchmarks', [])
             if b.get('reference_n') and b.get('query_relation') in
             ('within_range', 'above_range', 'below_range')]
    above = [b for b in marks if b['query_relation'] == 'above_range']
    return len(marks), above


def signals(data, assessment):
    """Category lights, each inherited from an assessment the screen already made."""
    sites = assessment['sites']
    cdr = [s for s in sites if str(s.get('region', '')).startswith('CDR')]
    framework = [s for s in sites if not str(s.get('region', '')).startswith('CDR')]
    patches = {p['kind']: p for p in assessment['patch_risk']['categories']}
    pm = data.get('propermab') or {}
    compared, above = _panel_position(pm)

    def patch_state(kind):
        entry = patches.get(kind)
        if not entry:
            return 'unknown', 'no structure'
        state = {'high': 'high', 'review': 'review', 'lower': 'good', 'unknown': 'unknown'}[entry['priority']]
        area = entry.get('largest_area_a2')
        return state, (f"{entry['count']} patches, max {area:,.0f} A^2" if area is not None else 'not assessed')

    high_cdr = sum(s.get('priority') == 'high' for s in cdr)
    cdr_state = 'high' if high_cdr else 'review' if cdr else 'good'
    cdr_detail = (f"{high_cdr} high, {len(cdr) - high_cdr} review" if cdr else 'none flagged')
    structure = data.get('structure_analysis')
    from .assessment import assessed
    lower_cdr = ([] if not structure else [
        s for s in structure.get('ptm_sites', [])
        if str(s.get('region', '')).startswith('CDR') and assessed(s)
        and s.get('priority') not in ('high', 'review')])
    if structure and not cdr:
        cdr_detail = 'No high/review-priority sites'
    if lower_cdr:
        cdr_detail += f"; {len(lower_cdr)} lower-priority motif sites"
    unassessed_cdr = ([] if not structure else [
        s for s in (structure.get('unassessed_ptm_sites', []) +
                    [s for s in structure.get('ptm_sites', []) if not assessed(s)])
        if str(s.get('region', '')).startswith('CDR') or s.get('region') in (None, '', 'Unassigned')])
    if not structure:
        cdr_state, cdr_detail = 'unknown', 'exposure not assessed'
    elif unassessed_cdr:
        cdr_state = cdr_state if cdr else 'unknown'
        cdr_detail += f"; {len(unassessed_cdr)} unassessed"

    hyd = patch_state('hydrophobic')
    neg = patch_state('negative')
    pos = patch_state('positive')
    # Concentration is a charge-surface proxy, not measured concentration
    # behaviour. Optional SCM availability must not look like the concern.
    scm = _descriptor(pm, 'scm')
    # A missing negative-charge result must not hide an observed positive patch.
    charge_name, charge = max((('negative', neg), ('positive', pos)),
                              key=lambda item: {'good': 0, 'unknown': 1, 'review': 2, 'high': 3}[item[1][0]])
    conc_state = 'review' if charge[0] in ('high', 'review') else charge[0]
    conc_detail = (f"{charge_name.capitalize()} patch concern" if conc_state == 'review'
                   else 'Charge-patch evidence incomplete' if conc_state == 'unknown'
                   else 'No charge-patch concern flagged')
    conc_detail += (f"; SCM {scm:,.0f}" if isinstance(scm, (int, float))
                    else '; optional SCM not assessed')

    # Framework deviation from germline, which the residue table already
    # computes per position. It belongs here rather than only on screen: if
    # it is worth flagging at all it is worth carrying onto the one-page
    # brief that gets forwarded. CDR mismatches are counted but never drive
    # the state -- a matured antibody is supposed to differ there, and
    # colouring that amber would bury the framework positions that matter.
    table = data.get('residue_table') or {}
    germline = table.get('germline_counts') or {}
    fr_off = germline.get('mismatch_framework') or 0
    cdr_off = germline.get('mismatch_cdr') or 0
    from .residue_table import _base, GERMLINE_MAX_IMGT
    expected_chains = set(data.get('chains') or {})
    v_framework = [r for r in table.get('rows', [])
                   if not str(r.get('imgt_region', '')).startswith('CDR')
                   and _base(r.get('imgt')) is not None
                   and _base(r.get('imgt')) <= GERMLINE_MAX_IMGT]
    covered_chains = {r.get('chain') for r in v_framework if r.get('germline_match') is not None}
    partial = (not v_framework or not expected_chains or not expected_chains <= covered_chains
               or any(r.get('germline_match') is None for r in v_framework))
    if not any(germline.values()):
        gl_state, gl_detail = 'unknown', 'no V allele assigned'
    elif fr_off:
        gl_state = 'review'
        gl_detail = f"{fr_off} V-FR differences" + ("; partial coverage" if partial else "")
    elif partial:
        gl_state, gl_detail = 'unknown', 'partial V-FR comparison'
    else:
        gl_state, gl_detail = 'good', 'Assessed V-FR match'

    return [
        ('CDR liabilities', cdr_state, cdr_detail),
        ('Aggregation', hyd[0], hyd[1]),
        ('Charge surface', charge[0], f"{charge_name}: {charge[1]}"),
        ('Concentration', conc_state, conc_detail),
        ('Germline FR', gl_state, gl_detail),
    ], framework, compared, above


def make_brief(data, title='Developability brief', detail_link=None,
               detail_label=None, structure_image=None):
    """One page. ``detail_link`` is a URL or a path relative to this file."""
    output = BytesIO()
    c = canvas.Canvas(output, pagesize=(612, 792))
    c.setTitle(str(title)[:100]); c.setAuthor('Antibody Atlas')
    chains = data.get('chains', {})
    assessment = build_assessment(data)
    provenance = data.get('provenance') or {}
    pm = data.get('propermab') or {}
    scheme = str(data.get('scheme', 'imgt')).upper()

    def text(value, x, y, width=516, height=40, size=10, color=NAVY, bold=False):
        style = ParagraphStyle('b', fontName='Helvetica-Bold' if bold else 'Helvetica',
                               fontSize=size, leading=size * 1.35, textColor=colors.HexColor(color))
        p = Paragraph(escape(str(value).replace('\n', ' ')), style)
        _, h = p.wrap(width, height)
        if h > height:
            lo, hi = 0, len(str(value))
            s = str(value)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if Paragraph(escape(s[:mid] + '...'), style).wrap(width, height)[1] <= height: lo = mid
                else: hi = mid - 1
            p = Paragraph(escape(s[:lo].rstrip() + '...'), style); _, h = p.wrap(width, height)
        p.drawOn(c, x, y - h)
        return h

    c.setFillColor(colors.HexColor(TEAL)); c.rect(0, 780, 612, 12, fill=1, stroke=0)
    text('ANTIBODY ATLAS  /  DEVELOPABILITY BRIEF', 48, 758, height=15, size=9, color=TEAL, bold=True)
    text(str(title)[:70], 48, 730, height=34, size=22)
    text(f"{datetime.now(timezone.utc).strftime('%Y-%m-%d')} | {scheme} | "
         f"pH {data.get('ph', '-')} | input {provenance.get('input_sha256', 'not recorded')[:12]}",
         48, 696, height=16, size=8.5, color=MUTED)

    lights, framework, compared, above = signals(data, assessment)

    # Verdict
    worst = 'high' if any(l[1] == 'high' for l in lights) else \
            'review' if any(l[1] == 'review' for l in lights) else \
            'unknown' if any(l[1] == 'unknown' for l in lights) else 'good'
    verdict = {'high': 'One or more categories need attention before this molecule advances.',
               'review': 'No category is flagged for attention; several warrant routine follow-up.',
               'good': 'No category is flagged in this screen.',
               'unknown': 'Assessment incomplete: one or more categories lack sufficient data.'}[worst]
    c.setFillColor(colors.HexColor(PALE)); c.roundRect(48, 636, 516, 46, 7, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(LIGHT[worst])); c.roundRect(48, 636, 5, 46, 2, fill=1, stroke=0)
    text(verdict, 64, 674, width=486, height=34, size=11, bold=True)

    # Category lights
    y = 606
    text('At a glance', 48, y, height=20, size=13, bold=True)
    y -= 28
    box_w = 516 / len(lights)
    for i, (name, state, detail) in enumerate(lights):
        x = 48 + i * box_w
        c.setFillColor(colors.HexColor('#f2f7f8')); c.rect(x, y - 74, box_w - 6, 74, fill=1, stroke=0)
        c.setFillColor(colors.HexColor(LIGHT[state])); c.rect(x, y - 74, box_w - 6, 4, fill=1, stroke=0)
        text(name, x + 9, y - 10, width=box_w - 22, height=24, size=9, bold=True)
        text(LABEL[state], x + 9, y - 32, width=box_w - 22, height=14, size=10,
             color=LIGHT[state], bold=True)
        # two short lines rather than one truncated one
        text(detail, x + 9, y - 46, width=box_w - 22, height=24, size=7, color=MUTED)
    y -= 90

    # Numbers that drive a ranking decision
    text('Key numbers', 48, y, height=20, size=13, bold=True); y -= 26
    cells = []
    for label, chain in list(chains.items())[:2]:
        m = chain.get('metrics', {})
        cells.append((f'pI {label}', f"{m.get('isoelectric_point','-')}"))
        cells.append((f'Charge {label}', f"{m.get('net_charge_at_ph','-')}"))
    scm = _descriptor(pm, 'scm')
    cells.append(('SCM', f"{scm:,.0f}" if isinstance(scm, (int, float)) else '-'))
    cells.append(('In panel range', f"{compared - len(above)}/{compared}" if compared else '-'))
    width = 516 / len(cells)
    for i, (k, v) in enumerate(cells):
        x = 48 + i * width
        text(k.upper(), x, y, width=width - 4, height=12, size=6.5, color=MUTED, bold=True)
        text(v, x, y - 13, width=width - 8, height=20, size=14, bold=True)
    y -= 46

    if above:
        text('Above the reference panel range: ' + ', '.join(b.get('label', b['key']) for b in above)
             + '. The panel range is an observed minimum to maximum, not an acceptance interval.',
             48, y, width=516, height=26, size=8.5, color=AMBER); y -= 32

    # What is actually flagged
    text('Flagged sites', 48, y, height=20, size=13, bold=True); y -= 24
    shown = assessment['sites'][:4]
    if not shown:
        y -= text('No site met the focused screen criteria. Deferred and unassessed sites are in the full report.',
                  48, y, width=516, height=26, size=9.5) + 8
    for s in shown:
        comp = s.get('comparators') or {}
        region = s.get('region', '')
        tag = ''
        if not str(region).startswith('CDR'):
            tag = (f" — carried by {comp.get('approved',0)} approved and {comp.get('phase_iii',0)} Phase-III drugs"
                   if comp.get('searched') and (comp.get('approved') or comp.get('phase_iii'))
                   else ' — no approved or Phase-III precedent found')
        kind = str(s.get('kind', '')).replace('Asparagine deamidation motif', 'Deamidation') \
            .replace('Aspartate isomerization motif', 'Isomerization') \
            .replace('Oxidation-susceptible residue', 'Oxidation') \
            .replace('N-linked glycosylation sequon', 'N-glycosylation')
        exposure = s.get('relative_sidechain_sasa')
        detail = (f"{s.get('chain','?')} {s.get('residue','')}{s.get('sequence_position','?')} "
                  f"({region} · {scheme} {imgt_position(s.get('number'))}) — {kind} {s.get('motif','')}, "
                  f"{s.get('priority','').upper()}"
                  + (f", {exposure:.0%} exposed" if isinstance(exposure, (int, float)) else '')
                  + tag)
        y -= text('• ' + detail, 48, y, width=516, height=26, size=9) + 5
    y -= 6

    if framework:
        y -= text(f"{len(framework)} framework site(s) are capped below high priority and checked against "
                  f"approved and Phase-III drugs; see the full report for the comparator list.",
                  48, y, width=516, height=24, size=8.5, color=MUTED) + 6

    # Next measurements
    text('Measure next', 48, y, height=20, size=13, bold=True); y -= 24
    actions = ['Peptide mapping for the flagged deamidation and isomerization sites; forced degradation for oxidation.',
               'HIC and self-interaction (AC-SINS or equivalent); viscosity at the intended formulation concentration.',
               'Binding and activity for any substitution before it is adopted.']
    for action in actions:
        y -= text('• ' + action, 48, y, width=516, height=24, size=9) + 4

    # Link to the detailed report
    if detail_link:
        c.setFillColor(colors.HexColor(PALE)); c.roundRect(48, 92, 516, 38, 6, fill=1, stroke=0)
        label = detail_label or 'Open the full four-page report'
        text(label, 62, 122, width=380, height=18, size=10.5, color=TEAL, bold=True)
        text('Evidence, precedents, biophysical table and method citations.',
             62, 106, width=380, height=14, size=8, color=MUTED)
        c.linkURL(detail_link, (48, 92, 564, 130), relative=0)

    text(method_line('ptm_priority', max_methods=3), 48, 80, width=516, height=22, size=7,
         color=MUTED)
    c.setStrokeColor(colors.HexColor('#dce7e9')); c.line(48, 43, 564, 43)
    text('Screening heuristics for ordering experiments, not predictions of behaviour. '
         'One page by design; the full report carries the evidence and limits.',
         48, 36, width=516, height=20, size=7.5, color=MUTED)
    c.showPage(); c.save()
    return output.getvalue()
