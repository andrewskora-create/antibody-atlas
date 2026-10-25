"""A self-contained HTML dashboard for one batch run.

No network, no CDN: the file opens from the results folder on any machine,
which matters because the point of a batch is to hand the folder to someone.
Rows link to each molecule's brief, full report, data and structure.

Colour follows the same rule as everywhere else in this project: green where
the screen found no qualifying signal, amber where it called review or a
descriptor sits outside the reference panel, red only where the screen already
rated something high priority. Nothing here invents a threshold.
"""
from __future__ import annotations

import html
import json
from urllib.parse import quote
from datetime import datetime, timezone

STYLE = """
:root{--bg:#fbfbfa;--panel:#fff;--ink:#1a1a18;--muted:#6b6b66;--line:#e4e4e0;--accent:#1f5f4f;
--high:#8c2f23;--high-bg:#fdf1ef;--mid:#8a5a12;--mid-bg:#fdf7ec;--low:#2c6248;--low-bg:#f0f7f2;--neutral:#f5f5f3}
@media (prefers-color-scheme:dark){:root{--bg:#17181a;--panel:#1f2023;--ink:#ecece8;--muted:#9a9a94;
--line:#323337;--accent:#5fae96;--high:#f0948a;--high-bg:#2b1d1b;--mid:#e3bd7c;--mid-bg:#2a241a;
--low:#88ccae;--low-bg:#1a2622;--neutral:#26272a}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1280px;margin:0 auto;padding:36px 18px 80px}
h1{font-size:30px;letter-spacing:-.02em;margin:.1em 0 .25em}
.sub{color:var(--muted);margin:0 0 22px}
.stats{display:flex;gap:1px;background:var(--line);border:1px solid var(--line);border-radius:10px;overflow:hidden;margin:0 0 22px}
.stat{background:var(--panel);padding:13px 18px;flex:1}
.stat .k{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);font-weight:600}
.stat .v{font-size:23px;font-weight:650;margin-top:2px;font-variant-numeric:tabular-nums}
table{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:10px;overflow:hidden;font-size:13.5px}
th{text-align:left;font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);
font-weight:600;padding:11px 10px;border-bottom:1px solid var(--line);cursor:pointer;user-select:none;white-space:nowrap}
th:hover{color:var(--ink)}
th[aria-sort]::after{content:' \\2195';opacity:.45}
th[aria-sort=ascending]::after{content:' \\2191';opacity:1}
th[aria-sort=descending]::after{content:' \\2193';opacity:1}
td{padding:10px;border-bottom:1px solid var(--line);vertical-align:middle;font-variant-numeric:tabular-nums}
tr:last-child td{border-bottom:0}
tbody tr:hover{background:var(--neutral)}
.name{font-weight:650;font-variant-numeric:normal}
.pill{display:inline-block;font-size:11px;font-weight:600;padding:2px 8px;border-radius:999px;white-space:nowrap}
.p-high{background:var(--high-bg);color:var(--high)}
.p-review{background:var(--mid-bg);color:var(--mid)}
.p-lower{background:var(--low-bg);color:var(--low)}
.p-none{background:var(--neutral);color:var(--muted)}
.links a{color:var(--accent);text-decoration:none;margin-right:9px;font-size:12.5px;white-space:nowrap}
.links a:hover{text-decoration:underline}
.failed{background:var(--high-bg);border:1px solid var(--high);border-radius:10px;padding:14px 16px;margin:22px 0 0}
.failed code{font-size:12px;color:var(--high)}
.note{color:var(--muted);font-size:12.5px;margin-top:18px;line-height:1.6}
"""

SCRIPT = """
document.querySelectorAll('th[data-key]').forEach((th, index) => th.onclick = () => {
  const table = th.closest('table'), body = table.tBodies[0];
  const ascending = th.getAttribute('aria-sort') !== 'ascending';
  table.querySelectorAll('th').forEach(h => h.removeAttribute('aria-sort'));
  th.setAttribute('aria-sort', ascending ? 'ascending' : 'descending');
  const rows = [...body.rows];
  rows.sort((a, b) => {
    const x = a.cells[index].dataset.sort ?? a.cells[index].textContent.trim();
    const y = b.cells[index].dataset.sort ?? b.cells[index].textContent.trim();
    const nx = parseFloat(x), ny = parseFloat(y);
    const both = !isNaN(nx) && !isNaN(ny);
    const result = both ? nx - ny : String(x).localeCompare(String(y));
    return ascending ? result : -result;
  });
  rows.forEach(r => body.appendChild(r));
});
"""


def _pill(value: str | None) -> str:
    mapping = {'high': ('p-high', 'high'), 'review': ('p-review', 'review'),
               'lower': ('p-lower', 'lower'), 'unknown': ('p-none', 'no data')}
    css, label = mapping.get(str(value), ('p-none', '–'))
    return f'<span class="pill {css}">{label}</span>'


def _cell(value, digits=2) -> str:
    if value is None:
        return '<span style="color:var(--muted)">–</span>'
    if isinstance(value, float):
        return f'{value:.{digits}f}'
    return html.escape(str(value))


def render(state: dict, panel) -> str:
    molecules = state.get('molecules', {})
    complete = {n: m for n, m in molecules.items() if m.get('status') == 'complete'}
    failed = {n: m for n, m in molecules.items() if m.get('status') != 'complete'}
    options = state.get('options', {})
    folded = options.get('fold')

    columns = [
        ('Molecule', 'name'), ('CDR high', 'cdr_sites_high'), ('CDR review', 'cdr_sites_review'),
        ('Framework', 'framework_sites'), ('Approved precedent', 'framework_approved_precedent'),
        ('pI H', 'pI_H'), ('pI L', 'pI_L'), ('Charge H', 'charge_H'), ('Charge L', 'charge_L'),
        ('Sequence vs panel', 'sequence_conclusion'),
    ]
    if folded:
        columns += [('Hydrophobic', 'hydrophobic_patch_priority'),
                    ('Negative', 'negative_patch_priority'),
                    ('Contacts', 'predicted_contacts')]

    head = ''.join(f'<th data-key="{html.escape(key)}">{html.escape(label)}</th>'
                   for label, key in columns) + '<th>Reports</th>'

    rows = []
    for name in sorted(complete):
        row = complete[name].get('row') or {}
        cells = []
        for label, key in columns:
            value = row.get(key)
            if key == 'name':
                cells.append(f'<td class="name" data-sort="{html.escape(name)}">{html.escape(name)}</td>')
            elif key.endswith('_patch_priority'):
                cells.append(f'<td data-sort="{html.escape(str(value))}">{_pill(value)}</td>')
            elif key == 'sequence_conclusion':
                css = {'within': 'p-lower', 'outside': 'p-review'}.get(str(value), 'p-none')
                text = {'within': 'typical', 'outside': 'outside range'}.get(str(value), '–')
                cells.append(f'<td data-sort="{html.escape(str(value))}"><span class="pill {css}">{text}</span></td>')
            elif key == 'cdr_sites_high':
                css = 'p-high' if (value or 0) else 'p-lower'
                cells.append(f'<td data-sort="{value or 0}"><span class="pill {css}">{value or 0}</span></td>')
            else:
                cells.append(f'<td data-sort="{value if value is not None else -1}">{_cell(value)}</td>')
        url_name = html.escape(quote(name, safe=''), quote=True)
        links = []
        if options.get('pdf', True):
            links += [f'<a href="{url_name}_brief.pdf">brief</a>',
                      f'<a href="{url_name}_report.pdf">report</a>']
        links.append(f'<a href="{url_name}.json">data</a>')
        if folded:
            links.append(f'<a href="{url_name}.pdb">structure</a>')
        cells.append('<td class="links">' + ''.join(links) + '</td>')
        rows.append('<tr>' + ''.join(cells) + '</tr>')

    seconds = [m.get('seconds') for m in complete.values() if isinstance(m.get('seconds'), (int, float))]
    stats = [
        ('Molecules', len(molecules)),
        ('Complete', len(complete)),
        ('Failed', len(failed)),
        ('Median time', f"{sorted(seconds)[len(seconds)//2]:.0f}s" if seconds else '–'),
        ('Structure', 'predicted' if folded else 'not run'),
    ]
    stat_html = ''.join(f'<div class="stat"><div class="k">{html.escape(k)}</div>'
                        f'<div class="v">{html.escape(str(v))}</div></div>' for k, v in stats)

    failed_html = ''
    if failed:
        items = ''.join(f'<li><b>{html.escape(n)}</b> — <code>{html.escape(str(m.get("error",""))[:220])}</code></li>'
                        for n, m in sorted(failed.items()))
        failed_html = (f'<div class="failed"><b>{len(failed)} molecule(s) did not complete.</b> '
                       f'Re-running the same command retries only these.<ul>{items}</ul></div>')

    warnings = ''
    if getattr(panel, 'warnings', None):
        items = ''.join(f'<li>{html.escape(w["name"])}: {html.escape(w["detail"])}</li>'
                        for w in panel.warnings)
        warnings = f'<div class="note"><b>Input warnings</b><ul>{items}</ul></div>'

    started = state.get('started_at_utc', '')[:19].replace('T', ' ')
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Developability batch</title><style>{STYLE}</style></head><body><div class="wrap">
<h1>Developability batch</h1>
<p class="sub">{len(molecules)} molecule(s) &middot; started {html.escape(started)} UTC &middot;
{html.escape(str(options.get('scheme','imgt')).upper())} &middot; pH {html.escape(str(options.get('ph','7.4')))}
&middot; click any column to sort</p>
<div class="stats">{stat_html}</div>
<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>
{failed_html}
{warnings}
<div class="note">
<p><b>How to read this.</b> Green means the screen found no qualifying signal in that category;
amber means it called review, or a sequence descriptor sits outside the range observed across the
approved and Phase-III reference panel; red means the screen already rated it high priority.
Framework PTM sites are capped below high and reported with the number of approved drugs carrying
the same motif at the same position. &ldquo;Approved precedent&rdquo; sums those counts across a
molecule&rsquo;s framework sites, so a high number means its framework liabilities are well
precedented.</p>
<p>Every value here is a screening heuristic for ordering experiments. None is a prediction of
aggregation, clearance, immunogenicity or viscosity, and a within-range result on any panel is not
a pass. Reference panels were not selected for favourable experimental behaviour, so their range
is an observed minimum to maximum, never an acceptance interval. Confirm liabilities analytically
before acting on them.</p>
</div>
</div><script>{SCRIPT}</script></body></html>"""
