"""Place a sequence descriptor against the clinical reference panel.

A number on its own does not tell a reader whether to worry. pI 8.59 means
nothing without knowing what approved antibodies look like. This compares each
sequence descriptor with the range observed across the reference panel, which
is the same basis the structure-derived descriptors already use.

The range is an observed minimum to maximum across antibodies that reached
approval or Phase III. It is NOT an acceptance interval: the panel was not
selected for favourable experimental behaviour, plenty of approved antibodies
sit at its edges, and sitting outside it is a reason to look, not a failure.
"""
from __future__ import annotations

import statistics
from bisect import bisect_right
from functools import lru_cache
from typing import Any

METRICS = {
    'isoelectric_point': 'Isoelectric point',
    'net_charge_at_ph': 'Net charge',
    'gravy': 'GRAVY',
    'hydrophobic_windows': 'Hydrophobic windows',
}
BASIS = ('Observed minimum to maximum across the therapeutic reference panel, computed '
         'under this protocol, with a percentile position where the panel is large enough '
         'to support one. Neither is an acceptance interval: a percentile says where a '
         'molecule sits among approved antibodies, not whether that position is acceptable.')


@lru_cache(maxsize=1)
def panel_ranges() -> dict[str, dict[str, dict[str, float]]]:
    """chain -> metric -> {min, median, max, n}, or {} when the panel is absent."""
    try:
        from .reference_data import load
        panel = load('profiler_comparators')
    except Exception:
        return {}
    collected: dict[str, dict[str, list[float]]] = {}
    for reference in panel.get('references', []):
        for chain, values in (reference.get('sequence_descriptors') or {}).items():
            for metric in METRICS:
                value = values.get(metric)
                if isinstance(value, (int, float)):
                    collected.setdefault(chain, {}).setdefault(metric, []).append(float(value))
    return {
        chain: {
            metric: {'min': min(values), 'max': max(values),
                     'median': statistics.median(values), 'n': len(values),
                     # Kept sorted so band() can place a query without re-sorting
                     # on every descriptor of every chain of every molecule in a
                     # batch.
                     'sorted': sorted(values)}
            for metric, values in metrics.items() if values
        }
        for chain, metrics in collected.items()
    }


def normalise_chain(chain: str) -> str:
    """ANARCII reports kappa as K and lambda as L; the panel is keyed H / L.

    Without this a kappa chain falls through to the heavy-chain ranges and is
    judged against the wrong distribution entirely.
    """
    return 'H' if str(chain).upper().startswith('H') else 'L'


#: See profiler_interpretation.MIN_REFERENCES_FOR_PERCENTILE for the reasoning;
#: the two must agree, or the same molecule gets a percentile in one section of
#: the report and not in another.
MIN_REFERENCES_FOR_PERCENTILE = 60


def band(chain: str, metric: str, value: Any) -> dict[str, Any]:
    """Where one descriptor sits relative to the panel, for this chain class."""
    ranges = panel_ranges().get(normalise_chain(chain), {})
    window = ranges.get(metric)
    if not isinstance(value, (int, float)) or not window:
        return {'state': 'unknown', 'text': 'no panel comparison', 'range': None}
    low, high = window['min'], window['max']
    if value < low:
        state, text = 'below', f'below panel ({low:g} to {high:g})'
    elif value > high:
        state, text = 'above', f'above panel ({low:g} to {high:g})'
    else:
        state, text = 'within', f'within panel ({low:g} to {high:g})'
    result = {'state': state, 'text': text, 'n': window['n'],
              'range': {'min': low, 'median': window['median'], 'max': high},
              'basis': BASIS}
    # 'Within range' is a weak statement: it is true of a molecule sitting with
    # the bulk and of one pressed against the wall. The percentile separates
    # those two, which is the whole reason the panel was expanded to the full
    # approved set.
    ordered = window.get('sorted') or []
    if len(ordered) >= MIN_REFERENCES_FOR_PERCENTILE:
        percentile = round(100.0 * bisect_right(ordered, value) / len(ordered), 1)
        result['percentile'] = percentile
        if state == 'within':
            if percentile >= 90 or percentile <= 10:
                result['text'] = (f'within panel but in the tail: {percentile:g}th percentile '
                                  f'of {len(ordered)}')
            else:
                result['text'] = f'{percentile:g}th percentile of {len(ordered)}'
    else:
        result['percentile'] = None
    return result


def annotate_chain(metrics: dict[str, Any], chain_type: str) -> dict[str, Any]:
    """Band every supported descriptor of one chain."""
    source = {
        'isoelectric_point': metrics.get('isoelectric_point'),
        'net_charge_at_ph': metrics.get('net_charge_at_ph'),
        'gravy': metrics.get('gravy'),
        'hydrophobic_windows': (metrics.get('aggregation_proxy') or {}).get('patch_count'),
    }
    return {metric: {'value': value, 'label': METRICS[metric],
                     'chain_class': normalise_chain(chain_type),
                     **band(chain_type, metric, value)}
            for metric, value in source.items()}


def conclusion(annotated_by_chain: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """A plain-language verdict over the sequence descriptors."""
    outside = [(chain, entry) for chain, metrics in annotated_by_chain.items()
               for entry in metrics.values() if entry.get('state') in ('above', 'below')]
    if not any(panel_ranges().values()):
        return {'state': 'unknown',
                'text': 'No reference panel is built on this machine, so these descriptors '
                        'are reported without comparison. Build it with '
                        'setup/build_reference_data.py.'}
    if not outside:
        n = next(iter(next(iter(panel_ranges().values())).values()))['n']
        return {'state': 'within',
                'text': f'Every sequence descriptor sits inside the range observed across the '
                        f'{n}-antibody approved and Phase-III panel. That is similarity on these '
                        f'descriptors, not evidence of good behaviour.'}
    named = ', '.join(f"{chain} {entry['label'].lower()} {entry['value']:g} "
                      f"({entry['state']} {entry['range']['min']:g}–{entry['range']['max']:g})"
                      for chain, entry in outside)
    return {'state': 'outside',
            'text': f"{len(outside)} descriptor{'s' if len(outside) > 1 else ''} "
                    f"{'sit' if len(outside) > 1 else 'sits'} outside the panel range: {named}. The panel range is an observed minimum to maximum, "
                    f"not an acceptance interval; this marks where to look, not a failure."}
