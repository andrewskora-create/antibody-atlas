#!/usr/bin/env python3
"""Build a public, paired approved-Fv comparator panel with one local protocol.

The only input is the official Thera-SAbDab sequence CSV. Intermediate results
are checkpointed outside data; the published panel is replaced atomically only
after at least 16 full-input structural PROPERMAB results have completed.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import contextlib
import io
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SOURCE = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv'
# Selection is programmatic, not a curated list.
#
# This panel used to be 36 named molecules. That was its single biggest
# weakness and its own first stated limitation: a hand-picked subset is not a
# population, so an "observed range" over it has no defined coverage and a
# percentile over it means nothing. The criteria below are exactly the ones
# the old list was validated against, applied to every row in the snapshot
# instead of to thirty-six chosen ones. The panel is now the complete set of
# molecules meeting its own definition, which is what makes a percentile
# interpretable.
#
# The criteria are deliberately unchanged, so the existing per-molecule
# checkpoints stay valid and the new members are scored under the identical
# protocol rather than a revised one.
STAGE_COLUMN = "Highest_Clin_Trial (Feb '25)"
GENETICS_COLUMN = 'Genetics (Bispecifics delimited with semicolon)'
ELIGIBLE_STAGES = ('Approved',)
ELIGIBLE_GENETICS = ('Humanised', 'Genetically human')


def eligible(row):
    """Approved, human or humanised, paired whole monoclonal antibody."""
    return (row[STAGE_COLUMN] in ELIGIBLE_STAGES
            and row['Format'] == 'Whole mAb'
            and bool(row['HeavySequence']) and bool(row['LightSequence'])
            and row[GENETICS_COLUMN] in ELIGIBLE_GENETICS)


def select(rows):
    """Every eligible molecule in the snapshot, in a stable order."""
    return sorted({r['Therapeutic'].lower() for r in rows if eligible(r)})



def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _offline(event, args):
    if event in {'socket.connect', 'socket.getaddrinfo', 'socket.connect_ex'}:
        raise RuntimeError('Reference inference must remain offline.')


def _score(args):
    name, row, work_dir, prior_antiberty = args
    os.environ.update(HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    sys.addaudithook(_offline)
    from ImmuneBuilder import ABodyBuilder2
    from app.main import AnalyzeRequest, analyze
    from app.propermab import run

    seqs = {'H': row['HeavySequence'], 'L': row['LightSequence']}
    ref = {
        'id': name, 'name': row['Therapeutic'], 'sequences': seqs,
        'sequence_source_url': SOURCE,
        'sequence_sha256': {k: hashlib.sha256(v.encode()).hexdigest() for k, v in seqs.items()},
        'format': row['Format'], 'isotype': row['CH1 Isotype'], 'light_subclass': row['VD LC'],
        'approval_source_url': SOURCE, 'genetics': row['Genetics (Bispecifics delimited with semicolon)'], 'source_estimated_status': row['Est. Status'],
        'approval_status': "Approved in Thera-SAbDab Highest_Clin_Trial (Feb '25) snapshot",
        'origin_class': 'humanized' if name.endswith('zumab') else 'human',
        'origin_evidence': 'WHO historical INN -zumab/-umab source substem indicates humanized/human origin.',
    }
    try:
        record = analyze(AnalyzeRequest(heavy=seqs['H'], light=seqs['L'], ph=7.4, scheme='imgt'))
        ref['chain_types'] = {k: v['numbering']['chain_type'] for k, v in record['chains'].items()}
        # Sequence-level descriptors, which app/sequence_panel.py bands against.
        #
        # These MUST be written here. They survive the sequence-stripping step
        # that prepares the distributable stub -- they are computed values, not
        # source sequence text -- and the stub is what a user's rebuild
        # rehydrates from. A panel built without them silently disables every
        # sequence descriptor comparison in the product: panel_ranges() returns
        # {} and every chain reports "no panel comparison" with no error
        # anywhere. That is exactly what happened once, because the field was
        # originally produced by a one-off script that was never committed.
        ref['sequence_descriptors'] = {
            label: {
                'isoelectric_point': chain['metrics']['isoelectric_point'],
                'net_charge_at_ph': chain['metrics']['net_charge_at_ph'],
                'gravy': chain['metrics']['gravy'],
                'hydrophobic_windows': (chain['metrics'].get('aggregation_proxy') or {}).get('patch_count'),
                'length': len(chain['sequence']),
            }
            for label, chain in record['chains'].items()
        }
        if prior_antiberty is not None:
            ref['antiberty'] = prior_antiberty
        else:
            try:
                from antiberty import AntiBERTyRunner
                from app.antiberty_analysis import analyze_antiberty
                runner = AntiBERTyRunner(device='cpu')
                ref['antiberty'] = {'status': 'complete', 'chains': analyze_antiberty(runner, seqs, scheme='imgt', batch_size=16)}
            except Exception as exc:
                ref['antiberty'] = {'status': 'unavailable', 'reason': str(exc)}
        pdb = Path(work_dir) / f'{name}.pdb'
        fold_metadata = pdb.with_suffix('.fold.json')
        if fold_metadata.exists() and json.loads(fold_metadata.read_text()).get('status') != 'complete':
            raise ValueError('Prior structure refinement failed; excluded without rerunning PROPERMAB.')
        if not pdb.exists():
            model = ABodyBuilder2(model_ids=[1], numbering_scheme='imgt').predict(seqs)
            with contextlib.redirect_stdout(io.StringIO()) as fold_log:
                model.save(str(pdb))
            if 'FAILED TO REFINE' in fold_log.getvalue():
                fold_metadata.write_text(json.dumps({'status':'failed_refinement'}))
                raise ValueError('Structure refinement failed; excluded from reference comparisons.')
            fold_metadata.write_text(json.dumps({'status':'complete'}))
        ref['structure'] = {'status': 'complete', 'pdb_sha256': digest(pdb), 'method': 'ABodyBuilder2 model1; IMGT'}
        result = run(record, pdb)
        result.pop('selected_structure_analysis', None)
        ref['propermab'] = {'status': 'complete' if result.get('prepared_pdb_sha256') and not result.get('structure_partial') else 'partial', 'results': result}
    except Exception as exc:
        ref.setdefault('structure', {'status': 'unavailable', 'reason': str(exc)})
        ref['propermab'] = {'status': 'unavailable', 'reason': str(exc)}
    return ref


def _qualified(ref):
    pm = ref.get('propermab', {})
    if pm.get('status') != 'complete' or ref.get('structure',{}).get('status') != 'complete':
        return False
    result = pm['results']
    if result.get('comparison_protocol', {}).get('preparation') != 'Atlas full-input Fv preparation v2: exact residue identity; terminal reconstruction when needed; CHARMM hydrogens.':
        return False
    features = {f['key']: f.get('value') for f in result['features']}
    return all(features.get(key) is not None for key in ('hyd_patch_area', 'pos_patch_area', 'neg_patch_area'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', type=Path, required=True)
    ap.add_argument('--work-dir', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--limit', type=int, default=0,
                    help='score only the first N eligible molecules (for a smoke test)')
    args = ap.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    with args.csv.open(encoding='utf-8-sig') as stream:
        all_rows = list(csv.DictReader(stream))
    names = select(all_rows)
    rows = {r['Therapeutic'].lower(): r for r in all_rows if r['Therapeutic'].lower() in set(names)}
    if args.limit:
        names = names[:args.limit]
        rows = {n: rows[n] for n in names}
    if not names:
        raise ValueError('No eligible molecules in the source snapshot; check the stage column.')
    # Re-assert the criteria on the selected set. select() already applied them,
    # so a failure here means the two drifted apart, which is worth a crash.
    for name, row in rows.items():
        if not eligible(row):
            raise ValueError(f'{name} does not meet the panel criteria in the source snapshot')
    print(f'Selected {len(names)} eligible molecules from {len(all_rows)} snapshot rows', flush=True)
    from app.propermab import comparison_protocol
    from app.antiberty_analysis import _checkpoint_provenance
    from antiberty import AntiBERTyRunner
    from app.ptm_precedents import PRECEDENTS
    old = json.loads((ROOT/'data/profiler_comparators.json').read_text())
    prior = {r['id']: r for r in old['references']}
    provenance = _checkpoint_provenance(AntiBERTyRunner(device='cpu'))
    protocol = comparison_protocol(7.4, structure=True)
    protocol.update(scheme='imgt', numbering_scheme='IMGT',
                    generator_script='scripts/build_profiler_comparators_v2.py',
                    generator_script_sha256=digest(__file__),
                    antiberty=provenance,
                    propermab_source_commit='c97e41b95cccc103143d3123db953e694310c76e',
                    propermab_adapter_version='Atlas local PROPERMAB adapter 1')
    panel = {'panel_id': 'approved-human-humanized-fv-all-v3',
             'generated_at_utc': datetime.now(timezone.utc).isoformat(),
             'description': 'Locally scored descriptors for 36 selected human or humanized approved whole-antibody Fv pairs from the Thera-SAbDab February 2025 snapshot. These therapeutic references were not selected for favorable experimental behavior.',
             'limitations': [
                 'Complete within its definition -- every approved human or humanised paired whole mAb in the snapshot -- but approved antibodies are themselves a survivorship-filtered population, so this is not a representative sample of antibodies in general.',
                 'Neither the observed range nor a percentile over it is a safety interval, a clinical acceptance limit, or a ranking of drug quality. A percentile states where a molecule sits among approved drugs; it does not state whether that position is acceptable.',
                 'Fv sequences and single predicted structures omit Fc, formulation, concentration, dynamics and measured assays.',
                 'Structural reconstruction is a model hypothesis. Rebuilt terminal coordinates are uncertain.',
                 'Single-model relaxation/hydrogen placement can vary between runs.',
                 'Upstream ANN and Ripley descriptors use unseeded default_rng Monte Carlo (1000 draws), so replicate values fluctuate.',
                 'Training-set overlap with these marketed antibodies is unknown; this is not an independent performance-validation set.',
                 'Approval labels are from the Thera-SAbDab February 2025 stage column, not a live regulatory-status check.',
             ], 'protocol': protocol,
             'sources': [{'title': 'Oxford Thera-SAbDab public sequence download', 'url': SOURCE,
                          'retrieved_date': '2026-10-04', 'sha256': digest(args.csv),
                          'clinical_status_note': 'Highest clinical trial stage column dated February 2025; approval labels are source snapshot values.'}],
             'references': []}
    checkpoint = args.work_dir / 'profiler_comparators_v2.checkpoint.json'
    done = json.loads(checkpoint.read_text()) if checkpoint.exists() else {}
    if done.get('source_sha256') != digest(args.csv):
        done = {'source_sha256': digest(args.csv), 'generator_script_sha256': digest(__file__), 'references': {}}
    jobs = []
    for name in names:
        cached = done['references'].get(name)
        if cached and cached.get('sequences') == {'H':rows[name]['HeavySequence'], 'L':rows[name]['LightSequence']} and cached.get('propermab',{}).get('results',{}).get('comparison_protocol') == comparison_protocol(7.4, structure=True) and (_qualified(cached) or cached.get('propermab',{}).get('status') == 'excluded') and cached.get('sequence_descriptors'):
            cached['genetics'] = rows[name]['Genetics (Bispecifics delimited with semicolon)']
            cached['source_estimated_status'] = rows[name]['Est. Status']
            continue
        existing = prior.get(name)
        sequence_hash = {k: hashlib.sha256(rows[name][field].encode()).hexdigest() for k, field in [('H', 'HeavySequence'), ('L', 'LightSequence')]}
        anti = existing.get('antiberty') if existing and existing['sequence_sha256'] == sequence_hash else None
        jobs.append((name, rows[name], str(args.work_dir), anti))
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_score, job): job[0] for job in jobs}
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            ref = future.result()
            done['references'][name] = ref
            checkpoint.write_text(json.dumps(done, indent=2, allow_nan=False) + '\n')
            print(f'{name}: {ref["propermab"]["status"]}; qualified={_qualified(ref)}; completed={len(done["references"])}/{len(names)}', flush=True)
    panel['references'] = [done['references'][name] for name in names if name in done['references']]
    panel['description'] = (f"Locally scored panel of {len(panel['references'])} human or humanized approved "
        f"whole-antibody Fv pairs: every molecule in the Thera-SAbDab February 2025 snapshot meeting those "
        f"criteria, not a selected subset. Members were not chosen for favorable experimental behavior, and "
        f"failed structural calculations are excluded from compatible ranges.")
    good = sum(_qualified(ref) for ref in panel['references'])
    panel['control_counts'] = {'selected': len(names), 'propermab_full_input_success': good,
                               'light_kappa': sum(r['chain_types']['L'] == 'K' for r in panel['references'] if 'chain_types' in r),
                               'light_lambda': sum(r['chain_types']['L'] == 'L' for r in panel['references'] if 'chain_types' in r)}
    if good < 30:
        raise RuntimeError(f'Only {good} fully scored controls; existing published panel preserved. See {checkpoint}')
    # Refuse to publish a panel that would silently switch off the sequence
    # descriptor comparison. Nothing downstream raises when this field is
    # absent; the product just stops comparing and says nothing.
    without = [r.get('id') for r in panel['references'] if not r.get('sequence_descriptors')]
    if without:
        raise RuntimeError(
            f'{len(without)} reference(s) carry no sequence_descriptors '
            f'(e.g. {without[:3]}); publishing would disable every sequence descriptor '
            f'comparison in the product. Existing panel preserved. See {checkpoint}')
    output = ROOT/'data/profiler_comparators.json'
    staged = output.with_suffix('.json.tmp')
    staged.write_text(json.dumps(panel, indent=2, allow_nan=False) + '\n')
    staged.replace(output)
    print('Published', output, 'qualified', good, 'selected', len(names), flush=True)


if __name__ == '__main__':
    main()
