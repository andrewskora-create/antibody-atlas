#!/usr/bin/env python3
"""Recompute a small public reference panel locally; never accepts query sequences.

Usage: python scripts/build_profiler_comparators.py --csv PATH --work-dir PATH
CSV must be the official Thera-SAbDab sequence download. It is fetched separately,
so the inference process can block network calls throughout model execution.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
NAMES = ['bevacizumab','infliximab','cetuximab','trastuzumab']
SOURCE = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def offline(event,args):
    if event in {'socket.connect','socket.getaddrinfo','socket.connect_ex'}:
        raise RuntimeError('Reference computation is offline; install model checkpoints first.')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--csv',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True);args=parser.parse_args()
    args.work_dir=args.work_dir.resolve();args.work_dir.mkdir(parents=True,exist_ok=True)
    with args.csv.open(encoding='utf-8-sig') as stream:
        rows={r['Therapeutic'].lower():r for r in csv.DictReader(stream) if r['Therapeutic'].lower() in NAMES}
    if set(rows)!=set(NAMES):raise ValueError('Expected exactly four named public reference records.')
    os.environ.update(HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1')
    sys.addaudithook(offline)
    from antiberty import AntiBERTyRunner
    from ImmuneBuilder import ABodyBuilder2
    import ImmuneBuilder
    from app.main import analyze, AnalyzeRequest
    from app.antiberty_analysis import analyze_antiberty,_checkpoint_provenance
    from app.propermab import run
    from app.ptm_precedents import PRECEDENTS
    runner=AntiBERTyRunner(device='cpu')
    predictor=ABodyBuilder2(model_ids=[1],numbering_scheme='imgt')
    versions={p:importlib.metadata.version(p) for p in ['antiberty','ImmuneBuilder','anarcii','freesasa','openmm','pdb2pqr','apbs-binary','numpy','scipy','numba']}
    checkpoint_dir=Path(ImmuneBuilder.__file__).parent/'trained_model'
    # Model1 is the only structure checkpoint used; record the filename from upstream.
    structure_files=[p for p in checkpoint_dir.iterdir() if p.is_file() and p.name=='antibody_model_1']
    panel={'panel_id':'approved-fv-illustrative-four-v1','generated_at_utc':datetime.now(timezone.utc).isoformat(),
           'description':'Locally computed descriptors for four selected approved antibody variable-domain pairs.',
           'limitations':['Illustrative, deliberately selected controls; n=4 is not a representative therapeutic population.',
                          'Observed min/max is not a safety interval, percentile calibration, clinical acceptance limit, or ranking of drug quality.',
                          'Fv sequences and single predicted structures omit Fc, formulation, concentration, dynamics and measured assays.',
                          'All four light chains are kappa. No lambda-specific reference range is established.',
                          'Single-model relaxation/hydrogen placement can vary between runs.',
                          'Upstream ANN and Ripley descriptors use unseeded default_rng Monte Carlo (1000 draws), so replicate values fluctuate.',
                          'Training-set overlap with these marketed antibodies is unknown; this is not an independent performance-validation set.'],
           'protocol':{'generator_script':'scripts/build_profiler_comparators.py','generator_script_sha256':digest(__file__),
                       'ph':7.4,'scheme':'imgt','numbering_scheme':'IMGT','structure_model':'ImmuneBuilder ABodyBuilder2 model 1',
                       'structure_checkpoint_artifacts':[{'filename':p.name,'sha256':digest(p)} for p in structure_files],
                       'preparation':'Exact input identity check, explicit IMGT numbering, local CHARMM hydrogen addition at pH7.4.',
                       'antiberty':_checkpoint_provenance(runner),'software_versions':versions,
                       'propermab_source_commit':'c97e41b95cccc103143d3123db953e694310c76e',
                       'propermab_adapter_version':'Atlas local PROPERMAB adapter 1'},
           'sources':[{'title':'Oxford Thera-SAbDab public sequence download','url':SOURCE,'retrieved_date':'2026-10-04',
                       'sha256':digest(args.csv),'clinical_status_note':'Download stage column dated Feb2025; separate approval citations retained.'}],
           'references':[]}
    output=ROOT/'data/profiler_comparators.json'
    for name in NAMES:
        print('Reference',name,flush=True)
        row=rows[name];seqs={'H':row['HeavySequence'],'L':row['LightSequence']}
        record=analyze(AnalyzeRequest(heavy=seqs['H'],light=seqs['L'],ph=7.4,scheme='imgt'))
        precedent=next(p for p in PRECEDENTS if p['antibody'].lower().startswith(name))
        ref={'id':name,'name':row['Therapeutic'],'sequences':seqs,'sequence_source_url':SOURCE,
             'sequence_sha256':{k:hashlib.sha256(v.encode()).hexdigest() for k,v in seqs.items()},
             'chain_types':{k:v['numbering']['chain_type'] for k,v in record['chains'].items()},
             'format':row['Format'],'isotype':row['CH1 Isotype'],'light_subclass':row['VD LC'],
             'approval_source_url':precedent['status_url'],'approval_status':precedent['drug_status']}
        try:
            results=analyze_antiberty(runner,seqs,scheme='imgt',batch_size=16)
            ref['antiberty']={'status':'complete','chains':results}
            print('AntiBERTy',[(k,round(v['pseudo_log_likelihood'],4)) for k,v in results.items()],flush=True)
        except Exception as exc:
            ref['antiberty']={'status':'unavailable','reason':str(exc)}
        pdb_path=args.work_dir/f'{name}.pdb'
        try:
            if not pdb_path.exists():
                model=predictor.predict(seqs);model.save(str(pdb_path))
            ref['structure']={'status':'complete','pdb_sha256':digest(pdb_path),'method':'ABodyBuilder2 model1; IMGT'}
            pm=run(record,pdb_path)
            ref['propermab']={'status':'complete','results':pm}
            print('PROPERMAB',len(pm['features']),pm['unavailable'],flush=True)
        except Exception as exc:
            ref['structure']={'status':'unavailable','reason':str(exc)}
            try:ref['propermab']={'status':'partial','results':run(record)}
            except Exception as inner:ref['propermab']={'status':'unavailable','reason':str(inner)}
        panel['references'].append(ref)
        output.write_text(json.dumps(panel,indent=2,allow_nan=False)+'\n')
    print('Saved',output,flush=True)

if __name__=='__main__':main()
