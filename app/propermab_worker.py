"""Subprocess adapter for the real vendored PROPERMAB descriptor implementation."""
from __future__ import annotations
import contextlib
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor' / 'propermab'))

def _offline(event, args):
    if event in {'socket.connect', 'socket.getaddrinfo', 'socket.connect_ex'}:
        raise RuntimeError('Network access is disabled for local PROPERMAB analysis.')


def compute(data):
    sys.addaudithook(_offline)
    import numpy as np
    from propermab.features.seq_featurizer import SeqFeaturizer
    np.random.seed(0)
    seq = SeqFeaturizer((data['heavy'], data['light']), is_fv=True,
                        isotype='igg1', lc_type=data['light_type'], pH=data['ph'])
    result = {'features': [], 'unavailable': [], 'numbering_scheme': 'IMGT',
              'ph': data['ph'], 'scope': 'Official PROPERMAB Fv descriptors',
              'source_commit': 'c97e41b95cccc103143d3123db953e694310c76e',
              'model_checkpoint': None, 'adapter_version': 'Atlas local PROPERMAB adapter 1',
              'interpretation': 'Descriptors, not a fitted prediction of PK, viscosity, solubility, clinical performance or CMC readiness. No fitted endpoint model is distributed in the upstream repository.'}
    def add(obj, name, label, unit='', **kwargs):
        try:
            value = float(getattr(obj, name)(**kwargs))
            if not math.isfinite(value):
                raise ValueError('Nonfinite descriptor')
            result['features'].append({'key': name + ('_' + kwargs['prop'] if 'prop' in kwargs else ''), 'label': label, 'value': value, 'unit': unit})
        except Exception as exc:
            result['unavailable'].append({'feature': label, 'reason': type(exc).__name__ + ': calculation unavailable for this structure'})
    for name,label,unit in [('n_charged_res_fv','Fv charged-residue count','residues'),('vh_charge','VH charge','e'),
                            ('vl_charge','VL charge','e'),('fv_charge','Fv charge','e'),('fv_csp','VH × VL charge separation','e²')]:
        add(seq,name,label,unit)
    if data.get('pdb_path'):
        from propermab.features.struct_featurizer import StructFeaturizer
        f = StructFeaturizer(data['pdb_path'])
        for name,label,unit in [('hyd_asa','Apolar accessible area','Å²'),('hph_asa','Polar accessible area','Å²'),
            ('aromatic_asa','Aromatic accessible area','Å²'),('cdr_length','CDR-H3 length (IMGT)','residues'),
            ('aromatic_cdr','CDR aromatic count (IMGT)','residues'),('exposed_aromatic','Exposed aromatic count (RSA >5%)','residues'),
            ('net_charge','Structure formal charge','e'),('exposed_net_charge','Exposed structure charge','e'),
            ('net_charge_cdr','CDR formal charge (IMGT)','e'),('fv_chml','VH minus VL formal charge','e'),
            ('exposed_fv_chml','Exposed VH minus VL charge','e'),('dipole_moment','Dipole moment','Debye')]:
            add(f,name,label,unit)
        try:
            from propermab.structure.scm import raw_pdb_scm_scoring
            _, scm = raw_pdb_scm_scoring(data['pdb_path'])
            if not math.isfinite(float(scm)):
                raise ValueError('Nonfinite SCM')
            result['features'].append({'key':'scm','label':'Spatial charge map score (SCM)','value':float(scm),'unit':''})
        except Exception as exc:
            result['unavailable'].append({'feature':'SCM','reason':type(exc).__name__ + ': upstream calculation unavailable'})
        # The upstream hyd_moment uses uncentered coordinates; omit rather than
        # presenting a translation-dependent value as a molecular descriptor.
        result['unavailable'].append({'feature':'Hydrophobic moment','reason':'Withheld: upstream implementation depends on coordinate origin.'})
        for prop,label in [('pos','Positive'),('neg','Negative'),('aro','Aromatic')]:
            add(f,'ann_index',label+' average nearest-neighbor index',prop=prop,n=1000)
            add(f,'ripley_k',label+' Ripley K (6 Å)',prop=prop,distance=6.0,n=1000)
    else:
        result['unavailable'].append({'feature':'Structure descriptors','reason':'Predict a paired Fv structure and rerun PROPERMAB.'})
    if data.get('pdb_path') and data.get('patch_enabled'):
        try:
            import apbs_binary
            from propermab import defaults
            from propermab.features.feature_utils import calculate_patch_features
            defaults.system_config.override_defaults({
                'nanoshaper_binary_path': data['nanoshaper_path'],
                'atom_radii_file': str(ROOT/'vendor/propermab/tools/amber.siz'),
                'pdb2pqr_path': str(Path(sys.executable).parent/'pdb2pqr30'),
                'apbs_binary_path': str(apbs_binary.APBS_BIN_PATH),
                'multivalue_binary_path': str(Path(apbs_binary.APBS_BIN_PATH).parent/'multivalue'),
                'apbs_ld_library_paths': [str(apbs_binary.LIB_DIR)]})
            # APBS loads its shared libraries from this path. macOS reads
            # DYLD_LIBRARY_PATH; Linux reads LD_LIBRARY_PATH. Set the one the
            # running platform uses, preserving any existing value, so the
            # patch pipeline works on both.
            _lib = str(apbs_binary.LIB_DIR)
            _var = 'DYLD_LIBRARY_PATH' if sys.platform == 'darwin' else 'LD_LIBRARY_PATH'
            _existing = os.environ.get(_var)
            os.environ[_var] = f'{_lib}:{_existing}' if _existing else _lib
            os.chdir(data['work_dir'])
            patches = calculate_patch_features(data['pdb_path'], data['work_dir'])
            labels = {'hyd_patch_area':'Hydrophobic patch area', 'pos_patch_area':'Positive patch area',
                      'neg_patch_area':'Negative patch area', 'hyd_patch_area_cdr':'Hydrophobic patch area near CDR',
                      'pos_patch_area_cdr':'Positive patch area near CDR', 'neg_patch_area_cdr':'Negative patch area near CDR',
                      'heiden_score':'Heiden hydrophobic potential score'}
            if any(not math.isfinite(float(v)) for v in patches.values()):
                raise ValueError('Nonfinite patch descriptor')
            for key,value in patches.items():
                result['features'].append({'key':key,'label':labels[key], 'value':float(value),
                                           'unit':'' if key=='heiden_score' else 'Å²'})
            result['patch_pipeline'] = 'Complete: NanoShaper 0.7.8 / APBS 3.4.1 / PDB2PQR 3.7.1 / multivalue'
            result['patch_conditions'] = 'Upstream CHARMM/PDB2PQR default protonation and APBS defaults; not titrated to the sequence-charge pH control. Runtime differs from upstream APBS 3.0.0. Research descriptors, not therapeutic calibrated cutoffs.'
        except Exception as exc:
            result['unavailable'].append({'feature':'APBS/NanoShaper surface patches','reason':type(exc).__name__ + ': local patch calculation unavailable'})
    else:
        result['unavailable'].append({'feature':'APBS/NanoShaper surface patches','reason':'Requires a matching local structure and configured native dependencies.'})
    return result

if __name__ == '__main__':
    payload = json.load(sys.stdin)
    # Keep library output away from the JSON protocol; no sequence logs are retained.
    with open('/dev/null','w') as sink, contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
        output = compute(payload)
    print(json.dumps(output, allow_nan=False))
