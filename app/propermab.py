"""Local, explicitly invoked PROPERMAB integration for academic research."""
from __future__ import annotations
import importlib.util
from importlib.metadata import version, PackageNotFoundError
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

def availability():
    missing = [name for name in ('numpy','scipy','anarci','freesasa','openmm') if importlib.util.find_spec(name) is None]
    if not (ROOT/'vendor/propermab/propermab/features/seq_featurizer.py').exists():
        missing.append('vendored PROPERMAB source')
    return {'available': not missing, 'missing': missing,
            'scope': 'Official PROPERMAB Fv descriptors', 'full_patch_pipeline': patch_available(),
            'source_commit': 'c97e41b95cccc103143d3123db953e694310c76e'}


def patch_available():
    configured = os.environ.get('PROPERMAB_NANOSHAPER')
    native = Path(configured) if configured else ROOT/'vendor/propermab/tools/NanoShaper'
    compatible = bool(configured) or (sys.platform == 'darwin' and platform.machine() == 'arm64')
    return compatible and native.is_file() and os.access(native, os.X_OK) and all(
        importlib.util.find_spec(name) is not None for name in ('apbs_binary','pdb2pqr','plyfile','meshio','matplotlib','sklearn'))


def comparison_protocol(ph=7.4, structure=False):
    """Capture installed runtime for the app's model-1 structure pipeline."""
    versions = {}
    for name in ('antiberty', 'ImmuneBuilder', 'anarcii', 'freesasa', 'openmm',
                 'pdb2pqr', 'apbs-binary', 'numpy', 'scipy', 'numba'):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    protocol = {'ph': ph, 'software_versions': versions}
    if os.environ.get('PROPERMAB_NANOSHAPER'):
        native=Path(os.environ['PROPERMAB_NANOSHAPER'])
        protocol['nanoshaper_override_sha256']=hashlib.sha256(native.read_bytes()).hexdigest() if native.is_file() else 'unavailable'
    if structure:
        spec = importlib.util.find_spec('ImmuneBuilder')
        checkpoint = Path(spec.origin).parent/'trained_model/antibody_model_1' if spec and spec.origin else None
        if checkpoint and checkpoint.is_file():
            protocol.update(structure_model='ImmuneBuilder ABodyBuilder2 model 1',
                            structure_checkpoint_artifacts=[{'filename': checkpoint.name,
                               'sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest()}],
                            preparation='Atlas full-input Fv preparation v2: exact residue identity; terminal reconstruction when needed; CHARMM hydrogens.')
    return protocol


def _run_single(record, pdb_path=None):
    chains = record['chains']
    if set(chains) != {'H','L'}:
        raise ValueError('PROPERMAB requires paired heavy and light variable domains.')
    ph = record['ph']
    if not 4 <= ph <= 10:
        raise ValueError('PROPERMAB supports pH 4–10.')
    data = {'heavy': chains['H']['sequence'], 'light': chains['L']['sequence'], 'ph': ph,
            'light_type': 'kappa' if chains['L']['numbering']['chain_type']=='K' else 'lambda'}
    data['patch_enabled'] = patch_available()
    data['nanoshaper_path'] = str(Path(os.environ.get('PROPERMAB_NANOSHAPER', ROOT/'vendor/propermab/tools/NanoShaper')).resolve())
    preparation_error = None
    coverage = {}
    with tempfile.TemporaryDirectory(prefix='propermab-') as temp:
        if pdb_path:
            try:
                data['pdb_path'] = str(_prepare_imgt_pdb(record, pdb_path, Path(temp), coverage))
            except Exception as exc:
                preparation_error = str(exc) if isinstance(exc, ValueError) else 'Structure preparation unavailable: ' + type(exc).__name__
        data['work_dir'] = temp
        env = dict(os.environ, MPLCONFIGDIR=temp, HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', PYTHONHASHSEED='0')
        completed = subprocess.run([sys.executable,'-m','app.propermab_worker'], input=json.dumps(data),
                    capture_output=True, text=True, cwd=ROOT, env=env, timeout=600)
        if completed.returncode:
            raise RuntimeError('Local PROPERMAB calculation failed; check local dependencies and structure completeness.')
        result = json.loads(completed.stdout)
        result['comparison_protocol'] = comparison_protocol(ph, structure=bool(data.get('pdb_path')))
        result['input_pdb_sha256'] = hashlib.sha256(Path(pdb_path).read_bytes()).hexdigest() if pdb_path else None
        result['prepared_pdb_sha256'] = hashlib.sha256(Path(data['pdb_path']).read_bytes()).hexdigest() if data.get('pdb_path') else None
        if data.get('pdb_path'):
            result['structure_coverage'] = coverage
            result['structure_partial'] = any(x['omitted_positions'] for x in coverage.values())
            if result['structure_partial']:
                result['comparison_protocol']['partial_structure'] = True
                omitted = '; '.join(f"{label}: {info['model_residues']}/{info['input_residues']} residues; omitted input positions {', '.join(map(str, info['omitted_positions'])) or 'none'}" for label, info in coverage.items())
                result['structure_note'] = ('Partial structure: ' + omitted + '. Every modeled residue matches the input exactly. Structural and surface descriptors describe the modeled fragment, with artificial chain termini; missing residues may change charge and surface patches. Full-input sequence descriptors are retained. Structural comparator calls are withheld.')
        if preparation_error:
            result['unavailable'] = [x for x in result['unavailable'] if x['feature'] not in {'Structure descriptors', 'APBS/NanoShaper surface patches'}]
            result['unavailable'].append({'feature':'APBS/NanoShaper surface patches','reason':'Not run because structure preparation failed: ' + preparation_error})
            result['unavailable'].append({'feature':'Structure descriptors','reason':preparation_error})
        return result


def run(record, pdb_path=None, repair_output_path=None):
    result=_run_single(record,pdb_path)
    if not result.get('structure_partial'):
        return result
    comparison={'status':'failed','omitted_positions':{c:v['omitted_positions'] for c,v in result['structure_coverage'].items()},
                'note':'The complete repaired model is selected for structural assessment when identity, geometry and exposure checks pass. The original-fragment measurements are retained below. Reconstruction remains a model hypothesis, not experimental evidence or automatic risk clearance.'}
    try:
        with tempfile.TemporaryDirectory(prefix='atlas-terminal-repair-') as temp:
            target=Path(temp)/'repaired.pdb'
            payload={'input_path':str(pdb_path),'output_path':str(target),'sequences':{k:v['sequence'] for k,v in record['chains'].items()},'ph':record['ph']}
            process=subprocess.run([sys.executable,'-m','app.terminal_repair'],input=json.dumps(payload),capture_output=True,text=True,cwd=ROOT,timeout=180)
            if process.returncode or not target.exists():raise ValueError('Terminal reconstruction or peptide geometry validation failed.')
            metadata=json.loads(process.stdout)
            repaired=_run_single(record,target)
            if not repaired.get('prepared_pdb_sha256'):raise ValueError('Repaired model could not be prepared for PROPERMAB.')
            comparison.update(status='complete',**metadata,repaired_pdb_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
            selected={'hyd_patch_area','hyd_patch_area_cdr','pos_patch_area','pos_patch_area_cdr','neg_patch_area','neg_patch_area_cdr','hyd_asa','heiden_score'}
            left={f['key']:f for f in result['features']};right={f['key']:f for f in repaired['features']}
            comparison['metrics']=[]
            for key in sorted(selected):
                original=left.get(key,{}).get('value');new=right.get(key,{}).get('value');feature=left.get(key) or right.get(key) or {'label':key,'unit':''}
                comparison['metrics'].append(dict(key=key,label=feature['label'],unit=feature['unit'],original=original,repaired=new,delta=round(new-original,4) if original is not None and new is not None else None))
            if repair_output_path:
                import shutil
                shutil.copyfile(target,repair_output_path)
            comparison['unavailable']=repaired.get('unavailable',[])
            comparison['note']+=' Changes below are sensitivity results, not calibrated pass/fail thresholds. Comparisons use controls processed with the same full-input preparation policy.'
            from .structure_analysis import analyze_structure
            selected_analysis=analyze_structure(target,payload['sequences'],record['ph'],record.get('scheme','imgt'))
            if any(value<100 for value in selected_analysis['mapping_coverage_percent'].values()):
                raise ValueError('Reconstructed model did not map completely during exposure analysis.')
            repaired['original_input_pdb_sha256']=result['input_pdb_sha256']
            repaired['original_structure_coverage']=result['structure_coverage']
            repaired['selected_structure_analysis']=selected_analysis
            repaired['structure_selection']='terminal_repaired'
            repaired['structure_note']='Full-input terminal-repaired model selected for PROPERMAB, surface assessment and matched-protocol comparisons. Rebuilt terminal coordinates remain uncertain; original-fragment measurements are retained as sensitivity evidence.'
            result=repaired
    except Exception as exc:
        comparison['status']='failed'
        comparison['note']='Repair selection failed; original-model results retained. '+comparison['note']
        comparison['error']=str(exc) if isinstance(exc,ValueError) else type(exc).__name__+': local repair comparison unavailable.'
    result['terminal_repair']=comparison
    return result


def _terminal_identity(expected, observed, label):
    """Accept exact whole input or a unique exact fragment missing <=2 terminal residues.

    This is an adapter policy for reporting partial models, not a validated
    adequacy threshold for their physical descriptors. No mutations/internal gaps.
    """
    if expected == observed:
        offset = 0
    else:
        offset = expected.find(observed) if observed else -1
        omitted = len(expected) - len(observed)
        if (offset < 0 or omitted not in (1, 2) or len(observed) / max(1, len(expected)) < .98
                or expected.find(observed, offset + 1) >= 0):
            raise ValueError(f'PROPERMAB {label} structure does not exactly match the input or a uniquely matching terminal-trimmed fragment (maximum 2 omitted residues, at least 98% coverage). Substitutions and internal gaps are not accepted.')
    positions = list(range(offset + 1, offset + len(observed) + 1))
    missing = list(range(1, offset + 1)) + list(range(offset + len(observed) + 1, len(expected) + 1))
    return dict(input_residues=len(expected), model_residues=len(observed),
                coverage_percent=round(100 * len(observed) / len(expected), 2),
                omitted_positions=missing, input_positions=positions)


def _prepare_imgt_pdb(record, pdb_path, directory, coverage=None):
    """Check identity then explicitly renumber each PDB residue to IMGT."""
    from Bio.PDB import PDBParser, PDBIO
    from Bio.SeqUtils import seq1
    from .numbering import number_chain
    structure = PDBParser(QUIET=True).get_structure('fv', str(pdb_path))
    if {chain.id for chain in structure[0]} != {'H','L'}:
        raise ValueError('PROPERMAB structure must contain both H and L chains.')
    for chain in structure[0]:
        if chain.id not in {'H','L'}:
            raise ValueError('PROPERMAB structure must contain only H and L chains.')
        residues = list(chain.get_residues())
        seq = ''.join(seq1(r.resname) for r in residues)
        expected = record['chains'][chain.id]['sequence']
        match = _terminal_identity(expected, seq, chain.id)
        if coverage is not None:
            coverage[chain.id] = {k: v for k, v in match.items() if k != 'input_positions'}
        # Keep input IMGT numbering rather than renumbering a shortened fragment.
        entries = {r['sequence_position']: r for r in number_chain(expected, 'imgt')['residues']}
        if entries:
            first,last=min(entries),max(entries)
            import re
            first_num=int(re.search(r'\d+',entries[first]['number'])[0]);last_num=int(re.search(r'\d+',entries[last]['number'])[0])
            for pos in match['input_positions']:
                if pos not in entries and pos<first and first-pos<=2:
                    entries[pos]={'number':str(first_num-(first-pos))}
                elif pos not in entries and pos>last and pos-last<=2:
                    entries[pos]={'number':str(last_num+(pos-last))}
        numbering = [entries[p] for p in match['input_positions'] if p in entries]
        if len(numbering) != len(residues):
            raise ValueError('Incomplete IMGT structure numbering.')
        for index,residue in enumerate(residues):
            residue.id = (' ', -1000-index, ' ')
        for residue, entry in zip(residues,numbering):
            number = str(entry['number'])
            import re
            match = re.fullmatch(r'[HKL]?(\d+)([A-Za-z]?)', number)
            if not match:
                raise ValueError('Unsupported IMGT residue identifier.')
            residue.id = (' ', int(match[1]), match[2] or ' ')
    target = directory/'fv-imgt.pdb'
    writer = PDBIO(); writer.set_structure(structure); writer.save(str(target))
    # Add hydrogens locally for the upstream CHARMM charge descriptors.
    from openmm.app import PDBFile, Modeller, ForceField
    pdb = PDBFile(str(target)); modeller = Modeller(pdb.topology,pdb.positions)
    modeller.addHydrogens(ForceField('charmm36.xml','charmm36/tip3p-pme-b.xml'), pH=record['ph'])
    with target.open('w') as handle:
        PDBFile.writeFile(modeller.topology,modeller.positions,handle,keepIds=True)
    return target
