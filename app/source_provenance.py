"""Versioned local source records; never sends analysis inputs to external services."""
from __future__ import annotations
from functools import lru_cache
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
import hashlib
import importlib.util
from .method_sources import METHOD_SOURCES

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = {'antiberty':'antiberty','immunebuilder':'ImmuneBuilder','anarcii':'anarcii',
            'freesasa':'freesasa','biopython':'biopython','reportlab':'reportlab'}

def package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return 'Not installed'

@lru_cache(maxsize=32)
def _hash_file(path, size, mtime):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def file_record(path):
    path = Path(path)
    if not path.is_file():
        return {'file':path.name, 'status':'Not installed'}
    stat = path.stat()
    return {'file':path.name,'status':'Available locally','sha256':_hash_file(str(path), stat.st_size, stat.st_mtime_ns)}

def _antiberty_artifacts():
    """Where AntiBERTy's checkpoint actually resolved, or why it did not.

    get_weights() is the package's own resolution order and the same one the
    Runner uses, so this reports the checkpoint that an analysis would load.
    download=False keeps a provenance read from reaching the network.
    """
    try:
        from antiberty import get_weights
        base = Path(get_weights(download=False))
    except Exception as exc:
        return [{'file': 'AntiBERTy checkpoint', 'status': 'Not installed',
                 'detail': str(exc)[:200]}]
    return [file_record(base / 'config.json'), file_record(base / 'model.safetensors')]


def source_records():
    sources = [dict(item) for item in METHOD_SOURCES]
    for source in sources:
        key=source['id']
        source['version'] = package_version(PACKAGES[key]) if key in PACKAGES else 'See pinned source/reference record'
        if key=='antiberty':
            source['implementation']='antiberty.AntiBERTyRunner; local classification and pseudo_log_likelihood'
            source['artifacts']=_antiberty_artifacts()
        elif key=='immunebuilder':
            source['implementation']='ImmuneBuilder ABodyBuilder2 for paired Fv / NanoBodyBuilder2 for heavy-only; model_ids=[1]'
            spec=importlib.util.find_spec('ImmuneBuilder')
            base=Path(spec.origin).parent/'trained_model' if spec and spec.origin else ROOT/'models/missing'
            source['artifacts']=[file_record(base/'antibody_model_1'),file_record(base/'nanobody_model_1')]
        elif key=='freesasa':
            source['implementation']='freesasa.calc; default 1.4 Å probe; residueAreas relativeSideChain'
        elif key=='propermab':
            source['version']='c97e41b95cccc103143d3123db953e694310c76e'
            source['implementation']='Vendored official descriptor code; see per-run completed and unavailable descriptor groups'
        elif key=='ogrdb':
            from .humanness import REFERENCE_PATH
            source['artifacts']=[file_record(REFERENCE_PATH)]
            source['version']='IGH release 10; IGK release 5; IGL release 4'
        elif key=='therapeutic_precedents':
            source['artifacts']=[file_record(ROOT/'data/clinical_frameworks.json'),file_record(ROOT/'data/therapeutic_ptm_evidence.json'),file_record(ROOT/'data/therapeutic_ptm_sequences.json')]
            source['version']='Clinical stage snapshot February 2025; local evidence expansion 2026-10-04'
        elif key=='ngl':
            source['artifacts']=[file_record(ROOT/'static/vendor/ngl/ngl.js')]
    return sources
