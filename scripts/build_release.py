#!/usr/bin/env python3
"""Build an allowlisted source/model release, never a copy of the working directory."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[1]
FILES=['TESTED_ENVIRONMENT.txt','run.py','run.sh','requirements.txt','requirements-propermab.txt','README.md','START_HERE.md','MODULE_GUIDE.md','THIRD_PARTY_NOTICES.md','Public-hydrophobic-summary.pdf','PTM_EXPOSURE_EVIDENCE.md']
# 'models/antiberty' was listed here and never existed: AntiBERTy ships its
# own checkpoint inside the installed package and this tree is not populated.
TREES=['app','data','models','static','vendor/propermab','tests','scripts']
BLOCKED={'__pycache__','.git','.DS_Store','.venv','runs'}

def normalize(value):
    if isinstance(value,dict):return {k:('<resolved at runtime>' if k=='checkpoint_directory' else normalize(v)) for k,v in value.items()}
    if isinstance(value,list):return [normalize(v) for v in value]
    return value

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nanoshaper-source',type=Path,required=True)
    p.add_argument('--cgal-source',type=Path,required=True)
    a=p.parse_args();payload={}
    for name in FILES:
        path=ROOT/name;payload[name]=path.read_bytes()
    for name in TREES:
        for path in sorted((ROOT/name).rglob('*')):
            rel=path.relative_to(ROOT)
            if not path.is_file() or path.is_symlink() or any(x in BLOCKED for x in rel.parts) or path.suffix in ('.pyc','.pyo'):continue
            if str(rel)=='vendor/propermab/tools/NanoShaper':continue
            payload[str(rel)]=path.read_bytes()
    panel='data/profiler_comparators.json'
    payload[panel]=(json.dumps(normalize(json.loads(payload[panel])),indent=2)+'\n').encode()
    provenance=[]
    for name,source in [('nanoshaper',a.nanoshaper_source),('cgal42',a.cgal_source)]:
        # Git's tracked-file inventory excludes all machine-local build outputs.
        files=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0')
        revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD']).decode().strip()
        diff=subprocess.check_output(['git','-C',str(source),'diff','--no-ext-diff','--binary','HEAD'])
        if name=='nanoshaper':payload['native-source/nanoshaper/CMakeLists.txt']=(source/'CMakeLists.txt').read_bytes()
        provenance.append(f'{name}: upstream {revision}; working-source modifications in {name}.patch')
        payload[f'native-source/{name}.patch']=diff
        for filename in files:
            path=source/filename
            if not filename or not path.is_file() or path.is_symlink():continue
            # NanoShaper repo contains binary artifacts/archives; only compile inputs + docs.
            if name=='nanoshaper' and not (filename.startswith(('src/','CGALPatch/')) or filename in ('CMakeLists.txt','license.txt','README.md','NanoShaper0.7_Guide.pdf')):continue
            payload[f'native-source/{name}/{filename}']=path.read_bytes()
    payload['native-source/PROVENANCE.md']=('# Native source provenance\n\n'+'\n'.join(provenance)+'\n\nPatched source is already applied. Do not apply patches twice. Build with scripts/build_native.sh. CGAL changes cover modern CMake/Boost compatibility, ARM endian detection, and the NanoShaper weighted-point patch. This source build has not been validated on every recipient platform.\n').encode()
    # Strip the developer-only build log path from the bundled provenance note.
    payload['vendor/propermab/tools/BUILD.txt']=b'See native-source/PROVENANCE.md and scripts/build_native.sh. The original Apple ARM64 executable is intentionally not included.\n'
    manifest={k:hashlib.sha256(v).hexdigest() for k,v in sorted(payload.items())}
    payload['MANIFEST.sha256']=''.join(f'{digest}  {name}\n' for name,digest in manifest.items()).encode()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name,content in sorted(payload.items()):
            info=zipfile.ZipInfo('antibody-platform/'+name)
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(0o100755 if name.endswith('.sh') else 0o100644)<<16
            archive.writestr(info,content)
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    a.output.with_suffix('.zip.sha256').write_text(f'{digest}  {a.output.name}\n')
    print(f'{a.output.name}: {len(payload)} files; {a.output.stat().st_size/1_000_000:.1f} MB; SHA256 {digest}')
if __name__=='__main__':main()
