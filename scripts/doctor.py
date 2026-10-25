#!/usr/bin/env python3
"""Local install diagnostics; no sequence, network request, or model inference."""
import importlib
import importlib.util
import platform
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
print(f'Python {platform.python_version()} | {platform.system()} {platform.machine()}')
failed=[]
for name in ('fastapi','Bio','torch','anarcii','anarci','ImmuneBuilder','openmm','pdbfixer','freesasa','antiberty','reportlab','PIL'):
    try: importlib.import_module(name); print(f'OK  {name}')
    except Exception as exc: failed.append(name);print(f'FAIL {name}: {type(exc).__name__}: {exc}')
HMMER_HELP = {
    'Darwin': 'brew install hmmer',
    'Windows': 'HMMER has no supported native Windows build, so structure prediction, '
               'exposure, surface patches and paratope prediction cannot run here. Use WSL2 '
               '(wsl --install, then: sudo apt install hmmer) and run ./run_local.sh inside it. '
               'Everything sequence-level runs natively without HMMER. See PLATFORMS.md',
    'Linux': 'apt install hmmer  (or dnf install hmmer). If no distribution package is '
             'available, build it: git clone --depth 1 --branch hmmer-3.4 '
             'https://github.com/EddyRivasLab/hmmer && cd hmmer && git clone --depth 1 '
             '--branch easel-0.49 https://github.com/EddyRivasLab/easel && autoconf && '
             '(cd easel && autoconf) && ./configure --prefix=$HOME/local && make -j4 && '
             'make install, then add $HOME/local/bin to PATH',
}
for name in ('hmmscan',):
    found = shutil.which(name)
    if found:
        print(f'OK  {name}: {found}')
    else:
        hint = HMMER_HELP.get(platform.system(), 'install HMMER and put it on PATH')
        print(f'FAIL  {name}: not on PATH. ImmuneBuilder numbers its input through legacy '
              f'ANARCI, which runs this binary, so structure prediction will fail without it. '
              f'Install with: {hint}')
        failed.append(name)
try:
    from app.propermab import availability
    print('PROPERMAB availability:',availability())
except Exception as exc:print('PROPERMAB diagnostic:',type(exc).__name__,str(exc))
# DeepSP is optional and never fails the check: the screen runs without it, and
# it is MIT, so unlike PROPERMAB there is no eligibility question to answer.
try:
    from app import deepsp
    if deepsp.available():
        has_reference = (ROOT/'data'/'deepsp_reference.json').is_file()
        print('OK  DeepSP: weights installed'
              + ('' if has_reference else '; percentile table missing, values will '
                                          'print without a reference position'))
    else:
        print('SETUP  DeepSP: not installed. MIT-licensed, usable by anyone, and it gives '
              'spatial aggregation and charge descriptors without a structure. '
              'Install with: pip install h5py && python setup/fetch_deepsp.py')
except Exception as exc:print('DeepSP diagnostic:',type(exc).__name__,str(exc))
spec=importlib.util.find_spec('ImmuneBuilder')
if spec:
    for name in ('antibody_model_1','nanobody_model_1'):
        ready=(Path(spec.origin).parent/'trained_model'/name).is_file()
        print(f'{"OK" if ready else "SETUP"}  {name}: {"cached" if ready else "run run.py --warmup --setup-only"}')
print('Import checks do not establish prediction or native-patch accuracy. Run the bundled public example next.')
raise SystemExit(1 if failed else 0)
