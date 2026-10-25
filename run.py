#!/usr/bin/env python3
"""Set up and launch the local Antibody Atlas server."""

from __future__ import annotations

import argparse
import os
import hashlib
import subprocess
import sys
import venv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv"
PYTHON = VENV / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
REQUIREMENTS = ROOT / "requirements.txt"
STAMP = VENV / "requirements.sha256"
REQUIRED_FILES = (
    ROOT / "app/main.py",
    ROOT / "static/index.html",
    ROOT / "static/start.html",
    ROOT / "static/vendor/ngl/ngl.js",
    ROOT / "data/human_v_ogrdb_numbered.json",
    ROOT / "data/therapeutic_ptm_evidence.json",
    REQUIREMENTS,
)
# Only files this package actually ships, and cannot run without, belong
# above. Two other classes do not, and listing either one turns a normal
# state into "the package is incomplete" and refuses to start the server:
#
#   Engine checkpoints. AntiBERTy carries its own inside the installed pip
#   package; ImmuneBuilder and DeepSP are fetched after install. `--doctor`
#   reports which engines are available.
#
#   The Thera-SAbDab-derived datasets in app.reference_data.OPTIONAL_DATASETS.
#   This project deliberately does not redistribute them -- the upstream
#   source states no licence grant -- so they are absent from every published
#   package by design, and are rebuilt locally by setup/build_reference_data.py.
#   ReferenceDataMissing already explains that at the point of use, with the
#   command to run; refusing to boot pre-empts a message that was written to
#   handle exactly this.
#
# The invariant is checked by tests/test_report_clarity.py::PackagingTests.



def ensure_linux_cpu_torch(python: Path, platform_name: str | None = None) -> None:
    """Choose the official CPU wheel on a fresh Linux install; preserve existing torch."""
    if (platform_name or sys.platform) != 'linux':
        return
    present = subprocess.run([str(python), '-c', "import importlib.util;raise SystemExit(importlib.util.find_spec('torch') is None)"], capture_output=True)
    if present.returncode == 0:
        return
    print('Installing CPU PyTorch from the official index (no CUDA runtime download).', flush=True)
    subprocess.run([str(python), '-m', 'pip', 'install', 'torch', '--index-url', 'https://download.pytorch.org/whl/cpu'], cwd=ROOT, check=True)


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Check package contents only')
    parser.add_argument('--setup-only',action='store_true',help='Install without starting the server')
    parser.add_argument('--with-propermab',action='store_true',help='Install optional PROPERMAB Python dependencies')
    parser.add_argument('--doctor',action='store_true',help='Check installed engines without installing or downloading')
    parser.add_argument('--stability',action='store_true',help='Install or verify the pinned MIT stability model locally')
    parser.add_argument('--references',action='store_true',help='Explicitly build clinical reference data from its upstream source')
    parser.add_argument('--warmup',action='store_true',help='Download/load structure model 1 checkpoints without input sequences')
    args=parser.parse_args()
    if sys.version_info[:2] not in {(3, 11), (3, 12)}:
        print("Python 3.11 or 3.12 is required.", file=sys.stderr)
        return 2
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED_FILES if not path.is_file()]
    if missing:
        print("The package is incomplete: " + ", ".join(missing), file=sys.stderr)
        return 2
    if args.check:
        print("Package files and Python version look ready. Installation was not started.")
        return 0
    if args.doctor:
        if not PYTHON.is_file():
            print('Run python3.12 run.py --setup-only first.'); return 2
        return subprocess.call([str(PYTHON),str(ROOT/'scripts/doctor.py')],cwd=ROOT)
    if not args.setup_only:
        from scripts.launch import port_number, reuse_or_error
        try:
            port = port_number(os.environ.get("PORT", "8000"))
        except ValueError as exc:
            print(str(exc), file=sys.stderr); return 2
        existing = reuse_or_error(port)
        if existing is not None: return existing
    if PYTHON.exists():
        checked = subprocess.run([str(PYTHON), '-c', 'import sys;raise SystemExit(sys.version_info[:2] not in ((3,11),(3,12)))'])
        if checked.returncode:
            print('Existing .venv is incompatible; rename it and rerun. No files were removed.', file=sys.stderr); return 2
    elif VENV.exists():
        print('Existing .venv is incomplete; rename it and rerun. No files were removed.', file=sys.stderr); return 2
    if not PYTHON.is_file():
        print("Creating a private Python environment in .venv …", flush=True)
        venv.EnvBuilder(with_pip=True).create(VENV)
    files=[REQUIREMENTS]+([ROOT/'requirements-propermab.txt'] if args.with_propermab else [])
    stamp=VENV/('requirements-propermab.sha256' if args.with_propermab else 'requirements.sha256')
    requirements_hash = hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
    if not stamp.is_file() or stamp.read_text().strip() != requirements_hash:
        if not shutil.which('git'):
            print('Git is required for pinned upstream dependencies. macOS: xcode-select --install; Ubuntu: sudo apt install git', file=sys.stderr); return 2
        ensure_linux_cpu_torch(PYTHON)
        print("Installing Python dependencies. This first setup needs internet access …", flush=True)
        subprocess.run([str(PYTHON), "-m", "pip", "install"]+[item for p in files for item in ("-r",str(p))],
                       cwd=ROOT, check=True)
        stamp.write_text(requirements_hash + "\n")
    if args.stability:
        result = subprocess.run([str(PYTHON), str(ROOT/'scripts/setup_stability.py')], cwd=ROOT)
        if result.returncode:
            print('Stability model setup failed; other features remain available. Retry .venv/bin/python scripts/setup_stability.py when online.', file=sys.stderr)
            if args.setup_only: return result.returncode
    if args.references:
        subprocess.run([str(PYTHON),str(ROOT/'setup/build_reference_data.py'),'--all'], cwd=ROOT, check=True)
    if args.warmup:
        subprocess.run([str(PYTHON),'-c',"from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2; ABodyBuilder2(model_ids=[1]); NanoBodyBuilder2(model_ids=[1]); print('Structure model 1 checkpoints loaded. No input sequence was used.')"],cwd=ROOT,check=True)
    if args.setup_only: return 0
    if not shutil.which('hmmscan'):
        print('HMMER is missing; structure prediction needs it. macOS: brew install hmmer; Ubuntu: sudo apt install hmmer', flush=True)
    from scripts.launch import start
    return start(PYTHON, port)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"Setup command failed (exit {exc.returncode}). Resolve the error above, then rerun; completed downloads are kept. macOS may need xcode-select --install; Ubuntu may need build-essential and python3-dev.", file=sys.stderr)
        raise SystemExit(exc.returncode)
