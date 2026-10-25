#!/usr/bin/env python3
"""Build checked source archives for a manual GitHub Release; never publishes."""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', required=True, help='Version such as 0.1.0 or 0.1.0-rc1')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9.]+)?', args.version):
        parser.error('Use a version such as 0.1.0 or 0.1.0-rc1')
    args.out.mkdir(parents=True, exist_ok=True)
    base = f'antibody-atlas-{args.version}'
    targets = [args.out / f'{base}.zip', args.out / f'{base}.tar.gz', args.out / f'{base}-SHA256SUMS.txt']
    if any(p.exists() for p in targets):
        parser.error('Release files already exist; use a new output directory or version.')
    with tempfile.TemporaryDirectory(prefix='atlas-release-') as scratch:
        tree = Path(scratch)/base
        subprocess.run([sys.executable, str(ROOT/'tools/make_oss_export.py'), '--out', str(tree), '--panel-descriptors', 'permissive-only'], check=True)
        for name in ('run_local.sh', 'run.sh', 'Start Antibody Atlas.command'):
            (tree/name).chmod(0o755)
        subprocess.run([sys.executable, str(tree/'run.py'), '--check'], cwd=tree, check=True)
        for fmt in ('zip', 'gztar'):
            shutil.make_archive(str(args.out.resolve()/base), fmt, root_dir=scratch, base_dir=base)
    targets[2].write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n' for path in targets[:2]))
    print('\nRelease files prepared. Nothing was uploaded:')
    for path in targets:
        print(path.resolve())


if __name__ == '__main__':
    main()
