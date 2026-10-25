#!/usr/bin/env python
"""Install the optional components this project deliberately does not ship.

One optional analysis path depends on PROPERMAB, which Regeneron licenses for
academic research use only. This project distributes none of it. If you are
eligible, this script shows you Regeneron's licence, asks you to affirm
eligibility, and then fetches it from Regeneron's own repository.

    python setup/fetch_optional.py --propermab
    python setup/fetch_optional.py --status
    python setup/fetch_optional.py --native-help

If you work at a commercial entity, do not run --propermab. Contact Regeneron
for terms. Everything else in this project is permissively licensed and needs
no permission; see USE_TIERS.md for what you still get without it.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENDOR = ROOT / "vendor" / "propermab"

PROPERMAB_REPO = "https://github.com/regeneron-mpds/propermab"
PROPERMAB_COMMIT = "c97e41b95cccc103143d3123db953e694310c76e"

ELIGIBILITY = """
================================================================================
PROPERMAB is NOT part of this project. It belongs to Regeneron Pharmaceuticals,
Inc. and is licensed on these terms, quoted from its LICENSE.md:

    "You" refers to an academic institution or academically employed full-time
    personnel only.

    Regeneron hereby grants You a right to use, reproduce, modify, or distribute
    the PROPERMAB source code ... for academic research purposes only.

    Prohibited Uses: The rights granted herein do not include any right to use
    by commercial entities or commercial use of any kind, including, without
    limitation, (1) any integration into other code or software that is used for
    further commercialization, (2) any reproduction, copy, modification or
    creation of a derivative work that is then incorporated into a commercial
    product or service or otherwise used for any commercial purpose,
    (3) distribution of the source code, in whole or in part, or any resulting
    executables, in any commercial product, or (4) use of the source code, in
    whole or in part, or any resulting executables, in any commercial online
    service.

Two things follow:

  * Both conditions must hold. You must be an academic institution or
    academically employed full-time personnel, AND the use must be for academic
    research purposes. Nonprofit intent alone does not establish eligibility.

  * The authors of THIS project make no representation about your eligibility
    and grant you no rights in Regeneron's code. If you proceed, your licence is
    with Regeneron, not with us. Read the full LICENSE.md after fetching, and
    direct any question about scope to Regeneron.

If you are at a company, stop here. USE_TIERS.md lists what this project does
without PROPERMAB, and the free Therapeutic Antibody Profiler is a reasonable
substitute for charge-patch context.
================================================================================
"""

NATIVE_HELP = """
Native dependencies for the PROPERMAB surface-patch path
========================================================

This project builds and distributes none of these. They are needed only by
PROPERMAB, so install them only if you have already established eligibility.

  NanoShaper   molecular surface construction      GPLv3
  CGAL         geometry library NanoShaper needs   GPLv3+ / LGPL, varies by package
  APBS         Poisson-Boltzmann electrostatics    see upstream
  PDB2PQR      protonation and charge assignment   see upstream

NanoShaper and CGAL are GPL-licensed. Building them for your own use imposes no
obligation on you beyond the GPL's own terms, but if you distribute a binary you
built, you must provide corresponding source. That is why this project ships
neither the source nor a binary: your GPL obligations stay yours, on your
machine.

Install APBS and PDB2PQR from conda-forge:

    conda install -c conda-forge apbs pdb2pqr

Build NanoShaper from its own repository, following its instructions:

    https://gitlab.iit.it/SDecherchi/nanoshaper

On Linux, CGAL headers and the GMP/MPFR/Boost dependencies come from
conda-forge: libboost-headers, libboost-devel, gmp, mpfr. LINUX_AND_HEADLESS.md
records a verified Linux build, including the one portability fix this project
applies at runtime: PROPERMAB's APBS wrapper sets DYLD_LIBRARY_PATH, which is
macOS-only, and app/propermab_worker.py sets LD_LIBRARY_PATH instead on Linux.

Once NanoShaper is on PATH, verify with:

    python setup/fetch_optional.py --status
"""


def confirm_eligibility(assume_yes: bool) -> bool:
    print(ELIGIBILITY)
    if assume_yes:
        print("--i-am-eligible given; proceeding on your affirmation.\n")
        return True
    if not sys.stdin.isatty():
        print("Not an interactive terminal and --i-am-eligible was not given. "
              "Refusing to fetch.\n")
        return False
    print("Type exactly: I am eligible under Regeneron's academic licence")
    answer = input("> ").strip().rstrip(".")
    if answer.lower() != "i am eligible under regeneron's academic licence".lower():
        print("\nNot affirmed. Nothing fetched. This is the right outcome if you "
              "are unsure.\n")
        return False
    print()
    return True


def fetch_propermab(assume_yes: bool) -> int:
    if VENDOR.is_dir() and any(VENDOR.iterdir()):
        print(f"PROPERMAB already present at {VENDOR.relative_to(ROOT)}. "
              f"Remove it to re-fetch.")
        return 0
    if not confirm_eligibility(assume_yes):
        return 1
    if not shutil.which("git"):
        print("git is required to fetch from Regeneron's repository.")
        return 1

    VENDOR.parent.mkdir(parents=True, exist_ok=True)
    print(f"Cloning {PROPERMAB_REPO}")
    commands = [
        ["git", "clone", "--quiet", PROPERMAB_REPO, str(VENDOR)],
        ["git", "-C", str(VENDOR), "checkout", "--quiet", PROPERMAB_COMMIT],
    ]
    for command in commands:
        result = subprocess.run(command)
        if result.returncode:
            print(f"\nFailed: {' '.join(command)}")
            print("Fetch it manually from Regeneron and place it at "
                  f"{VENDOR.relative_to(ROOT)}.")
            return result.returncode

    licence = VENDOR / "LICENSE.md"
    print(f"\nFetched at commit {PROPERMAB_COMMIT[:12]}.")
    print(f"Regeneron's full licence is at {licence.relative_to(ROOT)}. Read it.")
    print("\nThis project applies two runtime adaptations without changing any "
          "feature equation:")
    print("  * app/propermab_worker.py sets LD_LIBRARY_PATH on Linux, since "
          "PROPERMAB's APBS wrapper sets the macOS-only DYLD_LIBRARY_PATH.")
    print("  * the NanoShaper subprocess is invoked via sys.executable rather "
          "than an unqualified 'python'.")
    print("\nNext: python setup/fetch_optional.py --native-help")
    return 0


def status() -> int:
    print("Optional components\n")
    present = VENDOR.is_dir() and any(VENDOR.iterdir())
    print(f"  [{'present' if present else 'ABSENT '}] PROPERMAB  "
          f"{'(Regeneron, academic research use only)' if present else '- not installed; this is the default'}")

    for binary, note in (("NanoShaper", "GPLv3"), ("apbs", "APBS"), ("pdb2pqr30", "PDB2PQR")):
        found = shutil.which(binary)
        print(f"  [{'present' if found else 'ABSENT '}] {binary:<10} "
              f"{found or f'- not on PATH ({note})'}")

    if present:
        try:
            sys.path.insert(0, str(ROOT))
            from app.propermab import availability
            print(f"\n  adapter reports: {availability()}")
        except Exception as exc:
            print(f"\n  adapter check failed: {type(exc).__name__}: {exc}")
    else:
        print("\n  The permissively-licensed core needs none of the above. "
              "See USE_TIERS.md.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--propermab", action="store_true",
                        help="fetch PROPERMAB from Regeneron (academic eligibility required)")
    parser.add_argument("--i-am-eligible", action="store_true",
                        help="affirm academic eligibility non-interactively")
    parser.add_argument("--native-help", action="store_true",
                        help="show native dependency build guidance")
    parser.add_argument("--status", action="store_true",
                        help="report what is installed")
    args = parser.parse_args()

    if args.status:
        return status()
    if args.native_help:
        print(NATIVE_HELP)
        return 0
    if args.propermab:
        return fetch_propermab(args.i_am_eligible)
    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
