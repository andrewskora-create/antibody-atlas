# Running it: macOS, Linux, Windows

The full workflow targets macOS and Linux. Windows users should use WSL2.
The table describes intended engine availability, not proof of a clean installation on every platform.
The rc3 scientific workflow was exercised after a fresh dependency installation on Apple Silicon macOS with Python 3.12 and Chrome. Automated tests cover Linux Python 3.11/3.12 and macOS/Linux startup. End-to-end operation on Intel Macs and an actual Windows/WSL2 computer has not been verified. The rc4 update adds onboarding documentation and a guide link; it does not expand platform validation.
Native Windows is not a supported full-workflow release.

| | macOS | Linux | Windows (WSL2) | Windows (native) |
|---|---|---|---|---|
| Sequence screen, PTM detection | yes | yes | yes | yes |
| Numbering (IMGT / Kabat) | yes | yes | yes | yes |
| Charge, pI, hydrophobicity, aggregation proxies | yes | yes | yes | yes |
| Human framework identity | yes | yes | yes | yes |
| AntiBERTy scoring | yes | yes | yes | yes |
| Clinical comparator search | yes | yes | yes | yes |
| PTM fix proposals and ranking | yes | yes | yes | yes |
| **Structure prediction** | yes | yes | yes | **no** |
| **Exposure, surface patches** | yes | yes | yes | **no** |
| **Paratope prediction** | yes | yes | yes | **no** |
| PROPERMAB descriptors (academic only) | yes | yes | yes | no |

## Use Python 3.11 or 3.12

The launcher accepts Python 3.11 and 3.12. Use the same versions for installation and testing.

The Python supplied with your system or developer tools may be too old. The launcher searches for a usable one and refuses with
instructions rather than building a virtual environment that cannot install.
If you only have the system Python:

```bash
brew install python@3.12
```

or, with conda or micromamba already installed:

```bash
conda create -y -n atlas python=3.12 && conda activate atlas
PYTHON=python ./run_local.sh
```

Set `PYTHON` to pick a specific interpreter. A virtual environment that was
built with an unsupported Python is detected and left untouched. Rename `.venv` before recreating it.

## macOS and Linux

```bash
# In the extracted release folder:
bash run_local.sh
```

It creates the virtual environment, installs dependencies and the MIT stability
predictor, starts the server on <http://127.0.0.1:8000> and opens your
browser. Re-running it is cheap; it skips what is already done.

For the structure-dependent features, install HMMER first:

```bash
brew install hmmer          # macOS
sudo apt install hmmer      # Debian/Ubuntu
```

The script tells you if it is missing and runs without it.

## Windows

**WSL2 is the intended route for the complete workflow**, but has not been verified on an actual Windows machine. Inside a WSL2
Ubuntu shell:

```bash
sudo apt update
sudo apt install hmmer python3-venv python3-dev build-essential git
# In the extracted release folder, with Python 3.11 or 3.12:
bash run_local.sh
```

Then open <http://127.0.0.1:8000> in Windows — WSL2 forwards localhost, so the
browser on the Windows side reaches the server without extra configuration.

**Native Windows** has an older launcher for sequence-level use. It has not been validated as part of this release; use WSL2 for the supported workflow:

```powershell
.\run_local.ps1
```

### Why structure prediction does not run on native Windows

ImmuneBuilder validates and numbers its input through *legacy ANARCI*
(`ImmuneBuilder/sequence_checks.py` imports `anarci`), and legacy ANARCI shells
out to `hmmscan` from HMMER. HMMER is distributed for Unix-like systems and has
no supported native Windows build. Without `hmmscan` the structure step fails,
and everything downstream of coordinates — relative side-chain exposure,
surface patches, the structure figure, and paratope prediction — goes with it.

This is not a limitation of this project and cannot be worked around by
installing a Python package. The app detects the absence, says so, and returns
the sequence-level results rather than failing the run.

Note that the project's *own* numbering is ANARCII, a language model with no
HMMER dependency, which is why numbering, CDR boundaries and the comparator
search all work natively on Windows. Only the structure branch is affected.

## Checking an installation

```bash
python scripts/doctor.py
```

reports every engine separately with a platform-specific fix, rather than
failing later at analysis time.

## Notes

- **AntiBERTy weights** are resolved automatically: the bundled copy if
  present, otherwise the `antiberty` package's own, which it downloads on
  first use. Nothing to install by hand.
- **Reference data** (clinical precedent and the comparison panel) is built on
  your machine from Oxford's CSV because it is not included in the public archive; review upstream terms in DATA_PROVENANCE.md. Run
  `.venv/bin/python run.py --setup-only --references` explicitly to build it.
  It is not downloaded by the default launcher.
- **Bind to localhost only.** There is no authentication; this is a
  single-user local tool.
- **Apple Silicon** works natively; no Rosetta needed. PyTorch installs an
  arm64 wheel.
