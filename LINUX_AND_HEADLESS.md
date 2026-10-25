# Linux setup and headless panel runs

This file records a verified Linux installation and the headless runner, as an addition
to `START_HERE.md` (which covers the Apple Silicon path). The web application, the
numbering, humanness, PTM, structure, surface and AntiBERTy engines and the PDF report
all run on Linux. The PROPERMAB native surface-patch stack does not, and is reported as
unavailable rather than substituted or zero-filled.

## Verified environment

| | |
| --- | --- |
| OS | Ubuntu 24.04, x86_64 |
| Python | 3.12.3 |
| Installed from | `requirements.txt` plus a CPU PyTorch wheel |
| Structure prediction | ImmuneBuilder 1.2 / ABodyBuilder2 model 1, CPU, about 25 s per paired Fv |
| AntiBERTy | checkpoint resolved by the installed package; downloads on first use, or set `ANTIBERTY_WEIGHTS_DIR` |
| PROPERMAB | Python descriptors only; native patch pipeline unavailable on this platform |

## Install

```sh
python3.12 -m venv .venv
./.venv/bin/python -m pip install --upgrade pip setuptools wheel
./.venv/bin/python -m pip install torch          # CPU wheel is sufficient
./.venv/bin/python -m pip install -r requirements.txt
```

### HMMER is required for structure prediction

ImmuneBuilder numbers its input through legacy ANARCI, which shells out to the
`hmmscan` binary. Without it, folding fails with an opaque ANARCI error while every
sequence-level readout continues to work. ANARCII, used for the application's own
numbering, does not need it.

```sh
apt install hmmer          # or: dnf install hmmer
```

If no distribution package is available (as on the minimal image used here), build it:

```sh
git clone --depth 1 --branch hmmer-3.4 https://github.com/EddyRivasLab/hmmer.git
cd hmmer
git clone --depth 1 --branch easel-0.49 https://github.com/EddyRivasLab/easel.git easel
autoconf && (cd easel && autoconf)
./configure --prefix="$HOME/local"
make -j4 && make install
export PATH="$HOME/local/bin:$PATH"
```

The easel submodule must be present before `autoconf`; a `--depth 1` clone of a tag does
not reliably fetch it, which is why it is cloned explicitly.

### Confirm

```sh
./.venv/bin/python scripts/doctor.py
./.venv/bin/python run_panel.py --check
```

`run_panel.py --check` reports each engine and names the specific remedy for anything
missing. On this platform it reports PROPERMAB's patch pipeline as unavailable, which is
expected: the bundled NanoShaper binary is Apple ARM64 only. Set
`PROPERMAB_NANOSHAPER` to a compatible binary to enable it, or run PROPERMAB on the
machine where it was built. `app.propermab.patch_available()` already verifies that any
configured override exists and is executable before claiming support.

## Headless panel runs

`run_panel.py` drives the same engines as the web application without a browser, so a
parent antibody and a set of engineered variants can be analysed in one unattended pass.

```sh
./.venv/bin/python run_panel.py --panel panel.json --out results/ \
    --fold --antiberty --substitutions --pdf --trim-to-domain
```

The panel file is JSON:

```json
{
  "baseline": "PARENT",
  "constructs": {
    "PARENT":  {"H": "QVQL...", "L": "SYEL..."},
    "VARIANT": {"H": "QVQL...", "L": "SYEL...",
                "preserve": {"L": [50, 91]},
                "notes": "carried into the output"}
  }
}
```

`preserve` lists one-based input positions to keep where possible, normally known
paratope contacts. Candidates that edit a preserved position are ranked last rather
than hidden, so the trade-off stays visible. The module cannot identify contacts
itself; that knowledge has to be supplied.

### Outputs

| File | Contents |
| --- | --- |
| `<construct>.json` | the full analysis record, same shape as the web API |
| `<construct>.pdb` | predicted Fv coordinates (with `--fold`) |
| `<construct>.pdf` | the three-page report (with `--pdf`) |
| `panel_summary.csv` | one row per chain: pI, charge, GRAVY, framework identity, PTM counts, CDR site list |
| `substitution_candidates.csv` | every ranked substitution with its evidence columns |
| `variant_specs.txt` | top-ranked proposals as `name \| H:N55Q` lines, pasteable into the variant comparison |
| `comparison_to_baseline.json` | the platform's variant comparison against the named baseline |
| `engine_report.json`, `run_manifest.json` | what ran, what did not, and why |

### `--trim-to-domain`

Sequence descriptors (pI, net charge, GRAVY) are computed over the whole supplied
chain, not over the numbered variable domain. A single stray constant-region residue
therefore moves them. This is easy to hit with patent-derived sequences: a "variable
light chain coding sequence" often carries one extra codon for the first residue of the
constant region, and that single Arg shifts light-chain net charge by a full unit and pI
by several tenths.

The runner always reports residues outside the numbered domain, in each record's
`domain_boundary`, in the `residues_outside_v_domain` column, and on stderr.
`--trim-to-domain` additionally removes them before the descriptors are computed.
Trimming is opt-in so that nothing is altered silently.

## What Linux does not change

Numbering, CDR assignment, framework identity, PTM motif detection, the structural
surface and PTM screens, AntiBERTy and the PDF report are the same calculations as on
macOS. Reference comparisons that require matching runtime and native-binary provenance
will be withheld across platforms; an unavailable comparison is not a low-risk result.

## PROPERMAB native patch pipeline on Linux

The full NanoShaper/APBS surface-patch pipeline **does** build and run on Linux x86_64. The
release notes describe recipient builds on other operating systems as unvalidated; this
records a working one. Two things were needed.

### 1. Native dependencies

`scripts/build_native.sh` needs Boost, GMP and MPFR development headers. The CGAL 4.2 and
NanoShaper 0.7.8 sources ship in `native-source/`, already patched. On a minimal image with
no distribution packages for these, conda-forge archives supply them without creating a
conda environment — download the package, unzip it (a `.conda` file is a zip holding
`pkg-*.tar.zst`), and extract the `include/` and `lib/` trees into a plain prefix:

| package | build used |
| --- | --- |
| `gmp` | `gmp-6.3.0-hac33072_2.conda` |
| `mpfr` | `mpfr-4.2.2-he0a73b1_0.conda` |
| `libboost-headers` | `libboost-headers-1.85.0-ha770c72_4.conda` |
| `libboost` / `libboost-devel` | `1.85.0` builds `h0ccab89_4` / `h00ab1b0_4` |

Boost 1.85 matches the patched CGAL 4.2; the bundled patch covers modern CMake/Boost
compatibility. Then:

```sh
export CMAKE_PREFIX_PATH="$HOME/local"
bash scripts/build_native.sh
export PROPERMAB_NANOSHAPER="$PWD/native-build/nanoshaper/NanoShaper"
export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"
```

CGAL and NanoShaper compile with warnings only (deprecated `std::random_shuffle` in CGAL's
spatial sort). APBS 3.4.1, PDB2PQR 3.7.1 and `multivalue` come from
`requirements-propermab.txt` as the `apbs-binary` wheel; nothing else native is required.
After this, `app.propermab.patch_available()` returns `True` and
`availability()['full_patch_pipeline']` is `True`.

### 2. One code fix

`app/propermab_worker.py` exported only `DYLD_LIBRARY_PATH` before invoking APBS. That is the
macOS dynamic-loader variable; on Linux APBS needs `LD_LIBRARY_PATH`, so the patch step
failed there and was reported as unavailable. It now sets whichever variable the running
platform uses and preserves any existing value. With that change the worker reports
`patch_pipeline: Complete: NanoShaper 0.7.8 / APBS 3.4.1 / PDB2PQR 3.7.1 / multivalue`
and returns all 31 descriptors, the seven patch features included.

### Comparing across platforms

A Linux NanoShaper/APBS build is not the build that produced `data/profiler_comparators.json`.
The release does not claim numerical equivalence across APBS runtimes, and the platform's own
comparison rules require matching runtime and structure preparation, so it withholds the
formal comparison. Split the descriptors when reading cross-platform results:

- **Sequence-derived** (`vh_charge`, `vl_charge`, `fv_charge`, `fv_csp`, `n_charged_res_fv`,
  `net_charge`, `fv_chml`, `cdr_length`, `aromatic_cdr`) are platform-independent and compare
  directly.
- **Patch features** (`*_patch_area`, `*_patch_area_cdr`, `heiden_score`) run through
  NanoShaper and APBS. Treat percentiles against the shipped panel as indicative.

To make the comparison self-consistent, rebuild the panel under this protocol with
`scripts/build_profiler_comparators_v2.py`. It reuses completed matching calculations, blocks
network access during inference and publishes atomically.
