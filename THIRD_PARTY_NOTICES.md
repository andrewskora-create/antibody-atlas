# Third-party notices

This file lists third-party components, their licences, and the notices their
licences require. Verified 2026-10-10 against upstream sources; see
COMPLIANCE.md for the audit and DATA_PROVENANCE.md for dataset terms.

Nothing here grants any right over a third-party component, dataset, trademark
or model beyond what that component's own licence grants.

## Bundled in this repository

### NGL viewer — MIT
`static/vendor/ngl/`. © 2014-2017 Alexander S. Rose. Full licence retained at
`static/vendor/ngl/LICENSE`.

### OGRDB germline sets — CC0
`data/human_v_ogrdb.json`, `data/human_v_ogrdb_numbered.json`. Published by the
AIRR Community under CC0; copyright in individual contributions retained by
contributors. CC0 imposes no attribution condition.

## Fetched at setup, redistributable, permissive

### AntiBERTy — MIT
Code and checkpoint, © 2023 Jeffrey Ruffolo. Licence retained at
the installed `antiberty` package. The checkpoint is fetched rather than committed
only because of its size.

## Declared dependencies — permissive

Obtained from their own distribution channels at install time under their
upstream terms. None imposes a non-commercial or academic-only restriction.

| Component | Licence |
|---|---|
| ANARCI | BSD-3-Clause |
| ANARCII | BSD-3-Clause |
| ImmuneBuilder / ABodyBuilder2 | BSD-3-Clause |
| Paragraph | BSD-3-Clause, © 2021 Lewis Alexander Chinery |
| FreeSASA | MIT |
| Biopython | Biopython License Agreement / BSD-3-Clause |
| OpenMM, PDBFixer | MIT (with permissively-licensed components) |
| FastAPI, Uvicorn, Pydantic | MIT / BSD-3-Clause |
| ReportLab | BSD-3-Clause |
| Pillow | MIT-CMU |
| PyTorch | BSD-3-Clause |

ImmuneBuilder structure checkpoints are downloaded from ImmuneBuilder's own
configured public locations during warmup, under ImmuneBuilder's terms.

HMMER is required by legacy ANARCI for numbering and is installed separately
(BSD-3-Clause). Build guidance is in `LINUX_AND_HEADLESS.md`.

## NOT distributed — user-installed, restricted

### PROPERMAB — Regeneron Pharmaceuticals, Inc. — ACADEMIC RESEARCH USE ONLY

**No Regeneron code is distributed with this software.** `setup/fetch_optional.py`
fetches it from Regeneron's own repository after the user affirms eligibility,
and reproduces Regeneron's licence in full at that point.

The licence limits the licensee class to "an academic institution or
academically employed full-time personnel only" and prohibits commercial use,
integration into software used for further commercialization, and distribution
in any commercial product or online service. Anyone obtaining it must identify
Regeneron Pharmaceuticals, Inc. as the original creator, retain its proprietary
notices, and include a copy of its licence. Nonprofit intent alone does not
establish eligibility. Commercial users must contact Regeneron.

Numeric descriptor values computed by PROPERMAB and published in this
repository's reference panel are research outputs, not Regeneron's source code.
Provenance is recorded per value, and the export tooling can omit them.

### NanoShaper — GPLv3; CGAL — GPLv3+ / LGPL by package

Required only by the PROPERMAB surface-patch path, and not distributed here.
`setup/fetch_optional.py` documents the build. A user who builds them takes on
the corresponding GPL obligations, including providing source for any binary
they themselves distribute.

### APBS and PDB2PQR

Required only by the PROPERMAB path; installed by the user from their own
channels under their own terms.

## Reference corpora — not distributed

Datasets derived from Oxford Thera-SAbDab are not distributed, because that
source publishes no licence grant. They are rebuilt locally from Oxford's CSV.
Oxford requests citation of Raybould, Marks et al., *Nucleic Acids Res.* (2019)
gkz827. See DATA_PROVENANCE.md.

## Excluded from this repository by design

User inputs, session records, generated reports and PDFs, temporary structures,
job output directories, credentials, virtual environments, machine-specific
native builds, logs and caches. `.gitignore` enforces this so analysis output
cannot be committed by accident.

---

## UniProtKB — human immunoglobulin constant-region sequences

`data/ig_constant_regions.json` contains sequence data obtained from UniProtKB
and redistributed under the Creative Commons Attribution 4.0 International
licence (<https://creativecommons.org/licenses/by/4.0/>).

> UniProtKB data are distributed under the Creative Commons Attribution 4.0 International (CC BY 4.0) License. Copyright 2002-2026 UniProt Consortium. See https://www.uniprot.org/help/license

UniProt release 2026_03, retrieved 2026-10-11.
Accessions: A0M8Q6, P01834, P01857-1, P01859-1, P01860-1, P01861-1, P0CG04, P0DOY2, P0DOY3.
Rebuild or verify with `python setup/fetch_constant_regions.py [--verify]`.

---

## DeepSP — spatial aggregation and charge prediction

DeepSP is distributed under the MIT License.

> MIT License
>
> Copyright (c) 2024 Lailabcode
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

Source: <https://github.com/Lailabcode/DeepSP>

What is derived from it here, and how:

- `models/deepsp/*.npz` — the authors' own trained weights, mechanically
  converted from Keras HDF5 to numpy arrays plus an ordered execution plan, so
  inference runs without a TensorFlow dependency. Dropout layers are omitted
  because they are the identity at inference time. No weight value is altered.
  Fetched and converted by `setup/fetch_deepsp.py`, which records the sha256 of
  every upstream file it read.
- `data/deepsp_reference.json` — percentile breakpoints computed from DeepSP's
  published training labels (molecular-dynamics values for 20,530 natural human
  Fv pairs). These are summary statistics, not the labels themselves.
- `tests/fixtures/deepsp_reference_molecules.json` — the sixteen reference
  molecules DeepSP publishes with the model, carried verbatim so the conversion
  can be verified offline.

Cite: Kalejaye L, Wu I-E, Terry T, Lai P-K. *DeepSP: Deep learning-based
spatial properties to predict monoclonal antibody stability.* Computational and
Structural Biotechnology Journal, 2024. doi:10.1016/j.csbj.2024.05.029

Neither the weights nor the training labels are redistributed in this
repository; `setup/fetch_deepsp.py` fetches them from the upstream repository
on the user's machine. MIT would permit redistribution; they are fetched
because they are 17 MB, the same reason AntiBERTy's weights are fetched.

## ThermoMPNN — predicted single-substitution folding stability

Upstream: <https://github.com/Kuhlman-Lab/ThermoMPNN>, pinned commit
`370f76ec62bd929f7425e311d8df04a0d094990f`. MIT License, copyright (c) 2023
Kuhlman-Lab. The original license is retained by `scripts/setup_stability.py`
at `vendor/thermompnn/LICENSE`, together with a SHA-256 source/checkpoint manifest.
The installer fetches source and model files; inference is entirely local.

Cite: Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B (2024).
Transfer learning to leverage larger datasets for improved prediction of protein
stability changes. PNAS 121(6):e2314853121.
<https://doi.org/10.1073/pnas.2314853121>.

Atlas presents these as learned predictions of folding stability, not measured
free energies, binding free energies, or a guarantee of retained antibody function.

ThermoMPNN's inherited ProteinMPNN, graph-protein-design and LightAttention
components also use MIT licenses. Their full upstream notices are retained in
`vendor/thermompnn/licenses/` by the installer, with pinned revisions and checksums
in `scripts/setup_stability.py` and the installation manifest. These permissions
include commercial use and redistribution subject to retaining the notices;
no academic-only restriction or license fee is introduced by this stability engine.
