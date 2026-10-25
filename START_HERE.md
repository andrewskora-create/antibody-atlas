# Antibody Atlas — start here

Antibody Atlas is a local, single-user workbench for research triage of antibody variable domains. It accepts a heavy chain, light chain, or paired VH/VL sequence and combines sequence analysis, structure prediction, local model scoring, comparisons, and a concise PDF report.

## Requirements and installation

The current verified setup is an Apple Silicon Mac. macOS Intel, Linux, and Windows/WSL2 are untested. Windows users should use WSL2; native Windows is unsupported. Use Python 3.11 or 3.12. Internet access is needed for initial Python package installation and structure-model downloads.

On macOS, install Apple's command-line build tools and Homebrew (https://brew.sh) first. In Terminal:

```sh
xcode-select --install
brew install python@3.12 hmmer
```

Add Homebrew's Python 3.12 to your `PATH` if needed (Homebrew prints the platform-specific instructions). HMMER provides `hmmscan`, required by legacy ANARCI through ImmuneBuilder; it is not used by ANARCII numbering.

On Ubuntu 24.04:

```sh
sudo apt update
sudo apt install python3.12 python3.12-venv python3.12-dev build-essential hmmer
```

From Terminal, enter the folder you received and run:

```sh
cd antibody-platform
python3.12 run.py --setup-only
python3.12 run.py --doctor
python3.12 run.py --warmup --setup-only
./run.sh
```

The doctor checks prerequisites; it does not install packages. Run it after Python dependency setup. `--setup-only` creates `.venv` and installs the base requirements. `--warmup` downloads structure-model weights without needing an antibody sequence. Initial setup and warmup need internet access. Open <http://127.0.0.1:8000> in a browser (set `PORT` to use another). Stop the server with Ctrl+C in Terminal. If `run.sh` is not executable, run `chmod +x run.sh` once.

For optional PROPERMAB Python packages, add `--with-propermab` to `--setup-only`. Its native surface-patch tools have additional prerequisites and are described in [MODULE_GUIDE.md](MODULE_GUIDE.md). The executable from the development Mac is not included in the portable ZIP.

## First analysis

1. Use **Load example** first to see the workflow with a non-confidential sequence.
2. Enter a heavy and/or light variable-domain sequence. Use standard amino-acid letters.
3. Choose IMGT or Kabat numbering, then start analysis.
4. Folding, AntiBERTy, and PROPERMAB run as workflow steps when enabled/available. The overall assessment is produced after the selected steps reach terminal status; a manual retry is available for a failed or interrupted step.
5. Review numbered residues, CDRs, framework identity, sequence descriptors, PTM flags, predicted structure, and profiler results with their limitations.
6. Download a PDF report or PDB if you need a local copy.

A favorable computational result is not a laboratory assay result or clinical safety conclusion. The app does not replace expression, stability, aggregation, viscosity, HIC, or other measurements.

## Privacy

The server binds to `127.0.0.1` and does not submit analysis sequences to a remote prediction API. Sequence inputs and analysis records are held in server memory; they are not saved automatically, and session structures are removed on normal shutdown. Structure-model weights may be fetched, but the sequence is not sent with that download request. Browser downloads, backups, and cloud-synced download folders are outside the app's control.

The local server has no login and other processes on the same computer may be able to access it. Use a trusted computer and account. Pasting a sequence into Codex, ChatGPT, or another hosted assistant sends that content through that service. The Codex in-app browser may expose the page to Codex; use a separate browser and avoid asking Codex to inspect it when the sequence must stay outside the conversation.

## Herceptin comparison

The optional “Herceptin-like” label compares selected PROPERMAB surface descriptors against trastuzumab references under a user-selected tolerance (default ±10%). It is a computational screening rule, not an assay comparison, validated cutoff, or evidence of experimental equivalence. See [MODULE_GUIDE.md](MODULE_GUIDE.md).

## Sharing the package
Send the ZIP unchanged to an eligible academic recipient, along with its SHA-256 file if desired. Ask them to extract it to a local folder, read this guide, and test **Load example**. No Codex, ChatGPT, API key, or paid prediction service is required. Prefer a folder outside cloud sync for confidential work. Allow several GB of free disk space; 16 GB RAM is a practical starting point, not a measured minimum. CPU operation is supported; a GPU is not required.

## Troubleshooting
- `No module named ...`: run `python3.12 run.py --setup-only` again and preserve the installation error if it fails.
- Missing `hmmscan`: install HMMER and restart Terminal.
- Missing PROPERMAB patches: install the extras and native tool described in MODULE_GUIDE; the local surface heatmap is a separate calculation and remains available.
- Reference comparisons unavailable after installing on another machine: runtime, native-binary or model differences may invalidate compatibility. Do not interpret unavailable values as low risk. Rebuild public references under the new protocol if needed.
- Port 8000 already in use: stop the previous app Terminal with Ctrl+C, then start again.
- The initial installation or model download fails: check internet connectivity and available disk space; no confidential sequence is needed during setup.

Upstream installation references: [ImmuneBuilder](https://github.com/oxpig/ImmuneBuilder), [FreeSASA](https://github.com/freesasa/freesasa-python), [HMMER](https://hmmer.org/installation.html), [WSL](https://learn.microsoft.com/windows/wsl/install). Windows users should install Ubuntu 24.04 in WSL2 and follow the Linux steps in its terminal.

## Release verification
The ZIP was extracted and its per-file hashes were checked. Package checks and API/report tests passed from the extracted copy using the development Python environment. The included NanoShaper/CGAL build helper successfully rebuilt the native tool on Apple Silicon. A fresh dependency installation on another computer, Linux or WSL2 has not been validated. `TESTED_ENVIRONMENT.txt` records the development package versions for troubleshooting; it is not a universal lockfile.
