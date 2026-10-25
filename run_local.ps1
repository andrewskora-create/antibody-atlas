# Antibody Atlas - local setup for Windows (PowerShell). Safe to re-run.
#
# Native Windows runs the sequence-level screen: numbering, PTM detection,
# charge, hydrophobicity, humanness, language-model scoring and the clinical
# comparator search. Structure prediction needs HMMER, which has no supported
# native Windows build, so structure, exposure, surface patches and paratope
# prediction are unavailable here. For the full pipeline use WSL2 and
# ./run_local.sh -- see PLATFORMS.md.

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
$Port = if ($env:PORT) { $env:PORT } else { "8000" }

Write-Host "Antibody Atlas - local setup (Windows)`n"

$py = Get-Command python -ErrorAction SilentlyContinue
if (-not $py) { $py = Get-Command py -ErrorAction SilentlyContinue }
if (-not $py) { Write-Host "No Python found. Install Python 3.11+ from python.org and re-run."; exit 1 }
$ver = & $py.Source -c "import sys;print('%d.%d'%sys.version_info[:2])"
Write-Host "  python $ver"

if (-not (Test-Path ".venv")) { Write-Host "  creating .venv"; & $py.Source -m venv .venv }
$venvPy = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
& $venvPy -m pip install --quiet --upgrade pip

if (-not (Test-Path ".venv\.deps-installed")) {
  Write-Host "  installing dependencies (first run only; a few minutes, torch is large)"
  & $venvPy -m pip install --quiet -r requirements.txt
  New-Item -ItemType File -Path ".venv\.deps-installed" -Force | Out-Null
} else { Write-Host "  dependencies already installed" }

if (Get-Command hmmscan -ErrorAction SilentlyContinue) {
  Write-Host "  hmmscan found"
} else {
  Write-Host ""
  Write-Host "  ! hmmscan is not available on native Windows."
  Write-Host "    Structure prediction, exposure, surface patches and paratope prediction"
  Write-Host "    will be skipped. Everything sequence-level still runs."
  Write-Host "    For the full pipeline, use WSL2 and ./run_local.sh (see PLATFORMS.md)."
  Write-Host ""
}

if (-not (Test-Path "data\clinical_frameworks.json")) {
  Write-Host ""
  Write-Host "  Clinical precedent data is not built yet (not redistributed; see DATA_PROVENANCE.md)."
  $answer = Read-Host "  Build now? [f]ull (~10 min) / [q]uick panel only (~1 min) / [s]kip"
  switch ($answer) {
    "f" { & $venvPy setup\build_reference_data.py --all }
    "q" { & $venvPy setup\build_reference_data.py --panel --ptm-sequences }
    default { Write-Host "  skipped; precedent columns will read 'not built'." }
  }
}

Write-Host "`n  starting on http://127.0.0.1:$Port"
Write-Host "  (Ctrl-C to stop)`n"
Start-Job -ScriptBlock { Start-Sleep 4; Start-Process "http://127.0.0.1:$using:Port" } | Out-Null
& $venvPy -m uvicorn app.main:app --host 127.0.0.1 --port $Port
