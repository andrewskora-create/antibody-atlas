#!/usr/bin/env bash
# Start from a downloaded release or clone. Never remove an existing environment.
set -euo pipefail
cd "$(dirname "$0")"
export PORT="${PORT:-8000}"
# Finder's PATH may omit Homebrew even when Python/HMMER are installed.
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"
usable() {
  command -v "$1" >/dev/null 2>&1 && "$1" -c 'import sys; raise SystemExit(sys.version_info[:2] not in ((3,11),(3,12)))' 2>/dev/null
}
PY=""
if [ -n "${PYTHON:-}" ]; then
  if usable "$PYTHON"; then PY="$PYTHON"; else echo 'PYTHON must point to Python 3.11 or 3.12.' >&2; exit 2; fi
elif [ -e .venv/bin/python ]; then
  if usable .venv/bin/python; then PY=.venv/bin/python; else
    echo 'The existing .venv is incompatible. It was not removed. Rename it, then rerun with Python 3.11 or 3.12.' >&2; exit 2
  fi
else
  for candidate in python3.12 python3.11 python3; do
    if usable "$candidate"; then PY="$candidate"; break; fi
  done
fi
if [ -z "$PY" ]; then
  echo 'Install Python 3.11 or 3.12 first.' >&2
  echo 'macOS (Homebrew): brew install python@3.12 hmmer git' >&2
  echo 'Ubuntu 24.04: sudo apt install python3 python3-venv python3-dev build-essential hmmer git' >&2
  exit 2
fi
exec "$PY" run.py --stability "$@"
