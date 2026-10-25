#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if command -v python3.12 >/dev/null 2>&1; then
  exec python3.12 run.py "$@"
elif command -v python3.11 >/dev/null 2>&1; then
  exec python3.11 run.py "$@"
else
  echo "Install Python 3.11 or 3.12, then run this script again." >&2
  exit 2
fi
