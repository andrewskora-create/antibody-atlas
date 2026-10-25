#!/usr/bin/env python3
"""Local-only startup and readiness checks using Python's standard library."""
from __future__ import annotations
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
import webbrowser

ROOT = Path(__file__).resolve().parents[1]


def port_number(value: str) -> int:
    try:
        port = int(value)
    except (ValueError, TypeError):
        raise ValueError("PORT must be a number from 1024 to 65535.") from None
    if not 1024 <= port <= 65535:
        raise ValueError("PORT must be a number from 1024 to 65535.")
    return port


def server_state(port: int) -> str:
    """Distinguish an available port, Atlas, and a different occupied service."""
    with socket.socket() as sock:
        sock.settimeout(1)
        if sock.connect_ex(("127.0.0.1", port)) != 0:
            return "free"
    try:
        # Do not allow proxy settings to send this local check elsewhere.
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        with opener.open(f"http://127.0.0.1:{port}/openapi.json", timeout=2) as response:
            data = json.loads(response.read(2_000_000))
        if isinstance(data, dict) and isinstance(data.get("info"), dict) and data.get("info", {}).get("title") == "Antibody Atlas" and "/api/analyze" in data.get("paths", {}):
            return "atlas"
    except (OSError, ValueError, urllib.error.URLError):
        pass
    return "occupied"


def open_browser(port: int) -> None:
    if os.environ.get("ATLAS_NO_BROWSER") == "1":
        return
    try:
        webbrowser.open(f"http://127.0.0.1:{port}")
    except Exception:
        pass  # The printed URL remains usable in headless/WSL environments.


def reuse_or_error(port: int) -> int | None:
    state = server_state(port)
    if state == "atlas":
        print(f"Atlas is already running at http://127.0.0.1:{port}.\n"
              "This keeps the running session. After a code update, stop its terminal with Ctrl-C and start again.", flush=True)
        open_browser(port)
        return 0
    if state == "occupied":
        print(f"Port {port} is used by another service. No process was stopped.\n"
              "Choose a different port, for example: PORT=8010 ./run_local.sh", file=sys.stderr)
        return 2
    return None


def start(python: Path, port: int) -> int:
    existing = reuse_or_error(port)
    if existing is not None:
        return existing
    print(f"Starting Atlas at http://127.0.0.1:{port}\nKeep this terminal open. Press Ctrl-C to stop.", flush=True)
    process = subprocess.Popen([str(python), "-m", "uvicorn", "app.main:app", "--host", "127.0.0.1", "--port", str(port)], cwd=ROOT)
    try:
        deadline = time.monotonic() + 90
        while process.poll() is None and time.monotonic() < deadline:
            if server_state(port) == "atlas":
                open_browser(port)
                return process.wait()
            time.sleep(0.3)
        if process.poll() is not None:
            return process.returncode
        print("Atlas did not become ready within 90 seconds; stopping this startup. Read the error above and run .venv/bin/python scripts/doctor.py.", file=sys.stderr)
        process.terminate()
        return 2
    except KeyboardInterrupt:
        process.terminate()
        return 0
    finally:
        if process.poll() is None:
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
