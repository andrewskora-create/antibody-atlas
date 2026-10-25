"""Startup regressions: occupied ports, session reuse and download-free checks."""
import importlib.util
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import socket
import threading
from unittest import mock
import pytest
from scripts import launch


@pytest.mark.parametrize('value', ['abc', '0', '80', '65536', '-1', '8000;echo nope'])
def test_port_rejects_invalid_input(value):
    with pytest.raises(ValueError):
        launch.port_number(value)


def test_valid_port():
    assert launch.port_number('8010') == 8010


def test_free_port():
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        port = s.getsockname()[1]
    assert launch.server_state(port) == 'free'


@pytest.mark.parametrize('atlas', [False, True])
def test_identifies_server_without_stopping_it(atlas):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            data = {'info': {'title': 'Antibody Atlas' if atlas else 'Other service'}, 'paths': {'/api/analyze': {}}}
            self.wfile.write(json.dumps(data).encode())
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        assert launch.server_state(port) == ('atlas' if atlas else 'occupied')
        with mock.patch.object(launch, 'open_browser') as browser:
            assert launch.reuse_or_error(port) == (0 if atlas else 2)
            assert browser.call_count == (1 if atlas else 0)
        assert thread.is_alive()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_browser_can_be_disabled(monkeypatch):
    monkeypatch.setenv('ATLAS_NO_BROWSER', '1')
    with mock.patch.object(launch.webbrowser, 'open') as browser:
        launch.open_browser(8000)
        browser.assert_not_called()


def test_launch_binds_only_loopback_and_waits_for_ready():
    process = mock.Mock()
    process.poll.side_effect = [None, 0]
    process.wait.return_value = 0
    with mock.patch.object(launch, 'reuse_or_error', return_value=None), mock.patch.object(launch.subprocess, 'Popen', return_value=process) as popen, mock.patch.object(launch, 'server_state', return_value='atlas'), mock.patch.object(launch, 'open_browser') as browser:
        assert launch.start(Path('/python'), 8010) == 0
        assert popen.call_args.args[0][-4:] == ['--host', '127.0.0.1', '--port', '8010']
        browser.assert_called_once_with(8010)


def test_start_reuses_existing_without_subprocess():
    with mock.patch.object(launch, 'reuse_or_error', return_value=0), mock.patch.object(launch.subprocess, 'Popen') as popen:
        assert launch.start(Path('/python'), 8000) == 0
        popen.assert_not_called()


@pytest.mark.parametrize('platform_name,installed,calls', [('linux', False, 2), ('linux', True, 1), ('darwin', False, 0)])
def test_cpu_torch_installs_only_for_empty_linux_environment(platform_name, installed, calls):
    import run as atlas_run
    result = mock.Mock(returncode=0 if installed else 1)
    with mock.patch.object(atlas_run.subprocess, 'run', return_value=result) as command:
        atlas_run.ensure_linux_cpu_torch(Path('/env/python'), platform_name)
        assert command.call_count == calls
        if calls == 2:
            assert command.call_args.args[0] == ['/env/python', '-m', 'pip', 'install', 'torch', '--index-url', 'https://download.pytorch.org/whl/cpu']
            assert command.call_args.kwargs['check'] is True
