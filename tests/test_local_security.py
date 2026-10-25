"""The local API must not become readable through a browser's remote origin."""
import pytest
from fastapi.testclient import TestClient
from app.main import app


@pytest.mark.parametrize('authority', ['127.0.0.1:8000', 'localhost:8000', '[::1]:8000'])
def test_loopback_cli_and_same_origin_browser_requests_work(authority):
    client = TestClient(app, base_url='http://' + authority)
    assert client.get('/api/health').status_code == 200
    assert client.get('/api/health', headers={'Origin': 'http://' + authority,
                                            'Sec-Fetch-Site': 'same-origin'}).status_code == 200


@pytest.mark.parametrize('host', ['attacker.example:8000', '127.0.0.1.attacker.example:8000',
                                 'localhost.attacker.example', 'user@localhost',
                                 '@localhost', 'localhost:garbage', '[invalid', 'localhost/path'])
def test_untrusted_or_malformed_host_is_rejected(host):
    client = TestClient(app)
    response = client.get('/api/health', headers={'Host': host})
    assert response.status_code == 400
    assert response.headers['cache-control'] == 'no-store'


@pytest.mark.parametrize('origin', ['https://attacker.example', 'null',
                                   'http://127.0.0.1:9999', 'https://127.0.0.1:8000',
                                   'http://localhost:8000', 'http://user@127.0.0.1:8000',
                                   'http://127.0.0.1:8000/path', 'http://[broken'])
def test_cross_origin_and_opaque_requests_cannot_reach_handlers(origin):
    client = TestClient(app, base_url='http://127.0.0.1:8000')
    for method in ('get', 'post'):
        response = getattr(client, method)('/api/health', headers={'Origin': origin})
        assert response.status_code == 403


def test_cross_site_without_origin_is_rejected():
    client = TestClient(app, base_url='http://localhost:8000')
    assert client.get('/api/health', headers={'Sec-Fetch-Site': 'cross-site'}).status_code == 403


def test_browser_can_open_local_home_from_an_external_link():
    client = TestClient(app, base_url='http://127.0.0.1:8000')
    navigation = {'Sec-Fetch-Site': 'cross-site', 'Sec-Fetch-Mode': 'navigate',
                  'Sec-Fetch-Dest': 'document'}
    assert client.get('/', headers=navigation).status_code == 200
    # Opening the shell is safe; navigation must not bypass API or Host guards.
    assert client.get('/api/health', headers=navigation).status_code == 403
    assert client.post('/', headers=navigation).status_code == 403
    assert client.get('/', headers={**navigation, 'Host': 'attacker.example'}).status_code == 400
    assert client.get('/', headers={**navigation, 'Origin': 'https://attacker.example'}).status_code == 403
    assert client.get('/', headers={**navigation, 'Sec-Fetch-Dest': 'iframe'}).status_code == 403
    assert client.get('/', headers={**navigation, 'Sec-Fetch-Mode': 'cors'}).status_code == 403


def test_testserver_only_works_for_in_process_transport():
    assert TestClient(app).get('/api/health').status_code == 200
    network_client = TestClient(app, client=('127.0.0.1', 41234))
    assert network_client.get('/api/health').status_code == 400


def test_default_origin_port_is_equivalent():
    client = TestClient(app, base_url='http://localhost')
    assert client.get('/api/health', headers={'Origin': 'http://localhost:80'}).status_code == 200
