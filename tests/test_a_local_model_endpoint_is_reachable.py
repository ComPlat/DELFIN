"""A model endpoint the session uses is reachable on its own port.

The egress proxy checks the port before it considers the host, and the
default ports are 80 and 443. A local model server listens on neither
(Ollama 11434, vLLM 8000, llama.cpp 8080), so registering a model
endpoint's host alone left every local endpoint unreachable from
anything the agent starts -- a script, a test, a probe. The configured
remote endpoints are on 443 and were unaffected, so the gap appeared
only where the whole installation is local.

Input: ``egress_proxy.note_endpoint`` with a model endpoint URL, and
``EgressProxy.decide(host, port)``.
Output: the verdict for that pair, for its host on another port, and for
another host on the same port.

Semantics: the exemption is the exact ``(host, port)`` pair the session
is configured to use. It relaxes the port gate only; the administrator's
ceiling, the host allow-list and the upload rules judge the request
afterwards exactly as they judge any known host.
"""
from __future__ import annotations

import pytest

from delfin.agent import egress_proxy as EP


@pytest.fixture(autouse=True)
def _no_endpoint_leaks(monkeypatch):
    """A fresh endpoint registry per test.

    The registry is module level and a process registers into it on the
    first model call, so without this one test's endpoint decides the
    next one's verdict.
    """
    monkeypatch.setattr(EP, "_ENDPOINT_HOSTS", set())
    # raising=False so this file states the DEFECT rather than an
    # AttributeError where the registry does not exist yet: a control run
    # at the baseline must fail on the verdict, which is the behaviour.
    monkeypatch.setattr(EP, "_ENDPOINT_ADDRESSES", set(), raising=False)
    yield


def _proxy():
    return EP.EgressProxy(allowed_domains=[], allowed_ports=(80, 443))


def test_a_local_ollama_is_refused_until_it_is_the_configured_endpoint():
    """The control first: the same pair, before anything registers it."""
    verdict, reason, _ = _proxy().decide("127.0.0.1", 11434)
    assert verdict == "deny"
    assert "port 11434 is not allowed" in reason

    EP.note_endpoint("http://127.0.0.1:11434/v1")
    verdict, reason, _ = _proxy().decide("127.0.0.1", 11434)
    assert verdict != "deny", reason


@pytest.mark.parametrize("url,host,port", [
    ("http://localhost:11434/v1", "localhost", 11434),      # Ollama
    ("http://127.0.0.1:8000/v1", "127.0.0.1", 8000),        # vLLM
    ("http://localhost:8080/v1", "localhost", 8080),        # llama.cpp
    ("http://10.0.0.7:1234/v1", "10.0.0.7", 1234),          # a server on the LAN
])
def test_every_local_server_port_is_reachable_once_configured(url, host, port):
    EP.note_endpoint(url)
    verdict, reason, _ = _proxy().decide(host, port)
    assert verdict != "deny", f"{host}:{port}: {reason}"


def test_the_exemption_does_not_open_the_port_on_another_host():
    """A configured endpoint on 11434 is not a rule about 11434."""
    EP.note_endpoint("http://127.0.0.1:11434/v1")
    verdict, reason, _ = _proxy().decide("models.example.org", 11434)
    assert verdict == "deny"
    assert "port 11434 is not allowed" in reason


def test_the_exemption_does_not_open_another_port_on_the_same_host():
    """Nor is it a rule about the host: only the pair passes."""
    EP.note_endpoint("http://127.0.0.1:11434/v1")
    verdict, reason, _ = _proxy().decide("127.0.0.1", 22)
    assert verdict == "deny"
    assert "port 22 is not allowed" in reason


def test_an_https_endpoint_registers_its_implied_port():
    """A URL without a port states one through its scheme."""
    EP.note_endpoint("https://api.example.org/v1")
    assert ("api.example.org", 443) in EP._ENDPOINT_ADDRESSES


def test_a_url_that_names_nothing_registers_nothing():
    for bad in ("", "   ", "not a url", "http://"):
        EP.note_endpoint(bad)
    assert EP._ENDPOINT_ADDRESSES == set()
    assert EP._ENDPOINT_HOSTS == set()
