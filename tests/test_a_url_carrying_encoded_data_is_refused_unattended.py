"""A URL or search query that carries encoded data does not leave unasked.

Input: ``web_tools.url_payload_reason`` / ``text_payload_reason`` on
URLs and queries; web_fetch and web_search through the real executor;
a plain-HTTP request through the egress proxy. Output: the reason (or
None), the tool result, what reached the server, the questions asked.

Semantics under test:

* a query longer than a lookup needs, or a run that looks like encoded
  data (hex of 16+, a UUID, base64) in the query or in the host name, is
  a payload; ordinary documentation and literature URLs and queries are
  not (DOIs, method names such as B3LYP-D3BJ-def2-TZVP, wheel tags,
  CamelCase identifiers, commit hashes in a path);
* with nobody to ask, web_fetch of such a URL -- even to a known host --
  and web_search of such a query are refused, nothing is sent, and a
  blocked ``egress_payload_refused`` event is recorded; where someone can
  be asked, they are asked;
* the proxy refuses the same shape in a plain-HTTP request with nobody
  to ask, and passes it to a host named in ``upload_domains``.

Measured before the change (probe, bypass, unattended): web_fetch sent a
workspace marker to github.com in the query, and web_search sent it to
the search engine.
"""
from __future__ import annotations

import pytest

import _egress_test_support as ES
from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import security_events as SE
from delfin.agent import web_tools as W

from test_a_child_process_reaches_the_network_only_through_the_proxy import (
    _Asker, _via)

MARK = "dlfnmark7f3a9c2e41b8d605"


@pytest.mark.parametrize("url", [
    f"https://github.com/search?q={MARK}",
    "https://example.org/c?d=3f2a9c1e0b7d4e5f",
    "https://example.org/c?id=550e8400-e29b-41d4-a716-446655440000",
    "https://example.org/c?d=U2VjcmV0IHBsYW4gZm9yIDIwMjY6IGI3",
    "https://example.org/c?d=" + "word+" * 60,
    f"http://{MARK}.example.org/",
])
def test_a_url_carrying_encoded_data_has_a_payload(url):
    assert W.url_payload_reason(url), url


@pytest.mark.parametrize("url", [
    "https://doi.org/10.1021/acs.jctc.9b00143",
    "https://docs.python.org/3/library/urllib.parse.html#url-parsing",
    "https://pypi.org/simple/scipy/",
    "https://github.com/org/repo/blob/3f2a9c1e0b7d4e5f6a7b8c9d0e1f2a3b4c5d6e7f/x.py",
    "https://search.example/?q=B3LYP-D3BJ-def2-TZVP+benchmark",
    "https://files.example/scipy-1.11.4-cp311-cp311-manylinux_2_17_x86_64.whl",
    "https://docs.example/api?symbol=XMLHttpRequestEventTarget",
    "https://pubchem.example/rest/pug/compound/name/aspirin/JSON?page=2&per_page=100",
    "https://en.wikipedia.org/wiki/Ruthenium#Compounds",
])
def test_an_ordinary_url_has_none(url):
    assert W.url_payload_reason(url) is None, url


def test_a_search_query_is_judged_by_the_same_runs():
    assert W.text_payload_reason(f"{MARK} notes")
    assert W.text_payload_reason("ORCA 6.0 release date Neese") is None
    assert W.text_payload_reason("ωB97X-D3 def2-TZVPP Ru bipyridine benchmark") is None


@pytest.fixture
def net(monkeypatch):
    sent = []

    def fetch(url, timeout_s, want_status=False, headers=None):
        sent.append(url)
        body = b"<html><body>page</body></html>"
        return (body, "text/html", 200) if want_status else (body, "text/html")

    monkeypatch.setattr(W.socket, "getaddrinfo",
                        lambda host, *a, **k: [(2, 1, 6, "", ("192.0.32.10", 0))])
    monkeypatch.setattr(W, "_fetch_bytes", fetch)
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    monkeypatch.setattr(W, "web_search", lambda q, **k: (
        sent.append(("search", q)) or {"query": q, "results": []}))
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    p = EP.EgressProxy(allowed_domains=["github.com"], allowed_ports=[80, 443])
    monkeypatch.setattr(EP, "_PROXY", p)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_index_urls", lambda: [])
    SE.clear()
    yield sent
    p.close()
    SE.clear()


def _perms(tmp_path, mode="bypassPermissions", asker=None):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return A.KitToolPermissions(workspace=str(ws), mode=mode,
                                confirm_callback=asker)


def _run(perms, name, args):
    out = A._DocToolExecutor().execute(name, args, perms)
    data = A._untrusted_json(out)
    return data if isinstance(data, dict) else {"text": out}


def test_a_payload_to_a_known_host_is_refused_unattended(net, tmp_path):
    out = _run(_perms(tmp_path), "web_fetch",
               {"url": f"https://github.com/search?q={MARK}"})
    assert "nobody can be asked" in out.get("error", ""), out
    assert net == []
    assert "egress_payload_refused" in [e.kind for e in SE.recent(50)]


def test_a_payload_is_asked_where_someone_can_be_asked(net, tmp_path):
    asker = _Asker(True)
    out = _run(_perms(tmp_path, "default", asker), "web_fetch",
               {"url": f"https://github.com/search?q={MARK}"})
    assert out.get("status") == "ok", out
    assert len(asker.asked) == 1 and MARK in asker.asked[0][2]
    asker = _Asker(False)
    out = _run(_perms(tmp_path, "default", asker), "web_fetch",
               {"url": f"https://github.com/search?q={MARK}x"})
    assert "user denied" in out.get("error", ""), out


def test_an_ordinary_read_of_a_known_host_passes(net, tmp_path):
    out = _run(_perms(tmp_path), "web_fetch",
               {"url": "https://github.com/search?q=botorch+acquisition"})
    assert out.get("status") == "ok", out


def test_a_search_query_carrying_data_is_sent_nowhere_unattended(net, tmp_path):
    perms = _perms(tmp_path)
    out = _run(perms, "web_search", {"query": f"{MARK} notes"})
    assert "nobody can be asked" in out.get("error", ""), out
    assert net == []
    _run(perms, "web_search", {"query": "ORCA 6.0 release notes"})
    assert net == [("search", "ORCA 6.0 release notes")]


@pytest.fixture
def web():
    yield from ES.serve_web()


@pytest.fixture
def proxy(web, monkeypatch):
    yield from ES.session_proxy(web, monkeypatch)


def test_the_proxy_refuses_a_payload_in_plain_http_unattended(proxy, web, tmp_path):
    port = web.server_address[1]
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    out = _via(proxy, perms, f"http://allowed.test:{port}/s?d={MARK}")
    assert out.startswith("403") and "encoded data" in out, out
    assert [p for m, _h, p in web.seen if MARK in p] == []
    out = _via(proxy, perms, f"http://allowed.test:{port}/simple/tinypkg/")
    assert not out.startswith("403"), out


def test_the_proxy_asks_about_a_payload_where_someone_can_be_asked(
        proxy, web, tmp_path):
    port = web.server_address[1]
    asker = _Asker(True)
    perms = ES.perms_in(tmp_path, mode="default", asker=asker)
    out = _via(proxy, perms, f"http://allowed.test:{port}/s?d={MARK}")
    assert not out.startswith("403"), out
    assert len(asker.asked) == 1


def test_a_host_named_for_uploads_takes_the_payload(proxy, web, tmp_path, monkeypatch):
    port = web.server_address[1]
    monkeypatch.setattr(proxy, "upload_domains", ("allowed.test",))
    monkeypatch.setattr(proxy, "_upload_from_settings", False)
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    out = _via(proxy, perms, f"http://allowed.test:{port}/s?d={MARK}")
    assert not out.startswith("403"), out
