"""What a session sends to each host is counted and recorded.

Input: requests a session's commands send through the egress proxy, and
web_fetch calls. Output: ``egress_proxy.sent_by_host(perms)``, the
``egress_sent`` security events and the ``egress`` lines of the audit log.

Semantics under test:

* every byte a command sends through the proxy is counted for its
  session and host (plain HTTP: request and body; a tunnel: everything
  the command writes into it); a web_fetch counts its URL plus a fixed
  header allowance;
* the total is recorded when it crosses a step (the first byte, 64 KiB,
  1 MiB, ...): a non-blocking ``egress_sent`` event and an audit-log
  line with tool ``egress``, the host and the session's total -- once per
  step, not per request;
* sessions are counted apart.

Before: the proxy decided per request and kept no count; nothing said
after an unattended run how much went where.
"""
from __future__ import annotations

import json

import pytest

import _egress_test_support as ES
from delfin.agent import api_client as A
from delfin.agent import audit_log
from delfin.agent import egress_proxy as EP
from delfin.agent import security_events as SE
from delfin.agent import web_tools as W

from test_a_child_process_reaches_the_network_only_through_the_proxy import (
    _Asker, _via)


@pytest.fixture
def web():
    yield from ES.serve_web()


@pytest.fixture
def proxy(web, monkeypatch, tmp_path):
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: tmp_path / "audit.log")
    SE.clear()
    yield from ES.session_proxy(web, monkeypatch)
    SE.clear()


def _audit(tmp_path):
    path = tmp_path / "audit.log"
    if not path.exists():
        return []
    return [json.loads(x) for x in path.read_text().splitlines()
            if json.loads(x).get("tool") == "egress"]


def test_a_commands_requests_are_counted_per_session_and_host(
        proxy, web, tmp_path):
    port = web.server_address[1]
    perms = ES.perms_in(tmp_path, mode="default", asker=_Asker(True))
    other = ES.perms_in(tmp_path, mode="default", asker=_Asker(True))
    _via(proxy, perms, f"http://allowed.test:{port}/up", data=b"x" * 5000)
    _via(proxy, perms, f"http://allowed.test:{port}/up", data=b"y" * 5000)
    sent = EP.sent_by_host(perms)
    assert sent.get("allowed.test", 0) >= 10000, sent
    assert EP.sent_by_host(other) == {}
    # One record for the first byte, none for the second request.
    events = [e for e in SE.recent(200) if e.kind == "egress_sent"]
    assert len(events) == 1 and not events[0].blocked
    lines = _audit(tmp_path)
    assert len(lines) == 1 and lines[0]["host"] == "allowed.test"
    assert lines[0]["bytes_sent"] > 0


def test_crossing_a_step_is_recorded_again(proxy, web, tmp_path):
    port = web.server_address[1]
    perms = ES.perms_in(tmp_path, mode="default", asker=_Asker(True))
    _via(proxy, perms, f"http://allowed.test:{port}/up", data=b"x" * 100)
    _via(proxy, perms, f"http://allowed.test:{port}/up", data=b"x" * 70000)
    totals = [rec["bytes_sent"] for rec in _audit(tmp_path)]
    assert len(totals) == 2 and totals[1] >= 64 * 1024, totals


def test_a_web_fetch_is_counted(tmp_path, monkeypatch):
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: tmp_path / "audit.log")
    monkeypatch.setattr(W.socket, "getaddrinfo",
                        lambda host, *a, **k: [(2, 1, 6, "", ("192.0.32.10", 0))])
    monkeypatch.setattr(W, "_fetch_bytes", lambda url, t, **k: (
        b"<html>ok</html>", "text/html"))
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    p = EP.EgressProxy(allowed_domains=["docs.allowed.test"])
    monkeypatch.setattr(EP, "_PROXY", p)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    try:
        perms = ES.perms_in(tmp_path, mode="bypassPermissions")
        url = "https://docs.allowed.test/guide.html"
        A._DocToolExecutor().execute("web_fetch", {"url": url}, perms)
        assert EP.sent_by_host(perms) == {
            "docs.allowed.test": len(url) + A._WEB_FETCH_HEADER_BYTES}
        assert [r["via"] for r in _audit(tmp_path)] == ["web_fetch"]
    finally:
        p.close()


def test_a_new_session_never_inherits_a_gone_ones_count(tmp_path, monkeypatch):
    """Sessions without an id were keyed by id(perms); CPython reuses an
    id once the object is gone, so a new session started with the old
    one's total and missed its own first-byte record (CI, shard 2).
    Two live objects given the same id stand for that reuse here."""
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: tmp_path / "audit.log")
    monkeypatch.setattr(EP, "id", lambda obj: 4242, raising=False)
    first = ES.perms_in(tmp_path, mode="default")
    second = ES.perms_in(tmp_path, mode="default")
    EP.note_sent(first, "reuse.test", 50000, "proxy")
    assert EP.sent_by_host(second) == {}
