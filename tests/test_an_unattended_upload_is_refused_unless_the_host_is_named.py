"""With nobody to ask, data leaving through the egress proxy is refused.

Input: a request a sandboxed command sends through the session's proxy
token, in a session with no approval channel (plan, bypassPermissions,
or no confirm callback). Output: the proxy's answer (pass, or 403 with
the reason), what reached the server, and the security events.

Semantics under test:

* a read (GET, a small request inside a tunnel) to a known host passes
  unattended, as before -- downloads keep working;
* a write (POST/PUT/PATCH/DELETE in plain HTTP, a body above
  ``upload_ask_bytes``, or more than that many bytes sent inside a
  CONNECT tunnel) to a known host is refused unattended unless the host
  is named in ``agent.sandbox_network.upload_domains`` or is the model
  endpoint the session itself uses;
* the refusal names the host and the setting, and is recorded as a
  blocked ``egress_upload_refused`` event (which reaches the attention
  inbox);
* attended sessions are unchanged: a small POST to a known host passes
  without a question.

Every request goes to a listener on 127.0.0.1 through the proxy's
injected connector; host names are labels only.
"""
from __future__ import annotations

import base64
import socket

import pytest

import _egress_test_support as ES
from delfin.agent import egress_proxy as EP
from delfin.agent import security_events as SE

from test_a_child_process_reaches_the_network_only_through_the_proxy import (
    _Asker, _via)


@pytest.fixture
def web():
    yield from ES.serve_web()


@pytest.fixture
def proxy(web, monkeypatch):
    yield from ES.session_proxy(web, monkeypatch)


@pytest.fixture(autouse=True)
def _events():
    SE.clear()
    yield
    SE.clear()


def _kinds():
    return [e.kind for e in SE.recent(200)]


def _posts(web):
    return [p for m, _h, p in web.seen if m == "POST"]


@pytest.mark.parametrize("mode", ["bypassPermissions", "plan", "default"])
def test_a_large_upload_to_a_known_host_is_refused_with_nobody_to_ask(
        proxy, web, tmp_path, mode):
    port = web.server_address[1]
    out = _via(proxy, ES.perms_in(tmp_path, mode=mode),
               f"http://allowed.test:{port}/up", data=b"x" * (2 * 1024 * 1024))
    assert out.startswith("403"), out
    assert "allowed.test" in out and "upload_domains" in out, out
    assert _posts(web) == []
    assert "egress_upload_refused" in _kinds()
    assert "egress_upload" not in _kinds()


def test_a_small_write_to_a_known_host_is_refused_with_nobody_to_ask(
        proxy, web, tmp_path):
    port = web.server_address[1]
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    for method in ("POST", "PUT", "PATCH", "DELETE"):
        out = _via(proxy, perms, f"http://allowed.test:{port}/api",
                   data=b'{"a": 1}' if method != "DELETE" else None,
                   method=method)
        assert out.startswith("403") and "upload_domains" in out, (method, out)
    assert _posts(web) == []


def test_a_read_from_a_known_host_still_passes_unattended(proxy, web, tmp_path):
    out = _via(proxy, ES.perms_in(tmp_path, mode="bypassPermissions"),
               f"http://allowed.test:{web.server_address[1]}/simple/tinypkg/")
    assert "tinypkg-0.1" in out, out
    assert "egress_upload_refused" not in _kinds()


def test_a_host_named_for_uploads_receives_them_unattended(proxy, web, tmp_path):
    port = web.server_address[1]
    proxy.upload_domains = ("allowed.test",)
    out = _via(proxy, ES.perms_in(tmp_path, mode="bypassPermissions"),
               f"http://allowed.test:{port}/up", data=b"x" * (2 * 1024 * 1024))
    assert out == f"got {2 * 1024 * 1024}", out
    # Recorded, not blocked: the user named the host, and nobody saw it.
    events = [e for e in SE.recent(200) if e.kind == "egress_upload"]
    assert events and not events[0].blocked


def test_the_setting_names_hosts_for_uploads(monkeypatch, web):
    monkeypatch.setattr(EP, "_settings", lambda: {
        "mode": "proxy", "upload_domains": ["up.example.test", "*.lab.test"]})
    p = EP.EgressProxy(allowed_domains=["allowed.test"],
                       allowed_ports=[web.server_address[1]])
    try:
        assert p.upload_allowed("up.example.test")
        assert p.upload_allowed("api.lab.test")
        assert not p.upload_allowed("api.github.com")
    finally:
        p.close()


def test_the_model_endpoint_still_receives_its_requests_unattended(
        proxy, web, tmp_path):
    port = web.server_address[1]
    EP.note_endpoint(f"http://models.test:{port}/v1")
    out = _via(proxy, ES.perms_in(tmp_path, mode="bypassPermissions"),
               f"http://models.test:{port}/v1/chat", data=b"{}")
    assert out == "got 2", out


def test_an_attended_small_post_to_a_known_host_is_unchanged(proxy, web, tmp_path):
    asker = _Asker(True)
    out = _via(proxy, ES.perms_in(tmp_path, asker=asker),
               f"http://allowed.test:{web.server_address[1]}/v1/chat",
               data=b'{"messages": []}')
    assert out == "got 16", out
    assert asker.asked == []


def _tunnel(proxy, perms, host, port):
    token = proxy.token_for(perms)
    auth = base64.b64encode(f"delfin:{token}".encode()).decode()
    s = socket.socket(socket.AF_UNIX)
    s.connect(proxy.socket_path)
    s.settimeout(10)
    s.sendall(f"CONNECT {host}:{port} HTTP/1.1\r\n"
              f"Proxy-Authorization: Basic {auth}\r\n\r\n".encode())
    assert s.recv(4096).startswith(b"HTTP/1.1 200")
    return s


def test_a_tunnel_to_a_known_host_carries_a_read_unattended(proxy, web, tmp_path):
    port = web.server_address[1]
    s = _tunnel(proxy, ES.perms_in(tmp_path, mode="bypassPermissions"),
                "allowed.test", port)
    s.sendall(b"GET / HTTP/1.0\r\nHost: allowed.test\r\n\r\n")
    reply = b""
    while True:
        chunk = s.recv(4096)
        if not chunk:
            break
        reply += chunk
    s.close()
    assert reply.startswith(b"HTTP/1.0 200") and b"hello allowed.test" in reply


def test_a_tunnel_to_a_known_host_is_cut_past_the_upload_limit_unattended(
        proxy, web, tmp_path):
    port = web.server_address[1]
    s = _tunnel(proxy, ES.perms_in(tmp_path, mode="bypassPermissions"),
                "allowed.test", port)
    body = b"x" * (2 * 1024 * 1024)
    try:
        s.sendall(b"POST /up HTTP/1.0\r\nHost: allowed.test\r\n"
                  + f"Content-Length: {len(body)}\r\n\r\n".encode() + body)
        reply = s.recv(4096)
    except OSError:
        reply = b""
    s.close()
    assert not reply.startswith(b"HTTP/1.0 200"), reply
    assert _posts(web) == []
    assert "egress_upload_refused" in _kinds()


def test_the_new_event_kind_has_a_label():
    assert "egress_upload_refused" in SE.known_kinds()
