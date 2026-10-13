"""A program the bash tool starts reaches the network only through the proxy.

The domain rules and the egress scan read the command TEXT. A Python
script or a binary started by that command opens its own sockets, which
no reading of the text sees. The kernel layer (socket guard, a network
namespace, or Landlock's TCP restriction) makes the egress proxy the only
door; the proxy then applies the domain decisions web_fetch is subject to:
the user's allowed domains, the permission rules for
``web_fetch(domain:...)``, and a question for an unknown domain or an
upload when someone can be asked.

Every request in this file goes to a listener on 127.0.0.1; a host name
like ``allowed.test`` reaches it through the proxy's injected connector.

This file runs on every host. The tests that start a child under a real
kernel layer (socket guard, Landlock ABI 4, a bubblewrap network
namespace) need that layer and live in
test_a_sandboxed_command_reaches_the_network_through_its_proxy.py.
"""
from __future__ import annotations

import base64
import os
import socket
import sys
import threading
import types
import urllib.error
import urllib.request

import pytest

import _egress_test_support as ES
from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import landlock_exec as LL
from delfin.agent import socket_guard as SG

HELPER = ES.HELPER
_perms = ES.perms_in


@pytest.fixture
def web():
    yield from ES.serve_web()


@pytest.fixture
def proxy(web, monkeypatch):
    yield from ES.session_proxy(web, monkeypatch)


class _Asker:
    def __init__(self, answer=True):
        self.answer = answer
        self.asked = []

    def __call__(self, name, args, preview):
        self.asked.append((name, dict(args), preview))
        return self.answer


def _via(proxy, perms, url, data=None, method=None):
    """A request through the proxy's Unix socket with the session token,
    the way the sandbox's forwarder delivers it."""
    token = proxy.token_for(perms)
    fwd = socket.socket(socket.AF_INET)
    fwd.bind(("127.0.0.1", 0))
    fwd.listen(4)
    port = fwd.getsockname()[1]

    def pipe():
        import select
        conn, _ = fwd.accept()
        up = socket.socket(socket.AF_UNIX)
        up.connect(proxy.socket_path)
        while True:
            r, _, _ = select.select([conn, up], [], [], 10)
            if not r:
                return
            for s in r:
                d = s.recv(65536)
                if not d:
                    return
                (up if s is conn else conn).sendall(d)

    threading.Thread(target=pipe, daemon=True).start()
    opener = urllib.request.build_opener(urllib.request.ProxyHandler(
        {"http": f"http://delfin:{token}@127.0.0.1:{port}"}))
    req = urllib.request.Request(url, data=data, method=method)
    try:
        return opener.open(req, timeout=10).read().decode()
    except urllib.error.HTTPError as exc:
        return f"{exc.code} {exc.read().decode()}"


# -- the proxy's decisions --------------------------------------------------

def test_an_allowed_domain_passes_without_a_question(proxy, web, tmp_path):
    asker = _Asker()
    perms = _perms(tmp_path, asker=asker)
    out = _via(proxy, perms, f"http://allowed.test:{web.server_address[1]}/")
    assert out.startswith("hello allowed.test"), out
    assert asker.asked == []


def test_an_unknown_domain_is_asked_when_someone_can_be_asked(proxy, web, tmp_path):
    port = web.server_address[1]
    yes = _Asker(True)
    out = _via(proxy, _perms(tmp_path, asker=yes), f"http://unknown.test:{port}/")
    assert out.startswith("hello unknown.test"), out
    assert len(yes.asked) == 1 and "unknown.test" in yes.asked[0][2]

    no = _Asker(False)
    out = _via(proxy, _perms(tmp_path, asker=no), f"http://other.test:{port}/")
    assert out.startswith("403") and "other.test" in out, out
    assert len(no.asked) == 1


def test_an_unknown_domain_with_nobody_to_ask_is_refused_by_name(proxy, web, tmp_path):
    port = web.server_address[1]
    for mode in ("default", "bypassPermissions", "plan"):
        out = _via(proxy, _perms(tmp_path, mode=mode), f"http://evil.test:{port}/")
        assert out.startswith("403"), (mode, out)
        assert "evil.test" in out and "web_fetch(domain:evil.test)" in out, out
        assert "allowed_domains" in out, out
    assert proxy.refused


def test_plan_mode_refuses_an_unknown_domain_without_asking(proxy, web, tmp_path):
    asker = _Asker(True)
    out = _via(proxy, _perms(tmp_path, mode="plan", asker=asker),
               f"http://unknown.test:{web.server_address[1]}/")
    assert out.startswith("403"), out
    assert asker.asked == []


def test_an_upload_is_asked(proxy, web, tmp_path):
    port = web.server_address[1]
    big = b"x" * (2 * 1024 * 1024)
    no = _Asker(False)
    out = _via(proxy, _perms(tmp_path, asker=no),
               f"http://allowed.test:{port}/up", data=big)
    assert out.startswith("403"), out
    assert len(no.asked) == 1 and "upload" in no.asked[0][2].lower()
    assert not any(m == "POST" for m, _h, _p in web.seen)

    yes = _Asker(True)
    out = _via(proxy, _perms(tmp_path, asker=yes),
               f"http://allowed.test:{port}/up", data=big)
    assert out == f"got {len(big)}", out

    # A small POST to an unknown host is an upload as well: asked as one.
    yes2 = _Asker(True)
    out = _via(proxy, _perms(tmp_path, asker=yes2),
               f"http://unknown.test:{port}/up", data=b"secret=1")
    assert out == "got 8", out
    assert len(yes2.asked) == 1 and "upload" in yes2.asked[0][2].lower()


def test_a_small_post_to_an_allowed_host_passes_unattended(proxy, web, tmp_path):
    """The model call of an unattended run: no question, no refusal."""
    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions"),
               f"http://allowed.test:{web.server_address[1]}/v1/chat",
               data=b'{"messages": []}')
    assert out == "got 16", out


def test_the_endpoint_the_session_uses_is_allowed(proxy, web, tmp_path):
    port = web.server_address[1]
    EP.note_endpoint(f"http://models.test:{port}/v1")
    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions"),
               f"http://models.test:{port}/v1/chat", data=b"{}")
    assert out == "got 2", out


def test_the_permission_rules_decide_first(proxy, web, tmp_path, monkeypatch):
    """The hook for web_fetch(domain:...) rules (permission_rules)."""
    port = web.server_address[1]
    rules = types.SimpleNamespace(
        deny=("allowed.test",), ask=("asked.test",), allow=("granted.test",))

    def friction(rs, name, args, workspace):
        host = args["url"].split("//", 1)[1].split(":", 1)[0].split("/", 1)[0]
        assert name == "web_fetch"
        if host in rs.deny:
            return "deny", types.SimpleNamespace(text=f"web_fetch(domain:{host})")
        if host in rs.ask:
            return "ask", types.SimpleNamespace(text=f"web_fetch(domain:{host})")
        return None

    def allows(rs, site, name, args, workspace, target=None):
        host = args["url"].split("//", 1)[1].split(":", 1)[0].split("/", 1)[0]
        return (types.SimpleNamespace(text=f"web_fetch(domain:{host})")
                if host in rs.allow else None)

    fake = types.SimpleNamespace(load=lambda workspace=None: rules,
                                 friction=friction, allows=allows)
    monkeypatch.setitem(sys.modules, "delfin.agent.permission_rules", fake)
    # `from . import permission_rules` reads the package attribute first,
    # which the real module (on main since #204) has already set.
    import delfin.agent as _agent_pkg
    monkeypatch.setattr(_agent_pkg, "permission_rules", fake, raising=False)

    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions"),
               f"http://allowed.test:{port}/")
    assert out.startswith("403") and "web_fetch(domain:allowed.test)" in out, out
    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions"),
               f"http://granted.test:{port}/")
    assert out.startswith("hello granted.test"), out
    asker = _Asker(True)
    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions", asker=asker),
               f"http://asked.test:{port}/")
    assert out.startswith("hello asked.test") and len(asker.asked) == 1, out


def test_a_real_permission_rule_decides_in_the_proxy(proxy, web, tmp_path,
                                                     monkeypatch):
    """The real permission_rules module, read from a settings file: a deny
    on an allow-listed host refuses even under bypass."""
    import json
    from delfin.agent import permission_rules as PR
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"permissions": {
        "deny": ["web_fetch(domain:allowed.test)"]}}), encoding="utf-8")
    monkeypatch.setattr(PR, "user_settings_path", lambda: settings)
    port = web.server_address[1]
    out = _via(proxy, _perms(tmp_path, mode="bypassPermissions"),
               f"http://allowed.test:{port}/")
    assert out.startswith("403") and "web_fetch(domain:allowed.test)" in out, out


def test_a_token_without_a_session_keeps_the_old_answer(proxy, web):
    """The process-wide token (Seatbelt, older callers): allow-list only."""
    port = web.server_address[1]
    perms = None
    assert _via(proxy, perms, f"http://allowed.test:{port}/").startswith("hello")
    assert _via(proxy, perms, f"http://unknown.test:{port}/").startswith("403")


# -- the kernel layer: how the bubblewrap path wires the proxy in ----------
# The runs under the real layers are in the sandbox network test file.

def test_the_bubblewrap_path_without_the_guard_binds_the_proxy_in(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setattr(SG, "available", lambda: False)
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    perms = A.KitToolPermissions(workspace=str(tmp_path), lock_workspace=True)
    argv = A._bash_isolation_argv("ls", tmp_path, perms)
    assert "--unshare-net" in argv
    i = argv.index(HELPER)
    assert "--netns" in argv[i:] and "--proxy-socket" in argv[i:], argv
    proxy_dir = os.path.dirname(argv[argv.index("--proxy-socket") + 1])
    assert argv[argv.index(proxy_dir) - 1] == "--ro-bind"


# -- honesty: where no layer exists, nothing claims one ---------------------

def _no_layers(monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: False)
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    monkeypatch.setattr(SG, "available", lambda: False)
    monkeypatch.setattr(LL, "abi_version", lambda: 0)
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})


def test_with_no_layer_the_notice_says_the_network_is_not_held(monkeypatch):
    _no_layers(monkeypatch)
    events = []
    monkeypatch.setattr(A, "_record_security_event",
                        lambda kind, tool, detail, blocked=False: events.append((kind, detail)))
    monkeypatch.setattr(A, "_CHILD_NETWORK_ANNOUNCED", False)
    state = A.child_network_enforcement()
    assert state["enforced"] is False and state["layer"] == "none"
    assert "NOT" in state["notice"]
    A._announce_child_network()
    assert any("NOT" in d and "network" in d for _k, d in events), events

    from delfin.agent import doctor
    rows = doctor._check_child_network({})
    assert rows and rows[0]["status"] == "WARN" and "NOT" in rows[0]["detail"]


def test_with_the_guard_the_capabilities_name_it(monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(SG, "available", lambda: True)
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    state = A.child_network_enforcement()
    assert state["enforced"] is True and state["layer"] == "socket guard"
    from delfin.agent import doctor
    assert doctor._check_child_network({})[0]["status"] == "PASS"


def test_landlock_tcp_alone_is_reported_as_partial(monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(SG, "available", lambda: False)
    monkeypatch.setattr(LL, "abi_version", lambda: 4)
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    state = A.child_network_enforcement()
    assert state["layer"] == "Landlock TCP" and state["enforced"] is True
    assert "UDP" in state["notice"]


def test_an_upload_inside_a_tunnel_is_asked(proxy, web, tmp_path):
    """HTTPS hides the method; the bytes leaving an admitted tunnel to a
    host the user let in by answering are counted, and the limit asks."""
    port = web.server_address[1]

    class Split(_Asker):
        def __call__(self, name, args, preview):
            self.asked.append(preview)
            return "UPLOAD" not in preview

    asker = Split()
    perms = _perms(tmp_path, asker=asker)
    token = proxy.token_for(perms)
    auth = base64.b64encode(f"delfin:{token}".encode()).decode()
    s = socket.socket(socket.AF_UNIX)
    s.connect(proxy.socket_path)
    s.settimeout(10)
    s.sendall(f"CONNECT unknown.test:{port} HTTP/1.1\r\n"
              f"Proxy-Authorization: Basic {auth}\r\n\r\n".encode())
    assert s.recv(4096).startswith(b"HTTP/1.1 200"), asker.asked
    body = b"x" * (2 * 1024 * 1024)
    try:
        s.sendall(b"POST /up HTTP/1.0\r\nHost: unknown.test\r\n"
                  + f"Content-Length: {len(body)}\r\n\r\n".encode() + body)
        reply = s.recv(4096)
    except OSError:
        reply = b""
    s.close()
    assert not reply.startswith(b"HTTP/1.0 200"), reply
    # At most upload_ask_bytes leave before the question; the rest does not.
    assert len(asker.asked) == 2 and "UPLOAD" in asker.asked[1], asker.asked


# -- what the session is configured to reach stays reachable ---------------

def test_the_configured_hosts_and_direct_endpoints(monkeypatch, tmp_path):
    monkeypatch.setattr(EP, "_DIRECT_TCP", {})
    monkeypatch.setattr(EP, "_GIT_REMOTES", {str(tmp_path): [
        "https://git.example.test/team/repo.git",
        "git@ssh.example.test:team/repo.git",
        "ssh://git@other.example.test:2222/team/repo.git"]})
    monkeypatch.setattr(EP, "slurm_controllers", lambda: ["10.0.0.1:6817"])
    monkeypatch.setattr(EP, "_resolve", lambda host, port: [f"{host}:{port}"])
    for var in [v for v in os.environ if v.upper().endswith("LICENSE_FILE")]:
        monkeypatch.delenv(var)
    monkeypatch.setenv("LM_LICENSE_FILE", "27000@lic.example.test:/opt/lic.dat")
    monkeypatch.setenv("PIP_INDEX_URL", "https://index.example.test/simple")
    direct = EP.direct_tcp_endpoints(tmp_path)
    assert {"10.0.0.1:6817", "ssh.example.test:22", "other.example.test:2222",
            "lic.example.test:0"} <= set(direct), direct
    hosts = EP.configured_domains(tmp_path)
    assert "git.example.test" in hosts and "index.example.test" in hosts


def test_a_license_server_admits_every_port_of_its_host():
    import ipaddress
    sup = SG.Supervisor.__new__(SG.Supervisor)
    sup.net_mode = "proxy"
    sup.allowed_tcp = {(ipaddress.ip_address("10.1.2.3"), 0)}
    assert sup.inet_refusal(ipaddress.ip_address("10.1.2.3"), 41234) == ""
    assert sup.inet_refusal(ipaddress.ip_address("10.1.2.4"), 41234)


# -- the switch belongs to the user -----------------------------------------

def test_a_session_cannot_loosen_the_network_mode_it_started_with(monkeypatch):
    """The setting is read when the session starts; a later write that
    loosens it (an agent's script editing the settings file) waits for the
    next start. Tightening applies at once."""
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    current = {"mode": "proxy"}
    monkeypatch.setattr(EP, "_settings", lambda: dict(current))
    assert EP.network_mode() == "proxy"
    current["mode"] = "open"
    assert EP.network_mode() == "proxy"
    current["mode"] = "none"
    assert EP.network_mode() == "none"
