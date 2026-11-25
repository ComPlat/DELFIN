"""A refused loopback port says what route exists; the agent stops retrying.

Measured 2026-11-26 in the product's own cage (the agent's bash isolation
path, on the host the swarm runs on): a command's own server on the
loopback is ALREADY reachable with no setting at all -- a raw connect
CONNECTED, a GET and a POST through the proxy env answered 200, and a
server started by the session's background job was reachable from a later
command -- while the host's services on the same loopback (11434, 5432,
8088) answered EACCES, because on that host the cage shares the host's
network namespace and the socket guard judges every connect. The refusals
were right. What was missing was the sentence saying so: the port refusal
named only ``agent.sandbox_network.allowed_ports``, so an agent that
wanted its own server (already open) or a model endpoint (unreachable at
ANY port -- a command gets no provider key) read "add the port" and
retried or asked the user to widen it. In one day the proxy refused 3,500
commands on port 11434 alone.

Input: ``EgressProxy.decide`` / ``settle`` for a refused loopback port and
for an unknown domain with nobody to ask; and one real request driven
through the proxy's Unix socket, so the text is the one the command sees.
Output: the refusal keeps its verdict (deny stays deny, the own server
stays allowed) and names the route that exists: the session's own listener
is let through without a setting, a refused loopback port is somebody
else's service, and a script calling a model API cannot work because a
caged command carries no provider key.

Semantics: this changes no verdict -- every pair that was allowed stays
allowed and every refusal stays a refusal; only the reason text grows.
"""
from __future__ import annotations

import contextlib
import select
import socket
import threading

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP


class _Perms:
    task_session_id = "sess-refusal"
    mode = "bypassPermissions"
    confirm_callback = None


@pytest.fixture(autouse=True)
def _quiet_ownership(monkeypatch):
    """No session job listens on anything, and no model endpoint is
    registered: the decisions below are about the rules, not about
    whatever the machine running the suite happens to hold.

    The seam is the session's listener INODES (``own_listening_sockets``),
    which is what ``_own_server`` reads now that it judges the address and
    not the port number: patching a port set would no longer reach the rule
    and this file would pass on a tree where nothing was patched.
    """
    monkeypatch.setattr(BJ, "own_listening_sockets", lambda sid: set())
    monkeypatch.setattr(EP, "_ENDPOINT_HOSTS", set())
    monkeypatch.setattr(EP, "_ENDPOINT_ADDRESSES", set(), raising=False)


@contextlib.contextmanager
def _proxy():
    p = EP.EgressProxy(allowed_domains=[], allowed_ports=(80, 443))
    try:
        yield p
    finally:
        p.close()


# -- the port refusal on the loopback ----------------------------------------

def test_a_refused_loopback_port_names_the_route_that_exists():
    """The refusal an agent hits on 11434 (the largest single group of
    refusals) must tell it that its OWN server was never blocked and that
    a model call from a command cannot work at any port -- the two facts
    that end the retrying. Before, the text offered only
    ``allowed_ports``, which is the widening it must not chase."""
    with _proxy() as p:
        verdict, reason, _ = p.decide("127.0.0.1", 11434, perms=_Perms())
    assert verdict == "deny"                      # the verdict is unchanged
    assert "port 11434 is not allowed" in reason  # the pinned prefix stays
    assert "no setting at all" in reason          # own server: never blocked
    assert "no provider key" in reason            # a model POST cannot work
    assert "Do not retry" in reason


def test_a_refused_port_on_another_host_gets_no_loopback_lecture():
    """Control, and the wrong-input case for the assertion above: the note
    is about THIS machine's loopback. A remote host on a refused port
    keeps the plain allowed_ports reason -- an agent reading the loopback
    sentence there would be told the wrong thing."""
    with _proxy() as p:
        verdict, reason, _ = p.decide("models.example.org", 11434,
                                      perms=_Perms())
    assert verdict == "deny"
    assert "port 11434 is not allowed" in reason
    assert "no setting at all" not in reason


def test_a_session_that_listens_on_the_port_is_still_allowed(monkeypatch):
    """Control passing before and after: the ownership rule is what opens
    the loopback, and the text change must not disturb it. A listener one of
    the session's background jobs holds -- on the address under test -- is
    allowed, with a POST included.

    Patched at the seam the rule reads: the session's inodes and the kernel
    table naming them (``socket_guard.loopback_listeners``), because the
    rule judges the address, not the port number.
    """
    from delfin.agent import socket_guard as SG
    monkeypatch.setattr(BJ, "own_listening_sockets",
                        lambda sid: {"8393"} if sid == "sess-refusal" else set())
    monkeypatch.setattr(SG, "loopback_listeners",
                        lambda ip, port: {"8393"} if port == 8393 else set())
    with _proxy() as p:
        assert p.settle("127.0.0.1", 8393, method="POST", body_bytes=10,
                        perms=_Perms()) is None
        why = p.settle("127.0.0.1", 8394, method="POST", body_bytes=10,
                       perms=_Perms())
    assert why and "port 8394 is not allowed" in why


# -- the unknown-domain refusal ----------------------------------------------

def test_an_unknown_domain_says_a_command_cannot_call_a_model():
    """The second group (750 refusals of the model endpoint's domain): the
    route for a real model answer is DELFIN's own run, not the domain
    setting -- a command gets no provider key, so allowing the domain does
    not make the script work. Unattended (nobody to ask), so this is the
    text the agent actually reads."""
    with _proxy() as p:
        why = p.settle("api.modelhost.test", 443, perms=_Perms())
    assert why and "not an allowed domain" in why
    assert "no provider key" in why
    assert "Do not retry" in why


# -- the text a command really sees -------------------------------------------

def _through(p: EP.EgressProxy, url: str) -> str:
    """GET ``url`` the way a caged command does: a loopback forwarder into
    the proxy's Unix socket, token-authenticated (``tests`` of the proxy
    use the same shape); returns the status line and body head."""
    fwd = socket.socket(socket.AF_INET)
    fwd.bind(("127.0.0.1", 0))
    fwd.listen(4)
    port = fwd.getsockname()[1]

    def pipe():
        try:
            conn, _ = fwd.accept()
        except OSError:
            return
        up = socket.socket(socket.AF_UNIX)
        try:
            up.connect(p.socket_path)
        except OSError:
            conn.close()
            return
        pair = [conn, up]
        while True:
            r, _, _ = select.select(pair, [], [], 10)
            if not r:
                break
            for s in r:
                d = s.recv(65536)
                if not d:
                    break
                (up if s is conn else conn).sendall(d)
            if not r:
                break
        for s in pair:
            try:
                s.close()
            except OSError:
                pass

    t = threading.Thread(target=pipe, daemon=True)
    t.start()
    req = (f"GET {url} HTTP/1.1\r\nHost: x\r\n"
           f"Proxy-Authorization: Basic "
           + _basic(p.token)
           + "\r\nConnection: close\r\n\r\n")
    c = socket.create_connection(("127.0.0.1", port), timeout=15)
    c.sendall(req.encode())
    data = b""
    c.settimeout(15)
    try:
        while len(data) < 8192:
            chunk = c.recv(4096)
            if not chunk:
                break
            data += chunk
    except OSError:
        pass
    c.close()
    t.join(timeout=5)
    fwd.close()
    return data.decode("latin-1")


def _basic(password: str) -> str:
    import base64
    return base64.b64encode(f"delfin:{password}".encode()).decode()


def test_the_command_reads_the_route_in_its_403():
    """The full path, not the helper: a request driven through the proxy's
    own socket for a refused loopback port comes back 403 carrying the
    note. A refusal that never reaches the client says nothing."""
    with _proxy() as p:
        raw = _through(p, "http://127.0.0.1:11434/v1")
    assert raw.startswith("HTTP/1.1 403"), raw[:120]
    assert "no setting at all" in raw
    assert "no provider key" in raw
