"""A caged command reaches the server its own session started, and nothing
else on this machine.

Session 0db62cdf (2026-10-19): `curl http://127.0.0.1:8324/` from the cage
was refused ("port 8324 is not allowed") -- the agent could not look at the
game it was serving, and asked the user to test every change. The egress
proxy now admits, for a session's token, the loopback ADDRESS that session's
running background jobs listen on (the listener inodes from
bash_jobs.own_listening_sockets, resolved against the kernel's listening
sockets), and connects to the loopback itself. Other addresses, other
sessions' servers and the process token stay refused.

A port number alone is not an endpoint -- ``[::1]:P`` and ``127.0.0.1:P``
coexist in one network namespace -- which is why the rule reads the address.
"""

from __future__ import annotations

import socket
import sys
import time

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP



class _Perms:
    task_session_id = "sess-own"
    mode = "bypassPermissions"


@pytest.fixture
def proxy(monkeypatch):
    """The session's own listener, named by inode as the cage names it.

    The rule is ``socket_guard.owns_loopback`` at the address the request
    goes to, so the seam is the session's LISTENER INODES
    (``bash_jobs.own_listening_sockets``) -- what ``--own-listener``
    carries -- not a set of port numbers. Measured 2026-11-25: a port number
    is not an endpoint (``[::1]:P`` and ``127.0.0.1:P`` coexist), so a port
    set cannot say which address the session's server is on.
    """
    monkeypatch.setattr(BJ, "own_listening_sockets",
                        lambda sid: {"8393"} if sid == "sess-own" else set())
    from delfin.agent import socket_guard as SG
    # The kernel table this file reads: one listener on the address under
    # test, held by inode "8393" -- the session's own.
    monkeypatch.setattr(
        SG, "loopback_listeners",
        lambda ip, port: {"8393"} if port == 8393 else set())
    return EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))



@pytest.mark.parametrize("host", ["127.0.0.1", "localhost", "::1",
                                  "[::1]", "127.0.1.1"])
def test_the_sessions_own_port_on_the_loopback_is_allowed(proxy, host):
    """Every spelling of "this machine" is judged at the address the proxy
    connects to (``egress_proxy._loopback_addr``).

    The five spellings do NOT exercise five behaviours, and should not be
    read as five. ``_loopback_addr`` keeps a literal loopback ADDRESS as
    written and resolves only a NAME: ``localhost`` -> ``127.0.0.1``,
    ``::1``/``[::1]`` -> ``::1``, and ``127.0.0.1``/``127.0.1.1`` ->
    themselves. So this test covers two addresses (IPv4 and IPv6) and one
    name; which ADDRESS a request is judged at, and the fact that the same
    port number on the other address is somebody else's, is the subject of
    ``test_a_sessions_own_server_is_judged_by_its_address.py`` -- that is the
    file that fails if the address half breaks.
    """
    assert proxy.decide(host, 8393, perms=_Perms())[0] == "allow"


@pytest.mark.parametrize("host,port,perms", [
    ("127.0.0.1", 22, _Perms()),          # not the session's
    ("127.0.0.1", 8394, _Perms()),        # nothing of it listens there
    ("example.org", 8393, _Perms()),      # not the loopback
    ("10.0.0.5", 8393, _Perms()),
    ("127.0.0.1", 8393, None),            # the process token: no session
])
def test_everything_else_is_judged_as_before(proxy, host, port, perms):
    assert proxy.decide(host, port, perms=perms)[0] != "allow"


def test_the_ports_are_read_from_the_jobs_processes(tmp_path):
    """The product's own view: ``open_ports`` -- what the terminal's /ports
    lists (``repl_commands.py``; the dashboard's Background panel reaches the
    same data through ``job_listening_ports`` directly, ``background_view.py``)
    -- learns a job's listening port from /proc, gives it to that job's
    session only, and stops listing it once the job has ended.

    Asserted through ``open_ports``/``job_listening_ports`` and not through a
    session-wide port set, because the egress rule judges the listener's
    ADDRESS (``own_listening_sockets``), for which a port number is not an
    endpoint: ``[::1]:P`` and ``127.0.0.1:P`` coexist.
    """
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    code = (f"import socket,time; s=socket.socket(); "
            f"s.bind(('127.0.0.1',{port})); s.listen(); time.sleep(30)")
    job = BJ.get_registry().start(
        command=f"python -c listen-{port}", cwd=str(tmp_path),
        workspace=str(tmp_path), session_id="sess-proc",
        argv=[sys.executable, "-c", code])
    try:
        deadline = time.time() + 10
        while time.time() < deadline and not any(
                port in entry["ports"]
                for entry in BJ.open_ports("sess-proc")):
            time.sleep(0.1)
        mine = [e for e in BJ.open_ports("sess-proc")
                if e["job_id"] == job.job_id]
        assert mine, "the job that listens on the port is not listed"
        assert port in mine[0]["ports"]
        assert mine[0]["session_id"] == "sess-proc"
        # Another session's view, and the unfiltered one, hold a job only
        # under the session that started it.
        assert [e for e in BJ.open_ports("someone-else")
                if e["job_id"] == job.job_id] == []
    finally:
        BJ.get_registry().kill(job.job_id)
        deadline = time.time() + 10
        while job.poll() is None and time.time() < deadline:
            time.sleep(0.05)
    assert [e for e in BJ.open_ports("") if e["job_id"] == job.job_id] == [], \
        "a job that has ended is still listed as listening"


def test_no_session_id_owns_no_server(tmp_path):
    """A caller that names no session owns NO server on the loopback.

    The by-port rule asserted this (``own_listening_ports("") == set()``);
    the rewrite that made the rule address-aware replaced that line with a
    question about a *known* job id, which a broken empty-session guard
    would still answer correctly. So the property is asked directly, with a
    listener that really exists: if an empty session id matched every job,
    the session-less caller (the process token, another front end) would be
    handed this session's server and ``decide`` would allow it.
    """
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    code = (f"import socket,time; s=socket.socket(); "
            f"s.bind(('127.0.0.1',{port})); s.listen(); time.sleep(30)")
    job = BJ.get_registry().start(
        command=f"python -c listen-{port}", cwd=str(tmp_path),
        workspace=str(tmp_path), session_id="sess-guard",
        argv=[sys.executable, "-c", code])
    try:
        deadline = time.time() + 10
        while time.time() < deadline and not BJ.own_listening_sockets(
                "sess-guard"):
            time.sleep(0.1)
        # It is the session's while the session is named ...
        assert BJ.own_listening_sockets("sess-guard"), \
            "the running job's listener is not found at all"
        # ... and nobody's when no session is named.
        assert BJ.own_listening_sockets("") == set(), \
            "a caller naming no session was handed a session's listener"
        assert BJ.own_listening_sockets("someone-else") == set()
    finally:
        BJ.get_registry().kill(job.job_id)
        deadline = time.time() + 10
        while job.poll() is None and time.time() < deadline:
            time.sleep(0.05)
