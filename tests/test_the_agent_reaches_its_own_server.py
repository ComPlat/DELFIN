"""A caged command reaches the server its own session started, and nothing
else on this machine.

Session 0db62cdf (2026-10-19): `curl http://127.0.0.1:8324/` from the cage
was refused ("port 8324 is not allowed") -- the agent could not look at the
game it was serving, and asked the user to test every change. The egress
proxy now admits, for a session's token, the loopback ports that session's
running background jobs listen on (bash_jobs.own_listening_ports, read from
/proc), and connects to the loopback itself. Other ports, other sessions'
servers and the process token stay refused.
"""

from __future__ import annotations

import socket
import sys
import time

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP

linux = pytest.mark.skipif(not sys.platform.startswith("linux"),
                           reason="/proc is Linux")


class _Perms:
    task_session_id = "sess-own"
    mode = "bypassPermissions"


@pytest.fixture
def proxy(monkeypatch):
    monkeypatch.setattr(BJ, "own_listening_ports",
                        lambda sid: {8393} if sid == "sess-own" else set())
    return EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))


@pytest.mark.parametrize("host", ["127.0.0.1", "localhost", "::1",
                                  "[::1]", "127.0.1.1"])
def test_the_sessions_own_port_on_the_loopback_is_allowed(proxy, host):
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


@linux
def test_the_ports_are_read_from_the_jobs_processes(tmp_path):
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
        while time.time() < deadline and \
                port not in BJ.own_listening_ports("sess-proc"):
            time.sleep(0.1)
        assert port in BJ.own_listening_ports("sess-proc")
        assert port not in BJ.own_listening_ports("someone-else")
        assert BJ.own_listening_ports("") == set()
    finally:
        BJ.get_registry().kill(job.job_id)
        deadline = time.time() + 10
        while job.poll() is None and time.time() < deadline:
            time.sleep(0.05)
    assert port not in BJ.own_listening_ports("sess-proc")
