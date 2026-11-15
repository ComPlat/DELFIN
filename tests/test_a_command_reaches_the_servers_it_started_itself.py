"""A sandboxed command reaches the servers it started itself on the
loopback, and no other server there.

Measured 2026-11-15 in an unattended swarm run: 118 refusals of the egress
proxy ("a POST request to 127.0.0.1 is an upload, and nobody can be
asked"), each a test posting to the scripted model endpoint it served to
itself; the agent could not run the test its task asked for. The loopback
of a shared host also carries other users' services and the user's own
unsandboxed ones, so it is not opened as a whole: a connect passes when
every socket listening there is the command's own -- held by its process
tree, or a listener of the session's background jobs (``--own-listener``).
The socket guard judges a direct connect that way, the proxy a request by
the pid of the helper that forwarded it (``egress_proxy._own_server``).

The tests that start a child need the socket guard and Landlock; the
address parsing and the ownership rule are checked on every Linux host.
"""
from __future__ import annotations

import ipaddress
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

import _egress_test_support as ES
from conftest import child_env
from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP
from delfin.agent import socket_guard as SG

HELPER = ES.HELPER
linux = pytest.mark.skipif(not Path("/proc/net/tcp").exists(), reason="reads /proc/net")
guard = pytest.mark.skipif(not SG.available(), reason="needs seccomp user notification")

LOOPBACK = ipaddress.ip_address("127.0.0.1")


@pytest.fixture
def web():
    yield from ES.serve_web()


@pytest.fixture
def proxy(web, monkeypatch):
    yield from ES.session_proxy(web, monkeypatch)


# -- the rule, on every Linux host ------------------------------------------

def test_a_proc_net_address_is_read_in_host_order():
    if sys.byteorder != "little":
        pytest.skip("the fixture rows are those of a little-endian kernel")
    assert SG._proc_net_address("0100007F:1F90") == (LOOPBACK, 8080)
    assert SG._proc_net_address("00000000000000000000000001000000:0050") == (
        ipaddress.ip_address("::1"), 80)


@linux
def test_a_listener_of_the_callers_own_tree_is_its_own():
    with socket.socket() as srv:
        srv.bind(("127.0.0.1", 0))
        srv.listen()
        port = srv.getsockname()[1]
        inode = str(os.fstat(srv.fileno()).st_ino)
        assert SG.loopback_listeners(LOOPBACK, port) == {inode}
        assert SG.owns_loopback(LOOPBACK, port, (), os.getpid())
        # Named as a job's listener, without the tree.
        assert SG.owns_loopback(LOOPBACK, port, {inode})
        # Below a process that holds nothing, it is somebody else's.
        other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
        try:
            assert not SG.owns_loopback(LOOPBACK, port, (), other.pid)
            # Below the root only: the root's own sockets do not count.
            assert not SG.owns_loopback(LOOPBACK, port, (), os.getpid(),
                                        include_root=False)
        finally:
            other.kill()
            other.wait()
    # Nothing listens there any more.
    assert not SG.owns_loopback(LOOPBACK, port, (), os.getpid())
    assert not SG.owns_loopback(ipaddress.ip_address("10.0.0.5"), port, (), os.getpid())


class _Perms:
    task_session_id = "sess-own"
    mode = "bypassPermissions"
    confirm_callback = None


def test_a_post_to_a_background_jobs_server_is_not_an_upload(monkeypatch):
    """Before, ``decide`` let the session's own port through and ``settle``
    then refused every write to it with nobody to ask."""
    monkeypatch.setattr(BJ, "own_listening_ports",
                        lambda sid: {8393} if sid == "sess-own" else set())
    p = EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))
    try:
        assert p.settle("127.0.0.1", 8393, method="POST", body_bytes=10,
                        perms=_Perms()) is None
        assert p.settle("127.0.0.1", 8394, method="POST", body_bytes=10,
                        perms=_Perms())
        assert p.settle("127.0.0.1", 8393, method="POST", body_bytes=10)
    finally:
        p.close()


# -- a child under the real guard --------------------------------------------

_PROBE = r'''
import http.server, socket, sys, threading, urllib.error, urllib.request
foreign = int(sys.argv[1])
class H(http.server.BaseHTTPRequestHandler):
    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        self.rfile.read(n)
        body = f"got {n}".encode()
        self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def log_message(self, *a):
        pass
srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
threading.Thread(target=srv.serve_forever, daemon=True).start()
own = srv.server_address[1]
def tcp(p):
    s = socket.socket(); s.settimeout(5)
    try:
        s.connect(("127.0.0.1", p)); return "connected"
    except OSError:
        return "refused"
    finally:
        s.close()
def post(p):
    try:
        return urllib.request.urlopen(urllib.request.Request(
            f"http://127.0.0.1:{p}/v1/chat", data=b"x" * 10), timeout=20
        ).read().decode().split()[0]
    except urllib.error.HTTPError as e:
        return f"HTTP{e.code}"
    except Exception:
        return "error"
print(tcp(own), post(own), tcp(foreign), post(foreign))
'''


def _helper_argv(proxy, perms, extra=()):
    return [sys.executable, "-I", HELPER, "--fs", "0", "--net", "proxy",
            "--proxy-socket", proxy.socket_path,
            "--proxy-token-file", proxy.token_file_for(perms),
            "--allow-socket", "/", *extra, "--"]


def _run(argv, ws, tmp_path):
    return subprocess.run(argv, cwd=str(ws), env=child_env(tmp_path),
                          capture_output=True, text=True, timeout=240)


@guard
def test_an_unattended_command_reaches_its_own_server_and_no_other(proxy, web, tmp_path):
    """Its own server directly and through the proxy, a POST included; the
    test's server on the same loopback (not the command's) neither way."""
    if not A._landlock_functional():
        pytest.skip("needs Landlock")
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_PROBE)
    out = _run(_helper_argv(proxy, perms) + [
        sys.executable, str(ws / "probe.py"), str(web.server_address[1])], ws, tmp_path)
    assert out.stdout.split() == ["connected", "got", "refused", "HTTP403"], out
    assert not any(m == "POST" for m, _h, _p in web.seen), web.seen


@guard
def test_the_bash_tools_cage_reaches_the_commands_own_server(proxy, web, tmp_path,
                                                             monkeypatch):
    """The same through the argv the bash tool builds (bubblewrap with the
    helper inside where it works, the helper alone otherwise)."""
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    if not (A._bwrap_functional() or A._landlock_functional()):
        pytest.skip("needs bubblewrap or Landlock for the default isolation")
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_PROBE)
    argv = A._bash_isolation_argv(
        f"{sys.executable} {ws / 'probe.py'} {web.server_address[1]}", ws, perms)
    out = _run(argv, ws, tmp_path)
    assert out.stdout.split() == ["connected", "got", "refused", "HTTP403"], out


@guard
def test_a_background_jobs_listener_is_reached_only_when_named(proxy, tmp_path):
    if not A._landlock_functional():
        pytest.skip("needs Landlock")
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    with socket.socket() as srv:
        srv.bind(("127.0.0.1", 0))
        srv.listen(8)
        port = srv.getsockname()[1]
        inode = str(os.fstat(srv.fileno()).st_ino)
        probe = (f"import socket; s = socket.socket(); s.settimeout(5)\n"
                 f"try:\n    s.connect(('127.0.0.1', {port})); print('connected')\n"
                 f"except OSError:\n    print('refused')\n")
        named = _run(_helper_argv(proxy, perms, ["--own-listener", inode])
                     + [sys.executable, "-c", probe], ws, tmp_path)
        unnamed = _run(_helper_argv(proxy, perms) + [sys.executable, "-c", probe],
                       ws, tmp_path)
    assert named.stdout.split() == ["connected"], named
    assert unnamed.stdout.split() == ["refused"], unnamed


@linux
def test_the_bash_tool_names_the_sessions_job_listeners(monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(BJ, "own_listening_sockets",
                        lambda sid: {"4242"} if sid == "sess-own" else set())
    args, _dir = A._sandbox_network_args(_Perms(), allow_slurm=False)
    assert args[args.index("--own-listener") + 1] == "4242", args


@linux
def test_a_jobs_listener_inodes_are_read_from_its_processes():
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    code = (f"import socket,time; s=socket.socket(); "
            f"s.bind(('127.0.0.1',{port})); s.listen(); time.sleep(30)")
    child = subprocess.Popen([sys.executable, "-c", code])
    try:
        for _ in range(100):
            found = BJ._listening_of(BJ._descendants(child.pid))
            if port in found.values():
                break
            time.sleep(0.05)
        assert port in found.values(), found
        assert BJ._listening_ports_of(BJ._descendants(child.pid)) == {port}
    finally:
        child.kill()
        child.wait()
