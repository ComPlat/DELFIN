"""Command isolation holds on macOS, through Seatbelt.

macOS has neither bubblewrap nor Landlock nor seccomp, so a locked session
there refused every shell command and an unattended run had no isolation.
Seatbelt (sandbox-exec) now gives the same promise. The profile text is
checked everywhere; the behaviour is driven on a macOS runner
(.github/workflows/sandbox-macos.yml) and skipped elsewhere.
"""
import socket
import sys
import threading

import pytest

from delfin.agent import seatbelt as SB


def test_the_profile_says_what_it_allows():
    prof = SB.profile(write_roots=["/w/ork space"], hide=['/h/.ss"h'],
                      allow_sockets=["/w/ork space"], net_mode="proxy", proxy_port=4321)
    assert "(deny file-write*)" in prof
    assert '(allow file-write* (subpath "/w/ork space")' in prof
    assert '(subpath "/h/.ss\\"h")' in prof            # quotes escaped
    assert "(deny network-outbound (remote unix-socket))" in prof
    assert '(allow network-outbound (remote ip "localhost:4321"))' in prof
    attended = SB.profile(restrict_files=False, deny_sockets=["/tmp/tmux-501", "regex:^/x\\."])
    assert "file-write" not in attended and "(remote unix-socket (regex" in attended


macos = pytest.mark.skipif(sys.platform != "darwin" or not SB.available(),
                           reason="Seatbelt runs on macOS only")


@pytest.fixture
def scene(monkeypatch):
    # A short base: macOS allows 104 bytes for a Unix socket path, and
    # pytest's temp directory under /private/var/folders is longer.
    import pathlib
    import shutil
    import tempfile
    from delfin.agent import api_client as A
    tmp_path = pathlib.Path(tempfile.mkdtemp(prefix="dsb", dir="/tmp")).resolve()
    yield from _scene(A, tmp_path, monkeypatch)
    shutil.rmtree(tmp_path, ignore_errors=True)


def _scene(A, tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "id_ed25519").write_text("PRIVATE-KEY-PROBE\n")
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    yield A, tmp_path, home, ws


def _run(A, argv, ws):
    from delfin.agent import contained_run
    return contained_run.run(argv, cwd=str(ws), env=A._scrubbed_bash_env(), timeout=90)


@macos
def test_a_locked_session_is_confined(scene):
    A, tmp, home, ws = scene
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    cmd = (f"cat {home}/.ssh/id_ed25519; echo x > {tmp}/outside.txt; "
           f"echo inside > {ws}/in.txt; echo done")
    out = _run(A, A._bash_isolation_argv(cmd, ws, perms), ws)
    assert "PRIVATE-KEY-PROBE" not in out.stdout and "done" in out.stdout, out
    assert not (tmp / "outside.txt").exists()
    assert (ws / "in.txt").read_text() == "inside\n"


_SOCK_PROBE = r'''
import socket, sys
for p in sys.argv[1:]:
    s = socket.socket(socket.AF_UNIX)
    try:
        s.connect(p); print("connected")
    except OSError:
        print("refused")
'''


def _unix_server(path):
    srv = socket.socket(socket.AF_UNIX)
    srv.bind(str(path))
    srv.listen(8)
    srv.settimeout(60)
    threading.Thread(target=lambda: [srv.accept() for _ in range(8)], daemon=True).start()
    return srv


@macos
def test_sockets_outside_the_workspace_are_refused(scene):
    A, tmp, home, ws = scene
    outside = tmp / "door.sock"
    inside = ws / "own.sock"
    _unix_server(outside)
    _unix_server(inside)
    probe = tmp / "probe.py"
    probe.write_text(_SOCK_PROBE)
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    cmd = f"{sys.executable} {probe} {outside} {inside}"
    out = _run(A, A._bash_isolation_argv(cmd, ws, perms), ws)
    assert out.stdout.split() == ["refused", "connected"], out


@macos
def test_the_network_goes_through_the_proxy_only(scene, monkeypatch):
    A, tmp, home, ws = scene
    import http.server
    from delfin.agent import egress_proxy as EP

    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"hello")

        def log_message(self, *a):
            pass

    web = http.server.HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=web.serve_forever, daemon=True).start()
    port = web.server_address[1]
    proxy = EP.EgressProxy(allowed_domains=["allowed.test"], allowed_ports=[port],
                           connector=lambda h, p, t: socket.create_connection(("127.0.0.1", p), timeout=t))
    monkeypatch.setattr(EP, "_PROXY", proxy)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    probe = tmp / "net.py"
    probe.write_text(
        "import socket, sys, urllib.request, urllib.error\n"
        "port = int(sys.argv[1])\n"
        "def get(u):\n"
        "    try:\n        return urllib.request.urlopen(u, timeout=10).read().decode()\n"
        "    except urllib.error.HTTPError as e:\n        return f'HTTP{e.code}'\n"
        "    except Exception:\n        return 'error'\n"
        "s = socket.socket()\n"
        "try:\n    s.connect(('127.0.0.1', port)); direct = 'connected'\n"
        "except OSError:\n    direct = 'refused'\n"
        "print(get(f'http://allowed.test:{port}/'), get(f'http://evil.test:{port}/'), direct)\n")
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    out = _run(A, A._bash_isolation_argv(f"{sys.executable} {probe} {port}", ws, perms), ws)
    proxy.close()
    web.shutdown()
    assert out.stdout.split() == ["hello", "HTTP403", "refused"], out


@macos
def test_an_attended_command_cannot_reach_a_session_door(scene, monkeypatch):
    """The door is shut, and the workspace is still reachable.

    The second socket used to sit OUTSIDE the workspace and was expected
    to connect: "auto" isolated under bypass and for a locked scope and
    nowhere else, so an attended session ran under the socket-only guard,
    which denies the session doors and nothing more. b01f6bd9 (2026-09-19)
    made the posture the strongest the host can give in every mode, and
    on a Seatbelt host that is the full profile -- every unix socket
    outside the workspace roots is denied, the doors among them.

    The intent of this test is the door. What it must NOT become is
    "everything is refused", which would pass just as well if the sandbox
    denied the workspace too and the agent could do nothing at all. So
    the second socket moves INSIDE the workspace: the door stays shut,
    and the half of the promise that has to survive is pinned beside it.
    """
    A, tmp, home, ws = scene
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    door_dir = tmp / "sessions"
    door_dir.mkdir()
    door = door_dir / "default"
    other = ws / "other.sock"
    _unix_server(door)
    _unix_server(other)
    monkeypatch.setattr(SB, "session_doors", lambda: [str(door_dir)])
    probe = tmp / "probe.py"
    probe.write_text(_SOCK_PROBE)
    perms = A.KitToolPermissions(workspace=str(ws), mode="acceptEdits")
    argv = A._bash_isolation_argv(f"{sys.executable} {probe} {door} {other}", ws, perms, mode="auto")
    out = _run(A, argv, ws)
    assert out.stdout.split() == ["refused", "connected"], out
