"""A sandboxed command cannot reach the user's sessions through a socket.

Driven on 2026-09-16: under the Landlock helper a command connected to a
Unix socket outside its folders and sent ``send-keys``; on the same host
``tmux ls`` listed the user's real sessions from inside the sandbox.
Landlock does not govern connect(); bubblewrap hides those sockets, and a
host without bubblewrap had nothing. These tests drive the helper for real
and skip where the kernel lacks Landlock or seccomp user notification.
"""
import os
import socket
import subprocess
import sys
import threading

import pytest

from delfin.agent import api_client as A
from delfin.agent import socket_guard as SG

HELPER = str(__import__("pathlib").Path(A.__file__).with_name("landlock_exec.py"))
needs = pytest.mark.skipif(not (A._landlock_functional() and SG.available()),
                           reason="needs Landlock and seccomp user notification")


def _server(path):
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    srv = socket.socket(socket.AF_UNIX)
    srv.bind(str(path))
    srv.listen(64)
    srv.settimeout(60)
    seen = []

    def loop():
        while True:
            try:
                conn, _ = srv.accept()
            except OSError:
                return
            try:
                seen.append(conn.recv(64))
                conn.sendall(b"pong")
            except OSError:
                pass            # a client that closes at once
            finally:
                conn.close()

    threading.Thread(target=loop, daemon=True).start()
    return seen


_PROBE = r'''
import os, socket, sys
def attempt(p):
    s = socket.socket(socket.AF_UNIX)
    try:
        s.connect(p); s.send(b"send-keys"); return s.recv(8).decode()
    except OSError as e:
        return "refused"
for p in sys.argv[1:]:
    print(attempt(p))
try:
    socket.socket(socket.AF_UNIX).connect("\0abstract-door"); print("abstract-connected")
except OSError:
    print("abstract-refused")
try:
    socket.socket(socket.AF_UNIX, socket.SOCK_DGRAM); print("dgram-created")
except OSError:
    print("dgram-refused")
'''


@pytest.fixture
def scene(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "sessions"
    outside.mkdir()
    (ws / "link.sock").symlink_to(outside / "door.sock")
    probe = tmp_path / "probe.py"
    probe.write_text(_PROBE)
    return ws, outside, probe


def _run(ws, probe, *targets):
    return subprocess.run(
        [sys.executable, "-I", HELPER, "--write", str(ws), "--tmpdir", str(ws), "--",
         sys.executable, str(probe), *map(str, targets)],
        capture_output=True, text=True, timeout=90)


@needs
def test_a_socket_outside_the_folders_is_refused(scene):
    ws, outside, probe = scene
    door = _server(outside / "door.sock")
    own = _server(ws / "own.sock")
    out = _run(ws, probe, outside / "door.sock", ws / "link.sock", ws / "own.sock")
    assert out.stdout.split() == ["refused", "refused", "pong",
                                  "abstract-refused", "dgram-refused"], out
    assert door == [] and own == [b"send-keys"]


_RACE = r'''
import ctypes, socket, struct, sys, threading, time
libc = ctypes.CDLL(None, use_errno=True)
good, bad = sys.argv[1].encode(), sys.argv[2].encode()
n = max(len(good), len(bad))
pack = lambda p: struct.pack("<H", socket.AF_UNIX) + p.ljust(n, b"\0") + b"\0"
buf = ctypes.create_string_buffer(pack(good), len(pack(good)))
stop = []
def flip():
    a, b = pack(good), pack(bad)
    while not stop:
        ctypes.memmove(buf, b, len(b)); ctypes.memmove(buf, a, len(a))
threading.Thread(target=flip, daemon=True).start()
end = time.time() + 2
while time.time() < end:
    s = socket.socket(socket.AF_UNIX); libc.connect(s.fileno(), buf, len(pack(good))); s.close()
stop.append(1)
'''


@needs
def test_rewriting_the_address_after_the_check_changes_nothing(scene):
    ws, outside, _ = scene
    door = _server(outside / "door.sock")
    own = _server(ws / "own.sock")
    race = ws.parent / "race.py"
    race.write_text(_RACE)
    subprocess.run([sys.executable, "-I", HELPER, "--write", str(ws), "--",
                    sys.executable, str(race), str(ws / "own.sock"), str(outside / "door.sock")],
                   capture_output=True, text=True, timeout=90)
    assert len(own) > 0
    assert door == [], "a connect reached the forbidden socket"


@needs
def test_tcp_and_the_exit_code_are_untouched(scene):
    ws, _, _ = scene
    srv = socket.socket()
    srv.bind(("127.0.0.1", 0))
    srv.listen(1)
    port = srv.getsockname()[1]
    out = subprocess.run(
        [sys.executable, "-I", HELPER, "--write", str(ws), "--", "/bin/bash", "-c",
         f"exec 3<>/dev/tcp/127.0.0.1/{port} && echo tcp-ok; exit 7"],
        capture_output=True, text=True, timeout=60)
    srv.close()
    assert "tcp-ok" in out.stdout and out.returncode == 7


def test_the_locked_session_asks_for_the_strict_guard(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    locked = A.KitToolPermissions(workspace=str(tmp_path), lock_workspace=True)
    argv = A._bash_isolation_argv("ls", tmp_path, locked)
    assert argv[argv.index("--strict") + 1] == "1"
    bypass = A.KitToolPermissions(workspace=str(tmp_path), mode="bypassPermissions")
    argv = A._bash_isolation_argv("ls", tmp_path, bypass, mode="auto")
    assert argv[argv.index("--strict") + 1] == "0"


def test_the_filter_program_decides_as_intended():
    if SG._arch() is None:
        pytest.skip("architecture without a filter")
    prog = SG._program()
    audit, _sc, nr_connect, nr_socket, nr_uring, _po, _pg = SG._arch()

    def run(nr, arch=audit, a0=0, a1=0):
        data = {0: nr, 4: arch, 16: a0, 24: a1}
        acc, i = 0, 0
        while True:
            code, jt, jf, k = prog[i]
            if code == 0x20:
                acc, i = data.get(k, 0), i + 1
            elif code == 0x54:
                acc, i = acc & k, i + 1
            elif code == 0x15:
                i += 1 + (jt if acc == k else jf)
            elif code == 0x45:
                i += 1 + (jt if acc & k else jf)
            else:
                return k

    assert run(nr_connect) == SG._RET_USER_NOTIF
    assert run(0) == SG._RET_ALLOW
    assert run(nr_socket, a0=socket.AF_INET, a1=1) == SG._RET_ALLOW
    assert run(nr_socket, a0=socket.AF_UNIX, a1=1) == SG._RET_ALLOW
    assert run(nr_socket, a0=socket.AF_UNIX, a1=2) & 0xFFFF == 13
    assert run(nr_uring) & 0xFFFF == 38
    assert run(0, arch=0x40000003) & 0xFFFF == 1
