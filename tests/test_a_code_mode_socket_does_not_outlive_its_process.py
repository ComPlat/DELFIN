"""Code Mode sockets were left behind in the shared /tmp.

MEASURED 2026-11-23 on a multi-user host: /tmp/.code_mode_rpc held 18
socket files, some from two days earlier and some still group-writable
from before the socket was narrowed, and /tmp held 78 cm1r_ directories
from three days of test runs -- about forty a day.

Two causes:

  * the RPC server's serve loop is a DAEMON thread. It closes and unlinks
    the socket when its idle window ends -- but a process that exits
    before then kills the thread, close() never runs, and the file stays.
    Every test that starts a server, and every agent session that ends
    soon after a Code Mode command, left one. Sockets still bound when
    the interpreter exits are now unlinked at exit;
  * three test files worked around the 108-byte AF_UNIX limit by creating
    workspaces directly in /tmp, and one passed /tmp ITSELF as the
    workspace, so its socket went into the shared directory every user's
    runs use. The limit is gone (a deep path reaches its socket through
    /proc/self/fd), so the helpers remove what they create, and no test
    hands /tmp itself to Code Mode.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

from delfin.agent import code_mode

TREE = str(pathlib.Path(code_mode.__file__).resolve().parents[2])
TESTS = pathlib.Path(__file__).resolve().parent


def test_a_socket_is_unlinked_when_its_process_exits(tmp_path):
    """A real process binds a server and exits without closing it."""
    probe = (
        "import sys\n"
        f"sys.path.insert(0, {TREE!r})\n"
        "from delfin.agent import code_mode as cm\n"
        f"p = cm.rpc_socket_path({str(tmp_path)!r}, 'exit')\n"
        "cm.RpcServer(object()).bind(p)\n"
        "print(p)\n")
    out = subprocess.run([sys.executable, "-c", probe], capture_output=True,
                         text=True, timeout=60)
    path = out.stdout.strip().splitlines()[-1]
    assert path.endswith("rpc-exit.sock"), out.stderr[-400:]
    assert not pathlib.Path(path).exists(), (
        "the socket outlived the process that bound it")


def test_a_closed_server_is_not_unlinked_twice(tmp_path):
    """close() and the exit handler must not both act on one path: the
    path could by then belong to another server."""
    p = code_mode.rpc_socket_path(str(tmp_path), "twice")
    srv = code_mode.RpcServer(object())
    srv.bind(p)
    srv.close()
    assert p not in code_mode._LIVE_SOCKETS


def test_no_test_hands_tmp_itself_to_code_mode():
    offenders = []
    for f in sorted(TESTS.glob("test_*.py")):
        text = f.read_text(encoding="utf-8", errors="ignore")
        if "_code_mode_child_env" not in text and "rpc_socket_path" not in text:
            continue
        if re.search(r"""_permissions\(\s*["']/tmp["']\s*\)"""
                     r"""|_code_mode_child_env\([^)]*["']/tmp["']""",
                     text, re.S):
            offenders.append(f.name)
    assert not offenders, offenders
