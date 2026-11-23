"""The Code Mode socket was reachable by the group, and it can start calculations.

MEASURED 2026-11-23 on the shared host (umask 002, seven other real
users): the RPC directory was created ``drwxrwxr-x`` and the socket
``srwxrwxr-x``, and the server accepted any connection -- no peer
credential was checked and no token was carried. The only thing keeping
another user out was that the operator's home directory happens to be
``drwx------``. A workspace in a shared group directory would have put a
group-connectable socket there.

Before 2026-11-23 a program could reach read-only tools through it,
which already meant reading this user's files. Since the change that
let a program run ``run_application`` and ``submit_application``, it
means starting calculations under this user's identity. The stakes went
up and the door was never checked; this file checks it.

THREE LAYERS, each enough alone on a normal host, so that one mistake
does not open it:

  * the directory is ``0700`` and must be a real directory OWNED by this
    user -- a symlink or a foreign directory planted at that name in a
    shared workspace is refused rather than used;
  * the socket is ``0600``;
  * the server asks the kernel who connected (``SO_PEERCRED``) and
    answers only its own uid. Measured through a real bubblewrap cage
    before this was written: the caged child's uid as the server sees it
    is the server's own, so the check does not shut the cage out.

AND A PATH THAT WAS TOO LONG. ``AF_UNIX`` paths are limited to 108
bytes; measured, Code Mode failed with "AF_UNIX path too long" for any
workspace path of 78 characters or more. Both ends now reach the socket
through ``/proc/self/fd/<dir>/<name>`` when the full path is too long.
"""
from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
import sys
import threading

import pytest

from delfin.agent import code_mode as cm

_LINUX_PROC = os.path.isdir("/proc/self/fd")


class _Echo:
    """An executor that answers every call, so a test can see a call land."""

    def execute(self, name, arguments, permissions=None):
        return json.dumps(f"ran {name}")


#: A tool the Code Mode allow-list admits, so a control case is refused
#: only by what it tests and never by the allow-list. An invented name
#: made two control cases red before the fix for the wrong reason.
_ALLOWED = "list_docs"


def _serve_once(path):
    srv = cm.RpcServer(_Echo())
    srv.bind(path)
    t = threading.Thread(target=srv.serve_one, kwargs={"idle_timeout_s": 10},
                         daemon=True)
    t.start()
    return srv, t


@pytest.fixture()
def group_umask():
    """The host's own umask, which is what made the socket group-writable."""
    old = os.umask(0o002)
    try:
        yield
    finally:
        os.umask(old)


# --- the modes ------------------------------------------------------

def test_the_directory_and_the_socket_are_the_owners_alone(tmp_path,
                                                           group_umask):
    path = cm.rpc_socket_path(str(tmp_path), "probe")
    srv, t = _serve_once(path)
    try:
        d = os.stat(os.path.dirname(path)).st_mode
        s = os.stat(path).st_mode
        assert stat.S_IMODE(d) == 0o700, stat.filemode(d)
        assert stat.S_IMODE(s) & 0o077 == 0, stat.filemode(s)
    finally:
        srv.close()


def test_an_existing_group_directory_is_tightened(tmp_path, group_umask):
    """A directory left behind by an earlier run under the old code."""
    base = tmp_path / ".code_mode_rpc"
    base.mkdir()
    base.chmod(0o775)
    path = cm.rpc_socket_path(str(tmp_path), "probe")
    assert stat.S_IMODE(os.stat(os.path.dirname(path)).st_mode) == 0o700


def test_a_symlink_planted_at_the_name_is_refused(tmp_path):
    """In a shared directory another user can create the name first and
    point it somewhere they control. The server must not follow it."""
    elsewhere = tmp_path / "theirs"
    elsewhere.mkdir()
    (tmp_path / ".code_mode_rpc").symlink_to(elsewhere)
    with pytest.raises(PermissionError):
        cm.rpc_socket_path(str(tmp_path), "probe")


# --- who connected ----------------------------------------------------

def test_the_owner_is_answered(tmp_path):
    path = cm.rpc_socket_path(str(tmp_path), "probe")
    srv, t = _serve_once(path)
    try:
        assert cm.RpcClient(path).call(_ALLOWED) == f"ran {_ALLOWED}"
    finally:
        srv.close()


def test_another_uid_is_not_answered(tmp_path, monkeypatch):
    """The decision, driven through a real connection with the kernel's
    answer replaced -- a second real uid needs root, and a test that
    only passes as root measures nothing on a normal host."""
    monkeypatch.setattr(cm, "_peer_uid", lambda conn: os.getuid() + 1)
    path = cm.rpc_socket_path(str(tmp_path), "probe")
    executed = []

    class _Recording(_Echo):
        def execute(self, name, arguments, permissions=None):
            executed.append(name)
            return super().execute(name, arguments, permissions)

    srv = cm.RpcServer(_Recording())
    srv.bind(path)
    t = threading.Thread(target=srv.serve_one, kwargs={"idle_timeout_s": 10},
                         daemon=True)
    t.start()
    try:
        with pytest.raises(RuntimeError, match="not the owner"):
            cm.RpcClient(path).call(_ALLOWED)
        assert executed == [], "a foreign caller reached the executor"
    finally:
        srv.close()


def test_the_check_is_on_the_path_that_serves():
    """Not a helper that nothing calls (rule 2 of the house brief)."""
    import inspect
    assert "_peer_uid(" in inspect.getsource(cm.RpcServer.serve_one)


# --- the length that broke it ----------------------------------------

@pytest.mark.skipif(not _LINUX_PROC, reason="the long-path route uses "
                    "/proc/self/fd; elsewhere a long path fails as before")
def test_a_deep_workspace_still_connects(tmp_path):
    deep = tmp_path
    while len(str(deep)) < 140:
        deep = deep / "a_rather_long_directory_name"
    deep.mkdir(parents=True)
    path = cm.rpc_socket_path(str(deep), "probe")
    assert len(path.encode()) > 108, "the case under test is not long"
    srv, t = _serve_once(path)
    try:
        assert cm.RpcClient(path).call(_ALLOWED) == f"ran {_ALLOWED}"
    finally:
        srv.close()


# --- and through the real cage ---------------------------------------

@pytest.mark.skipif(shutil.which("bwrap") is None, reason="no bubblewrap")
def test_a_caged_child_is_still_answered(tmp_path):
    """The check must not shut out the program it exists to serve."""
    path = cm.rpc_socket_path(str(tmp_path), "probe")
    srv, t = _serve_once(path)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(cm.__file__)))
    repo = os.path.dirname(repo)
    code = ("import sys; sys.path.insert(0, %r); "
            "from delfin.agent import code_mode as cm; "
            "print(cm.RpcClient(%r).call(%r))" % (repo, path, _ALLOWED))
    try:
        out = subprocess.run(
            ["bwrap", "--ro-bind", "/", "/", "--dev", "/dev",
             "--proc", "/proc", "--unshare-pid",
             "--bind", str(tmp_path), str(tmp_path),
             sys.executable, "-c", code],
            capture_output=True, text=True, timeout=60)
        if "Operation not permitted" in out.stderr or \
                "setting up uid map" in out.stderr:
            pytest.skip("bubblewrap cannot create a namespace here")
        assert f"ran {_ALLOWED}" in out.stdout, out.stderr[-400:]
    finally:
        srv.close()
