"""The dashboard's kernels get private Unix sockets on any machine.

Input: the Jupyter runtime directory the launcher would use, and the
machine's ``XDG_RUNTIME_DIR`` and ``/tmp``. Before: when the runtime path
was too long for a Unix socket (a cluster home: 119 bytes against 107 on
the KIT cluster, 2026-10-22), the launcher kept loopback TCP and printed a
note. Over TCP any user of a shared login node can connect to the
kernel's ports and read what the kernel prints (measured 2026-10-20,
tests/test_a_dashboard_kernel_is_reachable_by_its_user_only.py). A
runtime directory that other users could read was used as it was,
although it holds the connection files with the kernels' keys.

After: the runtime directory is used only when it is short enough AND
private (a directory, not a symlink, owned by this user, no group or
other permissions). Otherwise the launcher moves it, in this order, to
- ``$XDG_RUNTIME_DIR/delfin-jupyter``, if ``$XDG_RUNTIME_DIR`` is itself
  this user's private directory;
- ``/tmp/delfin-k-<uid>``, made 0700 by this user and the same on every
  start, so a restart keeps the login store that lives in it;
- a new directory made by ``mkdtemp`` under /tmp (mode 0700, a name no
  one can take first), removed when the launcher exits, only when the
  fixed name is taken by something that is not this user's private
  directory;
and exports it as ``JUPYTER_RUNTIME_DIR`` for the server. TCP is left
only when none of them works, with the note.
"""
from __future__ import annotations

import os
import stat

import pytest

from delfin import cli_voila as V


def _private(path: str) -> bool:
    st = os.lstat(path)
    return (stat.S_ISDIR(st.st_mode) and st.st_uid == os.getuid()
            and not st.st_mode & 0o077)


@pytest.fixture
def posix(monkeypatch):
    """POSIX, with the launcher's /tmp moved into a short private directory
    of the test's own: the fixed /tmp/delfin-k-<uid> is the user's real
    runtime directory on a machine whose default one is too long."""
    import shutil
    import tempfile
    monkeypatch.setattr(V.os, "name", "posix")
    monkeypatch.delenv("DELFIN_KERNEL_TRANSPORT", raising=False)
    root = tempfile.mkdtemp(prefix="kt", dir="/tmp")
    monkeypatch.setattr(V, "_TMP_ROOT", root)
    yield root
    shutil.rmtree(root, ignore_errors=True)


@pytest.fixture
def deep(tmp_path, monkeypatch):
    d = tmp_path / ("d" * 90)
    d.mkdir(mode=0o700)
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", str(d))
    return d


def test_a_deep_runtime_moves_to_the_users_xdg_runtime_dir(
        posix, deep, monkeypatch):
    import tempfile
    xdg = tempfile.mkdtemp(prefix="x", dir="/tmp")          # 0700, ours
    monkeypatch.setenv("XDG_RUNTIME_DIR", xdg)
    try:
        assert V._kernel_transport_args() == ["--KernelManager.transport=ipc"]
        chosen = os.environ["JUPYTER_RUNTIME_DIR"]
        assert chosen == os.path.join(xdg, "delfin-jupyter")
        assert _private(chosen)
    finally:
        import shutil
        shutil.rmtree(xdg, ignore_errors=True)


def test_without_a_private_xdg_dir_it_uses_its_own_under_tmp(
        posix, deep, monkeypatch, tmp_path):
    shared = tmp_path / "shared"
    shared.mkdir()
    os.chmod(shared, 0o755)                     # readable by everyone
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(shared))
    assert V._kernel_transport_args() == ["--KernelManager.transport=ipc"]
    chosen = os.environ["JUPYTER_RUNTIME_DIR"]
    assert not chosen.startswith(str(shared))
    assert _private(chosen)
    assert len(os.fsencode(chosen)) < 60


def test_a_runtime_dir_others_can_read_is_not_used(posix, tmp_path,
                                                   monkeypatch):
    open_dir = tmp_path / "rt"
    open_dir.mkdir()
    os.chmod(open_dir, 0o755)
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", str(open_dir))
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    assert V._kernel_transport_args() == ["--KernelManager.transport=ipc"]
    assert os.environ["JUPYTER_RUNTIME_DIR"] != str(open_dir)
    assert _private(os.environ["JUPYTER_RUNTIME_DIR"])


def test_a_symlinked_xdg_subdirectory_is_refused(posix, deep, monkeypatch,
                                                 tmp_path):
    import tempfile
    xdg = tempfile.mkdtemp(prefix="x", dir="/tmp")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir(mode=0o700)
    os.symlink(elsewhere, os.path.join(xdg, "delfin-jupyter"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", xdg)
    try:
        assert V._kernel_transport_args() == ["--KernelManager.transport=ipc"]
        assert not os.environ["JUPYTER_RUNTIME_DIR"].startswith(xdg)
    finally:
        import shutil
        shutil.rmtree(xdg, ignore_errors=True)


def test_only_when_nothing_private_exists_tcp_is_kept_and_said(
        posix, deep, monkeypatch, capsys):
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)

    def no_tmp(*_a, **_k):
        raise OSError("no /tmp")

    monkeypatch.setattr(V.tempfile, "mkdtemp", no_tmp)
    # Nothing on this machine counts as private: not the runtime
    # directory, not /tmp/delfin-k-<uid>.
    monkeypatch.setattr(V, "_is_private_dir", lambda _p: False)
    assert V._kernel_transport_args() == []
    assert "TCP" in capsys.readouterr().err


def test_the_tmp_directory_is_the_same_on_every_start(posix, deep,
                                                     monkeypatch):
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    V._kernel_transport_args()
    first = os.environ["JUPYTER_RUNTIME_DIR"]
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", str(deep))
    V._kernel_transport_args()
    assert os.environ["JUPYTER_RUNTIME_DIR"] == first
    assert _private(first)


def test_a_fixed_name_someone_else_could_hold_is_not_trusted(
        posix, deep, monkeypatch):
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    real = V.os.lstat
    fixed = os.path.join(posix, f"delfin-k-{os.getuid()}")

    class Foreign:
        def __init__(self, st):
            self._st = st

        def __getattr__(self, name):
            return getattr(self._st, name)

        @property
        def st_uid(self):
            return self._st.st_uid + 1          # owned by another user

    monkeypatch.setattr(V.os, "lstat",
                        lambda p: Foreign(real(p)) if p == fixed else real(p))
    V._kernel_transport_args()
    chosen = os.environ["JUPYTER_RUNTIME_DIR"]
    assert chosen != fixed and _private(chosen)


def test_a_tmp_others_may_rename_in_is_not_used(posix, deep, monkeypatch,
                                                capsys):
    """Without the sticky bit any user can rename or replace a directory
    in it; a world-writable /tmp of another owner is no place for the
    kernels' keys, under the fixed name or a mkdtemp one."""
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    real = V.os.lstat

    class Open:
        def __init__(self, st):
            self._st = st

        def __getattr__(self, name):
            return getattr(self._st, name)

        @property
        def st_mode(self):
            return (self._st.st_mode | 0o777) & ~stat.S_ISVTX

        @property
        def st_uid(self):
            return self._st.st_uid + 1

    monkeypatch.setattr(V.os, "lstat",
                        lambda p: Open(real(p)) if p == posix else real(p))
    assert V._kernel_transport_args() == []
    assert "TCP" in capsys.readouterr().err
