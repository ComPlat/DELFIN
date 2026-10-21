"""A hidden path inside a write root must not take that root's writes away.

Landlock only grants; a rule covers everything beneath its path. To hide
a path inside a write root, ``landlock_exec._grant_except`` grants the
root entry by entry, so the root directory itself loses MAKE_REG and a
new top-level file is refused ("Permission denied"). The approval-code
directory (``approval_code.codes_dir``) is created eagerly and handed to
the cage as ``--hide``; wherever it lay inside a granted root (the suite
puts HOME under the workspace, a session may grant ``~/.delfin``) every
shell write of a new file at the top of that root failed on a host with
Landlock and no usable bubblewrap (GitHub's runners, many HPC systems).

Input: a Landlock cage built by ``api_client._landlock_argv``. Checked:
a new top-level file is written with the code directory inside the
workspace; a code file outside every write root stays unreadable to a
caged command whose path is built at run time (no text gate sees it),
with the uncaged control reading it.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent import approval_code

from conftest import child_env

needs_landlock = pytest.mark.skipif(
    not A._landlock_functional(),
    reason="the kernel applies no Landlock policy here")


def _perms(ws):
    return A.KitToolPermissions(workspace=ws.resolve(),
                                mode="bypassPermissions")


def _hides(argv):
    return [b for a, b in zip(argv, argv[1:]) if a == "--hide"]


@needs_landlock
def test_a_new_top_level_file_is_written_with_the_code_dir_inside(
        tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    codes = ws / "home" / ".delfin" / "approval_codes"
    monkeypatch.setattr(approval_code, "codes_dir", lambda: codes)
    ws.mkdir()
    new = ws / "new.txt"
    argv = A._landlock_argv(
        ["/bin/bash", "-c", f"echo fresh > {new} && cat {new}"], _perms(ws))
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=60,
                          env=child_env(tmp_path))
    assert proc.returncode == 0, proc.stderr
    assert "fresh" in proc.stdout
    assert codes.is_dir(), "the code directory was not created eagerly"


@needs_landlock
def test_a_code_file_outside_the_write_roots_stays_unreadable(
        tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    codes = tmp_path / "codes" / "approval_codes"
    monkeypatch.setattr(approval_code, "codes_dir", lambda: codes)
    assert approval_code.write_code_file("rid1", "ABCDEF") is not None
    argv = A._landlock_argv(
        ["/bin/bash", "-c", f"d=approval; cat {codes.parent}/${{d}}_codes/rid1.code"],
        _perms(ws))
    assert str(codes.resolve()) in _hides(argv), argv
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=60,
                          env=child_env(tmp_path))
    assert "ABCDEF" not in proc.stdout, proc.stdout
    # The control: uncaged, the same command reads it.
    plain = subprocess.run(["/bin/bash", "-c",
                            f"cat {codes}/rid1.code"],
                           capture_output=True, text=True, timeout=60)
    assert "ABCDEF" in plain.stdout, plain.stderr


def test_the_hide_is_dropped_only_for_the_code_dir_inside_a_root(
        tmp_path, monkeypatch):
    """Argv level, every host: the credential hides are kept as they were."""
    ws = tmp_path / "ws"
    codes = ws / "approval_codes"
    monkeypatch.setattr(approval_code, "codes_dir", lambda: codes)
    ws.mkdir()
    hides = A._landlock_hide_paths([str(ws)])
    assert str(codes.resolve()) not in hides
    assert str((A.Path.home() / ".ssh").resolve()) in hides
    assert str(codes.resolve()) in A._landlock_hide_paths(
        [str(tmp_path / "elsewhere")])


# -- where the code directory is ------------------------------------------

#: The shipped resolver, taken at import, before conftest redirects it.
_SHIPPED_CODES_DIR = approval_code.codes_dir


def test_the_code_dir_is_one_place_on_the_home(tmp_path, monkeypatch):
    """~/.delfin/approval_codes, also where a private runtime directory
    exists: that one is local to a node, and a person on another login node
    of a shared home must still be able to answer."""
    home = tmp_path / "home"
    home.mkdir()
    rt = tmp_path / "rt"
    rt.mkdir(mode=0o700)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(rt))
    assert _SHIPPED_CODES_DIR() == home / ".delfin" / "approval_codes"
    assert approval_code.home_codes_dir() == _SHIPPED_CODES_DIR()
