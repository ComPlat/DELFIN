"""A run can hide paths from its agents: DELFIN_CAGE_HIDE.

Input: ``DELFIN_CAGE_HIDE``, absolute paths joined by ``os.pathsep``, set
by whoever starts the session (a swarm launcher) -- e.g. a lockbox of
held-out data that an analysis agent must never see while it may read the
folders around it. Output: every existing path in it is hidden exactly
like the home secret directories -- a tmpfs (directory) or /dev/null
(file) under bubblewrap, a Landlock hide, a firejail blacklist -- and the
file tools refuse it even inside a readable root. Semantics: hiding only
ever removes access; relative, empty or non-existent entries are ignored;
the variable is read by the session that builds the cage, so a command
inside it cannot change what its own cage hides.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent import sandbox


@pytest.fixture
def lockbox(tmp_path, monkeypatch):
    # Not under /tmp: every cage mounts a fresh /tmp, so a lockbox there is
    # invisible with or without this feature and the real-run tests would
    # pass on code that hides nothing. /var/tmp is shown read-only like the
    # rest of the host, as a real lockbox in a home directory would be.
    base = Path(tempfile.mkdtemp(prefix="delfin-hide-", dir="/var/tmp"))
    box = base / "LOCKBOX2"
    box.mkdir()
    (box / "held_out.txt").write_text("secret-frame\n")
    single = base / "stage1.csv"
    single.write_text("refcode\n")
    monkeypatch.setenv("DELFIN_CAGE_HIDE", os.pathsep.join(
        [str(box), str(single), "relative/ignored", "", str(base / "missing")]))
    yield box, single
    shutil.rmtree(base, ignore_errors=True)


def test_the_variable_names_existing_absolute_paths_only(lockbox):
    box, single = lockbox
    assert sandbox.run_hidden_paths() == [str(box.resolve()), str(single.resolve())]


def test_unset_hides_nothing(monkeypatch):
    monkeypatch.delenv("DELFIN_CAGE_HIDE", raising=False)
    assert sandbox.run_hidden_paths() == []


def test_landlock_and_seatbelt_hide_them(lockbox):
    box, single = lockbox
    hides = A._home_secret_paths()
    assert str(box.resolve()) in hides
    assert str(single.resolve()) in hides


def test_both_bwrap_builders_mask_them(lockbox, tmp_path):
    box, single = lockbox
    ws = tmp_path / "ws"
    ws.mkdir()
    argv = sandbox._bwrap_argv("true", ws, allow_network=False)
    joined = list(zip(argv, argv[1:], argv[2:]))
    assert ("--tmpfs", str(box.resolve())) in list(zip(argv, argv[1:]))
    assert ("--ro-bind", "/dev/null", str(single.resolve())) in joined
    # The shell path builds the cage this host can run: bubblewrap where it
    # works, Landlock where it does not (a CI runner). Either way the
    # hidden paths are hidden, in that layer's own form.
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    argv2 = A._bash_isolation_argv("true", str(ws), perms, mode="bwrap")
    pairs = list(zip(argv2, argv2[1:]))
    if argv2 and Path(argv2[0]).name == "bwrap":
        assert ("--tmpfs", str(box.resolve())) in pairs
        assert ("--ro-bind", "/dev/null", str(single.resolve())) in list(
            zip(argv2, argv2[1:], argv2[2:]))
    else:
        assert ("--hide", str(box.resolve())) in pairs, argv2[:12]
        assert ("--hide", str(single.resolve())) in pairs, argv2[:12]


def test_firejail_blacklists_them(lockbox, tmp_path):
    box, single = lockbox
    argv = sandbox._firejail_argv("true", tmp_path, allow_network=False)
    assert f"--blacklist={box.resolve()}" in argv
    assert f"--blacklist={single.resolve()}" in argv


def test_the_file_tools_refuse_them_inside_a_readable_root(lockbox):
    box, _single = lockbox
    # The whole folder around the lockbox is the workspace: readable, yet
    # the hidden part stays refused.
    around = box.parent
    perms = A.KitToolPermissions(workspace=around, mode="default")
    assert A._is_secret_path(perms, box / "held_out.txt")
    err = A._DocToolExecutor()._check_read_access(perms, box / "held_out.txt")
    assert err is not None
    # A neighbour is still readable.
    (around / "audit.md").write_text("ok\n")
    assert A._DocToolExecutor()._check_read_access(perms, around / "audit.md") is None


@requires_a_working_sandbox
def test_a_caged_command_cannot_read_the_hidden_file(lockbox, tmp_path, monkeypatch):
    box, single = lockbox
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("DELFIN_AGENT_SANDBOX", "bwrap")
    ws = tmp_path / "ws"
    ws.mkdir()
    res = sandbox.run_agent_command(
        f"cat {box}/held_out.txt; cat {single}; ls -A {box} | wc -l", ws)
    assert "secret-frame" not in res.stdout
    assert "refcode" not in res.stdout
    assert res.stdout.strip().splitlines()[-1] == "0"


@requires_a_working_sandbox
def test_a_script_the_agent_writes_cannot_read_it_either(lockbox, tmp_path, monkeypatch):
    # The case the read scan cannot see: a program opening the file itself.
    box, _single = lockbox
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("DELFIN_AGENT_SANDBOX", "bwrap")
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "peek.py").write_text(
        f"import pathlib\np = pathlib.Path({str(box / 'held_out.txt')!r})\n"
        "print('READ:' + p.read_text() if p.exists() else 'ABSENT')\n")
    res = sandbox.run_agent_command("python peek.py", ws)
    assert "secret-frame" not in res.stdout
    assert "ABSENT" in res.stdout
