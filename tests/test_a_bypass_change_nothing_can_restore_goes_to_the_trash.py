"""In bypass, a command that would destroy what nothing can restore is
refused, and told to move the file to the workspace's trash instead.

Measured in this installation's undo journals (1185 sessions, 2026-11-13):
3362 shell changes to existing files had no undo -- 1676 because the file
was too large for a pre-image, 1686 because the session's 500-entry cap had
dropped it. Nothing asked about any of them in bypass. Restorable means: the
journal can hold it, git holds it (tracked, identical to HEAD), or this
session wrote its current bytes. The journal's per-file limit is lowered
here so that the files stay small.
"""

from __future__ import annotations

import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import shell_changes, unrestorable as U

BIG = b"x" * 400            # over the lowered per-file limit
SMALL = b"y" * 40


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(shell_changes, "MAX_PRE_IMAGE_BYTES", 100)
    monkeypatch.setattr(shell_changes, "MAX_HELD_BYTES", 1000)
    monkeypatch.setattr(A._DocToolExecutor, "_own_written_file",
                        staticmethod(lambda path, perms: None))
    w = tmp_path / "ws"
    w.mkdir()
    (w / "big.dat").write_bytes(BIG)
    (w / "small.txt").write_bytes(SMALL)
    return w


def _perms(ws, mode="bypassPermissions"):
    return A.KitToolPermissions(workspace=ws, mode=mode)


def _git(ws, *args):
    subprocess.run(["git", "-C", str(ws), *args], check=True, capture_output=True,
                   env={"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.org",
                        "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.org",
                        "PATH": "/usr/bin:/bin", "HOME": str(ws)})


@pytest.mark.parametrize("cmd", [
    "rm big.dat", "rm -f ./big.dat", "unlink big.dat", "shred -u big.dat",
    "truncate -s 0 big.dat", "echo new > big.dat", "python3 run.py >big.dat 2>&1",
    "cp small.txt big.dat", "mv small.txt big.dat", "dd if=/dev/zero of=big.dat count=1",
    "rm *.dat", "nice -n 5 rm big.dat", "cd . && rm big.dat",
])
def test_destroying_an_unrestorable_file_is_refused_with_the_trash(ws, cmd):
    why = U.refusal(cmd, ws, _perms(ws))
    assert why and "big.dat" in why
    assert "mv " in why and ".delfin-trash" in why
    assert (ws / "big.dat").read_bytes() == BIG
    assert (ws / ".delfin-trash").is_dir()


@pytest.mark.parametrize("cmd", [
    "rm small.txt", "echo x >> big.dat", "cp big.dat copy.dat", "cat big.dat > /dev/null",
    "ls -l big.dat", "rm missing.dat", "mv big.dat .delfin-trash/", "ln -s big.dat link",
    "find . -name '*.tmp' -delete", "echo x 2>&1 | tee -a log.txt",
])
def test_what_keeps_or_can_restore_the_file_runs(ws, cmd):
    assert U.refusal(cmd, ws, _perms(ws)) is None


def test_a_file_git_holds_may_go(ws):
    _git(ws, "init", "-q")
    _git(ws, "add", "big.dat")
    _git(ws, "commit", "-q", "-m", "c")
    assert U.refusal("rm big.dat", ws, _perms(ws)) is None
    (ws / "big.dat").write_bytes(BIG + b"changed")
    assert U.refusal("rm big.dat", ws, _perms(ws))


def test_the_trash_is_kept_out_of_git(ws):
    _git(ws, "init", "-q")
    U.refusal("rm big.dat", ws, _perms(ws))
    exclude = (ws / ".git" / "info" / "exclude").read_text()
    assert "/.delfin-trash/" in exclude.splitlines()
    U.refusal("rm big.dat", ws, _perms(ws))
    assert (ws / ".git" / "info" / "exclude").read_text().count("/.delfin-trash/") == 1


def test_the_session_s_own_output_may_go(ws, monkeypatch):
    monkeypatch.setattr(A._DocToolExecutor, "_own_written_file",
                        staticmethod(lambda path, perms: path))
    assert U.refusal("rm big.dat", ws, _perms(ws)) is None


def test_a_tree_holding_an_unrestorable_file_is_refused(ws):
    d = ws / "run1"
    d.mkdir()
    (d / "a.txt").write_bytes(SMALL)
    (d / "out.gbw").write_bytes(BIG)
    why = U.refusal("rm -r run1", ws, _perms(ws))
    assert why and "out.gbw" in why
    (d / "out.gbw").unlink()
    assert U.refusal("rm -r run1", ws, _perms(ws)) is None


def test_find_follows_its_name_filter(ws):
    (ws / "x.tmp").write_bytes(BIG)
    assert "x.tmp" in U.refusal("find . -name '*.tmp' -delete", ws, _perms(ws))
    (ws / "x.tmp").unlink()
    assert U.refusal("find . -name '*.tmp' -delete", ws, _perms(ws)) is None
    # A predicate this reading does not follow: the whole tree counts.
    assert U.refusal("find . -mtime +3 -delete", ws, _perms(ws))


def test_the_trash_protects_itself(ws):
    day = U.prepare_trash(ws)
    (ws / "big.dat").rename(day / "big.dat")
    assert U.refusal(f"rm {day}/big.dat", ws, _perms(ws))
    assert U.refusal("rm -r .delfin-trash", ws, _perms(ws))


def test_expired_trash_days_are_removed(ws):
    trash = ws / ".delfin-trash"
    old = (datetime.now(timezone.utc) - timedelta(days=U.TRASH_DAYS + 2)).strftime("%Y-%m-%d")
    (trash / old).mkdir(parents=True)
    (trash / old / "f").write_bytes(BIG)
    (trash / "notes").mkdir()
    U.prepare_trash(ws)
    assert not (trash / old).exists() and (trash / "notes").exists()


def test_the_bypass_gate_refuses_and_the_file_stays(ws, monkeypatch):
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    perms = _perms(ws)
    why = A._DocToolExecutor()._run_permission_gate("bash", {"command": "rm big.dat"}, perms)
    assert why and ".delfin-trash" in why
    assert (ws / "big.dat").exists()


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_the_modes_that_ask_are_not_changed(ws, monkeypatch, mode):
    """There the person reads the command; this rule is bypass's."""
    monkeypatch.setattr(U, "refusal", lambda *a, **k: pytest.fail("called outside bypass"))
    perms = _perms(ws, mode)
    perms.confirm_callback = lambda *a, **k: False
    A._DocToolExecutor()._run_permission_gate("bash", {"command": "rm big.dat"}, perms)
