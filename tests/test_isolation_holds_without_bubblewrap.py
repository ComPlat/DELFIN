"""Filesystem isolation holds on a host where bubblewrap cannot run.

The audit of 2026-09-16 found that every filesystem protection of the
agent's shell rested on bubblewrap: without it a locked session and an
unattended run fell back to plain bash, the first with a panel note, the
second silently. Landlock (Linux 5.13+) is the second way; where neither
exists a locked session refuses the command. These tests switch bubblewrap
off and drive the real tools; they skip where the kernel has no Landlock.
"""
import json

import pytest

from delfin.agent import api_client as A

landlock = pytest.mark.skipif(not A._landlock_functional(),
                              reason="this kernel offers no Landlock")


@pytest.fixture
def host_without_bwrap(tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".ssh").mkdir(parents=True)
    (home / ".ssh" / "id_ed25519").write_text("PRIVATE-KEY-PROBE\n")
    (home / "notes.txt").write_text("readable\n")
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    return tmp_path, home, ws


@landlock
@pytest.mark.parametrize("profile", [dict(lock_workspace=True, mode="default"),
                                     dict(mode="bypassPermissions")])
def test_the_shell_cannot_read_a_key_or_write_outside(host_without_bwrap, profile):
    tmp, home, ws = host_without_bwrap
    perms = A.KitToolPermissions(workspace=str(ws), **profile)
    cmd = (f"cat {home}/.ssh/id_ed25519; cat {home}/notes.txt; "
           f"echo x > {tmp}/outside.txt; echo inside > {ws}/in.txt; echo done")
    argv = A._bash_isolation_argv(cmd, ws, perms)
    from delfin.agent import contained_run
    out = contained_run.run(argv, cwd=str(ws), env=A._scrubbed_bash_env(), timeout=60)
    assert "PRIVATE-KEY-PROBE" not in out.stdout
    assert "readable" in out.stdout and "done" in out.stdout
    assert not (tmp / "outside.txt").exists()
    assert (ws / "in.txt").read_text() == "inside\n"


@landlock
def test_a_test_run_still_reports_under_landlock(host_without_bwrap):
    tmp, home, ws = host_without_bwrap
    (ws / "test_probe.py").write_text(
        "import os\n"
        "def test_key_hidden():\n"
        f"    try:\n        open({str(home / '.ssh' / 'id_ed25519')!r}).read()\n"
        "    except PermissionError:\n        return\n"
        "    raise AssertionError('key readable')\n")
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    out = json.loads(A._doc_executor._execute_run_tests(
        {"target": "test_probe.py", "timeout_s": 120}, perms))
    assert out["status"] == "ok", out
    assert out["summary"].get("passed") == 1, out


def test_a_locked_session_with_no_isolation_at_all_refuses(host_without_bwrap, monkeypatch):
    tmp, home, ws = host_without_bwrap
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setattr(A, "_ISOLATION_GAP_ANNOUNCED", True)
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    argv = A._bash_isolation_argv(f"echo x > {tmp}/outside.txt", ws, perms)
    from delfin.agent import contained_run
    out = contained_run.run(argv, cwd=str(ws), timeout=30)
    assert out.returncode == 126 and "refused" in out.stderr
    assert not (tmp / "outside.txt").exists()
