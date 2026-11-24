"""A command an agent runs that outgrows its memory is the one that dies,
not the agent (2026-11-24: an integrator's full test suite outgrew the
member's limit, and the kill ended the integrator)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from delfin.agent import landlock_exec
from delfin.agent.swarm import wave


def test_the_cage_marks_itself_the_first_victim(tmp_path):
    target = tmp_path / "oom_score_adj"
    target.write_text("0")
    assert landlock_exec.prefer_as_oom_victim(str(target)) is True
    assert target.read_text() == str(landlock_exec.OOM_SCORE_ADJ) == "1000"


def test_without_the_file_nothing_changes(tmp_path):
    assert landlock_exec.prefer_as_oom_victim(str(tmp_path / "no" / "such")) is False


def test_a_caged_command_carries_the_mark():
    """End to end through the cage's entry: the command it runs reads the
    value back. Where the kernel has no such file the cage still runs the
    command, and the command reports that."""
    probe = ("import os; p='/proc/self/oom_score_adj'; "
             "print(open(p).read().strip() if os.path.exists(p) else 'none')")
    out = subprocess.run([sys.executable, "-I", landlock_exec.__file__, "--fs", "0", "--guard", "0",
                          "--net", "open", "--", sys.executable, "-I", "-c", probe],
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    expected = "1000" if Path("/proc/self/oom_score_adj").exists() else "none"
    assert out.stdout.strip() == expected


def test_a_member_scope_lets_the_member_outlive_a_killed_command(monkeypatch):
    monkeypatch.setattr(wave.shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setenv("XDG_RUNTIME_DIR", "/run/user/0")
    argv = wave._scoped(["member"], dict(wave._DEFAULTS))
    assert argv[:3] == ["systemd-run", "--user", "--scope"]
    assert "OOMPolicy=continue" in argv
    assert "MemoryMax=12G" in argv          # room for a full test suite beside the agent
