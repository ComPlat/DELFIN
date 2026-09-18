"""An agent is not a stranger to its own writes.

``edit_file`` refuses a file whose mtime moved since the last
``read_file``: somebody changed it under you, re-read before editing.
The file tools keep that baseline current when they write. The shell did
not — so an agent that appends with

    cat >> delfin/agent/repl_box.py <<'EOF'

was told by its very next ``edit_file`` that the file "was modified since
last read_file". True, and it reads as a warning about somebody else's
change; the somebody was the agent, one call earlier, through a tool that
passed the same gates.

Measured 2026-09-18: three such refusals in one session, each directly
after its own heredoc.

What the guard still covers is the case it exists for: a change this
session did not make. That half is asserted below, because removing it
would be the worse bug.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture()
def bash(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    engine = _DocToolExecutor.__new__(_DocToolExecutor)
    engine._permissions = perms

    def _run(command: str):
        return json.loads(engine._execute_bash(
            {"command": command, "description": "t"}, perms))
    return perms, tmp_path, _run


def _read_baseline(perms, path):
    """Exactly what read_file leaves behind: the mtime as it is now.

    Not an older value -- a stale baseline is refused before the command
    runs, by the write gate, which is a different guard and not the one
    under test here.
    """
    perms.read_tracker[str(path.resolve())] = path.stat().st_mtime


def test_appending_through_the_shell_moves_the_baseline(bash):
    perms, ws, run = bash
    target = ws / "notes.py"
    target.write_text("one\n", encoding="utf-8")
    _read_baseline(perms, target)

    out = run(f"cat >> {target} <<'EOF'\ntwo\nEOF")
    assert out.get("exit_code") == 0

    assert perms.read_tracker[str(target.resolve())] >= target.stat().st_mtime, (
        "the next edit_file would refuse a file this agent just wrote")


def test_a_plain_redirect_moves_it_too(bash):
    perms, ws, run = bash
    target = ws / "log.txt"
    target.write_text("old\n", encoding="utf-8")
    _read_baseline(perms, target)

    run(f"echo new > {target}")

    assert perms.read_tracker[str(target.resolve())] >= target.stat().st_mtime


def test_a_file_the_command_did_not_write_keeps_its_baseline(bash):
    """The guard is for changes this session did not make. That is what it
    must keep doing."""
    perms, ws, run = bash
    other = ws / "untouched.py"
    other.write_text("x\n", encoding="utf-8")
    stale_at = other.stat().st_mtime - 60.0
    perms.read_tracker[str(other.resolve())] = stale_at

    run("echo hello > " + str(ws / "elsewhere.txt"))

    assert perms.read_tracker[str(other.resolve())] == stale_at


def test_a_command_that_changes_nothing_moves_nothing(bash):
    perms, ws, run = bash
    target = ws / "same.txt"
    target.write_text("same\n", encoding="utf-8")
    _read_baseline(perms, target)
    was = perms.read_tracker[str(target.resolve())]

    # Writes the identical bytes: no change, so nothing to re-baseline.
    run(f"printf 'same\\n' > {target}")

    tracked = perms.read_tracker[str(target.resolve())]
    assert tracked == was or tracked >= target.stat().st_mtime


def test_a_file_created_by_the_shell_can_then_be_edited(bash):
    """Created here, so there is nothing to re-read: an edit must not be
    refused for want of a baseline the agent could not have."""
    perms, ws, run = bash
    fresh = ws / "created.py"
    run(f"printf 'a\\n' > {fresh}")
    assert str(fresh.resolve()) in perms.read_tracker
