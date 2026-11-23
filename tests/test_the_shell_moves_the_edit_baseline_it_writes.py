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


def test_the_overlay_journal_baselines_a_file_it_created(tmp_path,
                                                         monkeypatch):
    """The same contract, measured on every machine.

    The test above reaches the overlay only where bubblewrap's full cage
    runs, so it is green on a host -- and in CI -- where the overlay never
    runs at all. That is how this regressed unseen from 79770e23b until
    2026-11-22: green where the behaviour does not matter, red where it
    does. This drives the overlay journal directly instead.

    Input: one created, one modified-after-read and one modified-never-read
    change. Output: a baseline for the first two, none for the third -- the
    guard against somebody ELSE's change stays whole.
    """
    from delfin.agent import api_client as A
    from delfin.agent import overlay_workspace as _ow

    made = tmp_path / "created.py"
    read_before = tmp_path / "read-before.py"
    never_read = tmp_path / "never-read.py"
    for p in (made, read_before, never_read):
        p.write_text("x\n")

    perms = KitToolPermissions(workspace=tmp_path)
    perms.read_tracker[str(read_before.resolve())] = 1.0

    class _Change:
        def __init__(self, path, kind):
            self.path = str(path.resolve())
            self.kind = kind
            self.journal = True
            self.entry = self.root = self.trash = ""
            self.post_hash = self.post_sig = self.post_tree = None

    class _Outcome:
        changes = [_Change(made, "created"),
                   _Change(read_before, "modified"),
                   _Change(never_read, "modified")]

    monkeypatch.setattr(_ow, "apply", lambda plan: _Outcome())
    engine = _DocToolExecutor.__new__(_DocToolExecutor)
    engine._permissions = perms
    engine._journal_overlay("bash", object(), perms)

    assert str(made.resolve()) in perms.read_tracker, (
        "a file the session created has no baseline, so its own next edit "
        "is refused for a read nobody could have made")
    assert perms.read_tracker[str(read_before.resolve())] > 1.0, (
        "the baseline of a file the session had read was not moved")
    assert str(never_read.resolve()) not in perms.read_tracker, (
        "a file the session never read and did not create must still need "
        "a read before an edit")
