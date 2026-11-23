"""The shell has two write paths, and only one moved the edit baseline.

``edit_file`` refuses a file whose mtime moved since the last
``read_file``: somebody changed it under you, re-read before editing.
``test_the_shell_moves_the_edit_baseline_it_writes`` establishes that a
shell write is not "somebody" -- it went through the same gates -- and
so the baseline moves with it.

That holds on the direct path, which captures a pre-image of each
redirect target and refreshes the baseline beside it. Over an overlay
the pre-image capture is skipped deliberately (``_ov.covers(path)``):
the layer is applied afterwards through the trash and is undoable at any
size, so no pre-image is needed. The baseline update lived in that same
skipped loop, and the overlay's own journalling refreshed a path only
when it was ALREADY tracked -- so a file the command created arrived at
the editor as a stranger.

Measured 2026-11-23 on a workstation: ``printf 'a\\n' > created.py``
through the bash tool left ``read_tracker`` empty, and the next
``edit_file`` on that file would be refused as "modified since last
read_file". The same test passes on a runner, where bubblewrap is
absent, the overlay is never planned and the direct path runs -- which
is why continuous integration never saw it.

Input: a journalled overlay outcome. Output: the edit baseline.

``overlay_workspace.apply`` is substituted rather than driven, so the
assertion is about the journalling step and holds on a machine with no
bubblewrap, no overlayfs and no privileges.
"""
from __future__ import annotations

import pathlib

import pytest

from delfin.agent import api_client as ac
from delfin.agent import overlay_workspace as ow
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture()
def journal(tmp_path, monkeypatch):
    """Drive ``_journal_overlay`` with an outcome we dictate."""
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    engine = _DocToolExecutor.__new__(_DocToolExecutor)
    engine._permissions = perms
    # The journal writes into the session's tool context, which only an
    # ``__init__``-built executor has, and the seq it records is not what
    # is under test here.
    monkeypatch.setattr(_DocToolExecutor, "_note_change_seq",
                        lambda self, rec, perms: None)
    monkeypatch.setattr(ac, "_report_file_changed", lambda path: None)

    def _run(*changes: ow.Applied) -> list[str]:
        monkeypatch.setattr(
            ow, "apply",
            lambda plan: ow.Outcome(changes=list(changes)))
        plan = ow.Plan(cmd_id="t", start_ns=0)
        return engine._journal_overlay("bash", plan, perms)

    return perms, tmp_path, _run


def _applied(path: pathlib.Path, kind: str) -> ow.Applied:
    return ow.Applied(path=str(path), kind=kind, entry="file",
                      root=str(path.parent))


def test_a_file_the_overlay_created_is_known_to_the_editor(journal):
    """The case that was broken: nothing to re-read, and an edit must
    not be refused for want of a baseline the agent could not have."""
    perms, ws, run = journal
    created = ws / "created.py"
    created.write_text("a\n", encoding="utf-8")

    run(_applied(created, "created"))

    assert str(created) in perms.read_tracker, (
        "a file the command created over an overlay reached the editor "
        "as a stranger")
    assert perms.read_tracker[str(created)] == created.stat().st_mtime


def test_a_file_the_overlay_modified_keeps_its_baseline_current(journal):
    """The case that already worked, kept: an append through the shell
    does not make the next edit look like somebody else's change."""
    perms, ws, run = journal
    target = ws / "notes.py"
    target.write_text("one\n", encoding="utf-8")
    perms.read_tracker[str(target)] = target.stat().st_mtime - 50

    run(_applied(target, "modified"))

    assert perms.read_tracker[str(target)] == target.stat().st_mtime


def test_a_deleted_file_is_not_given_a_baseline(journal):
    """A path that is gone has no mtime to record, and inventing one
    would tell the editor a file is current that no longer exists."""
    perms, ws, run = journal
    gone = ws / "removed.py"
    perms.read_tracker[str(gone)] = 1.0

    run(_applied(gone, "deleted"))

    assert perms.read_tracker[str(gone)] == 1.0


def test_a_change_the_overlay_does_not_journal_is_left_alone(journal):
    """``journal=False`` marks a change the session must not claim; it
    is not a write this agent is answerable for, so the baseline does
    not move for it."""
    perms, ws, run = journal
    other = ws / "untouched.py"
    other.write_text("x\n", encoding="utf-8")
    change = _applied(other, "created")
    change.journal = False

    run(change)

    assert str(other) not in perms.read_tracker


def test_a_path_that_vanished_before_the_journal_does_not_raise(journal):
    """Never raises is part of this method's contract: the notes are
    shown to the user, and an exception here would replace the whole
    report with one swallowed message."""
    perms, ws, run = journal
    missing = ws / "never-written.py"

    notes = run(_applied(missing, "created"))

    assert str(missing) not in perms.read_tracker
    assert not any("could not be applied" in n for n in notes)
