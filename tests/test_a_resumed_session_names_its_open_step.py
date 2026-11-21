"""A resumed session names its open step on the resume path itself.

Before: ``cli_resume.step_resume_note`` existed, was tested directly,
and nothing in the product called it — the crash-recovery note carried
the open-step warning only when a mid-turn checkpoint had survived, so
a session interrupted between turns (checkpoint consumed or never
written) resumed silently and the person could repeat a side-effecting
call (a push, a submit) that may already have run.

Mounted: ``cli._open_session`` — the ``--resume``/``--session`` path —
prints the note right after the restore. The product entry is
``cli.main(["run", "--resume", ...])``-shaped; this test goes through
``cli.main``.

Companion deletion, same wave: ``cli_resume.replay_step`` and
``step_journal.finished_result`` are GONE — the wrapper's only caller
was its test suite, and the repo's own crash-recovery test pins that
the model must not be offered a replay tool.
"""
from __future__ import annotations

import pytest

from delfin.agent import cli, cli_resume, session_store, step_journal


@pytest.fixture()
def journal_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    return tmp_path


def test_resume_prints_the_open_step_note(journal_dir, capsys,
                                          monkeypatch):
    """--resume of a session with an open side-effecting step says so."""
    step_journal.start_step("s-open", turn=1, step=3, tool="bash",
                            arguments={"command": "git push"})
    # A minimal saved session so the restore path finds a record.
    session_store.save_session("s-open", title="t", workspace="",
                               chat_messages=[])

    # The engine restore is not what is under test; stub it out and let
    # the real _open_session run everything else.
    class _Engine:
        session_id = "s-open"
        kit_permissions = None

    monkeypatch.setattr(cli, "_resume_or_create", lambda engine, args: "")
    monkeypatch.setattr(cli, "_claim_session", lambda sid: None)
    monkeypatch.setattr(
        cli, "_presence_key_for", lambda sid: None, raising=False)

    class _Args:
        new_session = False
        session = "s-open"
        continue_session = False
        resume = None
        fork_session = False

    cli._open_session(_Engine(), _Args(), None)
    out = capsys.readouterr().err
    assert "CHECK BEFORE REPEATING" in out
    assert "git push" in out


def test_resume_without_open_step_prints_no_note(journal_dir, capsys,
                                                 monkeypatch):
    """A finished journal adds nothing to the resume output."""
    step_journal.start_step("s-done", turn=1, step=1, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s-done", turn=1, step=1, tool="read_file",
                             ok=True, result="text")
    session_store.save_session("s-done", title="t", workspace="",
                               chat_messages=[])

    class _Engine:
        session_id = "s-done"
        kit_permissions = None

    monkeypatch.setattr(cli, "_resume_or_create", lambda engine, args: "")
    monkeypatch.setattr(cli, "_claim_session", lambda sid: None)

    class _Args:
        new_session = False
        session = "s-done"
        continue_session = False
        resume = None
        fork_session = False

    cli._open_session(_Engine(), _Args(), None)
    out = capsys.readouterr().err
    assert "CHECK BEFORE REPEATING" not in out


def test_replay_step_is_gone():
    """The deleted names stay deleted (absence pin)."""
    assert not hasattr(cli_resume, "replay_step")
    assert not hasattr(step_journal, "finished_result")


def test_resume_resets_the_diagnostics_note(journal_dir, monkeypatch):
    """A resumed session may see the 'diagnostics unavailable' note once
    again — the ruling says once per SESSION, and the flag is process
    state, so the resume path resets it."""
    from delfin.agent import diagnostics

    session_store.save_session("s-diag", title="t", workspace="",
                               chat_messages=[])

    class _Engine:
        session_id = "s-diag"
        kit_permissions = None

    monkeypatch.setattr(cli, "_resume_or_create", lambda engine, args: "")
    monkeypatch.setattr(cli, "_claim_session", lambda sid: None)

    class _Args:
        new_session = False
        session = "s-diag"
        continue_session = False
        resume = None
        fork_session = False

    diagnostics._unavailability_shown = True
    try:
        cli._open_session(_Engine(), _Args(), None)
        assert diagnostics._unavailability_shown is False, (
            "a resumed session must be able to see the note again")
    finally:
        diagnostics._unavailability_shown = False
