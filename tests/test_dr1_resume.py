"""DR1 phase 2 red tests: step-resume surfaces on cli_resume and
working_state.

The surfaces under test:

* ``cli_resume.step_resume_note(session_id)`` — the resume path's read
  of the step journal: the open-step note when a step is open, None
  otherwise, "" for an unknown session.
* ``cli_resume.replay_step(session_id, step)`` — a finished step's
  recorded result, so a resume can hand the recorded result back
  instead of re-running the tool.
* ``working_state.build_working_state_block(..., session_id=...)`` — a
  "Last recorded finished step" section built from the step journal, so
  compaction keeps the last finished step across restarts.

Storage is redirected: ``session_store._SESSIONS_DIR`` points at
``tmp_path`` for the whole module (per-call resolution makes that
stick), so no real home is touched.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import cli_resume, session_store, step_journal


@pytest.fixture()
def journal_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    return tmp_path


# -- cli_resume.step_resume_note ------------------------------------------

def test_step_resume_note_reports_open_side_effecting_step(journal_dir):
    step_journal.start_step("s1", turn=1, step=3, tool="bash",
                            arguments={"command": "git commit -F m.txt"},
                            )
    note = cli_resume.step_resume_note("s1")
    assert note is not None
    assert "step 3" in note
    assert "bash" in note
    assert "may or may not have run" in note
    assert "CHECK BEFORE REPEATING" in note


def test_step_resume_note_none_when_everything_is_finished(journal_dir):
    step_journal.start_step("s1", turn=1, step=1, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s1", turn=1, step=1, tool="bash", ok=True, result="content")
    assert cli_resume.step_resume_note("s1") is None


def test_step_resume_note_none_for_unknown_session(journal_dir):
    # No journal -> nothing to surface: None ("no open steps, no notes"),
    # not an error and not an empty note block.
    assert cli_resume.step_resume_note("no-such-session") is None


# -- cli_resume.replay_step -----------------------------------------------

def test_replay_step_gives_back_recorded_result(journal_dir):
    step_journal.start_step("s1", turn=1, step=2, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s1", turn=1, step=2, tool="read_file", ok=True,
                             result="# heading\ncontent")
    rec = cli_resume.replay_step("s1", 2)
    assert rec["ok"] is True
    # The record stores the scrubbed digest, never the raw text:
    assert rec["result_digest"]
    assert "content" not in rec["result_digest"]
    assert rec["tool"] == "read_file"

    # turn-narrowed replay matches the same record:
    assert cli_resume.replay_step("s1", 2, turn=1) is not None
    assert cli_resume.replay_step("s1", 2, turn=9) is None


def test_replay_step_unknown_step_is_none(journal_dir):
    step_journal.start_step("s1", turn=1, step=1, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s1", turn=1, step=1, tool="bash", ok=True, result="x")
    assert cli_resume.replay_step("s1", 99) is None
    assert cli_resume.replay_step("no-such", 1) is None


# -- working_state: the last finished step survives compaction -------------

def test_working_state_carries_last_finished_step(journal_dir):
    step_journal.start_step("s1", turn=1, step=1, tool="bash",
                            arguments={"command": "true"})
    step_journal.finish_step("s1", turn=1, step=1, tool="bash", ok=True, result="")
    step_journal.start_step("s1", turn=1, step=5, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s1", turn=1, step=5, tool="read_file", ok=True, result="text")

    from delfin.agent import working_state
    block = working_state.build_working_state_block([], session_id="s1")
    assert block, "a finished step must produce a non-empty block"
    step_line = [l for l in block.splitlines()
                 if "step" in l and "read_file" in l]
    assert step_line, f"no step line in block:\n{block}"
    assert "5" in step_line[0]
    assert "read_file" in step_line[0]
    assert "must not be repeated" in step_line[0]


def test_working_state_without_journal_has_no_step_section(journal_dir):
    from delfin.agent import working_state
    assert working_state.build_working_state_block([], session_id="nobody") == ""


def test_working_state_step_section_bounded(journal_dir):
    from delfin.agent import working_state
    step_journal.start_step("s1", turn=1, step=1, tool="x" * 400,
                            arguments={})
    step_journal.finish_step("s1", turn=1, step=1, tool="x" * 400, ok=True, result="r" * 900)
    block = working_state.build_working_state_block([], session_id="s1")
    assert len(block) <= working_state._MAX_BLOCK_CHARS + 40
