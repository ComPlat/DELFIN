"""After a crash, the recovery note names the step that may have run.

Input: the step journal (delfin/agent/step_journal.py) and the two places
that read it in a running product:
- session_store.consume_crash_recovery_note, which the session load calls
  when a mid-turn checkpoint survived a crash;
- working_state.build_working_state_block, which compaction calls.

Before: render_resume_note was reached only through
cli_resume.step_resume_note, which nothing in the product called, so a
crashed session never saw its open step. The compaction block carried the
last finished step of every session, a read_file too, as "must not be
repeated; result replayable via replay_step" -- a function the model has
no tool for.
"""
from __future__ import annotations

import time

import pytest

from delfin.agent import session_store, step_journal, working_state


@pytest.fixture()
def sessions(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    return tmp_path


def _crash_checkpoint(sid):
    session_store.save_turn_checkpoint(
        sid, {"user_message": "submit the job", "partial_response": "",
              "tool_calls": 1, "ts": time.time()})


def test_the_recovery_note_names_the_open_step(sessions):
    step_journal.auto_start("s1", tool="bash",
                            arguments={"command": "sbatch job.sh"})
    _crash_checkpoint("s1")
    note = session_store.consume_crash_recovery_note("s1")
    assert note and "[recovered]" in note
    assert "bash" in note and "sbatch job.sh" in note
    assert "CHECK BEFORE REPEATING" in note


def test_a_recovery_without_an_open_step_adds_nothing(sessions):
    rec = step_journal.auto_start("s1", tool="bash",
                                  arguments={"command": "true"})
    step_journal.auto_finish("s1", rec, ok=True, result="")
    _crash_checkpoint("s1")
    note = session_store.consume_crash_recovery_note("s1")
    assert note and "recovered steps" not in note


def test_compaction_names_the_last_side_effect_without_a_missing_tool(
        sessions):
    for tool, args in (("bash", {"command": "git commit -m x"}),
                       ("read_file", {"path": "a.md"})):
        rec = step_journal.auto_start("s1", tool=tool, arguments=args)
        step_journal.auto_finish("s1", rec, ok=True, result="ok")
    block = working_state.build_working_state_block([], session_id="s1")
    lines = [ln for ln in block.splitlines() if "step" in ln]
    assert any("bash" in ln for ln in lines), block
    assert not any("read_file" in ln for ln in lines), block
    assert "replay_step" not in block


def test_a_session_that_only_read_has_no_step_section(sessions):
    rec = step_journal.auto_start("s1", tool="read_file",
                                  arguments={"path": "a.md"})
    step_journal.auto_finish("s1", rec, ok=True, result="text")
    assert working_state.build_working_state_block([], session_id="s1") == ""


def test_a_long_journal_is_read_from_its_tail(sessions, monkeypatch):
    """Every tool call reads the journal; it must not read all of it."""
    path = step_journal.journal_path("s1")
    filler = ('{"kind": "finish", "ok": true, "result_digest": "' + "x" * 300
              + '", "step": 0, "tool": "bash", "turn": 1}\n')
    path.write_text(filler * 40000)          # ~14 MB
    whole = []
    real = type(path).read_text
    monkeypatch.setattr(type(path), "read_text",
                        lambda self, *a, **k: whole.append(self)
                        or real(self, *a, **k))
    rec = step_journal.auto_start("s1", tool="bash",
                                  arguments={"command": "true"})
    assert rec["recorded"] and rec["step"] == 1
    assert path not in whole
