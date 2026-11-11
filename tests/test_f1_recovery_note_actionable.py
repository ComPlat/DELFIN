"""F1 builder red-control: the crash-recovery note must be actionable.

Target 1 (reviewer s2, finding 5): for an OPEN side-effecting step, the
recovery note must carry a concrete, verifiable check instruction phrased
as a question that names the evidence deciding 'did it happen' — not a
passive hint. Written against the code as it stands (RED here); it turns
green once the renderer in delfin/agent/step_journal.py emits such a
per-step question (built by the operator from .gate/px copies).
"""

from __future__ import annotations

import pytest

from delfin.agent import session_store
from delfin.agent import step_journal
from delfin.agent.state_paths import RedirectedUserState


@pytest.fixture()
def crash_state(tmp_path, monkeypatch):
    """A redirected state tree with an open side-effecting step on disk."""
    with RedirectedUserState(tmp_path) as red:
        # A surviving mid-turn checkpoint is what signals "the process died
        # mid-turn" (session_store.consume_crash_recovery_note refuses to
        # build a note without one). Write it through the REAL checkpoint
        # writer, then open a side-effecting step through the REAL step
        # writer, so the note's step section comes from render_resume_note.
        session_store.save_turn_checkpoint("sess-act", {
            "ts": 1000.0, "user_message": "fix the bug",
            "tool_calls": 1, "partial_response": "",
        })
        assert step_journal.start_step(
            "sess-act", turn=7, step=2, tool="git_commit",
            arguments={"msg": "fix the bug"})
        yield red


def test_recovery_note_names_a_verifiable_check_for_open_side_effect(
        crash_state):
    """The resumed session must receive, per OPEN side-effecting step, a
    check-or-ask instruction phrased as a question that names the evidence
    deciding whether it happened."""
    note = session_store.consume_crash_recovery_note("sess-act")
    assert note is not None, "an open side-effecting step yields a note"
    # names the OPEN step and its tool
    assert "step 2" in note
    assert "git_commit" in note
    assert "CHECK BEFORE REPEATING" in note
    # a concrete, question-phrased check that names evidence
    assert any(marker in note for marker in (
        "how to check", "how to verify", "check whether",
        "did it happen", "did the commit", "did the push",
        "is the file", "did the job",
    )), f"note has no verifiable check question: {note!r}"
