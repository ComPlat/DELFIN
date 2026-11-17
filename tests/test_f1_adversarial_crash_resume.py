"""Adversarial tests for package F1 Target 1 (crash-resume, reviewer s2).

Checks the recovery surface the resumed session actually receives: does
it get an actionable check-or-ask instruction for an OPEN side-effecting
step, not just a passive hint?  Written against the code as it stands on
main (04bf9703b) — these are expected RED here and turn green with the
builder's Target-1 fix.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import session_store
from delfin.agent import step_journal
from delfin.agent.state_paths import RedirectedUserState


@pytest.fixture()
def crash_state(tmp_path, monkeypatch):
    """A redirected state tree with a mid-turn crash already on disk."""
    with RedirectedUserState(tmp_path) as red:
        # Redact-then-save through the REAL writer so the record is one the
        # product produces (wave-15b: real records, not typed dicts).
        payload = {
            "user_message": "commit the fix and submit the job",
            "partial_response": "partially written answer...",
            "tool_calls": 2,
            "ts": 1000.0,
        }
        p = session_store.save_turn_checkpoint("sess-crash", payload)
        assert p is not None and p.exists()
        # Real journal writers: step 1 finished, step 2 OPEN (side-effecting).
        assert step_journal.start_step(
            "sess-crash", turn=7, step=1, tool="read_file", arguments={"path": "a.py"})
        assert step_journal.finish_step(
            "sess-crash", turn=7, step=1, tool="read_file", ok=True, result="ok")
        assert step_journal.start_step(
            "sess-crash", turn=7, step=2, tool="git_commit", arguments={"msg": "fix"})
        yield red


class TestRecoveryNoteIsActionable:
    def test_open_side_effect_step_gets_a_verifiable_check_instruction(
            self, crash_state):
        """Target 1: the resumed session must receive, per OPEN
        side-effecting step, a concrete, verifiable check instruction
        (what evidence decides 'did it happen'), and the note must be
        marked as a question the session answers before anything else.
        """
        note = session_store.consume_crash_recovery_note("sess-crash")
        assert note is not None, "a mid-turn crash must produce a recovery note"
        # (1) names the OPEN side-effecting step
        assert "step 2" in note
        assert "git_commit" in note
        assert "CHECK BEFORE REPEATING" in note
        # (2) says HOW to verify — a concrete evidence target the session
        # can check, not only the hint to check. The commit case: the
        # git log / HEAD state decides whether the step ran.
        assert "verify the workspace state first" in note  # existing hint
        # NEW requirement (the defect): a per-step check-or-ask that names
        # the evidence that decides, and is phrased as a question.
        lowered = note.lower()
        assert any(
            marker in lowered
            for marker in ("how to check", "how to verify", "check whether",
                           "did it happen", "did the commit")
        ), (
            "recovery note must tell the resumed session HOW to verify "
            "whether the open side-effecting step already ran, not only "
            "that it should: got note:\n" + note
        )

    def test_a_finished_step_is_never_listed_for_recheck(self, crash_state):
        """The finished step (read_file step 1) must not be offered as an
        OPEN step, and must not be suggested for re-running."""
        note = session_store.consume_crash_recovery_note("sess-crash")
        assert note is not None
        assert "step 1" not in note, (
            "a finished step must not be offered for re-check; note:\n" + note)

    def test_stale_checkpoint_produces_no_note(self, tmp_path):
        """A checkpoint older than the last full save is stale — no note,
        and the checkpoint is cleared either way."""
        with RedirectedUserState(tmp_path):
            assert session_store.consume_crash_recovery_note("none") is None
            payload = {"user_message": "t", "ts": 10.0}
            p = session_store.save_turn_checkpoint("sess-stale", payload)
            assert p is not None
            # A full save newer than the checkpoint makes it stale.
            session_store.save_session("sess-stale", mode="quick",
                                       title="full save after the crash")
            note = session_store.consume_crash_recovery_note("sess-stale")
            assert note is None
            assert session_store.load_turn_checkpoint("sess-stale") is None
