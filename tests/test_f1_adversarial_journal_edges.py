"""Adversarial edge cases for the step journal (F1 reviewer s2).

Attacks the journal reading/writing surface with the shapes a real
crash produces: torn lines, unicode session ids, turn/step masking,
non-integer fields. Written against main 04bf9703b.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import step_journal
from delfin.agent.state_paths import RedirectedUserState


@pytest.fixture()
def jroot(tmp_path):
    """Redirected state; the journal resolves against the redirected
    _SESSIONS_DIR (journal_path reads it per call)."""
    with RedirectedUserState(tmp_path):
        yield tmp_path


class TestJournalEdgeCases:
    def test_unicode_session_id_maps_to_a_stable_path(self, jroot):
        """A unicode session id must produce one deterministic, safe
        journal path (write then read-back round-trip)."""
        p1 = step_journal.journal_path("Sitzung-ö → セッション")
        p2 = step_journal.journal_path("Sitzung-ö → セッション")
        assert p1 == p2
        assert step_journal.start_step(
            "Sitzung-ö → セッション", turn=1, step=1, tool="read_file")
        steps, notes = step_journal.read_steps("Sitzung-ö → セッション")
        assert len(steps) == 1, notes
        assert steps[0]["tool"] == "read_file"

    def test_open_step_survives_a_torn_line_after_it(self, jroot):
        """A crash mid-write leaves a torn last line. An OPEN step before
        the torn line must still be reported as open (the torn line is a
        note, never a finish)."""
        sid = "torn"
        assert step_journal.start_step(sid, turn=3, step=1,
                                       tool="git_commit", arguments={"m": "x"})
        jp = step_journal.journal_path(sid)
        with open(jp, "a", encoding="utf-8") as f:
            f.write('{"kind": "finish", "ts": 1.0, "turn": 3, "step": 1, ')
        open_out, notes = step_journal.open_steps(sid)
        assert len(open_out) == 1
        assert open_out[0]["step"] == 1
        assert any("torn" in n.lower() or "corrupt" in n.lower()
                   for n in notes), (
            "a torn line must be reported where the journal is read; "
            f"notes: {notes}")

    def test_stale_finish_of_another_turn_does_not_mask_an_open_start(
            self, jroot):
        """Step numbers restart per turn: a finish of step 1 in turn 7
        must not close the step-1 start of turn 8."""
        sid = "masking"
        assert step_journal.start_step(sid, turn=7, step=1, tool="read_file")
        assert step_journal.finish_step(sid, turn=7, step=1,
                                        tool="read_file", ok=True)
        assert step_journal.start_step(sid, turn=8, step=1,
                                       tool="git_commit", arguments={})
        open_out, _ = step_journal.open_steps(sid)
        assert [(e["turn"], e["step"]) for e in open_out] == [(8, 1)], (
            "a stale finish from turn 7 must not mask turn 8's open start")

    def test_resume_note_handles_a_huge_open_step_count(self, jroot):
        """max_open caps the lines and says how many more there are."""
        sid = "many"
        for s in range(1, 13):
            assert step_journal.start_step(sid, turn=9, step=s,
                                           tool="read_file")
        note = step_journal.render_resume_note(sid, max_open=5)
        assert note is not None
        step_lines = [ln for ln in note.splitlines()
                      if ln.startswith("  step ")]
        assert len(step_lines) == 5  # max_open caps the lines
        assert "… and 7 more open steps" in note

    def test_resume_note_none_when_everything_finished(self, jroot):
        sid = "clean"
        assert step_journal.start_step(sid, turn=2, step=1, tool="read_file")
        assert step_journal.finish_step(sid, turn=2, step=1,
                                        tool="read_file", ok=True)
        assert step_journal.render_resume_note(sid) is None
