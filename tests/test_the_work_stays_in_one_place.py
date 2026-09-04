"""A write that drifts into a second folder is said out loud, mid-turn.

From the recorded run, report 20260904-104748. Three files —
`extract_hyperpolarizability.py`, the results CSV and the results JSON —
were written into the data archive AND into the home directory above it.
Two roots, three files, and the prompt says in as many words: lock to
ONE place, never write the same project under two roots, throwaway
scripts belong in `~/agent_workspace/<task-slug>/`.

The pin that records the first write already existed and already reached
the prompt. It did not bind — the same finding as the language rule and
the formula rule before it. So it now also rides the rail that reaches
the model WHILE it is working.

It advises rather than refuses: a gate would have to decide that a
legitimate second write is a mistake, and it cannot know that.
"""

from __future__ import annotations

from delfin.agent.engine import AgentEngine


class _Client:
    def __init__(self):
        self.notes: list[str] = []

    def push_run_note(self, text):
        self.notes.append(text)


def _engine(pinned="/w/archive"):
    engine = AgentEngine.__new__(AgentEngine)
    engine.client = _Client()
    engine._project_dir = pinned
    engine._stray_write_noted = False
    return engine


def _write(path):
    return ("mcp__kit-coding__write_file", {"path": path})


def test_a_second_root_is_named_with_both_places():
    engine = _engine()
    engine._note_stray_write(*_write("/w/extract_hyperpolarizability.py"))
    assert len(engine.client.notes) == 1
    note = engine.client.notes[0]
    assert "/w/archive" in note and "/w/extract_hyperpolarizability.py" in note
    assert "agent_workspace" in note


def test_a_write_in_the_pinned_folder_says_nothing():
    engine = _engine()
    engine._note_stray_write(*_write("/w/archive/results.csv"))
    assert engine.client.notes == []


def test_a_subfolder_of_the_pinned_folder_is_the_same_place():
    engine = _engine()
    engine._note_stray_write(*_write("/w/archive/run1/results.csv"))
    assert engine.client.notes == []


def test_a_sibling_that_merely_starts_the_same_is_not_the_same_place():
    """`/w/archive-old` is not inside `/w/archive`."""
    engine = _engine()
    engine._note_stray_write(*_write("/w/archive-old/results.csv"))
    assert len(engine.client.notes) == 1


def test_it_is_said_once_per_turn():
    engine = _engine()
    for path in ("/w/a.py", "/w/b.csv", "/w/c.json"):
        engine._note_stray_write(*_write(path))
    assert len(engine.client.notes) == 1


def test_nothing_is_said_before_the_first_write_pinned_a_place():
    """With no pin there is no second root to speak of."""
    engine = _engine(pinned="")
    engine._note_stray_write(*_write("/w/anywhere.py"))
    assert engine.client.notes == []


def test_a_read_is_not_a_write():
    engine = _engine()
    engine._note_stray_write("mcp__kit-coding__read_file",
                             {"path": "/w/elsewhere/x.out"})
    assert engine.client.notes == []


def test_a_client_that_cannot_take_notes_does_not_cost_the_turn():
    engine = _engine()
    engine.client = object()
    engine._note_stray_write(*_write("/w/stray.py"))     # must not raise
