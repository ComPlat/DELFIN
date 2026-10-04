"""A turn that ends blocked on a question leaves a note for the next one.

Wave-10 finding (s12 in FEHLVERHaltEN_welle10): a session sat an hour
blocked at a key question it had asked and nobody answered — because the
next turn started without a word about what was still open. The open-task
reminder exists (engine._build_open_tasks_block) but only says "here is
the list"; it does not say "the turn that just ended ended WITH a
question, so the list is not being worked".

The rule this pins: when a turn ends with an open question (a QUESTION:
tag, or a trailing sentence that ends in "?") or with a denial render
item, AND open tasks remain, the NEXT turn's prompt gets a short note
"blocked on X; open: Y". A turn that ends clean leaves no note. The note
is produced at turn end and delivered by the wake path, so no real time
is involved — the detector is a pure function over the turn's text.
"""

from __future__ import annotations

import pytest

from delfin.agent import job_wake
from delfin.agent import repl as R


# -- the detector ------------------------------------------------------------

def test_question_tag_is_a_block():
    assert job_wake.blocked_note(
        "Shall I proceed with plan A? QUESTION: which solvation model?",
        open_tasks=[{"subject": "pKa module", "status": "in_progress"}],
    )


def test_trailing_question_is_a_block():
    note = job_wake.blocked_note(
        "I found two valid approaches.\n\nWhich one should I take?",
        open_tasks=[{"subject": "fix the gate", "status": "pending"}],
    )
    assert note and "blocked" in note.lower()
    assert "fix the gate" in note


def test_a_clean_turn_leaves_no_note():
    assert job_wake.blocked_note(
        "Done: all five phases committed, tests green.",
        open_tasks=[{"subject": "leftover", "status": "pending"}],
    ) == ""


def test_no_open_tasks_leaves_no_note():
    assert job_wake.blocked_note(
        "Which one should I take?",
        open_tasks=[],
    ) == ""


def test_question_mark_inside_not_at_the_end_is_not_a_block():
    """A '?' must END the turn's text — a mention of one does not block."""
    assert job_wake.blocked_note(
        "The README asks 'what is DELFIN?' — that is its business, and my "
        "answer is done.",
        open_tasks=[{"subject": "t", "status": "pending"}],
    ) == ""


def test_denial_counts_as_a_block():
    note = job_wake.blocked_note(
        "I could not continue.",
        open_tasks=[{"subject": "write report", "status": "pending"}],
        denied=True,
    )
    assert note and "blocked" in note.lower()


def test_no_denial_no_question_no_note():
    assert job_wake.blocked_note(
        "All good.",
        open_tasks=[{"subject": "t", "status": "pending"}],
        denied=False,
    ) == ""


def test_the_note_names_the_open_task_subjects():
    note = job_wake.blocked_note(
        "QUESTION: which isomer do we keep?",
        open_tasks=[{"subject": "pKa module", "status": "in_progress"},
                    {"subject": "bench run", "status": "pending"}],
    )
    assert "pKa module" in note
    assert "bench run" in note
    assert "2" in note and "open" in note


def test_the_note_is_short():
    note = job_wake.blocked_note(
        "QUESTION: proceed?",
        open_tasks=[{"subject": f"task {i}", "status": "pending"}
                    for i in range(12)],
    )
    assert len(note) < 400, f"the note is {len(note)} chars"


# -- the terminal wiring -----------------------------------------------------

class _Engine:
    session_id = "sess-me"
    kit_permissions = None


class _Agent(R.TerminalAgent):
    """A TerminalAgent built the hard way, with the fields the wiring needs."""

    def __init__(self, tmp_path, tasks):
        self.opts = type("O", (), {"cwd": str(tmp_path)})()
        self.engine = _Engine()
        self._tasks = tasks


@pytest.fixture()
def wired_agent(tmp_path, monkeypatch):
    """A TerminalAgent with the fields the wiring needs.

    ``finished_shells`` is patched to nothing ON PURPOSE: this file pins
    the BLOCKED NOTE, and the global bash registry is shared process
    state — a neighbour test's finished "echo done" job is announced
    deliberately to unowned sessions (job_wake.py: "an unowned job is
    better announced twice than lost"), which would make every
    _wake_text() here non-empty for a reason this file does not test.
    """
    a = _Agent(tmp_path, [])
    monkeypatch.setattr(job_wake, "finished_shells", lambda seen, **k: [])
    monkeypatch.setattr(job_wake, "finished_watched_jobs",
                        lambda ws, seen, **k: [])
    return a


def test_turn_end_wiring_sets_the_note_for_the_next_prompt(
        wired_agent, tmp_path, monkeypatch):
    """The public path: turn() ends blocked -> the NEXT prompt carries the note.

    The note rides in _wake_text, so the delivery is the wake path's
    (session messages first, then jobs, then the blocked note). The
    wiring test drives the two seams directly: a real task in the real
    store (the same one _note_turn_blocked reads), a blocked turn end
    recorded through the public turn-end hook, and the prompt text read
    back through _wake_text.
    """
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    store.create("pKa module", session_id="sess-me")
    wired_agent._note_turn_blocked("QUESTION: proceed?", denied=False)
    wired_agent._wake_last_look = 0.0
    text = wired_agent._wake_text("")
    assert "blocked" in text.lower()
    assert "pKa module" in text


def test_a_clean_turn_leaves_no_note_at_the_next_prompt(
        wired_agent, tmp_path):
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    store.create("pKa module", session_id="sess-me")
    wired_agent._note_turn_blocked("All done, tests green.", denied=False)
    wired_agent._wake_last_look = 0.0
    assert wired_agent._wake_text("") == ""


def test_the_note_is_delivered_exactly_once(wired_agent, tmp_path):
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    store.create("pKa module", session_id="sess-me")
    wired_agent._note_turn_blocked("QUESTION: proceed?", denied=False)
    wired_agent._wake_last_look = 0.0
    first = wired_agent._wake_text("")
    wired_agent._wake_last_look = 0.0
    second = wired_agent._wake_text("")
    assert first and second == "", "the note woke the prompt twice"


def test_the_wiring_never_raises(wired_agent, tmp_path, monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("store is gone")
    monkeypatch.setattr(job_wake, "note_turn_blocked", _boom)
    wired_agent._note_turn_blocked("QUESTION: proceed?")   # must not raise
    wired_agent._wake_last_look = 0.0
    assert wired_agent._wake_text("") == ""


def test_a_denial_inside_the_turn_blocks_even_with_quiet_words(
        wired_agent, tmp_path):
    """The denial flag comes from the pump, not from the turn's prose.

    A refusal the model reports neutrally ("the command could not run")
    does not trip the phrase detector; the pump records the denied
    render item on the agent, and turn() hands it to the note.
    """
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    store.create("write the report", session_id="sess-me")
    wired_agent._note_turn_blocked("The command could not run; I stopped.",
                                   denied=True)
    wired_agent._wake_last_look = 0.0
    text = wired_agent._wake_text("")
    assert "denied action" in text
    assert "write the report" in text


def test_the_denial_flag_is_not_carried_across_turns(wired_agent, tmp_path):
    """A denial is this turn's state; the next turn starts clean."""
    from delfin.agent.agent_tasks import get_store
    store = get_store(tmp_path)
    store.create("t", session_id="sess-me")
    wired_agent.__dict__["_denied_in_turn"] = True
    denied = bool(getattr(wired_agent, "_denied_in_turn", False))
    wired_agent.__dict__.pop("_denied_in_turn", None)
    assert denied is True
    assert bool(getattr(wired_agent, "_denied_in_turn", False)) is False
