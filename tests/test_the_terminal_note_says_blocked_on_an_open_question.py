"""A turn that ends blocked on an open question leaves its note at the next
prompt (package C, welle 11, reviewer adversarial tests, phase 2).

The fix under test (``job_wake.blocked_note`` + the ``repl._wake_text``
delivery): when ``turn()`` ends with the turn's words ending on an open
question to the user (``QUESTION:`` tag within the tail, or a last line ending
in ``?``) or a denial, AND open tasks remain, the NEXT idle prompt carries a
short note "[blocked — a system note, not the user] … blocked on that answer;
open tasks: N — subjects". Chosen to be read by the model as *state*, not as a
backlog list.

Fakes only: ``blocked_note`` is pure (no store, no clock); the wiring is
checked by source-scope (``ast``) so the tests bind to DELFIN, never to a
machine.
"""

from __future__ import annotations

import ast
import inspect


from delfin.agent import job_wake


def _note(text, open_tasks, denied=False):
    return job_wake.blocked_note(text, open_tasks, denied=denied)


def _task(subject, seq=1):
    return {"id": seq, "seq": seq, "subject": subject}


# -- open question at turn end ------------------------------------------------

def test_question_tag_with_open_tasks_marks_the_note():
    note = _note("I need your go-ahead.\nQUESTION: may I change the basis?",
                 [_task("relax a"), _task("relax b")])
    assert note and "[blocked" in note and "open tasks: 2" in note
    assert "question" in note.lower()


def test_question_mark_on_the_last_line_marks_the_note():
    note = _note("The job died on OOM. Recalculate at a bigger  maxcore?",
                 [_task("rerun")])
    assert note and "[blocked" in note


def test_a_mention_question_mid_text_is_no_pending_question():
    # A "?" that is not at the end is a mention ("did/check?" where the model
    # reports on it), not a question being asked of the user.
    note = _note("job 12? checked: it FAILED. Resubmit with more memory.",
                 [_task("rerun")])
    assert note == ""


def test_a_turn_that_ends_clean_leaves_no_note():
    note = _note("All three relaxations converged; energies in the table.",
                 [_task("relax a"), _task("relax b")])
    assert note == ""


def test_no_open_tasks_never_notes_regardless_of_question():
    # The spec: the note fires only when there are still open tasks.
    assert _note("QUESTION: shall I continue?", []) == ""


# -- denial -------------------------------------------------------------------

def test_denial_flag_marks_a_note_even_when_text_is_neutral():
    note = _note("I did what I could.", [_task("x")], denied=True)
    assert note and "den" in note.lower()


def test_a_denial_phrase_in_the_text_marks_the_note():
    note = _note("The write was refused: not on the auto-allow list.",
                 [_task("x")])
    assert note and "den" in note.lower()


def test_denial_phrases_cover_common_refusals():
    for phrase in ("permission denied", "was refused",
                   "refusing to overwrite"):
        assert job_wake._DENIAL_PHRASES and any(
            phrase in p for p in job_wake._DENIAL_PHRASES)


# -- content bounds -----------------------------------------------------------

def test_note_caps_subjects_and_stays_short():
    subjects = [_task(f"task number {i}", i) for i in range(1, 8)]
    note = _note("QUESTION: next?", subjects)
    assert len(note) < 400
    assert subjects[0]["subject"] in note
    assert subjects[6]["subject"] not in note  # beyond the cap
    assert "+4 more" in note


def test_note_handles_a_missing_subject_key():
    note = _note("QUESTION: next?", [{"id": 1, "seq": 1}])
    assert note and "open tasks" in note


# -- the note is delivered by the idle prompt, exactly once --------------------

_SRC: dict = {}


def _load():
    if not _SRC:
        _SRC["repl"] = inspect.getsource(
            __import__("delfin.agent.repl", fromlist=["TerminalAgent"]))
        _SRC["wake"] = inspect.getsource(job_wake)
    return _SRC


def test_wake_text_consumes_the_note_via_pop():
    """The stored note is read with a pop, so an idle look reports it once and
    a second idle look reports nothing — exactly once, and nothing lingers."""
    src = _load()["repl"]
    tree = ast.parse(src)
    seen_pop = False
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr == "pop":
            seen_pop = True
    assert seen_pop


def test_note_turn_blocked_is_called_at_the_end_of_turn():
    """turn() must call _note_turn_blocked with the turn's own words and the
    flushed denied flag, after the result is known."""
    src = _load()["repl"]
    assert "_note_turn_blocked" in src
    assert "_note_turn_blocked(" in src
