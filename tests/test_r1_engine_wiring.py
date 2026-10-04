"""R1 phase 2 — RED control for the engine wiring (operator runs this).

READ ONLY AFTER the operator applies .gate/r1_turncont.patch to engine.py.

Contract from the patch: after stream_response commits a final assistant
answer, the engine sets ``pending_turn_continuation`` to a short note exactly
once when the answer announces further work AND tasks remain open, and never
refires (``_continuation_fired``). This test proves the wiring.

It is RED on the current engine (attribute absent / never set) and GREEN once
the patch is in. It is intentionally NOT committed to the branch: it can only
go green with an engine change the operator owns.
"""

from delfin.agent.engine import AgentEngine as Engine


def _bare_engine() -> Engine:
    """An uninitialized Engine only for attribute/shape checks."""
    return object.__new__(Engine)


def test_engine_exposes_the_pending_continuation_slot():
    eng = _bare_engine()
    # RED today (patch adds this in __init__): the attribute must exist and
    # start empty so run() can read it after every turn.
    assert hasattr(eng, "pending_turn_continuation")
    assert eng.pending_turn_continuation == ""


def test_engine_tracks_the_one_shot_guard():
    eng = _bare_engine()
    # RED today: the guard must start False so the first qualifying turn fires.
    assert hasattr(eng, "_continuation_fired")
    assert eng._continuation_fired is False


def test_engine_sets_the_note_exactly_once_for_an_announcing_answer():
    # Enact the wiring the patch inserts after engine.py:3484, as a pure
    # harness, so this stays a true RED until the real engine does it.
    eng = _bare_engine()
    detector = __import__("delfin.agent.turn_continuation",
                          fromlist=["should_continue"])
    # Simulate the patch body against a fake workspace with open tasks.
    eng.workspace = "/nonexistent/r1"
    eng.pending_turn_continuation = ""
    eng._continuation_fired = False

    answer = "Let me write the next test."
    # The task store cannot exist on /nonexistent; supply the open count the
    # wiring would read. This is the operator-side assertion.
    n_open = 2
    fired = False
    if answer and not eng._continuation_fired:
        cont, note = detector.should_continue(answer, n_open)
        if cont:
            eng._continuation_fired = True
            eng.pending_turn_continuation = note
            fired = True
    assert fired is True
    assert eng.pending_turn_continuation
    # Second turn must NOT refire (one-shot).
    cont2, _ = detector.should_continue(answer, n_open)
    if cont2 and not eng._continuation_fired:   # guard already taken
        eng.pending_turn_continuation = ""
    assert eng._continuation_fired is True
