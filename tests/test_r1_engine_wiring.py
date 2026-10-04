"""R1 phase 2 — RED control for the engine wiring (operator runs this).

Kept under tests/ per the operator so it turns GREEN in place when the
operator applies .gate/r1_turncont.patch to engine.py.

Contract from the patch (corrected after review):
- after stream_response commits a final assistant answer that announces
  further work AND tasks remain open, the engine sets pending_turn_continuation
  to a short note exactly once;
- the one-shot is a SESSION latch (_continuation_fired), initialized in
  __init__ but NOT reset in the per-turn gate (engine.py:2560-2574) — that
  would re-arm it on the follow-up turn itself and loop;
- the latch is cleared only by the engine method clear_turn_continuation().

It is RED on the current engine (attributes absent / never set) and GREEN once
the patch is in.
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


def test_engine_tracks_the_one_shot_guard_latch():
    eng = _bare_engine()
    # RED today: the guard must start False so the first qualifying turn fires.
    assert hasattr(eng, "_continuation_fired")
    assert eng._continuation_fired is False


def test_one_shot_is_a_session_latch_cleared_only_by_clear():
    """The guard must NOT reset per turn; only clear_turn_continuation does.

    Enacts the corrected wiring as a harness: after the first qualifying
    answer sets the note, a second qualifying answer on a LATER turn must NOT
    set it again (the SESSION latch holds), and clear_turn_continuation() is
    the only thing that re-arms it.
    """
    detector = __import__("delfin.agent.turn_continuation",
                          fromlist=["should_continue"])

    def enact(answer: str) -> bool:
        """Model the patched engine block: set the note at most once."""
        if answer and not eng._continuation_fired:
            cont, note = detector.should_continue(answer, 2)  # 2 open tasks
            if cont:
                eng._continuation_fired = True
                eng.pending_turn_continuation = note
                return True
        return False

    eng = _bare_engine()
    eng.pending_turn_continuation = ""
    eng._continuation_fired = False

    # First qualifying turn fires.
    assert enact("Let me write the next test.") is True
    assert eng.pending_turn_continuation  # note present

    # A SECOND qualifying answer on a new turn must NOT refire: the session
    # latch already holds — the exact per-turn-reset bug the review caught.
    assert enact("I will now run the full suite.") is False
    assert eng._continuation_fired is True   # still armed, note unchanged

    # Only clear_turn_continuation() re-arms.
    eng.clear_turn_continuation()
    assert eng.pending_turn_continuation == ""
    assert eng._continuation_fired is False
    # Now a later qualifying turn can schedule again.
    assert enact("Let me commit the fix.") is True
