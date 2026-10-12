"""Adversarial tests for the pure stall guard (H1E phase 2, reviewer).

These probe edges the builder's own tests left unpinned. The contract under
test (delfin/agent/stall_guard.py @ 66396a5a): the same (kind, fingerprint)
N times in a row (repeats, default 3) ends the turn; a different notice or
explicit progress resets the streak; no timing, no wall clock.
"""

from delfin.agent.stall_guard import (
    CONTINUE,
    STOP,
    DEFAULT_REPEATS,
    StallGuard,
)


def test_repeats_1_stops_on_the_first_notice():
    # The defensive parser allows repeats==1 (only <1 is rejected). With
    # repeats=1 a single notice ends the turn: the boundary must be exact.
    g = StallGuard(repeats=1)
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.has_stalled()


def test_large_repeats_never_stops_below_the_threshold():
    # A huge repeats must NEVER stall on a sub-threshold run (off-by-one in
    # the streak comparison would be caught here). Twenty identical notices
    # are far below 1000.
    g = StallGuard(repeats=1000)
    for _ in range(20):
        assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()
    assert g.streak == 20


def test_progress_mid_streak_then_rebuild_requires_fresh_repeats():
    # Two identical notices build streak=2; progress resets it; three more
    # identical notices must STOP at the THIRD post-reset call (streak
    # rebuilt from 0), not carry the old two forward. This pins that
    # note_progress genuinely discards the counter rather than accumulating.
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.streak == 2
    g.note_progress()
    assert g.streak == 0
    assert g.advance("dsml_leak", "fp-a") == CONTINUE   # post-reset #1
    assert g.advance("dsml_leak", "fp-a") == CONTINUE   # post-reset #2
    assert g.advance("dsml_leak", "fp-a") == STOP       # post-reset #3


def test_note_progress_after_stop_clears_the_stall():
    # After STOP, an explicit progress note clears the stalled state. This is
    # an explicit design choice (note_progress resets _stalled); pin it so a
    # later refactor cannot silently make a stopped turn un-stoppable or a
    # progressed turn stay latched.
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")  # STOP
    assert g.has_stalled()
    g.note_progress()
    assert not g.has_stalled()
    assert g.status_line() == ""


def test_whitespace_variant_fingerprint_resets_by_design():
    # The fingerprint is the caller's digest of the model text. A text that
    # differs only in trailing whitespace maps to a DIFFERENT fingerprint, so
    # by design (different fingerprint in between resets) those two variants
    # are never counted as the same repeat. Alternating them must therefore
    # never stall — each resets the other's streak.
    g = StallGuard()
    for _ in range(8):
        assert g.advance("dsml_leak", "call-now ") == CONTINUE
        assert g.advance("dsml_leak", "call-now") == CONTINUE
    assert not g.has_stalled()


def test_stop_is_sticky_for_consecutive_same_notices_after_the_threshold():
    # Once the threshold is crossed the turn must keep reporting STOP for
    # further identical notices (the wiring ends on the first STOP, but a
    # guard that went back to CONTINUE after one STOP would silently un-stall
    # if the caller kept polling). Pin that STOP persists while the streak
    # stays at or above repeats.
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.advance("dsml_leak", "fp-a") == STOP


def test_status_line_is_consistent_with_the_stalled_kind():
    # The status must name the EXACT stalled (kind, fingerprint) and the
    # repeats count, so the model knows what it repeated and how often.
    text = "do-thing-do-thing-do-thing"
    g = StallGuard()
    g.advance("dsml_leak", text)
    g.advance("dsml_leak", text)
    g.advance("dsml_leak", text)
    line = g.status_line()
    assert "dsml_leak" in line
    assert "do-thing-do-thing-do-thing" in line
    assert "3" in line                       # the repeats count
    assert "call the tool" in line           # the actionable fix
    assert "not as text" in line             # the specific correction


def test_empty_kind_and_empty_fingerprint_do_not_throw():
    # Defensive path: a notice with an empty kind or an empty fingerprint is
    # still a notice and must be processable without raising (the caller is
    # inside the stream loop).
    g = StallGuard()
    assert g.advance("", "") == CONTINUE
    assert g.advance("", "fp") == CONTINUE
    assert not g.has_stalled()


def test_numeric_fingerprint_coerces_and_compares_stably():
    # A numeric fingerprint is not str, but advance() coerces it via str()
    # and must not throw. Three identical 5's stall; a distinct numeric (0)
    # is a different fingerprint and resets the streak.
    g = StallGuard()
    assert g.advance("kind", 5) == CONTINUE
    assert g.advance("kind", 5) == CONTINUE
    assert g.advance("kind", 5) == STOP                # 3x "5" -> stop
    assert g.advance("kind", 0) == CONTINUE            # 0 != "5" -> reset
    assert g.advance("kind", 0) == CONTINUE
    assert g.streak == 2


def test_missing_fingerprint_never_stalls_or_latches():
    # Finding B (fixed in 24b2218f): a missing fingerprint means the caller
    # could not prove text identity, so it must NEVER count toward a stop
    # nor use a stale latch from an earlier stop. The old module coerced
    # None -> "None" and could false-STOP on 3 distinct un-digestible leaks;
    # the fix treats None/empty as "cannot prove identity" -> never stalls.
    g = StallGuard()
    assert g.advance("dsml_leak", None) == CONTINUE
    assert g.advance("dsml_leak", None) == CONTINUE
    assert g.advance("dsml_leak", None) == CONTINUE   # never STOP
    assert not g.has_stalled()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE  # resumes clean
    assert g.advance("dsml_leak", "fp-a") == CONTINUE


def test_alternating_kind_never_stalls_across_many_rounds():
    # Alternating between TWO notice kinds is not the reported stall (the
    # model is producing DIFFERENT output). Long alternation must never stop.
    g = StallGuard()
    for _ in range(30):
        assert g.advance("dsml_leak", "fp-a") == CONTINUE
        assert g.advance("verify", "fp-b") == CONTINUE
    assert not g.has_stalled()
