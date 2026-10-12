"""Unit tests for the pure stall guard (H1E phase 2).

The guard is the pure decision helper: it turns the sequence of harness
notices of one turn (kind + fingerprint of the text that caused the
notice) into a decision, continue or end the turn. Everything here tests
the module directly; the wiring into api_client.py is a separate red test
and an operator build (phase 3). Both directions are covered: the guard
fires when it should and does not fire when it should not.
"""

from delfin.agent.stall_guard import (
    CONTINUE,
    STOP,
    DEFAULT_REPEATS,
    StallGuard,
)


def test_three_identical_notices_end_the_turn():
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.streak == 3
    assert g.has_stalled()


def test_two_identical_notices_do_not_end_the_turn():
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()


def test_a_different_text_between_two_notices_resets_the_streak():
    # Different fingerprint (different model text that caused the notice)
    # resets: four notices alternating between two texts never stall.
    g = StallGuard()
    for _ in range(5):
        assert g.advance("dsml_leak", "fp-a") == CONTINUE
        assert g.advance("dsml_leak", "fp-b") == CONTINUE
    assert not g.has_stalled()


def test_a_different_kind_between_notices_resets_the_streak():
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    # A different notice kind resets the streak, so the dsml_leak streak
    # must be rebuilt from one.
    assert g.advance("empty_turn", "fp-x") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()


def test_progress_between_notices_never_ends_a_turn():
    # A turn that makes real progress (a tool actually ran) between two
    # identical notices is not ended: note_progress resets the streak.
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    g.note_progress()                 # a real tool executed
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    g.note_progress()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()


def test_repeats_is_configurable():
    g = StallGuard(repeats=4)
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == STOP


def test_defensive_repeats_parsing():
    # A corrupt or absurd repeats value must never throw in the stream loop.
    g = StallGuard(repeats="not-a-number")
    assert g.repeats == DEFAULT_REPEATS
    g0 = StallGuard(repeats=0)
    assert g0.repeats == DEFAULT_REPEATS
    g_neg = StallGuard(repeats=-3)
    assert g_neg.repeats == DEFAULT_REPEATS


def test_status_line_names_the_reason_and_the_fix():
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    line = g.status_line()
    assert "dsml_leak" in line
    assert "fp-a" in line
    assert "call the tool" in line          # what to do differently
    assert "text" in line


def test_status_line_empty_before_any_stall():
    g = StallGuard()
    assert g.status_line() == ""


def test_streak_is_1_after_a_single_notice():
    g = StallGuard()
    assert g.streak == 0
    g.advance("dsml_leak", "fp-a")
    assert g.streak == 1


# --- Finding A (reviewer s25): a different notice must clear a stale stall ---

def test_a_different_notice_after_a_stall_clears_status():
    # A,A,A stalls; a DIFFERENT notice after it must return CONTINUE, clear
    # has_stalled() and the status line. A stale stall is a lie about the
    # last advance.
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.has_stalled()
    # A different fingerprint resets: live again, status gone.
    assert g.advance("dsml_leak", "fp-b") == CONTINUE
    assert not g.has_stalled()
    assert g.status_line() == ""


def test_a_different_kind_after_a_stall_clears_status():
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.has_stalled()
    # A different notice kind resets and clears the earlier stall.
    assert g.advance("empty_turn", "fp-x") == CONTINUE
    assert not g.has_stalled()
    assert g.status_line() == ""


def test_a_streak_can_rebuild_after_a_stall_then_reset():
    # The guard must be reusable after it stalls once: a reset followed by a
    # FRESH trio of an identical notice stalls again, cleanly.
    g = StallGuard()
    for _ in range(3):
        g.advance("dsml_leak", "fp-a")
    assert g.has_stalled()
    g.advance("dsml_leak", "fp-b")          # reset
    assert not g.has_stalled()
    g.advance("dsml_leak", "fp-c")          # new streak of 1
    g.advance("dsml_leak", "fp-c")          # 2
    assert g.advance("dsml_leak", "fp-c") == STOP   # 3 -> stalls again
    assert "fp-c" in g.status_line()
    assert "fp-a" not in g.status_line()


# --- Finding B (reviewer s25): missing/empty fingerprint never false-stalls ---

def test_missing_fingerprint_never_false_stalls():
    # The caller could not produce a digest (e.g. leak parse failed): the
    # notice is non-identifiable. Identity cannot be proven, so it must never
    # accumulate toward a stall. Five structurally-distinct leaks with a
    # missing digest must all CONTINUE.
    g = StallGuard()
    for _ in range(5):
        assert g.advance("dsml_leak", None) == CONTINUE
    assert not g.has_stalled()
    assert g.status_line() == ""


def test_empty_fingerprint_never_false_stalls():
    g = StallGuard()
    for _ in range(5):
        assert g.advance("dsml_leak", "") == CONTINUE
    assert not g.has_stalled()


def test_whitespace_fingerprint_is_treated_as_missing():
    g = StallGuard()
    for _ in range(5):
        assert g.advance("dsml_leak", "   ") == CONTINUE
    assert not g.has_stalled()


def test_non_identifiable_notice_breaks_consecutiveness():
    # A missing-fingerprint notice in the MIDDLE of two otherwise-identical
    # notices breaks their consecutiveness, so they never reach repeats.
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    g.advance("dsml_leak", None)            # non-identifiable, resets
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    g.advance("dsml_leak", "")              # again, resets
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()
