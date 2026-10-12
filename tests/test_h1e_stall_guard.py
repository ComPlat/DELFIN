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
