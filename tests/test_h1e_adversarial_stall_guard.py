"""Adversarial tests for the pure stall guard (H1E phase 2 review, nacht-s25).

These probe the builder's delfin/agent/stall_guard.py beyond its own unit
tests. Each test that fails is a finding to the builder with file:line.

Both directions of the contract are covered (fires when it should, does not
fire when it should not), plus the degenerate edges the wave-15 checklist
names: None / empty inputs, huge strings, and the recurrent-transient stall.
"""

from delfin.agent.stall_guard import (
    CONTINUE,
    STOP,
    DEFAULT_REPEATS,
    StallGuard,
)


def test_a_stall_clears_when_a_different_notice_follows():
    """After a turn stalls (3 identical), a DIFFERENT notice returns
    CONTINUE and must clear the stall: has_stalled() reflects only the LAST
    advance. Currently _stalled is set only on STOP (stall_guard.py:84) and
    cleared only by note_progress (:95), never on a reset, so it goes stale."""
    g = StallGuard()
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == STOP
    # The turn then made a genuine change: a DIFFERENT notice. The last
    # advance says CONTINUE, so has_stalled() must be False and the status
    # must no longer describe the old stall.
    assert g.advance("dsml_leak", "fp-b") == CONTINUE
    assert not g.has_stalled(), (
        "has_stalled() must reflect the last advance, not a stale earlier "
        "STOP (stall_guard.py:80-86 never clears _stalled on a reset)"
    )
    assert g.status_line() == ""
    assert g.streak == 1


def test_status_internal_consistency_after_stall_then_reset():
    """The decision returned by advance() and has_stalled()/status_line()
    must never disagree. After STOP then a different notice (CONTINUE),
    status_line() reports the OLD stall instead of being empty."""
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")   # STOP
    g.advance("empty_turn", "fp-z")  # a different kind -> reset, CONTINUE
    assert not g.has_stalled()
    assert g.status_line() == ""


def test_no_fingerprint_never_false_stalls():
    """When the caller cannot compute a fingerprint (None / empty), three
    STRUCTURALLY DIFFERENT notices must not be treated as identical: without
    a fingerprint the guard cannot prove the text is the same, so it must
    not STOP. Currently str(None)=='None' (:77) collapses all three."""
    g = StallGuard()
    # Three different real texts (different tool names) whose digest was
    # never computed — the wire's default None must not look like one text.
    assert g.advance("dsml_leak", None) == CONTINUE
    assert g.advance("dsml_leak", None) == CONTINUE
    result = g.advance("dsml_leak", None)
    assert result == CONTINUE, (
        "3 missing fingerprints collapsed to ('dsml_leak','None') and "
        "false-STOPped (stall_guard.py:77 str(fingerprint))"
    )
    assert not g.has_stalled()


def test_empty_fingerprint_is_not_a_false_stall():
    """An empty-string fingerprint means 'unknown text', the same argument
    as None: three of them must not end the turn."""
    g = StallGuard()
    for _ in range(3):
        assert g.advance("dsml_leak", "") == CONTINUE
    assert not g.has_stalled()


def test_huge_fingerprint_does_not_crash_nor_distort():
    """A very long model text must not break the guard and must behave like
    any other distinct fingerprint (a different huge fp resets)."""
    g = StallGuard()
    big = "x" * 100_000
    assert g.advance("dsml_leak", big) == CONTINUE
    assert g.advance("dsml_leak", big) == CONTINUE
    assert g.advance("dsml_leak", big) == STOP
    assert g.advance("dsml_leak", big + "y") == CONTINUE  # distinct -> reset
    assert g.streak == 1


def test_stall_does_not_permanently_poison_future_notices():
    """A stalled guest guard instance that keeps being fed (the wire did not
    end the turn immediately) must not end the turn on a SINGLE fresh notice
    after note_progress. note_progress clears the stall; a fresh identical
    pair must rebuild the streak from one."""
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    assert g.advance("dsml_leak", "fp-a") == STOP
    g.note_progress()           # a real tool ran afterwards
    assert not g.has_stalled()  # stall lifted
    # One fresh notice after progress must be CONTINUE, not STOP.
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert not g.has_stalled()


def test_repeats_1_vs_2_matches_spec():
    """repeats=2 must end on the 2nd identical notice. The contract line
    ('.same (kind,fingerprint) N times in a row') is inclusive."""
    g = StallGuard(repeats=2)
    assert g.advance("dsml_leak", "fp-a") == CONTINUE
    assert g.advance("dsml_leak", "fp-a") == STOP


def test_repeats_is_inclusive_at_boundary():
    """With default 3, exactly three identical notices STOP (inclusive),
    and a fourth different one is CONTINUE, not STOP."""
    g = StallGuard()
    g.advance("dsml_leak", "fp-a")
    g.advance("dsml_leak", "fp-a")
    assert g.advance("dsml_leak", "fp-a") == STOP
    assert g.advance("dsml_leak", "fp-b") == CONTINUE
    assert not g.has_stalled()
