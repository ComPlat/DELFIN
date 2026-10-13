"""StallGuard: STOP at ``repeats`` identical consecutive endings, else
CONTINUE. The engine-level behaviour is in
test_zz_h1e_a_session_that_repeats_a_dead_turn_stops_asking.py."""

import pytest

from delfin.agent.stall_guard import CONTINUE, DEFAULT_REPEATS, STOP, StallGuard


def _feed(guard, *pairs):
    return [guard.advance(k, f) for k, f in pairs]


def test_the_third_identical_ending_stops_and_a_fourth_stays_stopped():
    g = StallGuard()
    assert _feed(g, *[("dsml_leak", "write_file")] * 4) == [
        CONTINUE, CONTINUE, STOP, STOP]
    assert g.streak == 4


@pytest.mark.parametrize("between", [("dsml_leak", "bash"), ("empty", "write_file")])
def test_a_different_fingerprint_or_kind_restarts_the_streak(between):
    g = StallGuard()
    seq = [("dsml_leak", "write_file")] * 2 + [between] + [("dsml_leak", "write_file")] * 2
    assert STOP not in _feed(g, *seq)


def test_progress_restarts_the_streak_including_after_a_stop():
    g = StallGuard()
    _feed(g, *[("dsml_leak", "write_file")] * 3)
    g.note_progress()
    assert g.streak == 0
    assert _feed(g, *[("dsml_leak", "write_file")] * 3) == [CONTINUE, CONTINUE, STOP]


@pytest.mark.parametrize("fp", [None, "", "   "])
def test_an_unidentifiable_ending_never_counts_and_breaks_a_streak(fp):
    g = StallGuard()
    assert STOP not in _feed(g, *[("dsml_leak", fp)] * 5)
    seq = [("dsml_leak", "write_file")] * 2 + [("dsml_leak", fp)] + [("dsml_leak", "write_file")] * 2
    assert STOP not in _feed(StallGuard(), *seq)


@pytest.mark.parametrize("raw, parsed", [
    (1, 1), ("5", 5), (0, DEFAULT_REPEATS), (-2, DEFAULT_REPEATS),
    ("abc", DEFAULT_REPEATS), (None, DEFAULT_REPEATS)])
def test_repeats_is_parsed_defensively(raw, parsed):
    assert StallGuard(raw).repeats == parsed


def test_repeats_one_stops_on_the_first_ending():
    assert StallGuard(1).advance("dsml_leak", "write_file") == STOP
