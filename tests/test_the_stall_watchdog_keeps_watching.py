"""The stall watchdog looked once and then stopped watching.

`_arm_stale_watcher` started a single non-repeating
`Timer(threshold + 1.0, _check_stale)`, and `_check_stale` had no re-arm
branch. So the check happened exactly once per turn, at threshold+1
seconds.

A turn that streamed normally for twelve minutes and THEN went silent was
never flagged: the one check had already fired at t=601s, seen recent
activity, and returned. The spinner stayed normal however long the
silence lasted, and the user had no signal that the provider had died.

Its sibling `_check_kill` re-arms itself for exactly this reason -- the
mechanism was understood, it was simply never applied to the warning
half. This is the watchdog the user relies on to notice a dead endpoint,
so "fires at most once, early" is the same as "does not fire".
"""

from __future__ import annotations

import pathlib

import pytest


_TAB = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
        / "dashboard" / "tab_agent.py")


def _body(name: str) -> str:
    src = _TAB.read_text(encoding="utf-8")
    start = src.index(f"def {name}(")
    nxt = src.find("\n        def ", start + 1)
    if nxt < 0:
        nxt = src.find("\n    def ", start + 1)
    return src[start:nxt if nxt > 0 else start + 4000]


def test_the_warning_check_re_arms_itself():
    body = _body("_check_stale")
    assert "_threading.Timer" in body, (
        "the stall check still fires once and stops watching")


def test_it_stops_re_arming_once_it_has_warned():
    """Warning every few seconds forever would be worse than silence."""
    body = _body("_check_stale")
    i = body.index('state["_stale_seen"] = True')
    j = body.index("_threading.Timer")
    assert i < j, "the re-arm must sit after the warning branch"
    assert "return" in body[i:j], (
        "the warning branch must return instead of re-arming")


def test_it_stops_when_the_turn_ends():
    body = _body("_check_stale")
    assert 'state.get("streaming")' in body
    i = body.index('state.get("streaming")')
    assert "return" in body[i:i + 120]


def test_the_re_armed_timer_is_a_daemon():
    """A live non-daemon timer would keep the process alive at exit."""
    body = _body("_check_stale")
    i = body.index("_threading.Timer")
    assert "daemon = True" in body[i:i + 300]


def test_the_handle_is_stored_so_it_can_be_cancelled():
    body = _body("_check_stale")
    i = body.index("_threading.Timer")
    assert '"_stale_timer"' in body[i:i + 400]


def test_the_kill_watchdog_still_re_arms():
    """The sibling this was modelled on must not have regressed."""
    body = _body("_check_kill")
    assert "_threading.Timer" in body or "again" in body
