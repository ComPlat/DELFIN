"""A Stop belongs to the turn it interrupted.

Field report 20260805-105313: the user stopped a running turn, then asked
the same question three times. Every one came back "Agent returned no
output", and no turn-metrics row was written for any of them — the
request never reached the backend at all.

The flag had no owner. `request_stop` set it; the ONLY place that cleared
it was the top of `stream_response`. The dashboard worker gates its loop
on that same flag BEFORE calling `stream_response`, so once it was set
the code that resets it could never run again: the reset sat behind the
check the flag itself failed. Permanent, in every mode, until /reset —
and /retry, which the message recommended, funnels through the same gate
and cannot help.

Five other call sites set the same flag: the stall watchdog, the
permission-denial breaker, the conflict stop, the inspector stop and the
per-turn cost circuit-breaker. Each one muted the agent the same way.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent.engine import AgentEngine
from delfin.dashboard.tab_agent import _should_run_step


class _Engine:
    """Just the two attributes the gate reads."""

    def __init__(self, stopped=False, complete=False):
        self._stop_requested = stopped
        self.is_cycle_complete = complete


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------

def test_a_new_turn_runs_even_with_a_stop_left_over():
    """The regression itself: step 0 must not read a previous turn's stop."""
    assert _should_run_step(_Engine(stopped=True), 0) is True


def test_a_stop_still_ends_a_multi_role_route_at_the_next_handoff():
    """From step 1 on the flag is the real between-roles stop and must hold."""
    assert _should_run_step(_Engine(stopped=True), 1) is False
    assert _should_run_step(_Engine(stopped=True), 3) is False


def test_a_completed_cycle_stops_at_every_step():
    for step in (0, 1, 5):
        assert _should_run_step(_Engine(complete=True), step) is False


def test_an_ordinary_turn_runs():
    for step in (0, 1, 5):
        assert _should_run_step(_Engine(), step) is True


def test_a_missing_engine_does_not_raise():
    assert _should_run_step(None, 0) is False


def test_an_engine_without_the_attributes_is_not_blocked():
    """getattr defaults matter: a stale stop must never be inferred from an
    attribute that simply is not there."""
    class _Bare:
        pass

    assert _should_run_step(_Bare(), 0) is True


# ---------------------------------------------------------------------------
# The flag has an owner now
# ---------------------------------------------------------------------------

def test_clear_stop_exists_and_clears():
    eng = AgentEngine.__new__(AgentEngine)
    eng._stop_requested = True
    eng.clear_stop()
    assert eng._stop_requested is False


def test_request_stop_still_sets_it():
    eng = AgentEngine.__new__(AgentEngine)
    eng._stop_requested = False
    eng.request_stop()
    assert eng._stop_requested is True


def test_clearing_is_idempotent():
    eng = AgentEngine.__new__(AgentEngine)
    eng._stop_requested = False
    eng.clear_stop()
    eng.clear_stop()
    assert eng._stop_requested is False


# ---------------------------------------------------------------------------
# The wiring — the reset must be reachable
# ---------------------------------------------------------------------------

_TAB_SOURCE = pathlib.Path(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
).read_text(encoding="utf-8")


def test_the_worker_clears_the_stop_before_the_loop():
    """Clearing it anywhere after the gate would be no clearing at all."""
    clear_at = _TAB_SOURCE.index("engine.clear_stop()")
    loop_at = _TAB_SOURCE.index("for _step in range(max_auto_steps):")
    assert clear_at < loop_at, (
        "clear_stop runs after the gate again — the reset is unreachable, "
        "which is the original bug")


def test_the_gate_no_longer_reads_the_flag_inline():
    """One definition, testable. The inline form is what hid this for months."""
    assert "if engine._stop_requested or engine.is_cycle_complete:" not in _TAB_SOURCE


def test_the_engine_still_resets_on_its_own_turn_start():
    """Belt and braces: the dashboard is not the only caller of the engine."""
    source = pathlib.Path(
        __import__("delfin.agent.engine", fromlist=["x"]).__file__
    ).read_text(encoding="utf-8")
    head = source[source.index("def stream_response"):]
    head = head[:head.index("_turn_t0")]
    assert "self._stop_requested = False" in head, (
        "stream_response no longer clears the flag itself; a caller that "
        "forgets to would be silently muted")
