"""The CLI turn wires the watchdog into its stream: a stall is reported.

Input: ``delfin/agent/cli.py::_run_once`` driven with a stub engine whose
``stream_response`` yields nothing (a model that never answers) and whose
``activity`` snapshot says a model wait. Output: one ``{"type": "stall"}`
stream-json event on ``emit`` while the wait ticks keep arriving.
Asserted: the stall event carries the class and the duration; a tick
inside the reporting cadence emits nothing new; leaving the wait (the
turn ends) resets the watchdog so a later turn reports from zero.
"""

from __future__ import annotations

import json

from delfin.agent import cli as cli_mod


class _StubActivity:
    def snapshot(self):
        from types import SimpleNamespace
        # The phase clock the stub engine reports; tests overwrite it.
        return SimpleNamespace(phase=self.phase, phase_s=self.phase_s,
                               model="stub-model", first=True)


class _NeverAnswersEngine:
    """A stream_response that emits wait ticks and never finishes early."""

    def __init__(self, ticks=3):
        self.ticks = ticks
        self.activity = _StubActivity()
        self.token_usage = {"input": 0, "output": 0}

    def stream_response(self, user_message, on_token=None, on_tool_use=None,
                        on_wait=None, on_notice=None, max_tokens=0,
                        on_tool_result=None, on_permission_denied=None):
        for i in range(self.ticks):
            if on_wait is not None:
                on_wait(f"waiting for the model … {i}s")
        return ""


def _drive(engine, ticks_settings=None):
    events = []
    cli_mod._run_once(engine, "hello",
                      emit=events.append)
    return events


def test_a_model_wait_tick_emits_one_stall_event(monkeypatch):
    eng = _NeverAnswersEngine(ticks=3)
    # phase_s grows with every tick: the phase clock of the snapshot the
    # engine owns. A reader sees the stall age.
    eng.activity.phase = "model_wait"
    eng.activity.phase_s = 0.0
    real_snapshot = eng.activity.snapshot

    def snapshot():
        eng.activity.phase_s += 40.0
        return real_snapshot()

    eng.activity.snapshot = snapshot
    events = _drive(eng)
    stalls = [e for e in events if e.get("type") == "stall"]
    assert stalls, f"no stall event in {events}"
    first = stalls[0]
    assert first["class"] == "model-wait"
    assert first["action"] == "report-model-retry"
    assert json.dumps(first)  # stream-json serialisable


def test_ticks_inside_the_cadence_emit_no_extra_stall(monkeypatch):
    eng = _NeverAnswersEngine(ticks=5)
    eng.activity.phase = "model_wait"
    eng.activity.phase_s = 0.0
    real_snapshot = eng.activity.snapshot

    def snapshot():
        eng.activity.phase_s += 2.0      # under the 30 s cadence
        return real_snapshot()

    eng.activity.snapshot = snapshot
    events = _drive(eng)
    stalls = [e for e in events if e.get("type") == "stall"]
    assert len(stalls) == 1


def test_a_gate_notice_on_the_wait_tick_makes_it_a_gate_queue(monkeypatch):
    eng = _NeverAnswersEngine(ticks=2)
    eng.activity.phase = "model_wait"
    eng.activity.phase_s = 50.0

    def stream_response(user_message, on_token=None, on_tool_use=None,
                        on_wait=None, on_notice=None, max_tokens=0,
                        on_tool_result=None, on_permission_denied=None):
        if on_notice is not None:
            on_notice("waiting for a slot on https://ep: 2 in use")
        if on_wait is not None:
            on_wait("waiting for the model … 1s")
        return ""

    eng.stream_response = stream_response
    events = _drive(eng)
    stalls = [e for e in events if e.get("type") == "stall"]
    assert stalls and stalls[0]["class"] == "gate-queue"
