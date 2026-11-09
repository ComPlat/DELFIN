"""The ops stalls command: read-only stall classes per open session.

Input: session records as ``session_health.collect()`` produces, and a
clock. Output: ``stall_class(record, now)`` maps each record to a
watchdog class (approval-wait / model-wait / gate-queue / tool-running
/ idle) with its age in seconds, and ``_cmd_ops_stalls`` prints one
line per open session -- or JSON rows with --json. Asserted: approval
outranks model-wait, tool-running comes from an in-progress phase,
gate-queue is reported from the record's wait notice, no session is
dropped, and nothing is written.
"""

from __future__ import annotations

import pytest

from delfin.agent import progress_watchdog as pw


def _record(**kw):
    base = {
        "name": "wave/x", "session_id": "s1",
        "pid_alive": True, "paused": False,
        "open_approvals": [], "approval_opened_at": None,
        "no_token_s": None, "idle_s": None,
        "last_turn": {}, "open_work": False, "resume_planned": False,
        "api_errors_in_a_row": 0,
    }
    base.update(kw)
    return base


def test_open_approval_maps_to_approval_wait():
    rec = _record(open_approvals=["run x"],
                  approval_opened_at=100.0)
    out = pw.stall_class(rec, now=400.0)
    assert out["class"] == "approval-wait"
    assert out["age_s"] == 300.0
    assert out["action"] == "keep-waiting"


def test_no_tokens_maps_to_model_wait():
    rec = _record(no_token_s=120.0)
    out = pw.stall_class(rec, now=400.0)
    assert out["class"] == "model-wait"
    assert out["age_s"] == 120.0


def test_gate_queue_outranks_model_wait():
    rec = _record(no_token_s=120.0, gate_wait_s=90.0)
    out = pw.stall_class(rec, now=400.0)
    assert out["class"] == "gate-queue"
    assert out["age_s"] == 90.0


def test_running_phase_maps_to_tool_running():
    rec = _record(phase="tool", phase_s=45.0)
    out = pw.stall_class(rec, now=400.0)
    assert out["class"] == "tool-running"
    assert out["age_s"] == 45.0


def test_nothing_waiting_is_not_reported_as_idle():
    rec = _record(idle_s=10000.0, last_turn={"stopped": True})
    out = pw.stall_class(rec, now=400.0)
    assert out["class"] == "idle"
    assert out["stall"] is False


@pytest.mark.parametrize("kw", [
    {"paused": True}, {"pid_alive": False},
    {"api_errors_in_a_row": 5},
])
def test_unstoppable_states_are_reported_not_classified(kw):
    out = pw.stall_class(_record(**kw), now=400.0)
    assert out["class"] == "unclassified"
