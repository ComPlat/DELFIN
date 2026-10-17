"""UX2 red tests: session_health.classify is a pure verdict function.

A record in, one verdict + reason + suggested next step out. Thresholds come
from the settings dict the caller passes (never read by classify itself) and
every out-of-range value falls back to its default instead of turning the
check off. No disk access, no wall clock: ``now`` is a parameter.
"""

import pytest

from delfin.agent import session_health


def _rec(**kw):
    rec = {
        "name": "nacht-s9",
        "session_id": "abc12345",
        "paused": False,
        "pid_alive": True,
        "open_approvals": [],
        "approval_opened_at": None,
        "last_turn_end": None,
        "last_turn": {},
        "idle_s": 0,
        "api_errors_in_a_row": 0,
        "resume_planned": False,
        "title": "demo session",
    }
    rec.update(kw)
    return rec


def _s(**kw):
    return {"health": dict({"wait_on_model_min": 10,
                            "stalled_min": 20,
                            "waiting_on_person_min": 15,
                            "api_errors": 3}, **kw)}


def _age(now, seconds):
    return now - seconds


def test_working_session_is_working():
    now = 1_000_000.0
    r = _rec(last_turn_end=_age(now, 30), idle_s=30)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "working"


def test_paused_beats_everything():
    now = 1_000_000.0
    r = _rec(paused=True, open_approvals=["a1"], idle_s=500)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "paused"
    assert "resume" in v["next"]


def test_waiting_on_person_open_approval():
    now = 1_000_000.0
    r = _rec(open_approvals=["a1"],
             approval_opened_at=_age(now, 120), idle_s=120)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "waiting_on_person"
    assert "approval" in v["reason"]


def test_waiting_on_model_no_token_for_ten_min():
    now = 1_000_000.0
    r = _rec(idle_s=660, last_turn_end=None, no_token_s=660)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "waiting_on_model"


def test_stalled_turn_ended_open_work_idle_over_threshold():
    now = 1_000_000.0
    r = _rec(last_turn_end=_age(now, 1800), idle_s=1800,
             last_turn={"stopped": True}, open_work=True)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "stalled"


def test_stopped_last_turn_stopped_nothing_resumes():
    now = 1_000_000.0
    r = _rec(last_turn_end=_age(now, 1800), idle_s=1800,
             last_turn={"stopped": True}, open_work=False)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "stopped"


def test_gone_pid_dead_lock_left():
    now = 1_000_000.0
    r = _rec(pid_alive=False)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "gone"


def test_erroring_after_k_api_errors():
    now = 1_000_000.0
    r = _rec(api_errors_in_a_row=3, last_turn_end=_age(now, 30), idle_s=30)
    v = session_health.classify(r, now, _s())
    assert v["verdict"] == "erroring"


def test_threshold_zero_negative_text_nan_fall_back():
    now = 1_000_000.0
    r = _rec(last_turn_end=_age(now, 1800), idle_s=1800,
             last_turn={"stopped": True}, open_work=True)
    for bad in (0, -5, "twenty", float("nan"), None, float("inf")):
        v = session_health.classify(r, now, _s(stalled_min=bad))
        assert v["verdict"] == "stalled", bad
