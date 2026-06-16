"""bash_status(wait_seconds=…) — server-side blocking wait.

Regression test for bug 20260615-152119: the agent polled bash_status every
few seconds for a ~10-min job, exhausting the tool-round budget. With a
blocking wait, one tool round covers a long wait and returns the instant the
job ends.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ


class _FakeJob:
    def __init__(self, finish_at: float):
        self._t0 = time.monotonic()
        self._finish_at = finish_at

    def _done(self) -> bool:
        return (time.monotonic() - self._t0) >= self._finish_at

    def poll(self):
        return 0 if self._done() else None

    def status_dict(self) -> dict:
        done = self._done()
        return {"job_id": "j1", "running": not done,
                "exit_code": 0 if done else None}


class _FakeReg:
    def __init__(self, job):
        self._job = job

    def get(self, job_id):
        return self._job


def _patch(monkeypatch, job):
    monkeypatch.setattr(BJ, "get_registry", lambda: _FakeReg(job))


@pytest.fixture(autouse=True)
def _clear_poll_history():
    # The executor is a process singleton; clear the busy-poll history so the
    # detector never leaks state across tests (which reuse job_id "j1").
    A._doc_executor.__dict__.pop("_bash_poll_ts", None)
    yield
    A._doc_executor.__dict__.pop("_bash_poll_ts", None)


def _status(**args) -> dict:
    return json.loads(A._doc_executor._execute_bash_status(args))


def test_no_wait_returns_immediately(monkeypatch):
    _patch(monkeypatch, _FakeJob(finish_at=100))   # still running
    t = time.monotonic()
    out = _status(job_id="j1")
    assert out["running"] is True
    assert time.monotonic() - t < 0.5              # did not block


def test_wait_returns_immediately_when_already_finished(monkeypatch):
    _patch(monkeypatch, _FakeJob(finish_at=0))     # already done
    t = time.monotonic()
    out = _status(job_id="j1", wait_seconds=10)
    assert out["running"] is False
    assert time.monotonic() - t < 0.5              # no needless blocking


def test_wait_returns_early_when_job_finishes(monkeypatch):
    _patch(monkeypatch, _FakeJob(finish_at=0.3))
    t = time.monotonic()
    out = _status(job_id="j1", wait_seconds=30)
    dt = time.monotonic() - t
    assert out["running"] is False
    assert dt < 3.0                                # returned ~when job ended, not at 30s


def test_wait_is_capped_for_a_long_running_job(monkeypatch):
    monkeypatch.setattr(A, "_BASH_STATUS_WAIT_CAP_S", 0.4)
    _patch(monkeypatch, _FakeJob(finish_at=1000))  # never finishes in-test
    t = time.monotonic()
    out = _status(job_id="j1", wait_seconds=30)    # asks 30s, capped to 0.4s
    dt = time.monotonic() - t
    assert out["running"] is True
    assert dt < 3.0                                # blocked ~cap, not 30s


def test_unknown_job_id(monkeypatch):
    monkeypatch.setattr(BJ, "get_registry", lambda: _FakeReg(None))
    out = _status(job_id="nope", wait_seconds=10)
    assert "error" in out


def test_busy_poll_is_throttled_after_first_snapshot(monkeypatch):
    # First check = instant snapshot; a rapid re-check of the SAME running job
    # (no wait_seconds) is auto-throttled so a tight loop can't burn rounds.
    monkeypatch.setattr(A, "_BASH_STATUS_BUSY_POLL_WAIT_S", 0.4)
    _patch(monkeypatch, _FakeJob(finish_at=100))   # never finishes in-test

    t = time.monotonic()
    out1 = _status(job_id="j1")                     # first → instant
    assert out1["running"] is True
    assert time.monotonic() - t < 0.3

    t = time.monotonic()
    out2 = _status(job_id="j1")                     # rapid re-poll → throttled
    dt = time.monotonic() - t
    assert out2["running"] is True
    assert 0.3 <= dt < 3.0                          # blocked ~busy-poll window


def test_busy_poll_guard_bypassed_by_explicit_zero(monkeypatch):
    # An explicit wait_seconds=0 is a deliberate snapshot request and is NOT
    # throttled, even on a repeated poll of a running job.
    monkeypatch.setattr(A, "_BASH_STATUS_BUSY_POLL_WAIT_S", 5.0)
    _patch(monkeypatch, _FakeJob(finish_at=100))
    _status(job_id="j1")                            # seed history
    t = time.monotonic()
    out = _status(job_id="j1", wait_seconds=0)      # explicit snapshot
    assert out["running"] is True
    assert time.monotonic() - t < 0.5               # not throttled


def test_busy_poll_throttle_returns_early_when_job_finishes(monkeypatch):
    # Throttled wait still returns the instant the job ends, not at window end.
    monkeypatch.setattr(A, "_BASH_STATUS_BUSY_POLL_WAIT_S", 30.0)
    _patch(monkeypatch, _FakeJob(finish_at=0.3))
    _status(job_id="j1")                            # seed history (running)
    t = time.monotonic()
    out = _status(job_id="j1")                      # throttled, but job ends ~0.3s
    dt = time.monotonic() - t
    assert out["running"] is False
    assert dt < 3.0
