"""An orchestration stage waited on its delegates with no deadline at all.

The tool-call fan-out learned this already: it collects with
``as_completed`` under one deadline for the whole fan-out, and its
docstring says why -- waiting in submission order bought nothing, because
the pool runs them concurrently.

The orchestration stage kept the older shape and had no deadline of any
kind. One delegate that stalls holds every sibling's finished work, the
stage, and the turn, for as long as it stalls. The child's own wall-clock
guard fires per streamed event, so a fully SILENT stream never trips it
and there is nothing else to end the wait.

What is asserted here is the property, not the plumbing: the stage comes
back inside its bound, the delegates that finished are in it with their
results, and the one that did not is named as not having finished rather
than dropped or invented.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import subagents as sa


@pytest.fixture
def quick_deadline(monkeypatch):
    monkeypatch.setattr(sa, "_orchestration_stage_timeout", lambda: 0.6)


def _jobs(n: int) -> list[dict]:
    return [{"subagent_type": "explore", "prompt": f"p{i}",
             "description": f"job {i}"} for i in range(n)]


def _run(monkeypatch, jobs, behaviour):
    """Run the stage with ``behaviour(job) -> payload`` per delegate."""
    monkeypatch.setattr(sa, "_ORCH_MAX_WORKERS", 4, raising=False)
    calls = {"n": 0}

    def _fake_one(job):
        calls["n"] += 1
        return behaviour(job)

    # The stage helper is a closure over _one; drive it through the public
    # shape instead by patching what _one calls would be invasive, so the
    # collection function is exercised directly.
    return sa._collect_stage(jobs, _fake_one), calls


# ---------------------------------------------------------------------------

def test_a_stalled_delegate_does_not_hold_the_finished_ones(quick_deadline,
                                                            monkeypatch):
    def _behaviour(job):
        if job["description"] == "job 1":
            time.sleep(30)                     # the straggler
        return {"subagent_type": "explore", "description": job["description"],
                "result": f"done {job['description']}", "error": ""}

    started = time.monotonic()
    results, _calls = _run(monkeypatch, _jobs(3), _behaviour)
    elapsed = time.monotonic() - started

    assert elapsed < 10, f"the stage waited {elapsed:.1f}s for a straggler"
    assert len(results) == 3
    assert results[0]["result"] == "done job 0"
    assert results[2]["result"] == "done job 2"


def test_the_one_that_did_not_finish_is_named_as_such(quick_deadline,
                                                      monkeypatch):
    """Not dropped, and not invented: the caller has to be able to tell
    which of six calls never came back."""
    def _behaviour(job):
        if job["description"] == "job 1":
            time.sleep(30)
        return {"subagent_type": "explore", "description": job["description"],
                "result": "ok", "error": ""}

    results, _calls = _run(monkeypatch, _jobs(3), _behaviour)
    straggler = results[1]
    assert straggler["description"] == "job 1"
    assert straggler["result"] == ""
    assert "did not finish" in straggler["error"]
    assert "left running" in straggler["error"]


def test_a_delegate_that_raises_is_reported_in_its_own_slot(quick_deadline,
                                                            monkeypatch):
    def _behaviour(job):
        if job["description"] == "job 0":
            raise RuntimeError("delegate exploded")
        return {"subagent_type": "explore", "description": job["description"],
                "result": "ok", "error": ""}

    results, _calls = _run(monkeypatch, _jobs(2), _behaviour)
    assert "delegate exploded" in results[0]["error"]
    assert results[1]["result"] == "ok"


def test_the_order_of_the_results_is_the_order_of_the_jobs(quick_deadline,
                                                           monkeypatch):
    """Collected as they finish, returned as they were asked for --
    ``{{stage:NAME}}`` substitution and the caller both index by job."""
    def _behaviour(job):
        # Finish in reverse order of submission.
        time.sleep(0.05 * (3 - int(job["description"].split()[-1])))
        return {"subagent_type": "explore", "description": job["description"],
                "result": job["description"], "error": ""}

    results, _calls = _run(monkeypatch, _jobs(3), _behaviour)
    assert [r["result"] for r in results] == ["job 0", "job 1", "job 2"]


def test_a_single_job_still_runs_inline(monkeypatch):
    """No pool, no deadline machinery, for the common case."""
    def _behaviour(job):
        return {"subagent_type": "explore", "description": job["description"],
                "result": "alone", "error": ""}

    results, calls = _run(monkeypatch, _jobs(1), _behaviour)
    assert results == [{"subagent_type": "explore", "description": "job 0",
                        "result": "alone", "error": ""}]
    assert calls["n"] == 1
