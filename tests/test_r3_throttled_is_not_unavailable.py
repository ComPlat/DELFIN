"""A query the shared scheduler client throttled is not a degraded scheduler.

check_agent_jobs runs on every agent turn and every wake. Mapping the
client's closed gate to STATE_UNAVAILABLE made each poll inside the 25 s
interval report "the scheduler could not be asked" for every watched job.
"""
from __future__ import annotations

import json

from delfin import scheduler_client
from delfin.agent import job_monitor as jm


def _throttled(cmd):
    return scheduler_client.THROTTLED


def test_a_throttled_query_yields_throttled_not_unavailable():
    states = jm.query_job_states_detailed(["101", "102"], _throttled)
    assert states == {"101": jm.STATE_THROTTLED, "102": jm.STATE_THROTTLED}
    assert jm.query_job_states(["101"], _throttled) == {}


def test_a_throttled_squeue_still_lets_sacct_answer():
    def run(cmd):
        if cmd[0] == "squeue":
            return scheduler_client.THROTTLED
        return "101 COMPLETED\n"
    states = jm.query_job_states_detailed(["101", "102"], run)
    assert states["101"] == "COMPLETED"
    assert states["102"] == ""


def test_a_throttled_poll_reports_nothing_and_keeps_the_watch(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_emit_watch_attention", lambda *a, **k: None)
    path = jm._agent_watch_path(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"jobs": {"101": {
        "kind": "slurm", "description": "opt", "added_at": 9e18}}}))
    assert jm.check_agent_jobs(tmp_path, _throttled) == []
    assert "101" in jm.load_watched(path)["jobs"]


def test_the_default_runner_passes_the_gate_through(monkeypatch):
    client = scheduler_client.SchedulerClient(runner=lambda c: "1 RUNNING\n")
    monkeypatch.setattr(scheduler_client, "default_client", client)
    assert jm._default_run(["squeue", "-j", "1"]) == "1 RUNNING\n"
    assert jm._default_run(["squeue", "-j", "2"]) is scheduler_client.THROTTLED
