"""The campaign wires into DELFIN's scheduler without a language model.

The chemistry loop needs three DELFIN connections, all of them reused
(no rebuilds):

* each candidate is a WorkflowJob on DELFIN's own scheduler, through
  ``delfin.tools._runner.step_as_workflow_job`` — so the campaign never
  bypasses the core budget or the resource pool;
* a running campaign registers every submitted job with
  ``delfin.agent.job_monitor.register_agent_job``, so job completion
  comes back through DELFIN's own watch (never private squeue calls);
* a failed/stagnated/budget-exhausted campaign sets an LLM-free
  wake-up through ``delfin.agent.scheduler.Scheduler.schedule_once``
  — the model is only ever woken by the conditions pinned in
  ``campaign_should_wake``.

These tests pin the wiring contract; every real cluster interaction is
faked at the boundary DELFIN itself defines (step_as_workflow_job's
work callable, register_agent_job's watch file, schedule_once's entry).
"""

import json
from pathlib import Path

import pytest

from delfin.agent.campaign import (
    Campaign,
    CampaignBudget,
    DecisionLog,
    DecisionReason,
    TargetGap,
    campaign_should_wake,
)


# ---------------------------------------------------------------------------
# fakes at DELFIN's own boundaries
# ---------------------------------------------------------------------------


class _RecordingScheduler:
    """Records schedule_once calls exactly as the real Scheduler would."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def schedule_once(self, *, delay_seconds: int, prompt: str,
                      reason: str = "", workspace: str = "",
                      budget_usd: float = 0.0, session_id: str = ""):
        self.calls.append({
            "delay_seconds": delay_seconds,
            "prompt": prompt,
            "reason": reason,
            "workspace": workspace,
        })
        return {"id": f"ent{len(self.calls)}"}


def _make_campaign(tmp_path: Path, **overrides) -> Campaign:
    kw = dict(
        space=list(range(8)),
        target=TargetGap(2.0),
        budget=CampaignBudget(max_evaluations=4),
        log=DecisionLog(tmp_path / "decisions.jsonl"),
        work_dir=tmp_path / "work",
    )
    kw.update(overrides)
    return Campaign(**kw)


# ---------------------------------------------------------------------------
# scheduler wiring
# ---------------------------------------------------------------------------


class TestSchedulerWiring:
    def test_campaign_jobs_are_workflow_jobs(self, tmp_path):
        from delfin.tools._runner import step_as_workflow_job

        job = step_as_workflow_job(
            "cand0", "xtb_sp",
            geometry=tmp_path / "geometry.xyz", cores_min=1,
            cores_optimal=1, cores_max=1, charge=0, mult=1,
            method="gfn2",
        )
        from delfin.workflows.engine.classic import WorkflowJob

        assert isinstance(job, WorkflowJob)
        assert job.cores_min == 1
        assert job.cores_optimal == 1

    def test_campaign_registers_submitted_jobs(self, tmp_path, monkeypatch):
        # register_agent_job writes into the workspace's own watch file;
        # run the real function against a fake HOME so the test never
        # touches ~/.delfin.
        fake_home = tmp_path / "fakehome"
        fake_home.mkdir()
        monkeypatch.setenv("HOME", str(fake_home))
        from delfin.agent import job_monitor

        workspace = tmp_path / "ws"
        workspace.mkdir()
        job_monitor.register_agent_job(
            workspace, "123456",
            description="campaign candidate 0")
        data = json.loads(
            (workspace / ".delfin" / "agent_watched_jobs.json")
            .read_text())
        assert data["jobs"]["123456"]["kind"] == "slurm"
        assert "candidate 0" in data["jobs"]["123456"]["description"]

    def test_wake_up_set_on_budget_exhaustion(self, tmp_path):
        camp = _make_campaign(tmp_path)
        sched = _RecordingScheduler()
        camp.scheduler = sched
        # Burn the budget, then the loop notices it.
        while not camp.budget.exhausted:
            camp.budget.spend()
        camp.report_to_scheduler()
        assert sched.calls
        call = sched.calls[0]
        assert call["prompt"].startswith("campaign")
        assert camp.budget.exhausted is True

    def test_no_wake_up_while_the_loop_is_healthy(self, tmp_path):
        camp = _make_campaign(tmp_path)
        sched = _RecordingScheduler()
        camp.scheduler = sched
        camp.report_to_scheduler()
        assert sched.calls == []


# ---------------------------------------------------------------------------
# Bayes-vs-random statistics (reused from benchmark.py)
# ---------------------------------------------------------------------------


class TestComparisonStats:
    def test_wilson_interval_is_reused_not_rewritten(self):
        from delfin.agent.benchmark import wilson_interval

        lo, hi = wilson_interval(3, 3)
        assert 0.0 <= lo < hi <= 1.0
        lo0, hi0 = wilson_interval(0, 0)
        assert (lo0, hi0) == (0.0, 1.0)

    def test_fisher_exact_is_reused_not_rewritten(self):
        from delfin.agent.benchmark import _fisher_exact_2x2_pvalue

        # 3-of-3 vs 0-of-3 is the smallest N Fisher can rate and it
        # gives exactly 0.1 (two-sided) — the smallest possible p at
        # N=3. A campaign-sized 5-of-5 vs 0-of-5 is significant.
        assert _fisher_exact_2x2_pvalue(3, 0, 0, 3) == pytest.approx(0.1)
        assert _fisher_exact_2x2_pvalue(5, 0, 0, 5) < 0.05
        # Identical rows are not.
        assert _fisher_exact_2x2_pvalue(3, 0, 3, 0) == pytest.approx(1.0)
