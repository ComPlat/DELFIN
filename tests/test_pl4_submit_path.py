"""PL4 phase-3 red control: submit-path wiring (bench_nightly P3).

Pin (brief target 3): jobs go through DELFIN's own submit path with
node-local staging — the staging lines come from
delfin.cluster_bench.slurm_script._cbatch_staging_lines (reused, never
copied), and the record shows which path each instance took.
"""
import json
import tempfile

from delfin.agent.bench_nightly import run
from delfin.agent.benchmark import Task
from delfin.cluster_bench.slurm_script import _cbatch_staging_lines

import tests.test_pl4_budget_and_seed as base


def _config(tmp, tasks, **extra):
    cfg = base._config(tmp, tasks)
    cfg.update(extra)
    return cfg


def _tasks():
    return base._tasks()


def test_submit_path_carries_staging_lines():
    """With a tool_python and a submitter, the submitter receives the
    REAL staging lines from _cbatch_staging_lines (not a copy)."""
    submitted = {}

    def submitter(task, staging_lines, job_name):
        submitted["lines"] = staging_lines
        submitted["name"] = job_name
        return "job-123"

    tasks = _tasks()
    with tempfile.TemporaryDirectory() as tmp:
        rec = run(_config(tmp, tasks,
                          tool_python="/opt/venv/bin/python",
                          submitter=submitter))
    assert submitted["lines"] == _cbatch_staging_lines(
        "/opt/venv/bin/python")
    assert len(submitted["lines"]) > 0
    assert submitted["name"].startswith("delfin-bench-")


def test_submit_path_recorded_on_instance():
    tasks = _tasks()
    with tempfile.TemporaryDirectory() as tmp:
        rec = run(_config(tmp, tasks,
                          tool_python="/opt/venv/bin/python",
                          submitter=lambda task, lines, name: "job-42"))
    inst = rec.instances[0]
    assert inst.ran
    assert inst.job_id == "job-42"
    assert inst.submit_path == "slurm"


def test_no_submitter_runs_locally():
    """Without a submitter the runner runs in-process and the record
    says so (the default path all earlier phases pin)."""
    tasks = _tasks()
    with tempfile.TemporaryDirectory() as tmp:
        rec = run(_config(tmp, tasks))
    inst = rec.instances[0]
    assert inst.ran
    assert inst.submit_path == "local"
    assert inst.job_id is None


def test_no_tool_python_empty_staging_lines():
    """s19's edge: _cbatch_staging_lines returns [] without a tool
    interpreter — the submitter must still be callable and receive []."""
    submitted = {}

    def submitter(task, staging_lines, job_name):
        submitted["lines"] = staging_lines
        return "job-7"

    tasks = _tasks()
    with tempfile.TemporaryDirectory() as tmp:
        rec = run(_config(tmp, tasks, tool_python=None,
                          submitter=submitter))
    assert submitted["lines"] == []
    assert rec.instances[0].job_id == "job-7"
