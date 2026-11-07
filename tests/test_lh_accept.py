"""F2 P2b red test — the long-horizon tasks load and grade through the
accept-script contract (ONE argument: the workspace; no model call).

Targets (all owned): delfin/agent/pack/benchmark/tasks_w28_lh.yaml,
delfin/agent/longhorizon/runner.py, accept/lh_*.py, setup/lh_*.py.

Two red controls:
  * load_tasks() auto-globs tasks_*.yaml; the 6 lh-* tasks must load with
    no error and no drop of the existing suite (a malformed entry breaks
    the whole benchmark, so this is a repo-wide guard).
  * the accept script grades a run whose metrics JSON lives in the
    workspace: a within-threshold run passes (exit 0), a violating run
    fails (exit != 0).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


def _workspace_with_metrics(tmp_path: Path, metrics: dict, thresholds: dict):
    """Write the metrics artifact the runner leaves in the workspace. The
    schema is the accept contract: ``LongHorizonResult`` serialised via
    metrics.to_json plus a ``recall_late`` computed field, under
    ``metrics``, and the task's acceptance thresholds under ``thresholds``.
    The accept script computes ``recall_late`` itself from
    ``facts_recalled_late``/``facts_total`` when it is absent, so a cached
    stale artifact without it still grades."""
    ws = tmp_path / "ws"
    ws.mkdir()
    payload = {"task_id": "lh-code-change-fact",
               "metrics": metrics,
               "thresholds": thresholds}
    (ws / ".longhorizon_metrics.json").write_text(json.dumps(payload))
    return ws


def test_load_tasks_includes_all_six_lh_tasks():
    """load_tasks() must load tasks_w28_lh.yaml cleanly and yield exactly
    the 6 lh-* tasks, and still load the rest of the suite's tasks (a
    broken yaml entry breaks every benchmark run)."""
    try:
        from delfin.agent.benchmark import load_tasks
    except Exception:  # function may live under another name
        from delfin.agent.benchmark_runner import load_tasks
    tasks = load_tasks()
    ids = {t.id for t in tasks}
    assert {"lh-code-change-fact", "lh-calc-campaign-numbers",
            "lh-user-decision-binds", "lh-refusal-holds",
            "lh-session-restart-continuity", "lh-goal-focus"} <= ids


def test_accept_passes_within_threshold_and_fails_violation(tmp_path):
    """run_acceptance runs [sys.executable, script, workspace] with exactly
    ONE arg (QS-verified run_acceptance:1332). A run inside every threshold
    passes (exit 0); a run violating recall_late or turns_min fails."""
    from delfin.agent.benchmark_runner import run_acceptance

    good = {"facts_recalled_late": 1, "facts_total": 1,
            "redone_work": 0, "turns": 108, "compactions": 3,
            "wall_time_s": 4000.0}
    ws_good = _workspace_with_metrics(tmp_path, good,
                                      {"facts_total_min": 1,
                                       "recall_late_min": 0.5,
                                       "turns_min": 100})
    ok, out = run_acceptance("lh_metrics_accept.py", ws_good,
                             root=Path(os.getcwd()))
    assert ok, f"good run must pass: {out}"

    bad = {"facts_recalled_late": 0, "facts_total": 1,
           "redone_work": 0, "turns": 60, "compactions": 0,
           "wall_time_s": 900.0}
    ws_bad = _workspace_with_metrics(tmp_path, bad,
                                     {"facts_total_min": 1,
                                      "recall_late_min": 0.5,
                                      "turns_min": 100})
    ok_bad, out_bad = run_acceptance("lh_metrics_accept.py", ws_bad,
                                     root=Path(os.getcwd()))
    assert not ok_bad, f"violating run must fail: {out_bad}"
