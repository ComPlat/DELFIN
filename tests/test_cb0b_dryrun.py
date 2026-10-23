"""Red control for CB0B phase 3: the offline dry run.

The brief requires: "prove it with a dry run on two practice instances
that needs no cluster (fixture outputs)". The dry run exercises the
same public path a runner would: run the setup scripts, read the task
records from the packaged YAML through load_tasks(), check that each
instance's grading levels match the sampler's, and grade the fixture
OUTPUTS (synthetic, written by the test, shaped like real outputs)
against the task's grading contract — internal consistency only, no
cluster, no ORCA/xtb call.
"""
from __future__ import annotations

import pathlib

import pytest
import yaml

from delfin.agent import bench_holdout as bh
from delfin.agent.benchmark import load_tasks

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_W20 = _PACK / "tasks_w20_workflows.yaml"


def _w20_tasks():
    data = yaml.safe_load(_W20.read_text(encoding="utf-8"))
    return data.get("tasks", [])


def test_w20_tasks_load_through_the_public_loader():
    ids = {t.id for t in load_tasks() if t.id.startswith("w20_")}
    assert ids == {t["id"] for t in _w20_tasks()}


def _run_setup(script: str, ws: pathlib.Path) -> None:
    import subprocess
    import sys
    r = subprocess.run(
        [sys.executable, str(_PACK / "setup" / script), str(ws)],
        capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr


def _grade_consistency(task_dir: pathlib.Path, task: dict) -> list[str]:
    """Level (a) grading on fixture outputs: the numbers a run reports
    must be internally consistent with the contract the task carries.
    Returns a list of failures (empty = pass)."""
    failures = []
    trap = task.get("trap", "")
    if trap:
        # Generic trap grading: the fixture answer must show the agent
        # noticed THIS task's designed pitfall (per F3, every task that
        # carries the trap level has a presentable trap).
        answer = (task_dir / "answer.txt").read_text(encoding="utf-8")
        key = trap.split()[0].lower()  # 'odd', 'bright', 'spin', ...
        if key not in answer.lower():
            failures.append(f"trap not reported: {trap} unmentioned")
    return failures


def test_dry_run_two_instances_execute_and_grade(tmp_path):
    tasks = {t["id"]: t for t in _w20_tasks()}
    # Two practice instances, one with a trap to check, one without.
    chosen = [tasks["w20_redox_1"], tasks["w20_excited_states_1"]]
    for k, t in enumerate(chosen):
        # w20_rest.py builds several folders and refuses if ANY exists,
        # so each instance gets a fresh workspace root.
        ws = tmp_path / f"run{k}"
        ws.mkdir()
        _run_setup(t["setup"], ws)
        d = ws / t["id"]
        assert (d / "start.xyz").is_file()
        assert (d / "brief.txt").is_file()
        brief = (d / "brief.txt").read_text(encoding="utf-8")
        # The brief is the sampler's prompt, so the YAML cannot drift.
        assert t["prompt"].split(". ", 1)[-1] in brief
        # Grading contract matches the sampler's levels for this
        # INSTANCE (the sampler may add 'trap' to the family's set for
        # an open-shell/charged molecule).
        fid = "_".join(t["id"].split("_")[1:-1])
        inst = next(i for i in bh.sample(bh.PUBLIC_SEED, len(bh.FAMILIES))
                    if i.family_id == fid)
        assert t["grading_levels"] == list(inst.grading_levels)
        # Level (a) grading against fixture-shaped outputs: for every
        # trapped task, an answer that notices nothing must fail and an
        # answer naming the trap must pass (F3 made more tasks carry
        # the trap level, so this is per-task, not redox-specific).
        if t.get("trap"):
            (d / "answer.txt").write_text("no issue noticed\n", encoding="utf-8")
            assert _grade_consistency(d, t), "a noticed-nothing answer must fail"
            (d / "answer.txt").write_text(
                f"noticed: {t['trap']}\n", encoding="utf-8")
            assert not _grade_consistency(d, t)


def test_setup_refuses_to_overwrite(tmp_path):
    tasks = _w20_tasks()
    _run_setup(tasks[1]["setup"], tmp_path)
    import subprocess
    import sys
    r = subprocess.run(
        [sys.executable, str(_PACK / "setup" / tasks[1]["setup"]), str(tmp_path)],
        capture_output=True, text=True, timeout=60)
    assert r.returncode == 1
    assert "refusing" in r.stderr
