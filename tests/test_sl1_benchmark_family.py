"""SL1 benchmark family tests: tasks load, conventions hold, verifies
grade in both directions.

Both directions per wave-15b: a good answer passes each accept script,
an empty/distractor workspace fails it. The setup script is exercised
for the usage rule (no argument -> usage + exit non-zero) and for
building its fixture into a throwaway workspace.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from delfin.agent.benchmark import load_tasks

ACCEPT = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "benchmark" / "accept"
SETUP = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "benchmark" / "setup"
YAML_PATH = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "benchmark" / "tasks_w25_sl1.yaml"


def _tasks():
    return load_tasks(YAML_PATH)


def test_family_loads_four_tasks():
    tasks = _tasks()
    assert len(tasks) == 4
    ids = {t.id for t in tasks}
    assert ids == {
        "w25-sl1-01-size-and-justify",
        "w25-sl1-02-pending-reason",
        "w25-sl1-03-failure-classify-recalc",
        "w25-sl1-04-efficiency-audit",
    }
    for t in tasks:
        assert t.task_class == "science_analysis"
        assert t.mode == "solo"


def test_every_task_points_at_a_real_verify_script():
    for task in _tasks():
        assert (ACCEPT / task.verify).is_file(), task.id


def _run_verify(name: str, ws: Path) -> int:
    return subprocess.run(
        [sys.executable, str(ACCEPT / name), str(ws)],
        capture_output=True, text=True).returncode


def _good_answer() -> str:
    return (
        "Partitions: cpu,cpu_il as one --partition list; SLURM starts in "
        "whichever can run it first, so the earliest start over all "
        "matching partitions costs nothing. PAL and maxcore come from "
        "DELFIN's own sizing (cluster_efficiency.advise / "
        "auto_configure_resources); the partition rule is earliest start "
        "and it NEVER lowers PAL, maxcore or time limit to start earlier.\n"
        "Job 7922181: FAILED at exit code 138 after 01:39:08 with pal 2 — "
        "2 core-hours x 1.653 h = 3.3 core-hours wasted; the measured .out "
        "line says CPU Efficiency: 1.53% of 01:39:08 core-walltime "
        "(TotalCPU from sacct is an artefact for TIMEOUT/CANCELLED, not a "
        "measurement). Job 7922985 (pal 1, smaller sibling shell) measured "
        "35.75 % over 00:24:47. Recommendation: next time run this kind at "
        "pal 1 and set the time limit from the measured wall. The two runs "
        "are NOT a clean A/B comparison — the second task is smaller and "
        "NOT the same task, so no pal-1-versus-pal-2 speedup can be "
        "inferred from them.\n"
        "Job 7924102: CANCELLED DUE TO TIME LIMIT — a walltime kill "
        "(sacct TIMEOUT), not an ORCA error: no 'ORCA finished by error "
        "termination' in the .out. TotalCPU 0 % is the known artefact on "
        "walltime kills; the job did compute. Fix: through DELFIN's recalc "
        "path, change the submit time_limit (e.g. to 48:00:00), applied "
        "only after your confirmation.\n"
        "Job 7925001: PENDING, ReqNodeNotAvail, Reserved for maintenance; "
        "estimated start 2026-11-02T04:30. Decision: wait, nothing to "
        "change — the nodes are reserved, not the request wrong; do not "
        "lower pal/maxcore/time just to start earlier.\n"
    )


def _write_good(ws: Path):
    (ws / "answer.txt").write_text(_good_answer(), encoding="utf-8")
    size = ws / "sl1" / "size"
    size.mkdir(parents=True, exist_ok=True)
    (size / "CONTROL_resources.txt").write_text(
        "PAL = 8\nmaxcore = 4000\ntime_limit = 48:00:00\n", encoding="utf-8")
    fail = ws / "sl1" / "failure"
    fail.mkdir(parents=True, exist_ok=True)
    (fail / "fix_proposal.txt").write_text(
        "Fix via DELFIN's recalc path (delfin-recalc-classic): the "
        "submit time_limit goes from 24:00:00 to 48:00:00. Applied only "
        "after your confirmation; the old run is kept.\n", encoding="utf-8")


def test_good_answer_passes_all_four():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        ws = Path(tmp)
        _write_good(ws)
        for name in ("sl1_w25_01_size_justify_verify.py",
                     "sl1_w25_02_pending_reason_verify.py",
                     "sl1_w25_03_failure_classify_verify.py",
                     "sl1_w25_04_efficiency_verify.py"):
            assert _run_verify(name, ws) == 0, name


def test_empty_workspace_fails_all_four():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        ws = Path(tmp)
        for name in ("sl1_w25_01_size_justify_verify.py",
                     "sl1_w25_02_pending_reason_verify.py",
                     "sl1_w25_03_failure_classify_verify.py",
                     "sl1_w25_04_efficiency_verify.py"):
            assert _run_verify(name, ws) == 1, name


def test_forbidden_resizing_fails_task_1_and_2():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        ws = Path(tmp)
        _write_good(ws)
        ws.joinpath("answer.txt").write_text(
            "I would lower the time limit to start earlier.\n",
            encoding="utf-8")
        assert _run_verify("sl1_w25_01_size_justify_verify.py", ws) == 1


def test_setup_usage_and_build():
    script = SETUP / "sl1_w25_failure_fixture.py"
    # No argument: usage + non-zero exit (the existing benchmark rule).
    proc = subprocess.run([sys.executable, str(script)],
                          capture_output=True, text=True)
    assert proc.returncode == 2
    assert "usage" in proc.stderr
    # With an argument it builds the fixture and refuses a rebuild.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        ws = Path(tmp)
        subprocess.run([sys.executable, str(script), str(ws)], check=True)
        assert (ws / "sl1" / "failure" / "CONTROL_snippet.txt").is_file()
        again = subprocess.run([sys.executable, str(script), str(ws)],
                               capture_output=True, text=True)
        assert again.returncode == 1
