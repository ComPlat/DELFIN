"""Package SW1 phase 3b: benchmark fixtures + tasks + skill.

Control: tests/test_sw1_benchmark.py, 6 cases, all 6 red on the
previous commit (tasks_w27_sw1.yaml did not exist).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import yaml

BENCH = Path("delfin/agent/pack/benchmark")


def _run(script: Path, ws: Path) -> subprocess.CompletedProcess:
    # The repo root on PYTHONPATH: the scripts must exercise THIS
    # checkout's delfin (with the new swarm package), not an installed
    # delfin from a venv, which would silently grade an older tree.
    env = dict(os.environ)
    repo_root = str(Path.cwd().resolve())
    env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run([sys.executable, str(script), str(ws)],
                            capture_output=True, text=True, timeout=120,
                            env=env)
    return result


def test_task_yaml_loads_and_is_wellformed():
    data = yaml.safe_load((BENCH / "tasks_w27_sw1.yaml").read_text(
        encoding="utf-8"))
    assert data["package"] == "SW1"
    assert len(data["tasks"]) == 3
    for t in data["tasks"]:
        assert t["task_class"] == "other"
        assert t["setup"].startswith("sw1_")
        assert t["verify"].startswith("sw1_")


def test_setup_builds_fixtures_and_refuses_overwrite(tmp_path):
    setup = BENCH / "setup" / "sw1_spec_fixtures.py"
    assert _run(setup, tmp_path).returncode == 0
    root = tmp_path / "swarm_fixtures"
    assert (root / "valid_swarm.yaml").exists()
    assert (root / "triage_queue.json").exists()
    # Second run refuses (a setup without argument must print usage and
    # exit non-zero is covered by the usage check here and below).
    assert _run(setup, tmp_path).returncode == 1


def test_setup_without_argument_prints_usage_and_fails():
    result = subprocess.run([sys.executable,
                             str(BENCH / "setup" / "sw1_spec_fixtures.py")],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "usage" in result.stderr.lower()


def test_accept_all_three_tasks_pass(tmp_path):
    setup = BENCH / "setup" / "sw1_spec_fixtures.py"
    assert _run(setup, tmp_path).returncode == 0
    for name in ("sw1_check_overlapping_scopes.py", "sw1_check_dry_run.py",
                 "sw1_check_triage.py"):
        assert _run(BENCH / "accept" / name, tmp_path).returncode == 0


def test_accept_without_fixtures_fails(tmp_path):
    result = subprocess.run(
        [sys.executable,
         str(BENCH / "accept" / "sw1_check_overlapping_scopes.py"),
         str(tmp_path)],
        capture_output=True, text=True)
    assert result.returncode != 0


def test_skill_declares_domain_and_carries_no_permission_widening():
    skill = Path("delfin/agent/pack/skills/swarm-orchestrate.md").read_text(
        encoding="utf-8")
    assert "domains:" in skill.split("---")[1]
    lowered = skill.lower()
    for banned in ("~/.ssh", "api_key", "password", "skip the check",
                   "bypasspermissions"):
        assert banned not in lowered, banned
