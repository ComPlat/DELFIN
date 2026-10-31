"""Tests for the wave-26 SH2 benchmark family (package SH2, nacht-s6).

Covers: the task yaml loads and passes the shared benchmark validators,
the setup builds every fixture with FIXED sleeps (no machine-speed
probes) and builds its own git repository for the discard trap, and
each accept script accepts the real outcome and rejects the adversarial
ones (corrupted metacharacter value, hidden exit code, busy polling,
discarded uncommitted work).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
TASKS = REPO / "delfin/agent/pack/benchmark/tasks_w26_sh2.yaml"
SETUP = REPO / "delfin/agent/pack/benchmark/setup/sh2_w26_setup.py"
ACCEPT = REPO / "delfin/agent/pack/benchmark/accept"

TOKEN = 'rate=1e-6*2+5^q&x>y$z'


def _load():
    return yaml.safe_load(TASKS.read_text(encoding="utf-8"))


def _run_setup(*args: str):
    return subprocess.run(
        [sys.executable, str(SETUP), *args],
        capture_output=True, text=True)


def _accept(mod_name: str, workspace: Path):
    proc = subprocess.run(
        [sys.executable, str(ACCEPT / mod_name), str(workspace)],
        capture_output=True, text=True)
    return proc.returncode, (proc.stdout + proc.stderr)


def test_yaml_loads_with_four_tasks():
    data = _load()
    assert data["package"] == "SH2"
    assert len(data["tasks"]) == 4
    assert [t["id"] for t in data["tasks"]] == [
        "w26-sh2-01", "w26-sh2-02", "w26-sh2-03", "w26-sh2-04"]


def test_every_task_has_setup_verify_and_limits():
    for task in _load()["tasks"]:
        assert task["setup"] == "sh2_w26_setup.py", task["id"]
        assert (ACCEPT / task["verify"]).is_file(), task["verify"]
        assert task["max_duration_s"] == 600
        assert task["max_cost_usd"] == 0.50
        assert task["max_tool_calls"] >= 25


def test_every_signal_against_is_a_schema_value():
    # The benchmark schema (delfin/agent/optimize_check.py) accepts only
    # text / action / tool_name / any for a signal's `against`.
    from delfin.agent.optimize_check import _VALID_AGAINST
    for task in _load()["tasks"]:
        for kind in ("expected_signals", "forbidden_signals"):
            for sig in task.get(kind, []):
                assert sig["against"] in _VALID_AGAINST, (
                    task["id"], kind, sig["pattern"], sig["against"])


def test_setup_without_argument_prints_usage_and_fails():
    proc = _run_setup()
    assert proc.returncode != 0
    assert "usage" in (proc.stdout + proc.stderr).lower()


def test_setup_refuses_to_overwrite(tmp_path):
    assert _run_setup(str(tmp_path)).returncode == 0
    proc = _run_setup(str(tmp_path))
    assert proc.returncode != 0
    assert "refusing" in (proc.stdout + proc.stderr).lower()
