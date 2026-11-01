"""Tests for the wave-25 MN1 benchmark family (package MN1, nacht-s9).

Covers: the task yaml loads and passes the shared benchmark validators,
the setup builds every fixture with FIXED sleeps (no machine-speed
probes), and each accept script accepts the real outcome and rejects
the adversarial ones (missing job, missing answer, busy polling,
blind retry, missing transcript record).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
TASKS = REPO / "delfin/agent/pack/benchmark/tasks_w25_mn1.yaml"
SETUP = REPO / "delfin/agent/pack/benchmark/setup/mn1_w25_setup.py"
ACCEPT = REPO / "delfin/agent/pack/benchmark/accept"


def _load():
    return yaml.safe_load(TASKS.read_text(encoding="utf-8"))


def test_yaml_loads_with_four_tasks():
    data = _load()
    assert data["package"] == "MN1"
    assert len(data["tasks"]) == 4
    assert [t["id"] for t in data["tasks"]] == [
        "w25-mn1-01", "w25-mn1-02", "w25-mn1-03", "w25-mn1-04"]


def test_every_task_has_setup_verify_and_limits():
    for task in _load()["tasks"]:
        assert task["setup"] == "mn1_w25_setup.py"
        assert (ACCEPT / task["verify"]).is_file(), task["verify"]
        assert task["max_duration_s"] == 600
        assert task["max_cost_usd"] == 0.50
        assert task["max_tool_calls"] >= 25


def test_every_signal_against_is_a_schema_value():
    # The benchmark schema (delfin/agent/optimize_check.py:40) accepts
    # only text / action / tool_name / any for a signal's `against`.
    # The original family shipped `against: transcript` on every tool-call
    # signal; the real run red: test_optimize_check.py::test_real_repo_is_safe
    # reported 12 "invalid against='transcript'" errors. This pins all four
    # tasks to the schema.
    from delfin.agent.optimize_check import _VALID_AGAINST
    for task in _load()["tasks"]:
        for kind in ("expected_signals", "forbidden_signals"):
            for sig in task.get(kind, []):
                assert sig["against"] in _VALID_AGAINST, (
                    task["id"], kind, sig["pattern"], sig["against"])


def test_setup_without_argument_prints_usage_and_fails():
    proc = subprocess.run(
        [sys.executable, str(SETUP)], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "usage" in (proc.stdout + proc.stderr).lower()


def test_setup_builds_all_fixtures_with_fixed_sleeps(tmp_path):
    proc = subprocess.run(
        [sys.executable, str(SETUP), str(tmp_path)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    root = tmp_path / "mn1_w25"
    for name in ("wait_line/run_slow_job.sh", "appear_file/make_report.py",
                 "fail_job/flaky_job.sh",
                 "two_jobs/fast_job.sh", "two_jobs/slow_job.sh"):
        assert (root / name).is_file(), name
    # FIXED sleeps, never machine-speed probes.
    wait = (root / "wait_line/run_slow_job.sh").read_text()
    assert "sleep 7" in wait and "sleep 8" in wait
    flaky = (root / "fail_job/flaky_job.sh").read_text()
    assert "sleep 5" in flaky and "FATAL: disk quota exceeded" in flaky


# -- accept-script adversarial tests -------------------------------------

_SIGNALS = "bash_background\nbash_status\nbash_status\n"   # legal: 2 polls


def _accept(mod_name: str, workspace: Path):
    proc = subprocess.run(
        [sys.executable, str(ACCEPT / mod_name), str(workspace)],
        capture_output=True, text=True)
    return proc.returncode, (proc.stdout + proc.stderr)


def test_accept_01_accepts_the_real_outcome(tmp_path):
    ws = tmp_path
    folder = ws / "mn1_w25" / "wait_line"
    folder.mkdir(parents=True)
    (folder / "job.log").write_text(
        "PHASE ONE started at cycle 1\n"
        "PHASE TWO REACHED: the optimizer left the first basin\n")
    (ws / "answer.txt").write_text(
        "Phase two was reached because the optimizer left the first basin.")
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 0, out


def test_accept_01_rejects_missing_job_and_missing_signals(tmp_path):
    ws = tmp_path
    (ws / "mn1_w25" / "wait_line").mkdir(parents=True)
    (ws / "answer.txt").write_text("the optimizer left the first basin")
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 1 and "never ran" in out
    (ws / "mn1_w25" / "wait_line" / "job.log").write_text("PHASE TWO REACHED: x")
    (ws / "transcript_signals.txt").unlink()
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 1 and "unmeasured" in out


def test_accept_01_rejects_busy_polling(tmp_path):
    ws = tmp_path
    folder = ws / "mn1_w25" / "wait_line"
    folder.mkdir(parents=True)
    (folder / "job.log").write_text("PHASE TWO REACHED: reason here")
    (ws / "answer.txt").write_text("the optimizer left the first basin")
    (ws / "transcript_signals.txt").write_text(
        "bash_background\n" + "bash_status\n" * 4)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 1 and "busy polling" in out


def test_accept_02_rejects_missing_unit(tmp_path):
    ws = tmp_path
    folder = ws / "mn1_w25" / "appear_file"
    folder.mkdir(parents=True)
    (folder / "report.json").write_text(
        json.dumps({"n": 5, "mean_kj_per_mol": 41.73}))
    (ws / "answer.txt").write_text("n=5, mean 41.73")   # no unit
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", ws)
    assert rc == 1 and "unit" in out


def test_accept_03_rejects_blind_retry(tmp_path):
    ws = tmp_path
    folder = ws / "mn1_w25" / "fail_job"
    folder.mkdir(parents=True)
    (folder / "job.log").write_text(
        "FATAL: disk quota exceeded while writing checkpoint.bin\n")
    (ws / "answer.txt").write_text(
        "It failed: disk quota exceeded. I am retrying now.")
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_03_fail_report_verify.py", ws)
    assert rc == 1 and "re-run" in out


def test_accept_04_rejects_one_missing_log(tmp_path):
    ws = tmp_path
    folder = ws / "mn1_w25" / "two_jobs"
    folder.mkdir(parents=True)
    (folder / "fast.log").write_text("fast job done: 3 iterations, residual 1.2e-6")
    (ws / "answer.txt").write_text(
        "fast finished first (3 iterations, residual 1.2e-6); slow still running")
    (ws / "transcript_signals.txt").write_text(
        "bash_background\nbash_background\nbash_status\n")
    rc, out = _accept("mn1_w25_04_first_finish_verify.py", ws)
    assert rc == 1 and "slow.log" in out


# -- Reviewer nacht-s10 adversarial battery (F1-F4, commit d5fc80c1f) ----
# Five cases were RED against 3a6d7a66c (the reviewer's gate run); these
# pin the fixes.

def _ws_base(tmp_path, folder):
    ws = tmp_path
    (ws / "mn1_w25" / folder).mkdir(parents=True)
    (ws / "answer.txt").write_text("placeholder")
    return ws


def test_f1_accept_01_rejects_bash_output_polling(tmp_path):
    ws = _ws_base(tmp_path, "wait_line")
    (ws / "mn1_w25" / "wait_line" / "job.log").write_text(
        "PHASE TWO REACHED: the optimizer left the first basin")
    (ws / "answer.txt").write_text("the optimizer left the first basin")
    (ws / "transcript_signals.txt").write_text(
        "bash_background\n" + "bash_output\n" * 6)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 1 and "busy polling" in out


def test_f1b_accept_01_rejects_repeated_cat_of_the_log(tmp_path):
    ws = _ws_base(tmp_path, "wait_line")
    (ws / "mn1_w25" / "wait_line" / "job.log").write_text(
        "PHASE TWO REACHED: the optimizer left the first basin")
    (ws / "answer.txt").write_text("the optimizer left the first basin")
    (ws / "transcript_signals.txt").write_text(
        "bash_background\n" + "cat mn1_w25/wait_line/job.log\n" * 5)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", ws)
    assert rc == 1 and "busy polling" in out


def test_f2_accept_03_rejects_rerun_wording(tmp_path):
    ws = _ws_base(tmp_path, "fail_job")
    (ws / "mn1_w25" / "fail_job" / "job.log").write_text(
        "FATAL: disk quota exceeded while writing checkpoint.bin\n")
    (ws / "answer.txt").write_text(
        "It failed: disk quota exceeded. Next step: I will rerun the script.")
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_03_fail_report_verify.py", ws)
    assert rc == 1 and "re-run" in out


def test_f3_accept_02_refuses_a_json_list_cleanly(tmp_path):
    ws = _ws_base(tmp_path, "appear_file")
    (ws / "mn1_w25" / "appear_file" / "report.json").write_text("[1, 2, 3]")
    (ws / "transcript_signals.txt").write_text(_SIGNALS)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", ws)
    assert rc == 1 and "not a JSON object" in out


def test_f4_accept_02_rejects_directory_probe_evasion(tmp_path):
    ws = _ws_base(tmp_path, "appear_file")
    (ws / "mn1_w25" / "appear_file" / "report.json").write_text(
        '{"n": 5, "mean_kj_per_mol": 41.73}')
    (ws / "answer.txt").write_text("n=5, mean 41.73 kJ/mol")
    (ws / "transcript_signals.txt").write_text(
        "bash_background\n" + "ls mn1_w25/appear_file\n" * 4)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", ws)
    assert rc == 1 and "probes" in out
