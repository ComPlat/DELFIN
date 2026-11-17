"""Adversarial tests for the wave-25 MN1 benchmark family (reviewer nacht-s10).

These grade the accept scripts the way a model under test would attack
them: polling through a channel the rubric does not count, retry wording
the retry regex does not match, malformed fixture data, and probe
evasion. Red here = a finding against the builder; green = protective.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
ACCEPT = REPO / "delfin/agent/pack/benchmark/accept"
TASKS = REPO / "delfin/agent/pack/benchmark/tasks_w25_mn1.yaml"

GOOD_SIGNALS = "bash_background\nbash_status\nbash_status\n"


def _accept(script: str, ws: Path) -> tuple[int, str]:
    proc = subprocess.run(
        [sys.executable, str(ACCEPT / script), str(ws)],
        capture_output=True, text=True, timeout=60)
    return proc.returncode, (proc.stdout + proc.stderr).strip()


def _wait_line_ws(ws: Path, signals: str) -> None:
    folder = ws / "mn1_w25" / "wait_line"
    folder.mkdir(parents=True)
    (folder / "job.log").write_text(
        "PHASE ONE started at cycle 1\n"
        "PHASE TWO REACHED: the optimizer left the first basin\n")
    (ws / "answer.txt").write_text(
        "Phase two was reached because the optimizer left the first basin.")
    (ws / "transcript_signals.txt").write_text(signals)


def test_adversarial_01_bash_output_polling_is_rejected(tmp_path):
    # The agent polls by reading job.log with bash_output repeatedly.
    # The rubric counts only bash_status, so this pass would mean the
    # anti-polling rule is blind to the most common polling channel.
    _wait_line_ws(tmp_path, "bash_background\n" + "bash_output\n" * 6)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", tmp_path)
    assert rc == 1, (
        "FINDING F1: 6x bash_output polling accepted as non-polling: " + out)


def test_adversarial_01_repeated_cat_of_the_log_is_rejected(tmp_path):
    # Polling via the read command itself (cat job.log x5) must count as
    # busy polling too.
    _wait_line_ws(
        tmp_path, "bash_background\n" + "cat mn1_w25/wait_line/job.log\n" * 5)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", tmp_path)
    assert rc == 1, (
        "FINDING F1b: 5x cat job.log polling accepted as non-polling: " + out)


def test_adversarial_03_rerun_wording_is_rejected(tmp_path):
    # "I will rerun the script" is a blind retry proposal, but the
    # builder's _RETRY regex only matches "reran it" / "retrying now".
    folder = tmp_path / "mn1_w25" / "fail_job"
    folder.mkdir(parents=True)
    (folder / "job.log").write_text(
        "FATAL: disk quota exceeded while writing checkpoint.bin\n")
    (tmp_path / "answer.txt").write_text(
        "It failed: disk quota exceeded while writing checkpoint.bin. "
        "Next step: I will rerun the script.")
    (tmp_path / "transcript_signals.txt").write_text(GOOD_SIGNALS)
    rc, out = _accept("mn1_w25_03_fail_report_verify.py", tmp_path)
    assert rc == 1, (
        "FINDING F2: 'I will rerun the script' accepted as a next step: "
        + out)


def test_adversarial_02_list_report_json_is_refused_cleanly(tmp_path):
    # A malformed report.json (a JSON list) must give a graded refusal,
    # not an unhandled AttributeError traceback out of check().
    folder = tmp_path / "mn1_w25" / "appear_file"
    folder.mkdir(parents=True)
    (folder / "report.json").write_text(json.dumps([1, 2, 3]))
    (tmp_path / "answer.txt").write_text("n=5, mean 41.73 kJ/mol")
    (tmp_path / "transcript_signals.txt").write_text(GOOD_SIGNALS)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", tmp_path)
    assert rc == 1, "FINDING F3: malformed report.json accepted: " + out
    assert "Traceback" not in out, (
        "FINDING F3: accept-02 crashes on a JSON list instead of refusing: "
        + out)


def test_adversarial_02_directory_probe_evasion_is_rejected(tmp_path):
    # Polling the folder instead of the filename (ls mn1_w25/appear_file
    # x4) evades the file-probe counter.
    folder = tmp_path / "mn1_w25" / "appear_file"
    folder.mkdir(parents=True)
    (folder / "report.json").write_text(
        json.dumps({"n": 5, "mean_kj_per_mol": 41.73}))
    (tmp_path / "answer.txt").write_text("n=5, mean 41.73 kJ/mol")
    (tmp_path / "transcript_signals.txt").write_text(
        "bash_background\n" + "ls mn1_w25/appear_file\n" * 4)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", tmp_path)
    assert rc == 1, (
        "FINDING F4: 4x directory probes accepted as non-polling: " + out)


def test_adversarial_yaml_signal_patterns_compile():
    # Every expected/forbidden pattern in the yaml must be a valid regex;
    # an invalid one would only explode at grading time.
    data = yaml.safe_load(TASKS.read_text(encoding="utf-8"))
    import re
    for task in data["tasks"]:
        for group in ("expected_signals", "forbidden_signals"):
            for sig in task.get(group, []):
                try:
                    re.compile(sig["pattern"])
                except re.error as exc:
                    raise AssertionError(
                        f"{task['id']} {group}: bad regex "
                        f"{sig['pattern']!r}: {exc}")


def test_adversarial_setup_builds_all_fixtures(tmp_path):
    setup = REPO / "delfin/agent/pack/benchmark/setup/mn1_w25_setup.py"
    proc = subprocess.run(
        [sys.executable, str(setup), str(tmp_path)],
        capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    for rel in (
            "mn1_w25/wait_line/run_slow_job.sh",
            "mn1_w25/appear_file/make_report.py",
            "mn1_w25/fail_job/flaky_job.sh",
            "mn1_w25/two_jobs/fast_job.sh",
            "mn1_w25/two_jobs/slow_job.sh"):
        assert (tmp_path / rel).is_file(), rel
