"""The eval report says whether the second run of a learning task was better.

``learning_metric`` is a pure function and ``compare_learning`` pairs two run
files; neither changes anything a person sees until something calls them.
The caller is the daily eval report: next to the benchmark drift note it
states, for the pre-registered learning tasks only, what the second run did
against the first. These tests drive ``eval_loop.build_report`` with real run
files and read the report text.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import benchmark as bm
from delfin.agent import eval_loop
from conftest import child_env

LEARNING = "l5_read_the_value_from_file"
ORDINARY = "office_summarise_a_table"


def _row(task_id, *, success=True, tool_calls=6, input_tokens=300,
         output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=True,
         learning_state="s"):
    """A five-replicate row with identical replicates."""
    return {"task_id": task_id, "model": "m",
            "learning_state": learning_state,
            "success": success, "quality_0_100": 90,
            "tool_calls": tool_calls, "input_tokens": input_tokens,
            "output_tokens": output_tokens, "duration_s": duration_s,
            "cost_usd": cost_usd, "verify_ok": verify_ok, "n_samples": 5,
            "per_run_tool_calls": [tool_calls] * 5,
            "per_run_tokens_total": [input_tokens + output_tokens] * 5}


@pytest.fixture
def runs(tmp_path, monkeypatch):
    """A runs directory of its own; ``write(first_rows, second_rows)``."""
    d = tmp_path / "benchmark_runs"
    d.mkdir()
    monkeypatch.setattr(bm, "_DEFAULT_RUNS_DIR", d)

    def write(first, second):
        now = time.time()
        for name, rows, age in (("run_1.jsonl", first, 100),
                                ("run_2.jsonl", second, 0)):
            p = d / name
            p.write_text("\n".join(json.dumps(r) for r in rows) + "\n",
                         encoding="utf-8")
            os.utime(p, (now - age, now - age))
    return write


def _learning_section(report: str) -> str:
    if "## Learning" not in report:
        return ""
    return report.split("## Learning", 1)[1].split("\n## ", 1)[0]


def test_fewer_tool_calls_on_the_second_run_is_reported_as_learned(runs):
    runs([_row(LEARNING, tool_calls=9)], [_row(LEARNING, tool_calls=4)])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "1 learned" in section
    assert f"{LEARNING}: learned (median tool calls 9→4" in section
    assert "N=5/5" in section


def test_the_same_counted_work_is_no_learning_whatever_the_clock_says(runs):
    runs([_row(LEARNING, duration_s=41.7, cost_usd=0.31)],
         [_row(LEARNING, duration_s=12.3, cost_usd=0.11)])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "0 learned, 1 no learning" in section
    assert f"{LEARNING}: no learning" in section


def test_a_cheaper_second_run_that_is_wrong_is_reported_as_wrong(runs):
    runs([_row(LEARNING, tool_calls=9)],
         [_row(LEARNING, tool_calls=2, success=False)])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "0 learned" in section and "1 wrong" in section


def test_a_second_run_that_failed_verification_is_wrong(runs):
    runs([_row(LEARNING, tool_calls=9)],
         [_row(LEARNING, tool_calls=2, verify_ok=False)])
    assert "1 wrong" in _learning_section(eval_loop.build_report(outcomes=[]))


def test_an_ordinary_task_is_never_called_learned(runs):
    """Fewer tool calls on a task with no learning setup is run-to-run
    spread; the report must not dress it up as learning."""
    runs([_row(ORDINARY, tool_calls=9)], [_row(ORDINARY, tool_calls=4)])
    report = eval_loop.build_report(outcomes=[])
    assert "## Learning" not in report


def test_the_counts_are_the_counts_of_the_lines_shown(runs):
    runs([_row(LEARNING, tool_calls=9), _row(ORDINARY, tool_calls=9)],
         [_row(LEARNING, tool_calls=4), _row(ORDINARY, tool_calls=4)])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "1 learned" in section
    assert ORDINARY not in section


def test_a_record_with_a_missing_field_is_an_error_line_not_a_zero(runs):
    broken = _row(LEARNING)
    del broken["tool_calls"]
    runs([_row(LEARNING)], [broken])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "1 error" in section and "0 learned" in section
    assert f"{LEARNING}: error (" in section


def test_one_run_alone_gives_no_learning_section(tmp_path, monkeypatch):
    d = tmp_path / "benchmark_runs"
    d.mkdir()
    monkeypatch.setattr(bm, "_DEFAULT_RUNS_DIR", d)
    (d / "run_1.jsonl").write_text(json.dumps(_row(LEARNING)) + "\n",
                                   encoding="utf-8")
    assert "## Learning" not in eval_loop.build_report(outcomes=[])


def test_the_learning_tasks_are_the_pre_registered_four():
    ids = bm.learning_task_ids()
    assert len(ids) == 4 and all(i.startswith("l5_") for i in ids)
    assert LEARNING in ids


@pytest.mark.parametrize("task", sorted(bm.learning_task_ids()))
def test_every_learning_task_has_a_setup_that_runs_and_a_check_that_starts_red(
        task, tmp_path):
    """A task whose setup or acceptance script is missing can be scored but
    never run, and an acceptance script that is green on the untouched
    fixture would pass a run that did nothing. Every learning task answers
    in a file, so every one must start red."""
    t = next(x for x in bm._load_one(bm._LEARNING_TASKS_PATH) if x.id == task)
    base = Path(bm.__file__).parent / "pack" / "benchmark"
    setup, verify = base / "setup" / t.setup, base / "accept" / t.verify
    assert setup.is_file() and verify.is_file()
    env = child_env(tmp_path)
    made = subprocess.run([sys.executable, str(setup), str(tmp_path)],
                          cwd=tmp_path, env=env, capture_output=True,
                          text=True, timeout=60)
    assert made.returncode == 0, made.stderr
    again = subprocess.run([sys.executable, str(setup), str(tmp_path)],
                           cwd=tmp_path, env=env, capture_output=True,
                           text=True, timeout=60)
    assert again.returncode == 1, "a second setup must refuse to overwrite"
    untouched = subprocess.run([sys.executable, str(verify), str(tmp_path)],
                               cwd=tmp_path, env=env, capture_output=True,
                               text=True, timeout=60)
    assert untouched.returncode == 1, (
        "acceptance is green on the untouched fixture: "
        + untouched.stdout + untouched.stderr)


def _write_runs(d: Path, files: list[list[dict]]) -> None:
    now = time.time()
    for i, rows in enumerate(files):
        p = d / f"run_{i + 1}.jsonl"
        p.write_text("\n".join(json.dumps(r) for r in rows) + "\n",
                     encoding="utf-8")
        age = 100 * (len(files) - i)
        os.utime(p, (now - age, now - age))


@pytest.mark.parametrize("newest", [
    [_row(LEARNING, tool_calls=3, learning_state="")],
    [_row(ORDINARY, tool_calls=3, learning_state="")],
])
def test_a_later_ordinary_run_does_not_hide_the_learning_pair(
        tmp_path, monkeypatch, newest):
    """The latest two run files are not necessarily the two arms: an
    ordinary run made after them must neither replace the pair nor turn
    its lines into errors."""
    d = tmp_path / "benchmark_runs"
    d.mkdir()
    monkeypatch.setattr(bm, "_DEFAULT_RUNS_DIR", d)
    _write_runs(d, [[_row(LEARNING, tool_calls=9)],
                    [_row(LEARNING, tool_calls=4)],
                    newest])
    section = _learning_section(eval_loop.build_report(outcomes=[]))
    assert "1 learned" in section and "0 error" in section
    assert "run_1.jsonl" in section and "run_2.jsonl" in section


def test_two_runs_under_different_learning_homes_are_not_a_pair(
        tmp_path, monkeypatch):
    d = tmp_path / "benchmark_runs"
    d.mkdir()
    monkeypatch.setattr(bm, "_DEFAULT_RUNS_DIR", d)
    _write_runs(d, [[_row(LEARNING, tool_calls=9, learning_state="s1")],
                    [_row(LEARNING, tool_calls=4, learning_state="s2")]])
    assert "## Learning" not in eval_loop.build_report(outcomes=[])
