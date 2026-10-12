"""compare_learning pairs the rows of two run files and asks the metric.

Drives ``benchmark.compare_learning`` with two JSONL run files as
``write_run`` writes them; run_a is the first run in time, run_b the
second. Rows carry five identical replicates, the same model and the same
learning state, so what is under test is the pairing and the bookkeeping;
the decision rule is pinned in tests/test_the_learning_verdict.py and the
pairing refusals in
tests/test_a_learning_measurement_keeps_state_between_its_two_arms.py.
"""

import json

from delfin.agent.benchmark import compare_learning, read_run

N = 5


def _row(*, task_id="l5_read_the_value_from_file", success=True,
         quality_0_100=90, tool_calls=6, input_tokens=300,
         output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=None):
    return {
        "task_id": task_id,
        "model": "m",
        "learning_state": "s",
        "success": success,
        "quality_0_100": quality_0_100,
        "tool_calls": tool_calls,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "duration_s": duration_s,
        "cost_usd": cost_usd,
        "verify_ok": verify_ok,
        "n_samples": N,
        "per_run_tool_calls": [tool_calls] * N,
        "per_run_tokens_total": [input_tokens + output_tokens] * N,
    }


def _write_run_file(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text("\n".join(json.dumps(r) for r in rows) + "\n",
                 encoding="utf-8")
    return p


def _report(run_a_rows, run_b_rows, tmp_path):
    return compare_learning(
        _write_run_file(tmp_path, "run_a.jsonl", run_a_rows),
        _write_run_file(tmp_path, "run_b.jsonl", run_b_rows),
    )


def _verdicts(report):
    return {p["task_id"]: p["verdict"] for p in report["per_task"]}


def test_faster_but_wrong_is_wrong(tmp_path):
    a = _row(tool_calls=8, input_tokens=500, output_tokens=300)
    b = _row(success=False, tool_calls=3, input_tokens=90, output_tokens=40)
    rep = _report([a], [b], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "wrong"
    assert rep["summary"]["n_wrong"] == 1


def test_identical_runs_are_no_learning(tmp_path):
    rep = _report([_row()], [_row()], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "no_learning"


def test_missing_field_is_error(tmp_path):
    b = _row()
    b.pop("tool_calls")
    rep = _report([_row()], [b], tmp_path)
    entry = rep["per_task"][0]
    assert entry["verdict"] == "error"
    assert "tool_calls" in entry.get("error", "")
    assert rep["summary"]["n_error"] == 1


def test_faster_correct_is_learned(tmp_path):
    a = _row(tool_calls=6, input_tokens=300, output_tokens=150)
    b = _row(tool_calls=4, input_tokens=180, output_tokens=90)
    rep = _report([a], [b], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "learned"
    assert rep["summary"]["n_learned"] == 1


def test_each_task_gets_its_own_verdict(tmp_path):
    a1, b1 = _row(task_id="t1", tool_calls=6), _row(task_id="t1", tool_calls=4)
    a2 = _row(task_id="t2", tool_calls=6)
    b2 = _row(task_id="t2", success=False, tool_calls=2)
    v = _verdicts(_report([a1, a2], [b1, b2], tmp_path))
    assert v == {"t1": "learned", "t2": "wrong"}


def test_the_last_row_of_a_task_represents_its_run(tmp_path):
    a = [_row(tool_calls=2), _row(tool_calls=6)]
    b = [_row(tool_calls=4)]
    assert _verdicts(_report(a, b, tmp_path))[
        "l5_read_the_value_from_file"] == "learned"


def test_a_task_in_only_one_run_is_not_judged(tmp_path):
    rep = _report([_row(task_id="only_in_a")],
                  [_row(task_id="only_in_b", tool_calls=4)], tmp_path)
    assert rep["per_task"] == []
    assert rep["summary"]["only_in_a"] == ["only_in_a"]
    assert rep["summary"]["only_in_b"] == ["only_in_b"]
    assert rep["summary"]["n_learned"] == 0
    assert rep["summary"]["n_wrong"] == 0


def test_an_unknown_verdict_is_counted_as_an_error(tmp_path, monkeypatch):
    from delfin.agent.pack.benchmark import learning_metric as lm

    monkeypatch.setattr(
        lm, "learning_verdict",
        lambda first, second: {"verdict": "mystery", "evidence": {}})
    rep = _report([_row()], [_row(tool_calls=5)], tmp_path)
    assert rep["summary"]["n_error"] == 1
    assert rep["per_task"][0]["verdict"] == "mystery"


def test_the_run_file_round_trips(tmp_path):
    rows = [_row()]
    assert read_run(_write_run_file(tmp_path, "run_a.jsonl", rows)) == rows
