"""L5 learning hook: the runner-side caller of learning_verdict.

Drives ``delfin.agent.benchmark.compare_learning`` through its public path
with two real JSONL run files (one asdict(BenchmarkResult) row per line, as
``write_run`` writes them). run_a is the FIRST run in time, run_b the SECOND.
The hook compares the two runs OF THE SAME TASK and reports whether the
second run learned (fewer tool_calls/tokens/duration/cost and no correctness
regression), via learning_metric.learning_verdict (PREREG_W15.md).

The three cases the task names, driven through the hook with real files:
  * a second run that is faster but wrong fails (verdict wrong),
  * identical runs score no_learning, never learned,
  * a record with a missing field is an error, never a silent zero.
Plus the positive test (faster + correct -> learned) and the pairing rule
(a task seen only in one run is not judgeable).

This test is the RED control: ``compare_learning`` does not exist yet; the
hook is delivered as a .gate patch (benchmark.py is read-only for L5, placed
next to ab_compare per QS s29). The operator builds the patch and this test
then turns green.
"""

import json

from delfin.agent.benchmark import read_run, compare_learning


def _row(*, task_id="l5_read_the_value_from_file", success=True,
         quality_0_100=90, tool_calls=6, input_tokens=300,
         output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=None):
    """One valid asdict(BenchmarkResult) row, as write_run writes it."""
    return {
        "task_id": task_id,
        "success": success,
        "quality_0_100": quality_0_100,
        "tool_calls": tool_calls,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "duration_s": duration_s,
        "cost_usd": cost_usd,
        "verify_ok": verify_ok,
    }


def _write_run_file(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    return p


def _report(run_a_rows, run_b_rows, tmp_path):
    return compare_learning(
        _write_run_file(tmp_path, "run_a.jsonl", run_a_rows),
        _write_run_file(tmp_path, "run_b.jsonl", run_b_rows),
    )


def _verdicts(report):
    return {p["task_id"]: p["verdict"] for p in report["per_task"]}


# --- the three cases the task names, via the hook with real files -----------


def test_faster_but_wrong_fails(tmp_path):
    """Second run saves on every proxy but is wrong -> per-task wrong."""
    a = _row(success=True, tool_calls=8, input_tokens=500,
             output_tokens=300, duration_s=40.0, cost_usd=0.3)
    b = _row(success=False, tool_calls=3, input_tokens=90,
             output_tokens=40, duration_s=8.0, cost_usd=0.05)
    rep = _report([a], [b], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "wrong"


def test_identical_runs_are_no_learning(tmp_path):
    """Same task, identical economy and correctness -> no_learning."""
    row = _row()
    rep = _report([row], [row], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "no_learning"


def test_missing_field_is_error(tmp_path):
    """A record without tool_calls cannot be compared -> error, not zero."""
    a = _row()
    b = _row(tool_calls=99)
    b.pop("tool_calls")
    rep = _report([a], [b], tmp_path)
    entry = rep["per_task"][0]
    assert entry["verdict"] == "error"
    assert "tool_calls" in entry.get("error", "")


# --- positive case -----------------------------------------------------------


def test_faster_correct_is_learned(tmp_path):
    """Second run strictly better and still correct -> learned."""
    a = _row(tool_calls=6, input_tokens=300, output_tokens=150,
             duration_s=20.0, cost_usd=0.12)
    b = _row(tool_calls=4, input_tokens=180, output_tokens=90,
             duration_s=12.0, cost_usd=0.07)
    rep = _report([a], [b], tmp_path)
    assert _verdicts(rep)["l5_read_the_value_from_file"] == "learned"


def test_faster_wrong_over_faster_correct_pair(tmp_path):
    """Two tasks in one pair: one learns, one goes wrong -> both reported,
    the run must NOT be judged as a single pass/fail."""
    a1 = _row(task_id="t1", tool_calls=6, duration_s=20.0)
    b1 = _row(task_id="t1", tool_calls=4, duration_s=12.0)
    a2 = _row(task_id="t2", tool_calls=6, duration_s=20.0)
    b2 = _row(task_id="t2", success=False, tool_calls=2, duration_s=5.0)
    rep = _report([a1, a2], [b1, b2], tmp_path)
    v = _verdicts(rep)
    assert v["t1"] == "learned"
    assert v["t2"] == "wrong"


# --- pairing rule: a task must have a counterpart in both runs ---------------


def test_task_only_in_one_run_not_judged(tmp_path):
    """A task present in only one of the two runs has no first or second
    record and cannot be judged as learning: it is excluded from per_task
    and reported in the summary as only_in_a / only_in_b."""
    a = _row(task_id="only_in_a")
    b = _row(task_id="only_in_b", tool_calls=4, duration_s=12.0)
    rep = _report([a], [b], tmp_path)
    # no task has a counterpart -> per_task empty, no learning verdict at all
    assert rep["per_task"] == []
    assert rep["summary"]["only_in_a"] == ["only_in_a"]
    assert rep["summary"]["only_in_b"] == ["only_in_b"]
    assert rep["summary"]["n_learned"] == 0
    assert rep["summary"]["n_wrong"] == 0


# --- operator F2: an unknown verdict is an error count, never a KeyError -----


def test_unknown_verdict_counted_as_error_not_keyerror(tmp_path, monkeypatch):
    """Operator F2: compare_learning must survive a verdict it has no bucket
    for. If learning_verdict ever returned a value outside the five, the old
    stats[f"n_{verdict}"] += 1 raised KeyError; the fixed hook counts the
    unknown verdict against n_error and does not crash. The decision itself
    is still surfaced on the per_task row."""
    from delfin.agent.pack.benchmark import learning_metric as lm

    monkeypatch.setattr(
        lm, "learning_verdict",
        lambda first, second: {"verdict": "mystery", "evidence": {}})
    a = _row()
    b = _row(tool_calls=5, duration_s=18.0)
    rep = _report([a], [b], tmp_path)   # must not raise KeyError
    assert rep["summary"]["n_error"] == 1
    assert rep["per_task"][0]["verdict"] == "mystery"


# --- the hook reads what write_run wrote -------------------------------------


def test_hook_reads_the_record_shape(tmp_path):
    """compare_learning reads real JSONL via read_run; the file written here
    must round-trip so the hook sees exactly the rows write_run would."""
    rows = [_row()]
    p = _write_run_file(tmp_path, "run_a.jsonl", rows)
    assert read_run(p) == rows
