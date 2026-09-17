"""Tests for delfin.agent.report_flakes — fabricated runs only."""

from __future__ import annotations

import json
import time

from delfin.agent import report_flakes


def _run(target, ts, status="failed", failures=(), edits=()):
    return {
        "target": target, "ts": ts, "status": status,
        "failures": list(failures), "_edits": [{"ts": e} for e in edits],
    }


NODE = "tests.test_x::test_flaky"


def test_collect_empty_never_raises():
    data = report_flakes.collect(runs=[])
    assert data["flakes"] == []
    assert data["n_runs"] == 0
    assert "No test flipped colour" in report_flakes.format_text(data)


def test_red_then_green_no_edit_is_a_flake():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 10, status="ok"),
    ])
    assert len(data["flakes"]) == 1
    f = data["flakes"][0]
    assert f["node_id"] == NODE
    assert f["was"] == "red" and f["then"] == "green"


def test_green_then_red_is_also_a_flake():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, status="ok"),
        _run("tests/test_x.py", now - 10, failures=[NODE]),
    ])
    assert len(data["flakes"]) == 1
    assert data["flakes"][0]["was"] == "green"
    assert data["flakes"][0]["then"] == "red"


def test_flip_with_edit_between_is_not_a_flake():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 10, status="ok",
             edits=[now - 50]),
    ])
    assert data["flakes"] == []


def test_edit_after_both_runs_leaves_pair_valid():
    # Only edits BETWEEN the two timestamps disqualify a pair.
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 10, status="ok",
             edits=[now + 500]),
    ])
    assert len(data["flakes"]) == 1


def test_red_in_both_runs_is_not_a_flake():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 10, failures=[NODE]),
    ])
    assert data["flakes"] == []


def test_different_targets_are_never_paired():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_y.py", now - 10, status="ok"),
    ])
    assert data["flakes"] == []


def test_stopped_and_timeout_runs_are_ignored():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 60, status="stopped"),
        _run("tests/test_x.py", now - 10, status="ok"),
    ])
    # stopped is not a verdict: the red->green pair across it counts.
    assert len(data["flakes"]) == 1


def test_three_run_flip_chain():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 300, status="ok"),
        _run("tests/test_x.py", now - 200, failures=[NODE]),
        _run("tests/test_x.py", now - 100, status="ok"),
    ])
    # Two adjacent pairs, each flips once.
    assert len(data["flakes"]) == 2


def test_collect_from_trace_records(tmp_path, monkeypatch):
    """The real store: a fabricated tool_traces dir, read via root=."""
    now = time.time()
    sid = "sess1"
    entries = [
        # a mutating call between the runs
        {"ts": now - 50, "tool": "mcp__kit-coding__edit_file",
         "input": json.dumps({"path": "x.py"})},
        {"ts": now - 100, "tool": "mcp__kit-coding__run_tests",
         "input": json.dumps({"target": "tests/test_x.py"}),
         "output": json.dumps({
             "status": "failed", "framework": "pytest",
             "exit_code": 1,
             "summary": {"passed": 1, "failed": 1},
             "failures": [{"node_id": NODE, "message": "boom"}],
         })},
        {"ts": now - 10, "tool": "mcp__kit-coding__run_tests",
         "input": json.dumps({"target": "tests/test_x.py"}),
         "output": json.dumps({
             "status": "ok", "framework": "pytest", "exit_code": 0,
             "summary": {"passed": 2, "failed": 0}, "failures": [],
         })},
    ]
    d = tmp_path / "traces"
    d.mkdir()
    (d / f"{sid}.jsonl").write_text(
        "".join(json.dumps(e) + "\n" for e in entries), encoding="utf-8")
    data = report_flakes.collect(root=str(d))
    # The edit between the two runs disqualifies the pair.
    assert data["n_runs"] == 2
    assert data["flakes"] == []
    assert data["n_traces"] == 1

    # Same entries without the edit -> the flip is found.
    (d / f"{sid}.jsonl").write_text(
        "".join(json.dumps(e) + "\n" for e in entries if e is not entries[0]),
        encoding="utf-8")
    data = report_flakes.collect(root=str(d))
    assert len(data["flakes"]) == 1
    assert data["flakes"][0]["node_id"] == NODE


def test_truncated_output_counted_unreadable(tmp_path):
    now = time.time()
    d = tmp_path / "traces"
    d.mkdir()
    entries = [
        {"ts": now - 10, "tool": "mcp__kit-coding__run_tests",
         "input": json.dumps({"target": "tests/"}),
         "output": "…truncated garbage, not JSON…"},
    ]
    (d / "s.jsonl").write_text(
        "".join(json.dumps(e) + "\n" for e in entries), encoding="utf-8")
    data = report_flakes.collect(root=str(d))
    assert data["n_runs"] == 0
    assert data["n_unreadable_runs"] == 1


def test_format_text_lists_flakes():
    now = time.time()
    data = report_flakes.collect(runs=[
        _run("tests/test_x.py", now - 100, failures=[NODE]),
        _run("tests/test_x.py", now - 10, status="ok"),
    ])
    text = report_flakes.format_text(data)
    assert NODE in text
    assert "1 flip(s)" in text
    assert "red" in text and "green" in text


def test_main_smoke(capsys, monkeypatch, tmp_path):
    # Empty fabricated trace dir so the CLI runs without the real store.
    d = tmp_path / "traces"
    d.mkdir()
    monkeypatch.setattr(
        "delfin.agent.tool_trace._DIR", d)
    rc = report_flakes.main()
    out = capsys.readouterr().out
    assert rc == 0
    assert "No test flipped colour" in out
