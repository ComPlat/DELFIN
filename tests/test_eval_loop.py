"""Tests for the eval loop (outcome mining → draft tasks → report).

Contract: LLM-free, opt-in (default off), recurring-failure mining with
a threshold, drafts land in the shared review inbox and are valid YAML,
and the report reflects outcomes + suite integrity.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from delfin.agent import eval_loop as ev


def _o(verdict="FAIL", task_class="chemistry", error_type="timeout",
       mode="solo", task="optimize the complex"):
    return {"verdict": verdict, "task_class": task_class,
            "error_type": error_type, "mode": mode, "task": task}


# ---------------------------------------------------------------------------
# Opt-in
# ---------------------------------------------------------------------------

def test_default_settings_enabled():
    from delfin.user_settings import DEFAULT_SETTINGS
    cfg = DEFAULT_SETTINGS["agent"]["eval_loop"]
    assert cfg["enabled"] is True
    assert cfg["threshold"] >= 2


# ---------------------------------------------------------------------------
# Pattern mining
# ---------------------------------------------------------------------------

def test_recurring_failures_are_detected():
    outcomes = [_o() for _ in range(4)] + [_o(verdict="PASS")] * 10
    patterns = ev.analyze_outcomes(outcomes, threshold=3)
    assert len(patterns) == 1
    p = patterns[0]
    assert p.count == 4
    assert p.task_class == "chemistry"
    assert p.error_type == "timeout"
    assert p.examples and "optimize" in p.examples[0]


def test_below_threshold_is_not_a_pattern():
    outcomes = [_o(), _o()]                      # only 2 fails
    assert ev.analyze_outcomes(outcomes, threshold=3) == []


def test_distinct_fingerprints_bucket_separately():
    outcomes = ([_o(error_type="timeout")] * 3
                + [_o(error_type="crash")] * 3)
    patterns = ev.analyze_outcomes(outcomes, threshold=3)
    assert len(patterns) == 2
    assert {p.error_type for p in patterns} == {"timeout", "crash"}


def test_passes_never_count():
    outcomes = [_o(verdict="PASS")] * 10
    assert ev.analyze_outcomes(outcomes, threshold=1) == []


# ---------------------------------------------------------------------------
# Draft scaffolding
# ---------------------------------------------------------------------------

def test_drafts_are_valid_yaml_with_todo(tmp_path):
    patterns = ev.analyze_outcomes([_o()] * 3, threshold=3)
    drafts = ev.write_pattern_drafts(patterns, dest_dir=tmp_path)
    assert len(drafts) == 1
    parsed = yaml.safe_load(drafts[0].read_text())
    task = parsed["tasks"][0]
    assert task["id"].startswith("recur_")
    assert task["mode"] == "solo"
    assert any("TODO" in s["pattern"] for s in task["expected_signals"])
    assert "REVIEW before committing" in drafts[0].read_text()


# ---------------------------------------------------------------------------
# Report + full pass
# ---------------------------------------------------------------------------

def test_report_contains_outcomes_patterns_and_integrity():
    outcomes = [_o()] * 3 + [_o(verdict="PASS")] * 7
    report = ev.build_report(outcomes=outcomes, threshold=3)
    assert "failure rate: 3/10" in report
    assert "chemistry|timeout|solo" in report
    assert "Benchmark suite integrity" in report
    assert "OK — 0 errors" in report            # committed suite is clean


def _patch_telemetry_dirs(monkeypatch, tmp_path):
    """Point every telemetry reader the report consumes at tmp dirs."""
    from delfin.agent import benchmark as bm
    from delfin.agent import tool_trace as tt
    from delfin.agent import turn_metrics as tm
    monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turns")
    monkeypatch.setattr(bm, "_DEFAULT_RUNS_DIR", tmp_path / "runs")
    return tmp_path / "traces", tmp_path / "turns", tmp_path / "runs"


def test_report_contains_tool_and_turn_health(tmp_path, monkeypatch):
    import json
    import time
    traces, turns, _runs = _patch_telemetry_dirs(monkeypatch, tmp_path)
    traces.mkdir(parents=True)
    turns.mkdir(parents=True)
    now = time.time()
    (traces / "s.jsonl").write_text(
        json.dumps({"ts": now, "tool": "bash", "ok": False,
                    "duration_ms": 10, "error": "exit 1"}) + "\n"
        + json.dumps({"ts": now, "tool": "read_file", "ok": True,
                      "duration_ms": 5}) + "\n",
        encoding="utf-8")
    (turns / "s.jsonl").write_text(
        json.dumps({"ts": now, "model": "m", "total_ms": 90000,
                    "ttft_ms": 89000, "output_chars": 5, "tool_calls": 0,
                    "stopped": True}) + "\n",
        encoding="utf-8")
    report = ev.build_report(outcomes=[])
    assert "## Tool health" in report
    assert "**bash**: 100% errors (1/1 calls)" in report
    assert 'e.g. "exit 1"' in report
    assert "## Turn health" in report
    assert "turns: 1, stalls: 1, stopped: 1" in report
    assert "p90 89.0s" in report
    # fewer than two benchmark runs → no drift section
    assert "## Benchmark drift" not in report


def test_report_health_sections_present_without_telemetry(tmp_path, monkeypatch):
    _patch_telemetry_dirs(monkeypatch, tmp_path)        # dirs never created
    report = ev.build_report(outcomes=[])
    assert "## Tool health" in report
    assert "no tool calls recorded in window" in report
    assert "## Turn health" in report
    assert "no turns recorded in window" in report


def test_report_benchmark_drift_with_two_runs(tmp_path, monkeypatch):
    import os
    from delfin.agent import benchmark as bm
    _traces, _turns, runs = _patch_telemetry_dirs(monkeypatch, tmp_path)

    def _results(quality, success):
        return [bm.BenchmarkResult(task_id=f"t{i}", task_class="x",
                                   model="m", success=success,
                                   quality_0_100=quality, cost_usd=0.01)
                for i in range(3)]

    pa = bm.write_run(_results(50, False), model="m", runs_dir=runs,
                      run_id="a_base")
    pb = bm.write_run(_results(80, True), model="m", runs_dir=runs,
                      run_id="b_cand")
    os.utime(pa, (1000, 1000))                          # force mtime order
    os.utime(pb, (2000, 2000))
    report = ev.build_report(outcomes=[])
    assert "## Benchmark drift" in report
    assert "`a_base.jsonl` → `b_cand.jsonl`" in report
    assert "verdict: BETTER" in report
    assert "score delta: +30.0" in report


def test_run_eval_writes_report(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.outcome_tracker.load_outcomes",
        lambda max_entries=200: [_o()] * 3,
    )
    path = ev.run_eval(settings={"agent": {"eval_loop": {"enabled": True}}},
                       reports_dir=tmp_path / "reports",
                       drafts_dir=tmp_path / "drafts")
    assert path.is_file()
    body = path.read_text()
    assert "eval report" in body
    # the recurring pattern produced a draft, referenced in the report
    assert list((tmp_path / "drafts").glob("recur_*.yaml"))


def test_maybe_run_scheduled_disabled_returns_none(tmp_path):
    assert ev.maybe_run_scheduled(
        {"agent": {"eval_loop": {"enabled": False}}},
        reports_dir=tmp_path) is None


def test_maybe_run_scheduled_runs_once_per_interval(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        ev, "run_eval",
        lambda **kw: calls.append(1) or (tmp_path / "r.md"))
    cfg = {"agent": {"eval_loop": {"enabled": True}}}
    first = ev.maybe_run_scheduled(cfg, reports_dir=tmp_path)
    assert first is not None
    assert calls == [1]
    # Second call inside the interval is a no-op.
    assert ev.maybe_run_scheduled(cfg, reports_dir=tmp_path) is None
    assert calls == [1]
    # Expired stamp -> runs again.
    (tmp_path / ".last_run").write_text("0", encoding="utf-8")
    assert ev.maybe_run_scheduled(cfg, reports_dir=tmp_path) is not None
    assert calls == [1, 1]


def test_maybe_run_scheduled_survives_run_error(tmp_path, monkeypatch):
    def _boom(**kw):
        raise RuntimeError("mining failed")
    monkeypatch.setattr(ev, "run_eval", _boom)
    cfg = {"agent": {"eval_loop": {"enabled": True}}}
    assert ev.maybe_run_scheduled(cfg, reports_dir=tmp_path) is None
    # No stamp written on failure -> next call tries again.
    assert not (tmp_path / ".last_run").exists()
