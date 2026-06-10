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

def test_default_settings_disabled():
    from delfin.user_settings import DEFAULT_SETTINGS
    cfg = DEFAULT_SETTINGS["agent"]["eval_loop"]
    assert cfg["enabled"] is False
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
