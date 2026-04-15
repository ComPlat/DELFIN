"""Tests for delfin.agent.briefing — pre-task briefing generator."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.briefing import classify_task, generate_briefing, _analyse_outcomes
from delfin.agent.outcome_tracker import CycleOutcome


def _make_outcome(
    verdict: str = "PASS",
    task_class: str = "coding",
    provider: str = "claude",
    cost: float = 0.5,
    error_type: str | None = None,
    denied: list[str] | None = None,
    retries: int = 0,
    mode: str = "solo",
    task: str = "fix the bug",
) -> CycleOutcome:
    return CycleOutcome(
        task=task,
        provider=provider,
        model="opus",
        mode=mode,
        verdict=verdict,
        cost_usd=cost,
        retries=retries,
        denied_commands=denied or [],
        error_type=error_type,
        task_class=task_class,
        timestamp="2026-04-15T10:00:00",
    )


class TestClassifyTask:
    def test_chemistry_task(self):
        assert classify_task("optimize the DFT basis set for iron") == "chemistry"

    def test_dashboard_task(self):
        assert classify_task("fix the dashboard widget layout") == "dashboard"

    def test_coding_task(self):
        assert classify_task("refactor the pipeline module") == "coding"

    def test_default_to_coding(self):
        assert classify_task("hello world") == "coding"

    def test_empty_string(self):
        assert classify_task("") == "coding"


class TestAnalyseOutcomes:
    def test_empty_outcomes(self):
        assert _analyse_outcomes([], "coding") == {}

    def test_basic_stats(self):
        outcomes = [
            _make_outcome("PASS"),
            _make_outcome("PASS"),
            _make_outcome("FAIL", error_type="timeout"),
        ]
        stats = _analyse_outcomes(outcomes, "coding")
        assert stats["total"] == 3
        assert abs(stats["pass_rate"] - 2 / 3) < 0.01
        assert stats["fails"] == 1
        assert "timeout" in stats["error_counts"]

    def test_denied_commands_tracked(self):
        outcomes = [
            _make_outcome(denied=["git push"]),
            _make_outcome(denied=["git push", "rm -rf"]),
            _make_outcome(),
        ]
        stats = _analyse_outcomes(outcomes, "coding")
        assert "git" in stats["denied_counts"]

    def test_mode_stats(self):
        outcomes = [
            _make_outcome(mode="solo"),
            _make_outcome(mode="solo"),
            _make_outcome(mode="reviewed", verdict="FAIL"),
        ]
        stats = _analyse_outcomes(outcomes, "coding")
        assert stats["mode_stats"]["solo"]["pass"] == 2
        assert stats["mode_stats"]["reviewed"]["pass"] == 0


class TestGenerateBriefing:
    def test_empty_with_few_outcomes(self, tmp_path: Path):
        path = tmp_path / "outcomes.jsonl"
        path.write_text("")
        assert generate_briefing("claude", "fix bug", outcome_path=path) == ""

    def test_generates_briefing(self, tmp_path: Path):
        path = tmp_path / "outcomes.jsonl"
        outcomes = [
            _make_outcome("PASS", cost=0.3),
            _make_outcome("PASS", cost=0.4),
            _make_outcome("FAIL", error_type="timeout", cost=1.2),
            _make_outcome("PASS", cost=0.25),
        ]
        lines = [json.dumps(o.__dict__) for o in outcomes]
        path.write_text("\n".join(lines) + "\n")

        briefing = generate_briefing("claude", "fix the parser bug", outcome_path=path)
        assert "coding" in briefing
        assert "pass rate" in briefing.lower() or "%" in briefing

    def test_low_success_rate_warning(self, tmp_path: Path):
        path = tmp_path / "outcomes.jsonl"
        outcomes = [
            _make_outcome("FAIL", error_type="crash"),
            _make_outcome("FAIL", error_type="timeout"),
            _make_outcome("PASS"),
            _make_outcome("FAIL"),
            _make_outcome("FAIL"),
        ]
        lines = [json.dumps(o.__dict__) for o in outcomes]
        path.write_text("\n".join(lines) + "\n")

        briefing = generate_briefing("claude", "fix the bug", outcome_path=path)
        assert "WARNING" in briefing or "low" in briefing.lower()

    def test_briefing_length_capped(self, tmp_path: Path):
        path = tmp_path / "outcomes.jsonl"
        outcomes = [
            _make_outcome(
                "FAIL",
                error_type="crash",
                denied=["git push", "rm -rf"],
                task="a very long task " * 20,
            )
            for _ in range(20)
        ]
        lines = [json.dumps(o.__dict__) for o in outcomes]
        path.write_text("\n".join(lines) + "\n")

        briefing = generate_briefing("claude", "fix it", outcome_path=path)
        assert len(briefing) <= 1600

    def test_filters_by_provider(self, tmp_path: Path):
        path = tmp_path / "outcomes.jsonl"
        outcomes = [
            _make_outcome("PASS", provider="claude"),
            _make_outcome("PASS", provider="claude"),
            _make_outcome("PASS", provider="claude"),
            _make_outcome("FAIL", provider="openai"),
            _make_outcome("FAIL", provider="openai"),
            _make_outcome("FAIL", provider="openai"),
        ]
        lines = [json.dumps(o.__dict__) for o in outcomes]
        path.write_text("\n".join(lines) + "\n")

        briefing = generate_briefing("claude", "fix bug", outcome_path=path)
        assert "100%" in briefing or "pass rate" in briefing.lower()
