"""Tests for delfin.agent.context_tracker."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.context_tracker import (
    ContextUsageTracker,
    _detect_references,
)


class TestDetectReferences:
    def test_playbook_hit(self):
        hits = _detect_references(
            "Following the playbook step 1, I will grep first.",
            ["playbook", "repo_map"],
        )
        assert hits["playbook"] is True
        assert hits["repo_map"] is False

    def test_repo_map_hit(self):
        hits = _detect_references(
            "Looking at delfin/agent/engine.py and tests/test_agent.py",
            ["repo_map", "briefing"],
        )
        assert hits["repo_map"] is True
        assert hits["briefing"] is False

    def test_no_hits(self):
        hits = _detect_references(
            "I will fix the bug by changing the function.",
            ["playbook", "briefing", "profile"],
        )
        assert all(not v for v in hits.values())

    def test_empty_response(self):
        hits = _detect_references("", ["playbook"])
        assert hits["playbook"] is False


class TestContextUsageTracker:
    def test_record_and_get_rates(self, tmp_path: Path):
        tracker = ContextUsageTracker(path=tmp_path / "usage.jsonl")

        # Record 10 interactions where playbook is always referenced
        for _ in range(10):
            tracker.record_usage(
                sections_injected=["playbook", "repo_map"],
                response_text="Following the playbook step 1...",
                role_id="solo_agent",
                provider="claude",
            )

        rates = tracker.get_hit_rates()
        assert rates["playbook"] == 1.0
        assert rates["repo_map"] == 0.0

    def test_should_skip_low_hit_rate(self, tmp_path: Path):
        tracker = ContextUsageTracker(path=tmp_path / "usage.jsonl")

        # Record 10 interactions where briefing is never referenced
        for _ in range(10):
            tracker.record_usage(
                sections_injected=["briefing"],
                response_text="I will fix the bug directly.",
                role_id="solo_agent",
            )

        assert tracker.should_skip("briefing", role_id="solo_agent") is True

    def test_should_not_skip_insufficient_data(self, tmp_path: Path):
        tracker = ContextUsageTracker(path=tmp_path / "usage.jsonl")

        # Only 2 records — not enough to make a decision
        for _ in range(2):
            tracker.record_usage(
                sections_injected=["briefing"],
                response_text="Plain response.",
            )

        assert tracker.should_skip("briefing") is False

    def test_should_not_skip_good_hit_rate(self, tmp_path: Path):
        tracker = ContextUsageTracker(path=tmp_path / "usage.jsonl")

        for _ in range(10):
            tracker.record_usage(
                sections_injected=["playbook"],
                response_text="Following the playbook approach...",
            )

        assert tracker.should_skip("playbook") is False

    def test_filter_by_role(self, tmp_path: Path):
        tracker = ContextUsageTracker(path=tmp_path / "usage.jsonl")

        for _ in range(10):
            tracker.record_usage(
                sections_injected=["repo_map"],
                response_text="Looking at delfin/agent/engine.py",
                role_id="builder_agent",
            )
            tracker.record_usage(
                sections_injected=["repo_map"],
                response_text="The plan looks correct.",
                role_id="critic_agent",
            )

        rates_builder = tracker.get_hit_rates(role_id="builder_agent")
        rates_critic = tracker.get_hit_rates(role_id="critic_agent")
        assert rates_builder["repo_map"] > rates_critic["repo_map"]
