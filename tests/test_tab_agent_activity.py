"""Tests for the Agent Activity tab pure-data helpers.

ipywidgets construction (``create_tab``) is not exercised here — that path
needs a running notebook kernel.  All testable logic (filtering,
aggregation, HTML rendering) lives in module-level helpers.
"""
from __future__ import annotations

from delfin.agent.outcome_tracker import CycleOutcome
from delfin.dashboard import tab_agent_activity as tab


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_outcomes() -> list[CycleOutcome]:
    return [
        CycleOutcome(
            task="fix bug in foo", provider="claude", model="opus",
            mode="solo", verdict="PASS", cost_usd=0.12,
            duration_s=42.0, retries=0,
            denied_commands=[], error_type=None,
            task_class="coding", timestamp="2026-04-27T08:00:00",
        ),
        CycleOutcome(
            task="run NMR on Co complex", provider="claude", model="sonnet",
            mode="reviewed", verdict="FAIL", cost_usd=0.34,
            duration_s=180.0, retries=2,
            denied_commands=["rm -rf /"], error_type="permission",
            task_class="chemistry", timestamp="2026-04-27T09:30:00",
        ),
        CycleOutcome(
            task="dashboard refresh", provider="openai", model="gpt-5",
            mode="dashboard", verdict="PASS", cost_usd=0.05,
            duration_s=10.0,
            task_class="dashboard", timestamp="2026-04-27T10:00:00",
        ),
        CycleOutcome(
            task="partial work", provider="claude", model="haiku",
            mode="quick", verdict="PARTIAL", cost_usd=0.02,
            duration_s=30.0,
            task_class="coding", timestamp="2026-04-27T10:15:00",
        ),
    ]


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def test_filter_no_filters_returns_all():
    out = tab._filter_outcomes(_make_outcomes())
    assert len(out) == 4


def test_filter_by_provider():
    out = tab._filter_outcomes(_make_outcomes(), provider="claude")
    assert len(out) == 3
    assert all(o.provider == "claude" for o in out)


def test_filter_by_mode():
    out = tab._filter_outcomes(_make_outcomes(), mode="reviewed")
    assert len(out) == 1
    assert out[0].mode == "reviewed"


def test_filter_by_verdict():
    out = tab._filter_outcomes(_make_outcomes(), verdict="PASS")
    assert len(out) == 2
    assert all(o.verdict == "PASS" for o in out)


def test_filter_by_task_class():
    out = tab._filter_outcomes(_make_outcomes(), task_class="coding")
    assert len(out) == 2


def test_filter_combined():
    out = tab._filter_outcomes(
        _make_outcomes(), provider="claude", verdict="PASS",
    )
    assert len(out) == 1
    assert out[0].mode == "solo"


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def test_aggregate_empty():
    stats = tab._aggregate_stats([])
    assert stats["n"] == 0
    assert stats["passes"] == 0
    assert stats["success_rate"] == 0.0
    assert stats["total_cost"] == 0.0


def test_aggregate_typical():
    stats = tab._aggregate_stats(_make_outcomes())
    assert stats["n"] == 4
    assert stats["passes"] == 2
    assert stats["fails"] == 1
    assert stats["partials"] == 1
    # 2 passes / 4 rated outcomes = 0.5
    assert stats["success_rate"] == 0.5
    assert stats["total_cost"] == 0.12 + 0.34 + 0.05 + 0.02
    assert stats["avg_cost"] == stats["total_cost"] / 4
    # avg duration over outcomes with duration_s > 0
    assert stats["avg_duration"] == (42.0 + 180.0 + 10.0 + 30.0) / 4


def test_aggregate_unrated_excluded_from_rate():
    """Outcomes with empty verdict don't count toward success rate denominator."""
    raw = _make_outcomes() + [
        CycleOutcome(verdict="", task="unrated", provider="claude"),
    ]
    stats = tab._aggregate_stats(raw)
    assert stats["n"] == 5
    # Rate still computed over rated only (4 rated, 2 passes)
    assert stats["success_rate"] == 0.5


# ---------------------------------------------------------------------------
# Rendering — sanity checks on HTML output
# ---------------------------------------------------------------------------

def test_render_summary_empty_message():
    html = tab._render_summary(tab._aggregate_stats([]))
    assert "Keine Outcome-Einträge" in html


def test_render_summary_includes_metrics():
    stats = tab._aggregate_stats(_make_outcomes())
    html = tab._render_summary(stats)
    assert "Total runs" in html
    assert "Pass" in html
    assert "Success rate" in html
    assert "50%" in html  # 2 of 4 rated


def test_render_timeline_empty():
    html = tab._render_timeline([])
    assert "Keine Einträge" in html


def test_render_timeline_includes_data():
    html = tab._render_timeline(_make_outcomes())
    assert "PASS" in html
    assert "FAIL" in html
    assert "PARTIAL" in html
    assert "claude" in html
    assert "openai" in html


def test_render_timeline_escapes_user_text():
    """User-provided task strings must be HTML-escaped."""
    bad = [CycleOutcome(
        task="<script>alert('x')</script>", verdict="PASS",
        provider="claude", model="opus", mode="solo",
        timestamp="2026-04-27T08:00:00",
    )]
    html = tab._render_timeline(bad)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_render_timeline_handles_long_tasks():
    long_task = "x" * 500
    out = [CycleOutcome(
        task=long_task, verdict="PASS", provider="claude",
        mode="solo", timestamp="2026-04-27T08:00:00",
    )]
    html = tab._render_timeline(out)
    # Truncated to <80 chars
    assert "x" * 80 not in html


def test_render_timeline_shows_denied_commands_count():
    out = _make_outcomes()  # second one has 1 denied command
    html = tab._render_timeline(out)
    assert "⛔" in html  # warning marker present
