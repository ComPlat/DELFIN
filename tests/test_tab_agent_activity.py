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


# ---------------------------------------------------------------------------
# A6 — honest cost_usd_delta accounting
# ---------------------------------------------------------------------------

def test_outcome_cost_delta_uses_explicit_field():
    """When cost_usd_delta is set, _outcome_cost_delta returns it."""
    o = CycleOutcome(cost_usd=12.5, cost_usd_delta=0.75)
    assert tab._outcome_cost_delta(o) == 0.75


def test_outcome_cost_delta_zero_for_legacy():
    """Legacy entries (no cost_usd_delta) get 0 from the per-outcome helper."""
    o = CycleOutcome(cost_usd=12.5, cost_usd_delta=0.0)
    assert tab._outcome_cost_delta(o) == 0.0


def test_outcomes_total_delta_sums_explicit_deltas():
    """All-new entries: sum of deltas, NOT sum of cumulative."""
    outs = [
        CycleOutcome(provider="claude", mode="solo", cost_usd=0.10, cost_usd_delta=0.10,
                     timestamp="2026-04-28T10:00:00"),
        CycleOutcome(provider="claude", mode="solo", cost_usd=0.30, cost_usd_delta=0.20,
                     timestamp="2026-04-28T10:05:00"),
        CycleOutcome(provider="claude", mode="solo", cost_usd=0.45, cost_usd_delta=0.15,
                     timestamp="2026-04-28T10:10:00"),
    ]
    # Honest spend = 0.10 + 0.20 + 0.15 = 0.45 (matches the final cumulative)
    assert abs(tab._outcomes_total_delta(outs) - 0.45) < 1e-9


def test_outcomes_total_delta_legacy_session_aware_diff():
    """Legacy entries (delta=0) within a 4 h same-(provider,mode) window get
    diffed against their predecessor — that recovers the per-cycle cost."""
    outs = [
        CycleOutcome(provider="claude", mode="solo", cost_usd=5.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T10:00:00"),
        CycleOutcome(provider="claude", mode="solo", cost_usd=15.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T10:05:00"),
        CycleOutcome(provider="claude", mode="solo", cost_usd=22.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T10:10:00"),
    ]
    # Honest spend should be 5 (first) + 10 (15-5) + 7 (22-15) = 22
    # NOT 5 + 15 + 22 = 42 (the buggy sum-of-cumulative)
    assert abs(tab._outcomes_total_delta(outs) - 22.0) < 1e-9


def test_outcomes_total_delta_far_apart_treated_as_separate_sessions():
    """Outcomes > 4 h apart can't be diffed — each counts as its own first-of-session."""
    outs = [
        CycleOutcome(provider="claude", mode="solo", cost_usd=5.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T08:00:00"),
        CycleOutcome(provider="claude", mode="solo", cost_usd=12.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T20:00:00"),  # 12 h later
    ]
    # Two separate sessions → 5 + 12 = 17 (each cumulative treated as its own delta)
    assert abs(tab._outcomes_total_delta(outs) - 17.0) < 1e-9


def test_outcomes_total_delta_mixed_explicit_and_legacy():
    """A real-world history: some legacy + some new entries."""
    outs = [
        # Legacy first session
        CycleOutcome(provider="claude", mode="solo", cost_usd=5.0, cost_usd_delta=0.0,
                     timestamp="2026-04-28T08:00:00"),
        # New entry continues the same session — explicit delta
        CycleOutcome(provider="claude", mode="solo", cost_usd=8.0, cost_usd_delta=3.0,
                     timestamp="2026-04-28T08:30:00"),
        # New separate session, explicit delta
        CycleOutcome(provider="claude", mode="dashboard", cost_usd=0.5, cost_usd_delta=0.5,
                     timestamp="2026-04-28T15:00:00"),
    ]
    # 5 (legacy first) + 3 (explicit) + 0.5 (explicit) = 8.5
    assert abs(tab._outcomes_total_delta(outs) - 8.5) < 1e-9


def test_outcomes_total_delta_empty_returns_zero():
    assert tab._outcomes_total_delta([]) == 0.0


def test_aggregate_stats_uses_honest_total():
    """End-to-end: _aggregate_stats reports the corrected total_cost."""
    outs = [
        CycleOutcome(provider="x", mode="solo", verdict="PASS", cost_usd=2.0,
                     cost_usd_delta=2.0, timestamp="2026-04-28T08:00:00"),
        CycleOutcome(provider="x", mode="solo", verdict="PASS", cost_usd=5.0,
                     cost_usd_delta=3.0, timestamp="2026-04-28T08:05:00"),
    ]
    stats = tab._aggregate_stats(outs)
    # Was buggy 2+5=7, now honest 2+3=5
    assert stats["total_cost"] == 5.0
    assert stats["avg_cost"] == 2.5


# ---------------------------------------------------------------------------
# A1 — _format_live_state pure helper
# ---------------------------------------------------------------------------

def test_live_state_missing_engine_shows_placeholder():
    html = tab._format_live_state({})
    assert "Not started" in html
    assert "No active engine" in html


def test_live_state_idle_shows_idle_badge():
    html = tab._format_live_state({"engine_status": "idle"})
    assert ">Idle<" in html
    assert "(engine present but idle)" in html


def test_live_state_streaming_badge():
    html = tab._format_live_state({
        "engine_status": "streaming",
        "current_mode": "dashboard",
        "current_model": "haiku",
        "current_provider": "claude",
    })
    assert ">Streaming<" in html
    assert "dashboard" in html
    assert "claude/haiku" in html


def test_live_state_session_cost_and_turns():
    html = tab._format_live_state({
        "engine_status": "idle",
        "session_cost_usd": 0.123,
        "session_turns": 5,
    })
    assert "$0.123" in html or "$0.12¢" in html or "session" in html
    assert "5 turn" in html


def test_live_state_jobs_summary():
    html = tab._format_live_state({
        "engine_status": "idle",
        "active_jobs": {"RUNNING": 2, "PENDING": 1},
    })
    assert "RUNNING:" in html
    assert ">2<" in html
    assert "PENDING:" in html
    assert ">1<" in html


def test_live_state_jobs_zeroes_filtered():
    """Job statuses with count 0 must NOT show up (clutter)."""
    html = tab._format_live_state({
        "engine_status": "idle",
        "active_jobs": {"RUNNING": 1, "FAILED": 0},
    })
    assert "RUNNING:" in html
    assert "FAILED" not in html


def test_live_state_perm_profile_visible():
    html = tab._format_live_state({
        "engine_status": "idle",
        "perm_profile": "repo_free",
    })
    assert "repo_free" in html
    assert "perms" in html


def test_live_state_session_id_truncated():
    html = tab._format_live_state({
        "engine_status": "idle",
        "active_session_id": "abcdef1234567890" * 4,
    })
    assert "abcdef12" in html
    assert "abcdef1234567890abcdef1234567890" not in html  # truncated


# ---------------------------------------------------------------------------
# A3 — cost insights (period buckets + per-mode)
# ---------------------------------------------------------------------------

from datetime import datetime as _dt


def _outcome_at(ts: str, *, mode: str = "solo", delta: float = 0.5,
                cum: float | None = None, provider: str = "x") -> CycleOutcome:
    return CycleOutcome(
        provider=provider, mode=mode, verdict="PASS",
        cost_usd=cum if cum is not None else delta,
        cost_usd_delta=delta,
        timestamp=ts,
    )


def test_aggregate_costs_by_period_buckets_correctly():
    """Calendar boundaries: today / this-week / this-month / all."""
    now = _dt(2026, 4, 28, 12, 0, 0)  # Tuesday in week 18
    outcomes = [
        _outcome_at("2026-04-28T08:00:00", delta=0.10),  # today
        _outcome_at("2026-04-27T20:00:00", delta=0.20),  # this week
        _outcome_at("2026-04-15T08:00:00", delta=0.30),  # this month, NOT this week
        _outcome_at("2026-03-15T08:00:00", delta=1.00),  # all-time only
    ]
    p = tab._aggregate_costs_by_period(outcomes, now=now)
    assert p["today"] == 0.10
    assert abs(p["week"] - 0.30) < 1e-9   # 0.10 + 0.20
    assert abs(p["month"] - 0.60) < 1e-9  # 0.10 + 0.20 + 0.30
    assert abs(p["all"] - 1.60) < 1e-9    # all four


def test_aggregate_costs_by_period_empty():
    assert tab._aggregate_costs_by_period([]) == {
        "today": 0.0, "week": 0.0, "month": 0.0, "all": 0.0,
    }


def test_aggregate_costs_by_period_unparseable_ts_excluded():
    """An entry with a broken timestamp must NOT corrupt buckets."""
    now = _dt(2026, 4, 28, 12, 0, 0)
    outcomes = [
        _outcome_at("2026-04-28T08:00:00", delta=0.10),
        CycleOutcome(provider="x", mode="solo", verdict="PASS",
                     cost_usd=999.0, cost_usd_delta=999.0,
                     timestamp="not a real timestamp"),
    ]
    p = tab._aggregate_costs_by_period(outcomes, now=now)
    assert p["today"] == 0.10  # broken-ts entry excluded
    assert p["all"] == 0.10 + 999.0   # but still counted in all-time


def test_aggregate_costs_by_period_at_day_boundary():
    """An entry one second before midnight counts to that day, not the next."""
    now = _dt(2026, 4, 28, 0, 30, 0)  # right after midnight
    outcomes = [
        _outcome_at("2026-04-27T23:59:30", delta=0.50),  # yesterday
        _outcome_at("2026-04-28T00:00:30", delta=0.10),  # today
    ]
    p = tab._aggregate_costs_by_period(outcomes, now=now)
    assert p["today"] == 0.10
    # Both fall in the same iso week (Mon-Sun) → both in week
    assert abs(p["week"] - 0.60) < 1e-9


def test_aggregate_costs_by_mode_groups_correctly():
    outcomes = [
        _outcome_at("2026-04-28T08:00:00", mode="solo",      delta=0.10),
        _outcome_at("2026-04-28T08:01:00", mode="solo",      delta=0.20),
        _outcome_at("2026-04-28T08:02:00", mode="dashboard", delta=0.05),
        _outcome_at("2026-04-28T08:03:00", mode="quick",     delta=0.15),
    ]
    by_mode = tab._aggregate_costs_by_mode(outcomes)
    assert abs(by_mode["solo"] - 0.30) < 1e-9
    assert abs(by_mode["dashboard"] - 0.05) < 1e-9
    assert abs(by_mode["quick"] - 0.15) < 1e-9


def test_aggregate_costs_by_mode_empty_returns_empty_dict():
    assert tab._aggregate_costs_by_mode([]) == {}


def test_render_cost_insights_includes_all_period_cards():
    out = tab._render_cost_insights(
        {"today": 0.10, "week": 0.50, "month": 1.20, "all": 5.0},
        {"solo": 3.0, "dashboard": 2.0},
    )
    assert "Today" in out
    assert "This week" in out
    assert "This month" in out
    assert "All time" in out
    # Mode bars
    assert "Cost by mode" in out
    assert "solo" in out
    assert "dashboard" in out


def test_render_cost_insights_handles_empty_modes():
    out = tab._render_cost_insights(
        {"today": 0.0, "week": 0.0, "month": 0.0, "all": 0.0},
        {},
    )
    assert "Today" in out
    # No "Cost by mode" section when by_mode is empty
    assert "Cost by mode" not in out
