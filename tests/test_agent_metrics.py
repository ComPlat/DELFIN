"""Tests for the agent metrics recorder + aggregator.

Metrics drive the "got better?" verdict — we need them to round-trip
cleanly, aggregate correctly per model, and surface improvement /
regression deltas between windows.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from delfin.agent import agent_metrics as am


@pytest.fixture
def log_path(tmp_path, monkeypatch):
    p = tmp_path / "metrics.jsonl"
    monkeypatch.setattr(am, "_LOG_PATH", p)
    return p


# ---------------------------------------------------------------------------
# Recording
# ---------------------------------------------------------------------------


def test_record_turn_creates_log_file(log_path):
    am.record_turn(am.TurnMetrics(
        model="opus", duration_s=2.5, output_tokens=120,
        cost_usd=0.012, tool_calls=2,
    ), path=log_path)
    assert log_path.exists()
    text = log_path.read_text(encoding="utf-8")
    assert '"model": "opus"' in text
    assert '"cost_usd": 0.012' in text


def test_record_turn_appends_multiple_lines(log_path):
    for i in range(3):
        am.record_turn(am.TurnMetrics(model="opus", output_tokens=i),
                       path=log_path)
    lines = log_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3


def test_record_turn_never_raises_on_bad_path(tmp_path):
    """Best-effort: a non-writable path must NOT crash the agent path."""
    bad = tmp_path / "nope" / "deeper" / "x.jsonl"
    # mkdir succeeds, just confirm no exception
    am.record_turn(am.TurnMetrics(model="x"), path=bad)


def test_record_turn_auto_timestamps_when_ts_omitted(log_path):
    am.record_turn(am.TurnMetrics(model="x"), path=log_path)
    rec = am.read_turns(path=log_path)[0]
    assert rec["ts"] > 0
    assert abs(rec["ts"] - time.time()) < 5  # within 5 s


# ---------------------------------------------------------------------------
# Per-model aggregation
# ---------------------------------------------------------------------------


def test_aggregate_by_model_computes_correct_averages(log_path):
    am.record_turn(am.TurnMetrics(model="A", duration_s=2.0, cost_usd=0.01,
                                  output_tokens=100, tool_calls=2,
                                  tool_errors=0), path=log_path)
    am.record_turn(am.TurnMetrics(model="A", duration_s=4.0, cost_usd=0.03,
                                  output_tokens=300, tool_calls=4,
                                  tool_errors=1), path=log_path)
    am.record_turn(am.TurnMetrics(model="B", duration_s=10.0, cost_usd=0.10,
                                  output_tokens=500, tool_calls=2,
                                  tool_errors=2), path=log_path)
    out = am.aggregate_by_model(path=log_path)
    a = out["A"]
    assert a["n_turns"] == 2
    assert a["avg_duration_s"] == 3.0
    assert a["avg_cost_usd"] == pytest.approx(0.02)
    assert a["total_cost_usd"] == pytest.approx(0.04)
    assert a["avg_tokens_out"] == 200
    # 1 error out of 6 tool calls
    assert a["tool_error_rate"] == pytest.approx(1 / 6)
    b = out["B"]
    # 2 errors out of 2 tool calls
    assert b["tool_error_rate"] == pytest.approx(1.0)


def test_aggregate_filters_by_since(log_path):
    """``since_s`` should drop older records — useful for 'last hour'
    or 'last day' rollups."""
    am.record_turn(am.TurnMetrics(model="A", ts=time.time() - 7200,
                                  cost_usd=0.10), path=log_path)
    am.record_turn(am.TurnMetrics(model="A", ts=time.time() - 60,
                                  cost_usd=0.01), path=log_path)
    out = am.aggregate_by_model(since_s=3600, path=log_path)
    assert out["A"]["n_turns"] == 1
    assert out["A"]["total_cost_usd"] == pytest.approx(0.01)


def test_silent_exit_rate_tracked(log_path):
    am.record_turn(am.TurnMetrics(model="A", silent_exit=True),
                   path=log_path)
    am.record_turn(am.TurnMetrics(model="A", silent_exit=False),
                   path=log_path)
    am.record_turn(am.TurnMetrics(model="A", silent_exit=True),
                   path=log_path)
    out = am.aggregate_by_model(path=log_path)
    assert out["A"]["silent_exit_rate"] == pytest.approx(2 / 3)


def test_cooperative_stop_rate_tracked(log_path):
    am.record_turn(am.TurnMetrics(model="A", cooperative_stop=True),
                   path=log_path)
    am.record_turn(am.TurnMetrics(model="A", cooperative_stop=False),
                   path=log_path)
    out = am.aggregate_by_model(path=log_path)
    assert out["A"]["cooperative_stop_rate"] == 0.5


# ---------------------------------------------------------------------------
# Before/after window comparison — "did we improve?"
# ---------------------------------------------------------------------------


def test_compare_windows_signals_improvement_on_cost_drop(log_path):
    """Older window expensive (avg $0.10), newer window cheap ($0.01).
    Comparison must report ``improved=True`` for avg_cost_usd."""
    now = time.time()
    # 7 expensive turns in the "older" window (3 days ago)
    for _ in range(7):
        am.record_turn(am.TurnMetrics(
            model="A", ts=now - 3 * 86400, cost_usd=0.10,
            duration_s=10.0, output_tokens=500,
        ), path=log_path)
    # 7 cheap turns in the "newer" window (10 minutes ago)
    for _ in range(7):
        am.record_turn(am.TurnMetrics(
            model="A", ts=now - 600, cost_usd=0.01,
            duration_s=2.0, output_tokens=100,
        ), path=log_path)
    cmp = am.compare_windows(model="A", path=log_path)
    assert cmp is not None
    assert cmp["avg_cost_usd"]["improved"] is True
    assert cmp["avg_cost_usd"]["new"] < cmp["avg_cost_usd"]["old"]
    assert cmp["avg_duration_s"]["improved"] is True


def test_compare_windows_returns_none_on_thin_data(log_path):
    """Fewer than 5 turns in either window — refuse to claim a trend."""
    now = time.time()
    am.record_turn(am.TurnMetrics(model="A", ts=now - 600, cost_usd=0.01),
                   path=log_path)
    assert am.compare_windows(model="A", path=log_path) is None


def test_compare_windows_signals_regression_on_cost_rise(log_path):
    now = time.time()
    for _ in range(7):
        am.record_turn(am.TurnMetrics(
            model="A", ts=now - 3 * 86400, cost_usd=0.01,
        ), path=log_path)
    for _ in range(7):
        am.record_turn(am.TurnMetrics(
            model="A", ts=now - 600, cost_usd=0.10,
        ), path=log_path)
    cmp = am.compare_windows(model="A", path=log_path)
    assert cmp is not None
    assert cmp["avg_cost_usd"]["improved"] is False


# ---------------------------------------------------------------------------
# read_turns
# ---------------------------------------------------------------------------


def test_read_turns_handles_missing_file(tmp_path, monkeypatch):
    monkeypatch.setattr(am, "_LOG_PATH", tmp_path / "nonexistent.jsonl")
    assert am.read_turns() == []


def test_read_turns_skips_corrupt_lines(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        '{"model": "A", "cost_usd": 0.01}\n'
        "this is not json\n"
        '{"model": "B", "cost_usd": 0.02}\n',
        encoding="utf-8",
    )
    out = am.read_turns(path=log_path)
    assert len(out) == 2
    assert out[0]["model"] == "A"
    assert out[1]["model"] == "B"


def test_read_turns_last_n_filter(log_path):
    for i in range(10):
        am.record_turn(am.TurnMetrics(model="A", output_tokens=i),
                       path=log_path)
    last3 = am.read_turns(last_n=3, path=log_path)
    assert [r["output_tokens"] for r in last3] == [7, 8, 9]
