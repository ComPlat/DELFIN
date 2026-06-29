"""Cross-session per-tool aggregate for `/agents tools`.

Turns the per-call tool traces into usage/error/latency stats so the next
optimisation is data-driven (find the error + latency hotspots), not guessed.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import tool_trace as tt


def _write(d, session, entries):
    (d / f"{session}.jsonl").write_text(
        "\n".join(json.dumps(e) for e in entries) + "\n", encoding="utf-8")


@pytest.fixture
def traces(tmp_path):
    _write(tmp_path, "s1", [
        {"tool": "read_file", "ok": True, "duration_ms": 100},
        {"tool": "read_file", "ok": True, "duration_ms": 300},
        {"tool": "read_file", "ok": False, "duration_ms": 50, "error": "no"},
        {"tool": "bash", "ok": True, "duration_ms": 40},
    ])
    _write(tmp_path, "s2", [
        {"tool": "search_docs", "ok": True, "duration_ms": 900},
        {"tool": "search_docs", "ok": True, "duration_ms": 14000},
        {"tool": "read_file", "ok": True, "duration_ms": 120},
        {"tool": "view_image", "ok": False, "duration_ms": 300, "error": "x"},
    ])
    return tmp_path


def test_aggregates_across_sessions_sorted_by_calls(traces):
    rows = tt.aggregate_tools(dir_path=traces)
    # most-used first
    assert [r["tool"] for r in rows][:2] == ["read_file", "search_docs"]
    by = {r["tool"]: r for r in rows}
    # read_file spans both sessions: 4 calls, 1 error
    assert by["read_file"]["calls"] == 4
    assert by["read_file"]["errors"] == 1
    assert by["read_file"]["error_rate"] == 0.25


def test_latency_percentiles_surface_a_slow_tool(traces):
    by = {r["tool"]: r for r in tt.aggregate_tools(dir_path=traces)}
    # the 14s search_docs call shows up at p95 / max, not hidden by the fast one
    assert by["search_docs"]["p95_ms"] == 14000
    assert by["search_docs"]["max_ms"] == 14000
    assert by["search_docs"]["p50_ms"] in (900, 14000)


def test_error_hotspot_is_visible(traces):
    by = {r["tool"]: r for r in tt.aggregate_tools(dir_path=traces)}
    assert by["view_image"]["error_rate"] == 1.0   # 100% — a hotspot


def test_format_table_has_header_and_rows(traces):
    out = tt.format_tool_stats(tt.aggregate_tools(dir_path=traces))
    assert "8 calls" in out and "2 errors" in out
    assert "err%" in out and "p95ms" in out
    assert "read_file" in out and "search_docs" in out


def test_empty_and_corrupt_are_safe(tmp_path):
    # no traces at all
    assert tt.aggregate_tools(dir_path=tmp_path / "nope") == []
    assert "no tool traces" in tt.format_tool_stats([])
    # a corrupt line is skipped, valid ones still counted
    (tmp_path / "c.jsonl").write_text(
        '{"tool": "bash", "ok": true, "duration_ms": 5}\n'
        "not json at all\n", encoding="utf-8")
    rows = tt.aggregate_tools(dir_path=tmp_path)
    assert len(rows) == 1 and rows[0]["tool"] == "bash" and rows[0]["calls"] == 1
