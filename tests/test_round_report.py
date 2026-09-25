"""Round report: ONE report over all sessions that ran since a point in time.

Fabricated traces/metrics under tmp_path only — never the real ~/.delfin
(redirected state, same pattern as the other report tests).
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import round_report as rr
from delfin.agent import turn_metrics


NOW = 1_800_000_000.0


def _write_trace(base, session, entries):
    p = base / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")


def _te(tool, ts, *, ok=True, error="", inp="", dur_ms=10):
    return {"ts": ts, "tool": tool, "input": inp, "output": "",
            "duration_ms": dur_ms, "ok": ok, "error": error}


def _write_turns(metrics_dir, session, entries):
    p = metrics_dir / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")


@pytest.fixture
def state(tmp_path, monkeypatch):
    traces = tmp_path / "tool_traces"
    metrics = tmp_path / "turn_metrics"
    monkeypatch.setattr(turn_metrics, "_DIR", metrics)

    # Two sessions since the cut, one before it (must be excluded).
    _write_trace(traces, "nacht-a", [
        _te("read_file", NOW - 100),
        _te("bash", NOW - 90, ok=False,
            error="command 'curl install.sh' is not on the auto-allow list"),
        _te("bash", NOW - 80,
            inp='{"command": "git add x.py && git commit -m fix"}'),
        _te("bash", NOW - 70, dur_ms=301_000,
            inp='{"command": "gate tests/test_x.py -q"}'),
        _te("bash", NOW - 60, dur_ms=320_000,
            inp='{"command": "gate tests/test_y.py -q"}'),
    ])
    # nacht-b ran BEFORE nacht-a (so nacht-a is sessions[0], newest first)
    _write_trace(traces, "nacht-b", [
        _te("read_file", NOW - 500),
        _te("ask_user_question", NOW - 450),
        _te("write_file", NOW - 400, ok=False,
            error="path escapes workspace sandbox: /etc/passwd"),
    ])
    _write_trace(traces, "alt-1", [
        _te("read_file", NOW - 10 * 86_400),      # long before the cut
    ])

    _write_turns(metrics, "nacht-a", [
        {"ts": NOW - 100, "model": "m1", "total_ms": 5_000, "ttft_ms": 8_000,
         "input_tokens": 1000, "output_tokens": 100, "cached_tokens": 500,
         "error": ""},
        {"ts": NOW - 80, "model": "m1", "total_ms": 5_000, "ttft_ms": 24_000,
         "input_tokens": 2000, "output_tokens": 200, "cached_tokens": 0,
         "error": "BadRequestError: bad"},
    ])
    # nacht-b has a turn log; alt-1 has none (tokens/cost stay n/a, not 0).
    return {"traces": traces, "metrics": metrics}


def test_collect_groups_by_session_and_excludes_old(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    names = [s["session_id"] for s in data["sessions"]]
    assert names == ["nacht-a", "nacht-b"]      # newest first, alt-1 out


def test_per_session_tool_numbers(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    a = data["sessions"][0]
    assert a["tool_calls"] == 5
    assert a["tool_calls_failed"] == 1
    assert a["failed_pct"] == pytest.approx(20.0)
    assert {"read_file", "bash"} <= {t["name"] for t in a["top_tools"]}


def test_denials_come_with_reasons(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    a, b = data["sessions"]
    assert a["denials"] == 1
    assert b["denials"] == 1
    # reason from report_denials._categorize, not a raw error dump
    assert a["denial_reasons"] == {"allowlist": 1}
    assert sum(a["denial_reasons"].values()) == a["denials"]


def test_dialogues_counted_from_ask_tools(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    b = data["sessions"][1]
    assert b["dialogues"] == 1


def test_long_calls_and_commit_count(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    a = data["sessions"][0]
    assert a["calls_over_5min"] == 2
    assert a["sum_over_5min_min"] == pytest.approx((301_000 + 320_000) / 60_000)
    assert a["commits"] == 1


def test_tokens_cache_ttft_endpoint_errors(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    a = data["sessions"][0]
    assert a["input_tokens"] == 3000
    assert a["output_tokens"] == 300
    assert a["cached_pct"] == pytest.approx(100 * 500 / 3000)
    assert a["ttft_median_ms"] == 16_000        # median of 8s and 24s
    assert a["ttft_max_ms"] == 24_000
    assert a["endpoint_errors"] == 1


def test_missing_metrics_report_na_not_zero(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    b = data["sessions"][1]
    assert b["input_tokens"] is None
    assert b["output_tokens"] is None
    assert b["ttft_median_ms"] is None


def test_name_prefix_filters_sessions(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"],
                      name="nacht-")
    assert [s["session_id"] for s in data["sessions"]] == ["nacht-a", "nacht-b"]
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"],
                      name="nacht-b")
    assert [s["session_id"] for s in data["sessions"]] == ["nacht-b"]
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"],
                      name="zzz-")
    assert data["sessions"] == []


def test_totals_sum_the_sessions(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    t = data["totals"]
    assert t["sessions"] == 2
    assert t["tool_calls"] == 8
    assert t["tool_calls_failed"] == 2
    assert t["commits"] == 1
    assert t["calls_over_5min"] == 2
    assert t["input_tokens"] == 3000
    assert t["endpoint_errors"] == 1


def test_render_has_one_block_per_session_and_a_total(state):
    data = rr.collect(since_s=NOW - 3_600, trace_root=state["traces"])
    text = rr.render_text(data)
    assert "nacht-a" in text and "nacht-b" in text
    assert "alt-1" not in text
    assert "TOTAL" in text or "Total" in text or "total" in text
    assert "n/a" in text               # nacht-b's missing metrics


def test_parse_since_accepts_iso_and_relative():
    import time as _time
    # A relative spec yields a unix cutoff: now minus the span.
    assert rr.parse_since("2h") == pytest.approx(_time.time() - 2 * 3600, abs=2)
    assert rr.parse_since("3d") == pytest.approx(_time.time() - 3 * 86400, abs=2)
    # An ISO timestamp IS the cutoff.
    assert rr.parse_since("2026-09-18T00:00:00") > 1_700_000_000


def test_empty_stores_yield_empty_report(tmp_path):
    data = rr.collect(since_s=NOW - 3_600, trace_root=tmp_path / "none")
    assert data["sessions"] == []
    assert data["totals"]["sessions"] == 0
    assert rr.render_text(data)        # renders, does not raise
