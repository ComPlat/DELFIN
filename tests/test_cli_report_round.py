"""CLI wiring: `delfin-agent report --since <t> [--name PREFIX]` runs the
round report over all matching sessions instead of the single-session one.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import round_report as rr
from delfin.agent import tool_trace, turn_metrics


def _write_trace(base, session, entries):
    p = base / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")


NOW = 1_800_000_000.0


@pytest.fixture
def state(tmp_path, monkeypatch):
    traces = tmp_path / "tool_traces"
    metrics = tmp_path / "turn_metrics"
    monkeypatch.setattr(turn_metrics, "_DIR", metrics)
    _write_trace(traces, "nacht-a", [
        {"ts": NOW - 100, "tool": "read_file", "input": "", "output": "",
         "duration_ms": 5, "ok": True, "error": ""},
        {"ts": NOW - 90, "tool": "bash",
         "input": '{"command": "git add x && git commit -m t"}',
         "output": "", "duration_ms": 5, "ok": True, "error": ""},
    ])
    _write_trace(traces, "alt-1", [
        {"ts": NOW - 10 * 86_400, "tool": "read_file", "input": "",
         "output": "", "duration_ms": 5, "ok": True, "error": ""},
    ])
    return {"traces": traces}


def test_since_routes_to_round_report(state, capsys, monkeypatch):
    # The cutoff is wall-clock based; freeze it so NOW-100 is inside.
    monkeypatch.setattr(rr.time, "time", lambda: NOW)
    # Redirect the trace store the CLI-side collect() defaults to.
    monkeypatch.setattr(tool_trace, "_DIR", state["traces"])
    args = agent_cli.build_parser().parse_args(
        ["report", "--since", "1h", "--name", "nacht-"])
    rc = agent_cli.cmd_report(args)
    assert rc == 0
    out = capsys.readouterr().out
    assert "nacht-a" in out
    assert "alt-1" not in out
    assert "TOTAL" in out


def test_since_with_session_still_single_report(state, capsys, monkeypatch):
    # Explicit --session keeps the old single-session behaviour.
    from delfin.agent import session_report as sr
    args = agent_cli.build_parser().parse_args(
        ["report", "--session", "nacht-a"])
    monkeypatch.setattr(sr, "collect_session_report",
                        lambda sid: sr.SessionReport(session_id=sid))
    rc = agent_cli.cmd_report(args)
    assert rc == 0
    assert "nacht-a" in capsys.readouterr().out


def test_bad_since_is_rejected_not_guessed(capsys, monkeypatch):
    # An unparseable --since must fail loudly, not fall back to some
    # default window nobody asked for.
    args = agent_cli.build_parser().parse_args(
        ["report", "--since", "not-a-time"])
    rc = agent_cli.cmd_report(args)
    assert rc == 2
    err = capsys.readouterr().err
    assert "cannot parse" in err or "not-a-time" in err
