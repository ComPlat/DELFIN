"""Focused unit tests for delfin.agent.session_report.

Fabricates one session across all four real sources (tool_trace, change_journal,
security_events, agent_metrics) by monkeypatching the functions session_report
imports, then checks every SessionReport field. Also checks the never-raise
contract for missing/misbehaving sources.
"""

from __future__ import annotations

import dataclasses

import pytest

from delfin.agent import session_report as sr


SESSION = "test-session-abc"


def patch_sources(monkeypatch, *, trace=None, changes=None, events=None, turns=None):
    """Point session_report's source helpers at fabricated data."""
    import delfin.agent.tool_trace as tt
    import delfin.agent.change_journal as cj
    import delfin.agent.security_events as se
    import delfin.agent.agent_metrics as am

    monkeypatch.setattr(tt, "read", lambda s, **k: trace if s == SESSION else [])
    monkeypatch.setattr(cj, "list_changes", lambda s, **k: changes if s == SESSION else [])
    monkeypatch.setattr(se, "recent", lambda limit=20: events or [])
    monkeypatch.setattr(am, "read_turns", lambda **k: turns or [])


def sample_trace():
    return [
        {"ts": 100.0, "tool": "bash", "input": "ls -la\ncd /tmp", "output": "ok", "ok": True},
        {"ts": 105.0, "tool": "bash", "input": "pytest -q", "output": "3 passed in 1s", "ok": True},
        {"ts": 110.0, "tool": "run_tests", "input": '{"target": "tests/test_x.py", "pytest_args": ["-q"]}',
         "output": "5 passed, 1 failed in 2s", "ok": True},
        {"ts": 120.0, "tool": "read_file", "input": "foo.py", "output": "...", "ok": True},
        {"ts": 130.0, "tool": "read_file", "input": "bar.py", "output": "err", "ok": False, "error": "boom"},
    ]


def sample_changes():
    return [
        {"seq": 1, "ts": "2026-01-01T00:00:00", "tool": "write_file", "path": "a.py",
         "created": True, "extra": {}},
        {"seq": 2, "ts": "2026-01-01T00:00:01", "tool": "edit_file", "path": "a.py",
         "created": False, "extra": {}},
        {"seq": 3, "ts": "2026-01-01T00:00:02", "tool": "bash", "path": "b.py",
         "created": False, "extra": {"deleted": True}},
        # malformed entries must be skipped, not crash
        {"seq": 4, "no_path": True},
        "garbage",
    ]


def sample_events():
    ev1 = se_event("deny_pattern", "tool X blocked")
    ev2 = se_event("notice", "just a notice", blocked=False)
    return [ev1, ev2]


def se_event(kind, detail, blocked=True):
    from delfin.agent.security_events import SecurityEvent

    return SecurityEvent(seq=1, kind=kind, tool="bash", detail=detail, blocked=blocked)


def sample_turns():
    return [
        {"session_id": SESSION, "model": "kit.glm-5.3", "ts": 100.0,
         "cost_usd": 1.5, "input_tokens": 100, "output_tokens": 50},
        {"session_id": SESSION, "model": "", "ts": 120.0,
         "cost_usd": 0.5, "input_tokens": 30, "output_tokens": 20},
        {"session_id": "other-session", "model": "zzz", "ts": 999.0,
         "cost_usd": 99.0, "input_tokens": 9999, "output_tokens": 9999},
    ]


def test_full_report(monkeypatch):
    patch_sources(
        monkeypatch,
        trace=sample_trace(),
        changes=sample_changes(),
        events=sample_events(),
        turns=sample_turns(),
    )
    r = sr.collect_session_report(SESSION)

    assert r.session_id == SESSION
    assert r.model == "kit.glm-5.3"
    assert r.started_at == 100.0
    assert r.ended_at == 130.0
    # numbers only from THIS session's turns
    assert r.cost_usd == pytest.approx(2.0)
    assert r.input_tokens == 130
    assert r.output_tokens == 70

    tools = {t["name"]: t for t in r.tool_calls}
    assert tools["bash"] == {"name": "bash", "count": 2, "ok": 2, "failed": 0}
    assert tools["read_file"] == {"name": "read_file", "count": 2, "ok": 1, "failed": 1}

    # one entry per unique path, latest change wins, order oldest-first
    assert r.files_changed == [
        {"path": "a.py", "change": "modified"},   # created then modified
        {"path": "b.py", "change": "deleted"},
    ]

    # commands_run: first line of every bash input
    assert r.commands_run == ["ls -la", "pytest -q"]

    # tests_run: pytest-flavoured tools only, counts parsed from output
    assert r.tests_run == [
        {"target": "tests/test_x.py", "status": "failed", "passed": 5, "failed": 1},
    ]

    # denials: blocked events only
    assert r.denials == [{"kind": "deny_pattern", "detail": "tool X blocked"}]


def test_missing_sources_yield_empty(monkeypatch):
    patch_sources(monkeypatch)  # everything returns []/None
    r = sr.collect_session_report("nonexistent-session")
    assert r.session_id == "nonexistent-session"
    assert r.model == ""
    assert r.started_at == 0.0 and r.ended_at == 0.0
    assert r.tool_calls == [] and r.files_changed == []
    assert r.commands_run == [] and r.tests_run == []
    assert r.denials == []
    assert r.cost_usd == 0.0 and r.input_tokens == 0 and r.output_tokens == 0


def test_misbehaving_source_never_raises(monkeypatch):
    import delfin.agent.tool_trace as tt

    def boom(*a, **k):
        raise RuntimeError("source exploded")

    monkeypatch.setattr(tt, "read", boom)
    # ensure the other imports still succeed with empty data
    patch_sources(monkeypatch, changes=[], events=[], turns=[])
    r = sr.collect_session_report(SESSION)  # must not raise
    assert r.tool_calls == []
    assert isinstance(r, sr.SessionReport)


def test_timestamps_fall_back_to_turns(monkeypatch):
    # no tool trace at all -> started/ended from turn rows
    patch_sources(monkeypatch, trace=[], turns=sample_turns())
    r = sr.collect_session_report(SESSION)
    assert r.started_at == 100.0
    assert r.ended_at == 120.0


def test_session_report_contract_fields():
    """The shared contract: exact field names, order, and defaults."""
    names = [f.name for f in dataclasses.fields(sr.SessionReport)]
    assert names == [
        "session_id", "model", "started_at", "ended_at",
        "tool_calls", "files_changed", "commands_run", "tests_run",
        "denials", "cost_usd", "input_tokens", "output_tokens",
    ]
