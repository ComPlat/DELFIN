"""A report falls back to the record the store already keeps.

Supervision measured on 2026-09-19: ``delfin-agent report`` on a session
that had run 13 minutes said ``duration 1m 4s``, ``model -`` and
``tokens in 0 out 0``. The reason is structural, not a fluke: the report
reads its numbers from ``agent_metrics.jsonl``, and only the dashboard
ever writes rows there (``record_turn`` has no other caller) — a solo
session contributes none, so model, tokens and the timestamps come back
empty. The stored session record (``session_store.save_session``) has
carried model, token_usage, cost_usd, created_at and updated_at all
along; the report simply never looked at it.

These tests pin the fallback: metrics rows keep precedence when they
exist, the record fills the gaps when they do not, and a record without
the fields (the 52 pre-2026-09-19 ones) degrades to "unknown" instead of
breaking anything.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_report as sr

SESSION = "solo-session-no-metrics"


def patch_sources(monkeypatch, *, trace=None, turns=None, record=None):
    """Point session_report's sources at fabricated data.

    ``record`` is what ``session_store.load_session`` returns — the
    fallback under test. Trace and metrics default to empty, which is
    what a solo session leaves behind.
    """
    import delfin.agent.tool_trace as tt
    import delfin.agent.agent_metrics as am
    import delfin.agent.session_store as ss

    monkeypatch.setattr(tt, "read", lambda s, **k: trace if s == SESSION else [])
    monkeypatch.setattr(am, "read_turns", lambda **k: turns or [])
    monkeypatch.setattr(ss, "load_session", lambda s: record if s == SESSION else None)


def solo_record(**extra):
    """A stored solo-session record in the shape save_session writes."""
    base = {
        "session_id": SESSION,
        "model": "kit.glm-5.3",
        "provider": "kit",
        "token_usage": {"input": 458689, "output": 2204, "cached": 351104},
        "cost_usd": 0.0,
        "created_at": 1000.0,
        "updated_at": 1780.0,   # 13 minutes later
    }
    base.update(extra)
    return base


def test_model_comes_from_the_record_when_no_metrics_exist(monkeypatch):
    patch_sources(monkeypatch, record=solo_record())
    report = sr.collect_session_report(SESSION)
    assert report.model == "kit.glm-5.3"


def test_tokens_come_from_the_record(monkeypatch):
    patch_sources(monkeypatch, record=solo_record())
    report = sr.collect_session_report(SESSION)
    assert report.input_tokens == 458689
    assert report.output_tokens == 2204


def test_duration_spans_the_record_when_no_trace_exists(monkeypatch):
    patch_sources(monkeypatch, record=solo_record())
    report = sr.collect_session_report(SESSION)
    assert report.started_at == pytest.approx(1000.0)
    assert report.ended_at == pytest.approx(1780.0)
    # 780 s of session must not render as "1m 4s"
    rendered = sr._fmt_duration(report.started_at, report.ended_at)
    assert "13m" in rendered


def test_a_once_saved_session_takes_its_lifetime_from_the_run_clock(monkeypatch):
    # A solo session saves once at the end: created == updated, the span
    # is zero, and only run_elapsed_s — the clock the engine persists for
    # the resume budget — says how long it ran.
    rec = solo_record(updated_at=1000.0, run_elapsed_s=780.0)
    patch_sources(monkeypatch, record=rec)
    report = sr.collect_session_report(SESSION)
    assert report.started_at == pytest.approx(1000.0)
    assert report.ended_at == pytest.approx(1780.0)
    assert "13m" in sr._fmt_duration(report.started_at, report.ended_at)


def test_a_short_trace_does_not_shrink_the_sessions_lifetime(monkeypatch):
    # Real shape (2026-09-20): a solo session ran 300 s but its trace has
    # two calls 30 ms apart — both near the END of the run. The trace
    # spans the tool calls, not the session; the record's created_at +
    # run clock say how long it lived. The answer is the union, never
    # the shorter of the two.
    trace = [
        {"ts": 1299.97, "tool": "bash", "input": "ls", "output": "", "ok": True},
        {"ts": 1300.00, "tool": "read_file", "input": "x", "output": "", "ok": True},
    ]
    rec = solo_record(created_at=1000.0, updated_at=1000.0, run_elapsed_s=300.0)
    patch_sources(monkeypatch, trace=trace, record=rec)
    report = sr.collect_session_report(SESSION)
    assert report.started_at == pytest.approx(1000.0)
    assert report.ended_at == pytest.approx(1300.0)
    assert "5m" in sr._fmt_duration(report.started_at, report.ended_at)


def test_metrics_rows_still_win_when_they_exist(monkeypatch):
    patch_sources(
        monkeypatch,
        turns=[{"session_id": SESSION, "model": "other-model", "ts": 1100.0,
                "cost_usd": 0.25, "input_tokens": 10, "output_tokens": 5}],
        record=solo_record(),
    )
    report = sr.collect_session_report(SESSION)
    assert report.model == "other-model"
    assert report.input_tokens == 10
    assert report.started_at == pytest.approx(1100.0)


def test_a_record_without_the_fields_stays_unknown(monkeypatch):
    # The 52 pre-field records: no model, no token_usage, no timestamps.
    patch_sources(monkeypatch, record={"session_id": SESSION})
    report = sr.collect_session_report(SESSION)
    assert report.model == ""
    assert report.input_tokens == 0
    assert report.output_tokens == 0
    assert report.started_at == 0.0
    assert report.ended_at == 0.0


def test_no_record_at_all_keeps_the_old_answer(monkeypatch):
    patch_sources(monkeypatch, record=None)
    report = sr.collect_session_report(SESSION)
    assert report.model == ""
    assert report.input_tokens == 0
