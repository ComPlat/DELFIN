"""DELFIN /loop: interval parsing + scheduler-backed recurring runs."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import scheduler as sch


def test_parse_interval_seconds():
    assert sch.parse_interval_seconds("5m") == 300
    assert sch.parse_interval_seconds("2h") == 7200
    assert sch.parse_interval_seconds("1d") == 86400
    assert sch.parse_interval_seconds("90s") == 90
    assert sch.parse_interval_seconds("30s") == 60      # clamped to 60 min
    assert sch.parse_interval_seconds("10") is None     # no unit
    assert sch.parse_interval_seconds("abc") is None
    assert sch.parse_interval_seconds("") is None


def test_loop_registers_interval_and_fires(tmp_path):
    s = sch.Scheduler(path=tmp_path / "cron.json")
    fired = []
    s.set_fire_callback(lambda e: fired.append(e.prompt))
    ent = s.schedule_interval(every_seconds=60, prompt="/check",
                              reason="loop: /check")
    assert ent.kind == "interval"
    # Force-due and tick → the callback runs the prompt (the loop mechanism).
    ent.next_fire_at = 0
    assert s.tick() == 1
    assert fired == ["/check"]
    # Interval entries persist (reschedule), unlike "once".
    assert any(e.id == ent.id for e in s.list_entries())


def test_loop_wired_into_dashboard():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    # routed + handled + listed in autocomplete
    assert '"/loop"' in src and '"/trace"' in src           # in _SLASH_PREFIXES
    assert 'if cmd.startswith("/loop"):' in src
    assert "parse_interval_seconds" in src
