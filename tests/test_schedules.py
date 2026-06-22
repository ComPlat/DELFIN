"""Tests for delfin.dashboard.schedules."""
from __future__ import annotations

import json
from pathlib import Path

from delfin.dashboard import schedules


# ---------------------------------------------------------------------------
# discover_schedules
# ---------------------------------------------------------------------------

def test_discover_missing_file_returns_empty(tmp_path):
    assert schedules.discover_schedules(tmp_path / "nope.json") == []


def test_discover_corrupt_returns_empty(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{not json")
    assert schedules.discover_schedules(p) == []


def test_discover_full_entry(tmp_path):
    p = tmp_path / "schedules.json"
    p.write_text(json.dumps({
        "schedules": [{
            "name": "weekly-recalc",
            "cron": "0 9 * * MON",
            "command": "/skill recalc-failed",
            "enabled": True,
            "note": "Triage failed jobs every Monday at 09:00.",
        }]
    }))
    items = schedules.discover_schedules(p)
    assert len(items) == 1
    s = items[0]
    assert s.name == "weekly-recalc"
    assert s.cron == "0 9 * * MON"
    assert s.command == "/skill recalc-failed"
    assert s.enabled is True
    assert s.note.startswith("Triage")
    assert s.source_path == str(p)


def test_discover_skips_non_dict_entries(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps({
        "schedules": [
            "not a dict",
            {"name": "ok", "cron": "* * * * *", "command": "x"},
        ]
    }))
    items = schedules.discover_schedules(p)
    assert len(items) == 1
    assert items[0].name == "ok"


def test_discover_defaults_for_missing_keys(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps({"schedules": [{}]}))
    items = schedules.discover_schedules(p)
    assert len(items) == 1
    s = items[0]
    assert s.name == ""
    assert s.cron == ""
    assert s.command == ""
    # default enabled is True so a freshly-defined entry runs
    assert s.enabled is True
    assert s.note == ""


def test_discover_disabled_schedule_kept(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps({
        "schedules": [{"name": "x", "enabled": False}],
    }))
    items = schedules.discover_schedules(p)
    assert items[0].enabled is False


# ---------------------------------------------------------------------------
# render_schedules_html
# ---------------------------------------------------------------------------

def test_render_empty_shows_helpful_hint():
    html = schedules.render_schedules_html([])
    assert "schedules.json" in html


def test_render_enabled_shows_on_badge():
    items = [schedules.Schedule(
        name="weekly", cron="0 9 * * MON", command="/skill x",
        enabled=True, note="", source_path="/tmp/x",
    )]
    html = schedules.render_schedules_html(items)
    assert "ON" in html
    assert "weekly" in html
    assert "0 9 * * MON" in html
    assert "/skill x" in html


def test_render_disabled_shows_off_badge():
    items = [schedules.Schedule(
        name="paused", cron="* * * * *", command="x",
        enabled=False, note="", source_path="",
    )]
    html = schedules.render_schedules_html(items)
    assert "OFF" in html


def test_render_includes_note_when_present():
    items = [schedules.Schedule(
        name="x", cron="* * * * *", command="x",
        enabled=True, note="My note here", source_path="",
    )]
    html = schedules.render_schedules_html(items)
    assert "My note here" in html


def test_render_escapes_user_text():
    items = [schedules.Schedule(
        name="<script>", cron="<cron>", command="<cmd>",
        enabled=True, note="<note>", source_path="",
    )]
    html = schedules.render_schedules_html(items)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "&lt;note&gt;" in html


def test_render_truncates_long_command():
    long_cmd = "x" * 1000
    items = [schedules.Schedule(
        name="x", cron="*", command=long_cmd, enabled=True,
        note="", source_path="",
    )]
    html = schedules.render_schedules_html(items)
    assert "x" * 250 not in html  # capped at 200
