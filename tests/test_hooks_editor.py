"""Tests for hooks_editor — the backend of the /hooks slash commands."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import hooks_editor as he


@pytest.fixture
def settings_file(tmp_path):
    return tmp_path / "settings.json"


def test_add_hook_writes_settings_json(settings_file):
    rec = he.add_hook(
        "PreToolUse", "Edit|Write", "ruff check ${file}",
        settings_path=settings_file,
    )
    assert rec["event"] == "PreToolUse"
    assert settings_file.exists()
    data = json.loads(settings_file.read_text(encoding="utf-8"))
    bucket = data["hooks"]["PreToolUse"]
    assert len(bucket) == 1
    assert bucket[0]["matcher"] == "Edit|Write"
    assert bucket[0]["hooks"][0]["command"] == "ruff check ${file}"


def test_add_hook_rejects_unknown_event(settings_file):
    with pytest.raises(ValueError):
        he.add_hook("PostFlux", "", "do-it", settings_path=settings_file)


def test_add_hook_rejects_empty_command(settings_file):
    with pytest.raises(ValueError):
        he.add_hook("PreToolUse", "", "  ", settings_path=settings_file)


def test_add_then_remove_round_trip(settings_file):
    he.add_hook("PreToolUse", "", "a", settings_path=settings_file)
    he.add_hook("PreToolUse", "", "b", settings_path=settings_file)
    he.add_hook("PostToolUse", "", "c", settings_path=settings_file)
    removed = he.remove_hook("PreToolUse", 0, settings_path=settings_file)
    assert removed is not None
    data = json.loads(settings_file.read_text(encoding="utf-8"))
    assert len(data["hooks"]["PreToolUse"]) == 1
    assert data["hooks"]["PreToolUse"][0]["hooks"][0]["command"] == "b"


def test_remove_out_of_range_returns_none(settings_file):
    he.add_hook("PreToolUse", "", "a", settings_path=settings_file)
    assert he.remove_hook("PreToolUse", 99, settings_path=settings_file) is None


def test_remove_empties_bucket_then_drops_event_key(settings_file):
    he.add_hook("Stop", "", "only", settings_path=settings_file)
    he.remove_hook("Stop", 0, settings_path=settings_file)
    data = json.loads(settings_file.read_text(encoding="utf-8"))
    assert "Stop" not in (data.get("hooks") or {})


def test_list_hooks_picks_up_added_hooks(tmp_path, monkeypatch):
    """list_hooks reads via hooks.load_hooks, which checks Path.home();
    monkey-patch home to the tmp dir so we only see what we added."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    user_file = tmp_path / ".delfin" / "settings.json"
    he.add_hook("PreToolUse", "Edit", "echo hi", settings_path=user_file)
    rows = he.list_hooks(workspace=None)
    assert any(
        r["event"] == "PreToolUse" and r["command"] == "echo hi"
        for r in rows
    )


def test_dry_run_returns_empty_when_no_hooks(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    rows = he.dry_run_hook("PreToolUse", tool_name="Edit", workspace=tmp_path)
    assert rows == []


def test_dry_run_executes_registered_hook(tmp_path, monkeypatch):
    """End-to-end-ish: register a no-op `true` hook + dry-run it. The
    matcher we use empty so the hook ALWAYS matches; the command is a
    POSIX-safe `true` that exits 0."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    user_file = tmp_path / ".delfin" / "settings.json"
    he.add_hook("PreToolUse", "", "true", settings_path=user_file)
    rows = he.dry_run_hook(
        "PreToolUse", tool_name="Edit", workspace=tmp_path,
    )
    assert rows
    assert rows[0]["exit_code"] == 0
    assert rows[0]["matched"] is True
