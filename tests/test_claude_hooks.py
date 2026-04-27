"""Tests for delfin.dashboard.claude_hooks helpers."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.dashboard import claude_hooks


# ---------------------------------------------------------------------------
# discover_hooks
# ---------------------------------------------------------------------------

def test_discover_no_paths_returns_empty(tmp_path):
    """Both user and project paths missing → no hooks."""
    user = tmp_path / "user-settings.json"
    project = tmp_path / "project-settings.json"
    assert claude_hooks.discover_hooks(
        user_path=user, project_path=project,
    ) == []


def test_discover_user_hook(tmp_path):
    user = tmp_path / "user.json"
    user.write_text(json.dumps({
        "hooks": {
            "PostToolUse": [
                {"matcher": "Edit|Write",
                 "hooks": [{"type": "command", "command": "echo touched"}]}
            ]
        }
    }))
    hooks = claude_hooks.discover_hooks(
        user_path=user, project_path=tmp_path / "no.json",
    )
    assert len(hooks) == 1
    h = hooks[0]
    assert h.event == "PostToolUse"
    assert h.matcher == "Edit|Write"
    assert h.command == "echo touched"
    assert h.hook_type == "command"
    assert h.source_path == str(user)


def test_discover_project_then_user(tmp_path):
    user = tmp_path / "user.json"
    project = tmp_path / "proj.json"
    user.write_text(json.dumps({
        "hooks": {"Stop": [{"matcher": "*",
                            "hooks": [{"type": "command", "command": "user-cmd"}]}]}
    }))
    project.write_text(json.dumps({
        "hooks": {"Stop": [{"matcher": "*",
                            "hooks": [{"type": "command", "command": "project-cmd"}]}]}
    }))
    hooks = claude_hooks.discover_hooks(user_path=user, project_path=project)
    # project comes first (priority order)
    assert hooks[0].command == "project-cmd"
    assert hooks[1].command == "user-cmd"


def test_discover_skips_non_dict_entries(tmp_path):
    user = tmp_path / "user.json"
    user.write_text(json.dumps({
        "hooks": {
            "PostToolUse": [
                "broken",
                {"matcher": "Edit",
                 "hooks": [{"type": "command", "command": "ok"}]}
            ]
        }
    }))
    hooks = claude_hooks.discover_hooks(
        user_path=user, project_path=tmp_path / "missing.json",
    )
    assert len(hooks) == 1
    assert hooks[0].command == "ok"


def test_discover_corrupt_settings_returns_empty(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{ corrupt")
    hooks = claude_hooks.discover_hooks(
        user_path=bad, project_path=tmp_path / "missing.json",
    )
    assert hooks == []


def test_discover_skips_non_list_event_value(tmp_path):
    """An event whose value isn't a list shouldn't crash discovery."""
    user = tmp_path / "user.json"
    user.write_text(json.dumps({
        "hooks": {"PostToolUse": "not a list"}
    }))
    hooks = claude_hooks.discover_hooks(
        user_path=user, project_path=tmp_path / "x.json",
    )
    assert hooks == []


# ---------------------------------------------------------------------------
# add_hook / remove_hook
# ---------------------------------------------------------------------------

def test_add_hook_creates_file(tmp_path):
    settings = tmp_path / "fresh" / "settings.json"
    added = claude_hooks.add_hook(
        "PostToolUse", "Edit|Write", "echo touched",
        settings_path=settings,
    )
    assert added is True
    data = json.loads(settings.read_text())
    assert data["hooks"]["PostToolUse"][0]["matcher"] == "Edit|Write"
    assert data["hooks"]["PostToolUse"][0]["hooks"][0]["command"] == "echo touched"


def test_add_hook_reuses_matcher_group(tmp_path):
    settings = tmp_path / "settings.json"
    claude_hooks.add_hook("PostToolUse", "Edit", "cmd-1", settings_path=settings)
    claude_hooks.add_hook("PostToolUse", "Edit", "cmd-2", settings_path=settings)
    data = json.loads(settings.read_text())
    groups = data["hooks"]["PostToolUse"]
    # Single matcher group with two commands
    assert len(groups) == 1
    cmds = [h["command"] for h in groups[0]["hooks"]]
    assert cmds == ["cmd-1", "cmd-2"]


def test_add_hook_idempotent(tmp_path):
    settings = tmp_path / "settings.json"
    a = claude_hooks.add_hook("Stop", "*", "x", settings_path=settings)
    b = claude_hooks.add_hook("Stop", "*", "x", settings_path=settings)
    assert a is True
    assert b is False  # already present


def test_add_hook_rejects_unsupported_event(tmp_path):
    settings = tmp_path / "settings.json"
    with pytest.raises(ValueError):
        claude_hooks.add_hook(
            "Bogus", "*", "x", settings_path=settings,
        )


def test_remove_hook(tmp_path):
    settings = tmp_path / "settings.json"
    claude_hooks.add_hook("PostToolUse", "Edit", "first", settings_path=settings)
    claude_hooks.add_hook("PostToolUse", "Edit", "second", settings_path=settings)
    removed = claude_hooks.remove_hook(
        "PostToolUse", "Edit", "first", settings_path=settings,
    )
    assert removed is True
    data = json.loads(settings.read_text())
    cmds = [h["command"] for h in data["hooks"]["PostToolUse"][0]["hooks"]]
    assert cmds == ["second"]


def test_remove_hook_prunes_empty_groups(tmp_path):
    settings = tmp_path / "settings.json"
    claude_hooks.add_hook("Stop", "*", "only", settings_path=settings)
    claude_hooks.remove_hook("Stop", "*", "only", settings_path=settings)
    data = json.loads(settings.read_text())
    # The Stop event (and the whole hooks block) was pruned to keep the
    # file clean.
    assert "hooks" not in data


def test_remove_hook_unknown_returns_false(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"hooks": {}}))
    assert claude_hooks.remove_hook(
        "Stop", "*", "nope", settings_path=settings,
    ) is False


# ---------------------------------------------------------------------------
# render_hooks_html
# ---------------------------------------------------------------------------

def test_render_empty_shows_helpful_hint():
    html = claude_hooks.render_hooks_html([])
    assert "No hooks configured" in html


def test_render_groups_by_event():
    hooks = [
        claude_hooks.HookEntry("PostToolUse", "Edit", "command", "a", "/x.json"),
        claude_hooks.HookEntry("PostToolUse", "Write", "command", "b", "/x.json"),
        claude_hooks.HookEntry("Stop", "*", "command", "c", "/y.json"),
    ]
    html = claude_hooks.render_hooks_html(hooks)
    assert "POSTTOOLUSE (2)" in html.upper()
    assert "STOP (1)" in html.upper()


def test_render_escapes_user_text():
    hooks = [
        claude_hooks.HookEntry(
            "Stop", "<script>", "command",
            "<bad-cmd>", "/<path>",
        ),
    ]
    html = claude_hooks.render_hooks_html(hooks)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "&lt;bad-cmd&gt;" in html


def test_render_truncates_long_command():
    long_cmd = "x" * 1_000
    hooks = [claude_hooks.HookEntry("Stop", "*", "command", long_cmd, "/x")]
    html = claude_hooks.render_hooks_html(hooks)
    # Renderer caps display at 200 chars
    assert "x" * 250 not in html


def test_project_settings_path_helper(tmp_path):
    p = claude_hooks.project_settings_path(tmp_path)
    assert p == tmp_path / ".claude" / "settings.json"
