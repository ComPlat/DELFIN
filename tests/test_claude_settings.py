"""Tests for delfin.dashboard.claude_settings (read-only inspector)."""
from __future__ import annotations

import json
from pathlib import Path

from delfin.dashboard import claude_settings


# ---------------------------------------------------------------------------
# discover_settings
# ---------------------------------------------------------------------------

def test_discover_no_files(tmp_path):
    views = claude_settings.discover_settings(
        user_path=tmp_path / "u.json",
        project_path=tmp_path / "p.json",
    )
    assert views == []


def test_discover_user_view(tmp_path):
    user = tmp_path / "user.json"
    user.write_text(json.dumps({
        "model": "claude-opus-4-7",
        "outputStyle": "minimal",
        "env": {"DEBUG": "1"},
        "permissions": {"allow": ["Bash(git *)"], "deny": ["Bash(rm *)"]},
    }))
    views = claude_settings.discover_settings(
        user_path=user, project_path=tmp_path / "missing.json",
    )
    assert len(views) == 1
    v = views[0]
    assert v.model == "claude-opus-4-7"
    assert v.output_style == "minimal"
    assert v.env == {"DEBUG": "1"}
    assert v.allow == ["Bash(git *)"]
    assert v.deny == ["Bash(rm *)"]


def test_discover_project_first(tmp_path):
    user = tmp_path / "u.json"
    project = tmp_path / "p.json"
    user.write_text(json.dumps({"model": "user-model"}))
    project.write_text(json.dumps({"model": "proj-model"}))
    views = claude_settings.discover_settings(user_path=user, project_path=project)
    assert views[0].model == "proj-model"
    assert views[1].model == "user-model"


def test_discover_skips_missing_files(tmp_path):
    user = tmp_path / "u.json"
    user.write_text(json.dumps({"model": "x"}))
    views = claude_settings.discover_settings(
        user_path=user, project_path=tmp_path / "nope.json",
    )
    assert len(views) == 1
    assert views[0].model == "x"


def test_discover_corrupt_file_skipped(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    views = claude_settings.discover_settings(
        user_path=bad, project_path=tmp_path / "missing.json",
    )
    assert views == []


def test_discover_records_other_keys(tmp_path):
    user = tmp_path / "u.json"
    user.write_text(json.dumps({
        "model": "x", "extraField": 1, "anotherOne": True,
    }))
    views = claude_settings.discover_settings(
        user_path=user, project_path=tmp_path / "x.json",
    )
    # other_keys lists keys we don't render in detail
    assert "extraField" in views[0].other_keys
    assert "anotherOne" in views[0].other_keys


def test_discover_handles_missing_permissions_gracefully(tmp_path):
    user = tmp_path / "u.json"
    user.write_text(json.dumps({"model": "x"}))
    v = claude_settings.discover_settings(
        user_path=user, project_path=tmp_path / "nope.json",
    )[0]
    assert v.allow == []
    assert v.deny == []
    assert v.env == {}


# ---------------------------------------------------------------------------
# render_settings_html
# ---------------------------------------------------------------------------

def test_render_empty_shows_helpful_hint():
    html = claude_settings.render_settings_html([])
    assert "No <code>settings.json</code> found" in html


def test_render_includes_model_and_permissions(tmp_path):
    user = tmp_path / "u.json"
    user.write_text(json.dumps({
        "model": "claude-opus-4-7",
        "permissions": {"allow": ["Bash(pytest*)"]}
    }))
    views = claude_settings.discover_settings(
        user_path=user, project_path=tmp_path / "x.json",
    )
    html = claude_settings.render_settings_html(views)
    assert "claude-opus-4-7" in html
    assert "Bash(pytest*)" in html


def test_render_escapes_user_text():
    view = claude_settings.SettingsView(
        source_path="/<path>",
        model="<bad>",
        env={"<k>": "<v>"},
        allow=["<allow>"],
        deny=["<deny>"],
        output_style="<style>",
        other_keys=["<extra>"],
    )
    html = claude_settings.render_settings_html([view])
    assert "<bad>" not in html
    assert "&lt;bad&gt;" in html
    assert "&lt;allow&gt;" in html
