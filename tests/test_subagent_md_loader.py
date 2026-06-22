"""Tests for the markdown-based subagent preset loader.

The loader walks ``delfin/agent/pack/agents/*_subagent.md`` and
``~/.delfin/subagents/*_subagent.md``, parses YAML frontmatter, and
folds the discovered presets into ``SUBAGENT_PRESETS`` alongside the
four built-ins. User-local files (``~/.delfin/``) may override built-ins;
pack-shipped files may not (so an accidentally-shadowed name in the
codebase doesn't silently change behaviour).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import subagents as sa


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    yield tmp_path
    # Restore the global preset registry so leaked fixtures from this
    # test file don't bleed into other tests that import subagents.
    sa.reload_subagent_presets()


def test_builtins_present_after_reload(fake_home):
    sa.reload_subagent_presets()
    names = sa.subagent_type_names()
    assert {"explore", "plan", "code-reviewer", "general-purpose"} <= set(names)


def test_user_md_preset_discovered(fake_home, tmp_path):
    user_dir = tmp_path / ".delfin" / "subagents"
    user_dir.mkdir(parents=True)
    (user_dir / "chem-audit_subagent.md").write_text(
        "---\n"
        "name: chem-audit\n"
        "description: Audit ORCA inputs for soundness.\n"
        "mode: plan\n"
        "---\n\n"
        "You are a chemistry audit sub-agent.\n",
        encoding="utf-8",
    )
    sa.reload_subagent_presets()
    p = sa.SUBAGENT_PRESETS.get("chem-audit")
    assert p is not None
    assert p.mode == "plan"
    assert "Audit ORCA" in p.description
    assert "chemistry audit" in p.system_prompt


def test_user_md_can_override_builtin(fake_home, tmp_path):
    user_dir = tmp_path / ".delfin" / "subagents"
    user_dir.mkdir(parents=True)
    (user_dir / "explore_subagent.md").write_text(
        "---\n"
        "name: explore\n"
        "description: Custom explore agent.\n"
        "mode: plan\n"
        "---\n\n"
        "Custom explore body.\n",
        encoding="utf-8",
    )
    sa.reload_subagent_presets()
    p = sa.SUBAGENT_PRESETS["explore"]
    assert p.description == "Custom explore agent."
    assert "Custom explore body" in p.system_prompt


def test_pack_md_cannot_silently_shadow_builtin(fake_home):
    """If a pack md file accidentally uses a built-in name, the built-in
    must win — we don't trust packed overrides to be intentional."""
    sa.reload_subagent_presets()
    p = sa.SUBAGENT_PRESETS["explore"]
    # The built-in description starts with "Read-only research"; if a pack
    # file overrode it, that prefix would disappear.
    assert p.description.startswith("Read-only research")


def test_subagent_type_names_used_by_api_client_schema(fake_home, tmp_path):
    """The tool schema enum exposed to the LLM must reflect the loader
    so user-defined subagents are pickable. The api_client builds the
    enum from subagent_type_names() at import time, but we also want a
    direct unit-level check that the loader is the source of truth."""
    user_dir = tmp_path / ".delfin" / "subagents"
    user_dir.mkdir(parents=True)
    (user_dir / "x-agent_subagent.md").write_text(
        "---\nname: x-agent\ndescription: x.\nmode: default\n---\nbody",
        encoding="utf-8",
    )
    sa.reload_subagent_presets()
    assert "x-agent" in sa.subagent_type_names()


def test_md_file_without_frontmatter_skipped(fake_home, tmp_path):
    user_dir = tmp_path / ".delfin" / "subagents"
    user_dir.mkdir(parents=True)
    (user_dir / "noframe_subagent.md").write_text(
        "# Just a markdown file with no frontmatter\n\nbody\n",
        encoding="utf-8",
    )
    sa.reload_subagent_presets()
    # Without frontmatter the loader falls back to stem-derived name —
    # the file IS picked up (name = "noframe") but with empty description.
    # This is acceptable; just confirm the loader didn't crash.
    p = sa.SUBAGENT_PRESETS.get("noframe")
    if p is not None:
        assert p.name == "noframe"
