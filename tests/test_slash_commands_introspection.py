"""Tests for the three introspection slash commands: /agents /skills /context.

These commands surface DELFIN's existing subsystems (subagents.py, skills.py,
engine.py compaction state) into the chat as one-line summaries the user can
glance at without leaving the agent tab.
"""

from __future__ import annotations


def test_slash_palette_lists_the_core_commands():
    """The browse-commands palette (_SLASH_COMMANDS) must surface the real
    user commands — guards against a handler-only command being unbrowsable
    (2026-06-26: /batch was handled but missing from the palette)."""
    from delfin.dashboard.tab_agent import _SLASH_COMMANDS

    cmds = {row[1] for row in _SLASH_COMMANDS}
    expected = {
        "/loop", "/batch", "/submit", "/check", "/trace", "/review",
        "/explore", "/delegate", "/watch", "/jobs", "/usage", "/help",
    }
    missing = expected - cmds
    assert not missing, f"palette missing commands: {sorted(missing)}"
    # Every palette entry carries a non-empty description (col index 2).
    for row in _SLASH_COMMANDS:
        assert row[2], f"palette entry without description: {row}"


def test_list_subagents_returns_all_four_presets():
    from delfin.agent.subagents import list_subagents

    presets = list_subagents()
    names = {p.get("subagent_type") or p.get("name") for p in presets}
    assert {"explore", "plan", "code-reviewer", "general-purpose"} <= names, names
    # Every preset advertises a non-empty description so the slash output is useful
    for p in presets:
        assert p.get("description"), f"missing description: {p}"


def test_discover_skills_picks_up_project_skill(tmp_path):
    """When a `.delfin/skills/<name>/SKILL.md` exists in the workspace,
    /skills must include it. This is the path the slash command relies on."""
    from delfin.agent.skills import discover_skills

    skill_dir = tmp_path / ".delfin" / "skills" / "fix-bug"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        "---\n"
        "name: fix-bug\n"
        "description: Reproduce and patch a reported bug.\n"
        "---\n\n"
        "# Fix bug\n\nDo X then Y.\n",
        encoding="utf-8",
    )

    skills = discover_skills(tmp_path)
    names = {s.name for s in skills}
    assert "fix-bug" in names
    fix = next(s for s in skills if s.name == "fix-bug")
    assert "Reproduce" in (fix.description or "")


def test_context_status_reads_engine_estimate_tokens(monkeypatch):
    """The /context handler relies on engine._estimate_context_tokens(); the
    method must exist on AgentEngine and return a non-negative int for an
    empty conversation."""
    from delfin.agent.engine import AgentEngine

    # Build a minimal engine without spinning up any backend.
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None

    # The real method walks self.messages; with [] it must return 0 or close.
    n = AgentEngine._estimate_context_tokens(eng)
    assert isinstance(n, int)
    assert n >= 0
