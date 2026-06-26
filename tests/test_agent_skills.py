"""Tests for delfin.agent.skills (skill discovery + expansion)."""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import skills


# ---------------------------------------------------------------------------
# parse_skill_text
# ---------------------------------------------------------------------------

def test_parse_skill_minimal():
    text = "# Tune CONTROL\n\nbody line"
    skill = skills.parse_skill_text("tune", text)
    assert skill is not None
    assert skill.name == "tune"
    assert skill.title == "Tune CONTROL"
    assert skill.description == ""
    assert skill.body == "body line"


def test_parse_skill_with_description():
    text = "# Tune CONTROL\n> One-line description.\n\nMain body here."
    skill = skills.parse_skill_text("tune", text)
    assert skill is not None
    assert skill.title == "Tune CONTROL"
    assert skill.description == "One-line description."
    assert skill.body == "Main body here."


def test_parse_skill_strips_hash_chars():
    """Multiple '#' characters in the title are stripped."""
    skill = skills.parse_skill_text("a", "## Heading\nbody")
    assert skill is not None
    assert skill.title == "Heading"


def test_parse_skill_empty_returns_none():
    assert skills.parse_skill_text("a", "") is None
    assert skills.parse_skill_text("a", "   \n\n  ") is None


def test_parse_skill_no_title_returns_none():
    """File without a leading '# Title' is rejected."""
    assert skills.parse_skill_text("a", "no heading here\njust body") is None


def test_parse_skill_handles_blank_lines_before_title():
    skill = skills.parse_skill_text("a", "\n\n# Title\n\nbody")
    assert skill is not None
    assert skill.title == "Title"
    assert skill.body == "body"


def test_parse_skill_records_source_path():
    skill = skills.parse_skill_text(
        "a", "# Title\nbody", source_path="/tmp/a.md",
    )
    assert skill is not None
    assert skill.source_path == "/tmp/a.md"


def test_skill_slash_command_property():
    skill = skills.Skill(
        name="tune-control", title="x", description="", body="b",
        source_path="",
    )
    assert skill.slash_command == "/skill tune-control"


# ---------------------------------------------------------------------------
# discover_skills
# ---------------------------------------------------------------------------

def test_discover_skills_empty_dir(tmp_path):
    assert skills.discover_skills(tmp_path) == []


def test_discover_skills_missing_dir(tmp_path):
    assert skills.discover_skills(tmp_path / "nope") == []


def test_discover_skills_reads_md_files(tmp_path):
    (tmp_path / "a.md").write_text("# A\nbody a")
    (tmp_path / "b.md").write_text("# B\n> desc\n\nbody b")
    found = skills.discover_skills(tmp_path)
    names = [s.name for s in found]
    assert names == ["a", "b"]  # sorted alphabetically


def test_discover_skills_skips_invalid_files(tmp_path):
    (tmp_path / "good.md").write_text("# Good\nbody")
    (tmp_path / "bad.md").write_text("no heading")
    found = skills.discover_skills(tmp_path)
    assert [s.name for s in found] == ["good"]


def test_discover_skills_ignores_non_md(tmp_path):
    (tmp_path / "a.md").write_text("# A\nbody")
    (tmp_path / "b.txt").write_text("# Should be ignored\nbody")
    found = skills.discover_skills(tmp_path)
    assert [s.name for s in found] == ["a"]


def test_discover_packaged_skills_present():
    """The shipped skills directory has at least one valid skill."""
    found = skills.discover_skills()
    assert len(found) >= 1, "expected at least one packaged skill"
    names = {s.name for s in found}
    assert "tune-control" in names


# ---------------------------------------------------------------------------
# find_skill
# ---------------------------------------------------------------------------

def test_find_skill_by_name(tmp_path):
    (tmp_path / "tune.md").write_text("# Tune\nbody")
    skill = skills.find_skill("tune", skills_dir=tmp_path)
    assert skill is not None
    assert skill.name == "tune"


def test_find_skill_case_insensitive(tmp_path):
    (tmp_path / "Tune.md").write_text("# Tune\nbody")
    assert skills.find_skill("TUNE", skills_dir=tmp_path) is not None
    assert skills.find_skill("tune", skills_dir=tmp_path) is not None


def test_find_skill_unknown_returns_none(tmp_path):
    (tmp_path / "x.md").write_text("# X\nbody")
    assert skills.find_skill("y", skills_dir=tmp_path) is None


def test_find_skill_empty_name(tmp_path):
    (tmp_path / "x.md").write_text("# X\nbody")
    assert skills.find_skill("", skills_dir=tmp_path) is None
    assert skills.find_skill("   ", skills_dir=tmp_path) is None


# ---------------------------------------------------------------------------
# format_skill_message
# ---------------------------------------------------------------------------

def test_format_skill_message_includes_header():
    skill = skills.Skill(
        name="x", title="My Title", description="", body="Do the thing.",
        source_path="",
    )
    msg = skills.format_skill_message(skill)
    assert msg.startswith("[Skill: My Title]")
    assert "Do the thing." in msg


def test_format_skill_message_strips_whitespace():
    skill = skills.Skill(
        name="x", title="T", description="", body="\n\n  body  \n\n",
        source_path="",
    )
    msg = skills.format_skill_message(skill)
    # No leading/trailing blank lines
    assert msg == msg.strip()
    assert "body" in msg
