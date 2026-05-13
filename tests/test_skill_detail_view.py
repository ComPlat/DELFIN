"""Tests for /skills <name> body view + palette autocomplete bugfix.

Earlier code referenced ``Skill.title`` and ``Skill.slash_command`` which
do not exist on the dataclass — the palette crashed with AttributeError
the moment any skill was discovered. The bugfix derives the missing
metadata from the actual fields (``name``, ``description``, ``body``).
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _write_pack_skill(tmp_pack: Path, name: str, title: str, desc: str = "") -> Path:
    p = tmp_pack / f"{name}.md"
    front = f"# {title}\n\n"
    if desc:
        front = f"---\ndescription: {desc}\n---\n\n# {title}\n\n"
    p.write_text(front + "Body of the skill.\n", encoding="utf-8")
    return p


def test_skill_dataclass_does_not_carry_title_attr():
    """Regression guard: an instantiated Skill must NOT silently grow a
    `title` field. If someone re-adds one we want to know so the palette
    fallback can be removed cleanly."""
    from delfin.agent.skills import Skill

    sk = Skill(name="x", description="", body="b", source=Path("/tmp/x.md"))
    assert not hasattr(sk, "title")
    assert not hasattr(sk, "slash_command")


def test_discover_skills_picks_up_local_skill(tmp_path, monkeypatch):
    """Sanity check that the test fixture actually puts a skill on disk
    where ``discover_skills`` will find it — anchors the assertions below."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    sk_dir = tmp_path / ".delfin" / "skills" / "demo"
    sk_dir.mkdir(parents=True)
    (sk_dir / "SKILL.md").write_text(
        "---\nname: demo\ndescription: A demo skill\n---\n\n# Demo\nbody",
        encoding="utf-8",
    )
    from delfin.agent.skills import discover_skills, get_skill
    skills = discover_skills(tmp_path)
    assert any(s.name == "demo" for s in skills)
    sk = get_skill("demo", tmp_path)
    assert sk is not None
    assert "A demo skill" in (sk.description or "")
    assert "body" in sk.body


def test_palette_skill_filter_does_not_reference_missing_attrs():
    """The palette source code must no longer access ``skill.title`` or
    ``skill.slash_command`` (both crash because the dataclass has neither).
    We scan only NON-comment lines so the "why" comment can mention them."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    start = text.find("Append discovered skills, filtered by the same query")
    assert start > 0
    snippet = text[start: start + 2000]
    code_lines = [
        ln for ln in snippet.splitlines()
        if not ln.lstrip().startswith("#")
    ]
    code_only = "\n".join(code_lines)
    assert "skill.title" not in code_only, "palette still uses skill.title"
    assert "skill.slash_command" not in code_only, "palette still uses skill.slash_command"
