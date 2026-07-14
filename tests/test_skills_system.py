"""Tests for delfin.agent.skills loader + executor."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent import skills as S
from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _write_skill(base: Path, name: str, body: str = "do x", desc: str = "") -> Path:
    skill_dir = base / ".delfin" / "skills" / name
    skill_dir.mkdir(parents=True, exist_ok=True)
    fm = f"---\nname: {name}\n"
    if desc:
        fm += f"description: {desc}\n"
    fm += "---\n\n"
    p = skill_dir / "SKILL.md"
    p.write_text(fm + body, encoding="utf-8")
    return p


def _write_flat_skill(base: Path, name: str, body: str = "flat") -> Path:
    skills_dir = base / ".delfin" / "skills"
    skills_dir.mkdir(parents=True, exist_ok=True)
    p = skills_dir / f"{name}.md"
    p.write_text(f"# {name.title()}\n\n{body}", encoding="utf-8")
    return p


def test_discover_finds_project_skills():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_skill(ws, "myskill", body="step 1", desc="A test skill")
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            skills = S.discover_skills(ws)
        assert len(skills) == 1
        assert skills[0].name == "myskill"
        assert skills[0].description == "A test skill"
        assert "step 1" in skills[0].body


def test_discover_skips_invalid_names():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        skill_dir = ws / ".delfin" / "skills" / "../danger"
        # `..` traversal is sanitised by _SKILL_NAME_RE; we just write
        # a normal skill folder named with a valid identifier and a
        # second one named with garbage to confirm garbage is dropped.
        good = _write_skill(ws, "good")
        bad_dir = ws / ".delfin" / "skills" / "@bad"
        bad_dir.mkdir(parents=True)
        (bad_dir / "SKILL.md").write_text(
            "---\nname: @bad\n---\nx", encoding="utf-8",
        )
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            skills = S.discover_skills(ws)
        names = {s.name for s in skills}
        assert "good" in names
        assert "@bad" not in names


def test_discover_flat_md_files():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_flat_skill(ws, "review", body="Run tests, check diff")
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            skills = S.discover_skills(ws)
        names = [s.name for s in skills]
        assert "review" in names
        sk = next(s for s in skills if s.name == "review")
        assert "Run tests" in sk.body


def test_project_overrides_user():
    with tempfile.TemporaryDirectory() as user_d, \
         tempfile.TemporaryDirectory() as proj_d:
        user_ws = Path(user_d)
        proj_ws = Path(proj_d)
        _write_skill(user_ws, "shared", body="USER VERSION")
        _write_skill(proj_ws, "shared", body="PROJECT VERSION")
        with patch.object(S, "_skill_dirs", return_value=[
            user_ws / ".delfin" / "skills",
            proj_ws / ".delfin" / "skills",
        ]):
            sk = S.get_skill("shared", proj_ws)
        assert sk is not None
        assert "PROJECT VERSION" in sk.body
        assert "USER VERSION" not in sk.body


def test_get_skill_returns_none_when_missing():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            assert S.get_skill("nope", ws) is None


def test_render_includes_args():
    sk = S.Skill(name="x", description="d", body="DO IT", source=Path("/x"))
    out = S.render_skill_invocation(sk, args="arg1 arg2")
    assert "Skill: x" in out
    assert "arg1 arg2" in out
    assert "DO IT" in out
    assert "end skill" in out


def test_render_no_args():
    sk = S.Skill(name="x", description="d", body="BODY", source=Path("/x"))
    out = S.render_skill_invocation(sk, args="")
    assert "Arguments:" not in out
    assert "BODY" in out


def test_skill_tool_dispatch_returns_body():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_skill(ws, "mytool", body="MAGIC", desc="Magic skill")
        perms = KitToolPermissions(workspace=ws, mode="default")
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            out = _DocToolExecutor().execute(
                "skill", {"name": "mytool"}, permissions=perms,
            )
        payload = json.loads(out)
        assert payload["status"] == "ok"
        assert payload["skill"] == "mytool"
        assert "MAGIC" in payload["content"]


def test_skill_tool_unknown_lists_available():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_skill(ws, "alpha")
        _write_skill(ws, "beta")
        perms = KitToolPermissions(workspace=ws, mode="default")
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            out = _DocToolExecutor().execute(
                "skill", {"name": "missing"}, permissions=perms,
            )
        payload = json.loads(out)
        assert "error" in payload
        assert sorted(payload["available"]) == ["alpha", "beta"]


def test_skill_tool_empty_name_rejected():
    perms = KitToolPermissions(workspace=Path("/tmp"), mode="default")
    out = _DocToolExecutor().execute(
        "skill", {"name": ""}, permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload


def test_first_heading_used_when_no_description():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        skill_dir = ws / ".delfin" / "skills" / "noheader"
        skill_dir.mkdir(parents=True)
        # No frontmatter, just a # heading
        (skill_dir / "SKILL.md").write_text(
            "# Auto-Header\n\nbody text", encoding="utf-8",
        )
        with patch.object(S, "_skill_dirs",
                          return_value=[ws / ".delfin" / "skills"]):
            sk = S.get_skill("noheader", ws)
        assert sk is not None
        assert sk.description == "Auto-Header"


def test_builtin_pack_skills_are_discoverable():
    """Regression: the 9 shipped 'curated set' pack skills must be reachable
    via the skill tool (they were previously undiscovered because _skill_dirs
    omitted the pack dir, so get_skill returned None for all of them)."""
    names = {s.name for s in S.discover_skills()}
    for expected in ("diagnose-failed-run", "casscf-setup",
                     "tddft-excited-states"):
        assert expected in names, f"pack skill {expected!r} not discovered"
        assert S.get_skill(expected) is not None


def test_project_skill_overrides_pack_skill(tmp_path):
    """A project-scoped skill wins over a same-named pack skill."""
    d = tmp_path / ".delfin" / "skills"
    d.mkdir(parents=True)
    (d / "casscf-setup.md").write_text(
        "---\nname: casscf-setup\ndescription: PROJECT_OVERRIDE\n---\nbody",
        encoding="utf-8",
    )
    sk = S.get_skill("casscf-setup", tmp_path)
    assert sk is not None and sk.description == "PROJECT_OVERRIDE"
