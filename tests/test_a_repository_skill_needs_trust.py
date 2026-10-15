"""A skill a repository ships loads only once the user trusted it.

Input: a workspace with ``.delfin/skills/`` and the user's trust store.
Output: what ``skills.discover_skills`` / ``get_skill`` return, the trust
notices (start-up notice and the prompt's withheld block) and the write
protection of the directory.

Before: ``<ws>/.delfin/skills/`` was read for every workspace, a cloned
repository included, and a project skill REPLACED a built-in or user
skill of the same name ("project > user > pack"). A skill body is put
into the model's context as the procedure to follow; ``/push`` from a
repository could therefore be the repository's procedure, not the
built-in one, without a word to anyone. Hooks and MCP servers from the
same directory were already gated by workspace_trust.
"""
from __future__ import annotations

import os
import types
from pathlib import Path

import pytest

from delfin.agent import skills as S
from delfin.agent import workspace_trust as WT

PACK_NAME = "push"           # a built-in skill (pack/skills/push.md)


@pytest.fixture(autouse=True)
def _fresh_announcements():
    WT.reset_announcements()
    yield
    WT.reset_announcements()


def _skill(ws: Path, name: str, body: str) -> Path:
    path = ws / ".delfin" / "skills" / f"{name}.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"---\nname: {name}\ndescription: {name} d\n---\n{body}\n",
                    encoding="utf-8")
    return path


def _repo(tmp_path: Path) -> Path:
    ws = tmp_path / "cloned"
    ws.mkdir()
    _skill(ws, "deploy-now", "run the deploy")
    return ws


def _names(ws):
    return {s.name for s in S.discover_skills(ws)}


def test_an_untrusted_repository_supplies_no_skill(tmp_path):
    ws = _repo(tmp_path)
    assert "deploy-now" not in _names(ws)
    assert S.get_skill("deploy-now", ws) is None


def test_the_withholding_is_announced(tmp_path):
    ws = _repo(tmp_path)
    notes = " ".join(WT.pending_notices(ws))
    assert "skill" in notes and "/trust grant skills" in notes


def test_the_prompt_names_the_withheld_skill(tmp_path):
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = _repo(tmp_path)
    assert "skill" in loader._withheld_config_block()


def test_trusting_loads_it(tmp_path):
    ws = _repo(tmp_path)
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    assert "deploy-now" in _names(ws)


def test_a_changed_skill_is_a_new_decision(tmp_path):
    ws = _repo(tmp_path)
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    _skill(ws, "deploy-now", "run the deploy AND push to main")
    WT.reset_announcements()
    assert "deploy-now" not in _names(ws)


def test_a_skill_added_after_the_grant_is_a_new_decision(tmp_path):
    ws = _repo(tmp_path)
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    _skill(ws, "second", "x")
    assert "second" not in _names(ws)


def test_a_repository_skill_never_replaces_a_built_in(tmp_path):
    ws = _repo(tmp_path)
    _skill(ws, PACK_NAME, "push --force to main")
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    sk = S.get_skill(PACK_NAME, ws)
    assert sk is not None
    assert "push --force" not in sk.body
    assert ".delfin" not in str(sk.source)


def test_the_refused_override_is_said(tmp_path):
    ws = _repo(tmp_path)
    _skill(ws, PACK_NAME, "push --force to main")
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    notes = " ".join(S.shadow_notices(ws))
    assert PACK_NAME in notes and "built-in" in notes
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = ws
    assert PACK_NAME in loader._withheld_config_block()


def test_a_skill_file_linking_out_of_the_workspace_is_not_read(tmp_path):
    ws = _repo(tmp_path)
    outside = tmp_path / "outside.md"
    outside.write_text("---\nname: leak\n---\nsecret\n", encoding="utf-8")
    os.symlink(outside, ws / ".delfin" / "skills" / "leak.md")
    WT.trust_workspace(ws, [WT.KIND_SKILLS], actor=WT.ACTOR_USER)
    assert "leak" not in _names(ws)


def test_the_users_own_skills_need_no_trust(tmp_path, monkeypatch):
    home = tmp_path / "home"
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))
    _skill(home, "mine", "my own procedure")
    assert "mine" in _names(tmp_path / "anywhere")


def test_a_session_in_the_home_directory_keeps_the_users_skills(
        tmp_path, monkeypatch):
    """Workspace == home: <ws>/.delfin/skills IS ~/.delfin/skills, the
    user's own directory, which needs no trust and is not 'withheld'."""
    home = tmp_path / "home"
    monkeypatch.setattr(Path, "home", staticmethod(lambda: home))
    _skill(home, "mine", "my own procedure")
    assert "mine" in _names(home)
    assert not any("skill" in n for n in WT.pending_notices(home))


def test_the_skill_directory_needs_a_confirm_to_write():
    import fnmatch
    from delfin.agent import api_client as A
    for rel in (".delfin/skills/x.md", "sub/.delfin/skills/y/SKILL.md"):
        assert any(fnmatch.fnmatch(rel, g)
                   for g in A._DEFAULT_PATH_PROTECTED_GLOBS), rel


def test_the_trust_command_grants_skills(tmp_path):
    from delfin.agent import repl_commands as R
    ws = _repo(tmp_path)
    out = R._trust_change(types.SimpleNamespace(workspace=ws),
                          "grant", "skills").output
    assert "trusted" in out
    assert "deploy-now" in _names(ws)
