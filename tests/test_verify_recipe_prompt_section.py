"""Controls for the verify-recipe prompt section (prompt_loader wiring).

Written BEFORE the wiring exists: on the unmodified tree no
``verify_recipe`` section is composed, so these fail.
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent.prompt_loader import PromptLoader


def _project(tmp_path: Path) -> Path:
    ci = tmp_path / ".github/workflows/ci.yml"
    ci.parent.mkdir(parents=True)
    ci.write_text(
        "on: [push]\njobs:\n  test:\n    steps:\n"
        "      - name: Run tests\n        run: pytest -q\n"
        "  lint:\n    steps:\n"
        "      - name: Lint\n        run: ruff check .\n",
        encoding="utf-8")
    return tmp_path


def _sections(loader: PromptLoader):
    return loader.compose_sections(
        role_id="solo_agent", mode_id="solo", task_text="fix a bug",
        session_key="vr-1")


def test_workspace_with_ci_gets_a_verify_recipe_section(tmp_path):
    loader = PromptLoader()
    loader.workspace_root = _project(tmp_path)
    names = [s.name for s in _sections(loader)]
    assert "verify_recipe" in names


def test_recipe_section_lists_the_discovered_commands(tmp_path):
    loader = PromptLoader()
    loader.workspace_root = _project(tmp_path)
    sec = next(s for s in _sections(loader) if s.name == "verify_recipe")
    assert "pytest -q" in sec.content
    assert "ruff check ." in sec.content
    assert "To check your work" in sec.content


def test_recipe_section_is_volatile_and_before_session_env(tmp_path):
    # The recipe is read from disk at build time; it belongs with the
    # per-turn material and must not sit in the cacheable head.
    loader = PromptLoader()
    loader.workspace_root = _project(tmp_path)
    sections = _sections(loader)
    idx = {s.name: i for i, s in enumerate(sections)}
    assert idx["verify_recipe"] > idx["role_prompt"]
    assert idx["verify_recipe"] < idx["session_env"]


def test_workspace_without_checks_gets_no_recipe_section(tmp_path):
    # Empty recipe renders to "" and the section is dropped entirely —
    # no guessing, no boilerplate.
    loader = PromptLoader()
    loader.workspace_root = tmp_path  # nothing in it
    names = [s.name for s in _sections(loader)]
    assert "verify_recipe" not in names
