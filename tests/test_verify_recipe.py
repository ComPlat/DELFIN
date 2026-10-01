"""Controls for delfin/agent/verify_recipe.py.

Written BEFORE the module exists: every test here must fail on the
unmodified tree (ModuleNotFoundError at best), so a pass afterwards
measures the module, not the absence of the defect.
"""
from __future__ import annotations

import sys
from pathlib import Path

import delfin.agent.verify_recipe as vr

REPO_ROOT = Path(__file__).resolve().parents[1]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# ---- fixtures: explicit project file ---------------------------------------

def test_explicit_file_wins_over_conflicting_ci(tmp_path):
    _write(tmp_path / ".delfin/verify.toml",
           '[[step]]\nkind = "test"\ncommand = "pytest -q tests/unit"\n'
           '[[step]]\nkind = "lint"\ncommand = "ruff check src"\n')
    _write(tmp_path / ".github/workflows/ci.yml",
           "on: [push]\njobs:\n  t:\n    steps:\n"
           "      - run: pytest tests/\n      - run: ruff check .\n")
    recipe = vr.discover(tmp_path)
    assert [s.command for s in recipe.steps] == [
        "pytest -q tests/unit", "ruff check src"]
    assert all(s.origin == ".delfin/verify.toml" for s in recipe.steps)
    assert [s.kind for s in recipe.steps] == ["test", "lint"]


def test_explicit_file_unknown_kind_is_not_invented_into_a_check(tmp_path):
    # kind must be one of test/lint/typecheck; anything else is refused
    # whole-file rather than silently dropped (a typo in the recipe should
    # be loud, not a quiet missing check).
    _write(tmp_path / ".delfin/verify.toml",
           '[[step]]\nkind = "tset"\ncommand = "pytest"\n')
    recipe = vr.discover(tmp_path)
    assert recipe.steps == []


# ---- fixtures: CI definition only ------------------------------------------

def test_ci_only_project(tmp_path):
    _write(tmp_path / ".github/workflows/main.yml",
           "on: [push]\njobs:\n"
           "  test:\n    steps:\n"
           "      - name: Run tests\n"
           "        run: pytest -q\n"
           "  lint:\n    steps:\n"
           "      - name: Lint\n"
           "        run: |\n"
           "          pip install ruff\n"
           "          ruff check src/ \\\n"
           "            --select F821\n")
    recipe = vr.discover(tmp_path)
    kinds = {s.kind: s for s in recipe.steps}
    assert set(kinds) == {"test", "lint"}
    assert kinds["test"].command == "pytest -q"
    assert "Run tests" in kinds["test"].origin
    assert "main.yml" in kinds["test"].origin
    # install lines are skipped; continuation backslashes joined
    assert kinds["lint"].command == "ruff check src/ --select F821"
    assert "Lint" in kinds["lint"].origin


def test_ci_job_without_a_check_step_contributes_nothing(tmp_path):
    _write(tmp_path / ".github/workflows/ci.yml",
           "on: [push]\njobs:\n  deploy:\n    steps:\n"
           "      - run: docker push example/app\n")
    recipe = vr.discover(tmp_path)
    assert recipe.steps == []
    assert vr.render(recipe) == ""


# ---- fixtures: pyproject / Makefile ----------------------------------------

def test_pyproject_only(tmp_path):
    _write(tmp_path / "pyproject.toml",
           "[tool.pytest.ini_options]\naddopts = \"-q\"\n")
    recipe = vr.discover(tmp_path)
    assert [s.kind for s in recipe.steps] == ["test"]
    assert recipe.steps[0].command == "pytest"
    assert "pyproject.toml" in recipe.steps[0].origin


def test_pyproject_without_ruff_section_gives_no_lint(tmp_path):
    # [tool.ruff] absent -> no lint step invented.
    _write(tmp_path / "pyproject.toml",
           "[tool.pytest.ini_options]\naddopts = \"-q\"\n")
    recipe = vr.discover(tmp_path)
    assert all(s.kind != "lint" for s in recipe.steps)


def test_makefile_only(tmp_path):
    _write(tmp_path / "Makefile",
           "test:\n\tpytest -q\n\nlint:\n\truff check .\n\n.PHONY: test lint\n")
    recipe = vr.discover(tmp_path)
    kinds = {s.kind: s for s in recipe.steps}
    assert kinds["test"].command == "make test"
    assert "Makefile" in kinds["test"].origin
    assert kinds["lint"].command == "make lint"


# ---- fixtures: nothing at all ----------------------------------------------

def test_empty_project_yields_empty_recipe(tmp_path):
    recipe = vr.discover(tmp_path)
    assert recipe.steps == []
    assert vr.render(recipe) == ""


# ---- DELFIN's own repo -------------------------------------------------------

def test_delfin_repo_recipe_comes_from_real_ci():
    recipe = vr.discover(REPO_ROOT)
    tests = [s for s in recipe.steps if s.kind == "test"]
    lints = [s for s in recipe.steps if s.kind == "lint"]
    assert tests and all(s.command.startswith("pytest") for s in tests)
    assert lints and all("ruff" in s.command for s in lints)
    assert all("ci.yml" in s.origin or ".yml" in s.origin
               for s in recipe.steps)
