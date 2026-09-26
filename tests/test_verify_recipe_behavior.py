"""Controls part 2: for_files, render, and the read-only safety contract."""
from __future__ import annotations

import ast
from pathlib import Path

import delfin.agent.verify_recipe as vr


def _step(kind="test", command="pytest", origin="x"):
    return vr.Step(kind=kind, command=command, origin=origin)


def _proj(tmp_path, *test_files):
    for name in test_files:
        (tmp_path / "tests").mkdir(exist_ok=True)
        (tmp_path / "tests" / name).write_text("", encoding="utf-8")
    return tmp_path


def test_for_files_maps_source_file_to_its_test_file(tmp_path):
    _proj(tmp_path, "test_agent_prompt_loader.py", "test_api_client_cli.py")
    recipe = vr.Recipe(steps=[_step(command="pytest {paths}")])
    out = vr.for_files(recipe, ["delfin/agent/prompt_loader.py",
                                "delfin/agent/api_client.py"],
                       workspace=tmp_path)
    assert out == [vr.Step(kind="test", command="pytest "
                           "tests/test_agent_prompt_loader.py "
                           "tests/test_api_client_cli.py",
                           origin="x (narrowed to changed files)")]


def test_for_files_passes_test_files_through_unchanged(tmp_path):
    _proj(tmp_path, "test_verify_recipe.py")
    recipe = vr.Recipe(steps=[_step(command="pytest {paths}")])
    out = vr.for_files(recipe, ["tests/test_verify_recipe.py"],
                       workspace=tmp_path)
    assert out[0].command == "pytest tests/test_verify_recipe.py"


def test_for_files_keeps_steps_without_placeholder_intact(tmp_path):
    _proj(tmp_path, "test_agent_prompt_loader.py")
    recipe = vr.Recipe(steps=[_step(kind="lint", command="ruff check .")])
    out = vr.for_files(recipe, ["delfin/agent/prompt_loader.py"],
                       workspace=tmp_path)
    assert out[0].command == "ruff check ."


def test_for_files_drops_test_step_when_no_test_file_exists(tmp_path):
    recipe = vr.Recipe(steps=[_step(command="pytest {paths}")])
    out = vr.for_files(recipe, ["docs/readme.md"], workspace=tmp_path)
    assert out == []


def test_for_files_keeps_shared_prefix_once(tmp_path):
    _proj(tmp_path, "test_foo.py")
    recipe = vr.Recipe(steps=[_step(command="pytest {paths}")])
    out = vr.for_files(recipe, ["delfin/foo.py", "delfin/bar.py"],
                       workspace=tmp_path)
    # bar has no test file; foo maps once and is not duplicated
    assert out[0].command == "pytest tests/test_foo.py"


def test_render_lists_every_step_with_kind_and_origin():
    recipe = vr.Recipe(steps=[
        _step(kind="test", command="pytest -q",
              origin=".github/workflows/ci.yml, step 'Run tests'"),
        _step(kind="lint", command="ruff check .",
              origin="pyproject.toml"),
    ])
    text = vr.render(recipe)
    assert "To check your work" in text
    assert "pytest -q" in text and "ruff check ." in text
    assert "ci.yml" in text and "pyproject.toml" in text
    assert "Tests" in text and "Lint" in text


def test_render_empty_recipe_is_empty():
    assert vr.render(vr.Recipe(steps=[])) == ""


# ---- safety: read and show, never execute ---------------------------------

def test_module_imports_nothing_that_runs_commands():
    # The module must never import subprocess, os.system, pty, asyncio
    # subprocess machinery or anything else that can execute a command.
    src = Path(vr.__file__).read_text(encoding="utf-8")
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        else:
            continue
        for name in names:
            assert not name.startswith(("subprocess", "os.system", "pty",
                                        "shutil")), name


def test_docstring_states_read_only_contract():
    doc = vr.__doc__ or ""
    assert "read" in doc.lower()
    assert "never" in doc.lower() and "execut" in doc.lower()
    assert "grant" in doc.lower() or "permission" in doc.lower()
