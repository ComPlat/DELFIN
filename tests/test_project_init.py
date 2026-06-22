"""Tests for delfin.agent.project_init — the /init backend."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import project_init as pi


# --- detectors --------------------------------------------------------------

def test_detect_python_pyproject(tmp_path):
    (tmp_path / "pyproject.toml").write_text(
        '[project]\nname = "wibble"\ndescription = "a thing"\n',
        encoding="utf-8",
    )
    p = pi.detect_project(tmp_path)
    assert p.language == "Python"
    assert p.name == "wibble"
    assert p.description == "a thing"
    assert "pytest" in p.test_command


def test_detect_python_picks_uv_when_lockfile_present(tmp_path):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")
    (tmp_path / "uv.lock").write_text("# stub", encoding="utf-8")
    assert pi.detect_project(tmp_path).package_manager == "uv"


def test_detect_node_package_json(tmp_path):
    (tmp_path / "package.json").write_text(
        json.dumps({
            "name": "frontend",
            "description": "ui",
            "scripts": {"test": "jest", "build": "webpack", "lint": "eslint"},
        }), encoding="utf-8",
    )
    p = pi.detect_project(tmp_path)
    assert p.language.startswith("JavaScript")
    assert p.test_command.endswith("test")
    assert p.lint_command.endswith("lint")


def test_detect_rust_cargo_toml(tmp_path):
    (tmp_path / "Cargo.toml").write_text(
        '[package]\nname = "krate"\ndescription = "rusty thing"\n',
        encoding="utf-8",
    )
    p = pi.detect_project(tmp_path)
    assert p.language == "Rust"
    assert p.name == "krate"
    assert p.test_command == "cargo test"


def test_detect_go_module(tmp_path):
    (tmp_path / "go.mod").write_text("module example.com/foo\n", encoding="utf-8")
    p = pi.detect_project(tmp_path)
    assert p.language == "Go"
    assert p.test_command.startswith("go test")


def test_detect_unknown_falls_back(tmp_path):
    p = pi.detect_project(tmp_path)
    assert p.language == "(undetected)"
    assert p.name == tmp_path.name


def test_makefile_targets_surface_as_extras(tmp_path):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")
    (tmp_path / "Makefile").write_text(
        "test:\n\tpytest\n\nfmt:\n\truff format\n\nrelease:\n\tbuild-and-push\n",
        encoding="utf-8",
    )
    p = pi.detect_project(tmp_path)
    assert "make test" in (p.extras.get("test") or "")
    assert "make fmt" in (p.extras.get("fmt") or "")
    assert "make release" in (p.extras.get("release") or "")


# --- renderer ---------------------------------------------------------------

def test_render_agents_md_has_required_sections(tmp_path):
    p = pi.ProjectProfile(
        name="my-app", description="cool app",
        language="Python", package_manager="uv",
        test_command="pytest", build_command="uv build",
        lint_command="ruff check .", run_command="python -m my_app",
    )
    body = pi.render_agents_md(p, tmp_path)
    for needle in ("# AGENTS.md", "Stack", "Commands",
                   "Workflow", "Don't touch"):
        assert needle in body
    assert "uv build" in body
    assert "ruff check" in body


def test_render_agents_md_handles_empty_profile(tmp_path):
    p = pi.ProjectProfile(name="empty", language="(undetected)")
    body = pi.render_agents_md(p, tmp_path)
    assert "(set me)" in body or "(undetected)" in body
    assert "fill in" in body


# --- end-to-end init --------------------------------------------------------

def test_init_writes_files(tmp_path):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")
    out = pi.init_project(tmp_path)
    assert (tmp_path / "AGENTS.md").exists()
    assert (tmp_path / ".delfin" / "settings.json").exists()
    assert len(out["files"]) == 2
    assert out["profile"].language == "Python"


def test_init_skips_existing_unless_force(tmp_path):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")
    (tmp_path / "AGENTS.md").write_text("# preserved", encoding="utf-8")
    out = pi.init_project(tmp_path)
    assert (tmp_path / "AGENTS.md").read_text(encoding="utf-8") == "# preserved"
    assert any("AGENTS.md" in s for s in out["skipped"])


def test_init_force_overwrites(tmp_path):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")
    (tmp_path / "AGENTS.md").write_text("# old", encoding="utf-8")
    pi.init_project(tmp_path, overwrite=True)
    body = (tmp_path / "AGENTS.md").read_text(encoding="utf-8")
    assert "# AGENTS.md" in body  # the new header, not the old "# old"


def test_init_settings_json_is_valid(tmp_path):
    pi.init_project(tmp_path)
    settings = tmp_path / ".delfin" / "settings.json"
    data = json.loads(settings.read_text(encoding="utf-8"))
    assert "hooks" in data
    assert data["_meta"]["created_by"].endswith("/init")
