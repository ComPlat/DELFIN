"""Tests for delfin.agent.project_introspect."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from delfin.agent import project_introspect as PI
from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _write_pyproject(d: Path, deps: list[str] | None = None,
                     pep621: bool = True) -> None:
    deps = deps or []
    text = ""
    if pep621:
        text += "[project]\nname = \"x\"\nversion = \"0.1\"\n"
        text += "dependencies = [\n"
        for d_ in deps:
            text += f'  "{d_}",\n'
        text += "]\n"
    (d / "pyproject.toml").write_text(text, encoding="utf-8")


def _make_fake_venv(base: Path, name: str = ".venv") -> Path:
    """Create a venv pointing at the live interpreter — fast, no install."""
    venv_dir = base / name
    subprocess.run(
        [sys.executable, "-m", "venv", "--without-pip", str(venv_dir)],
        check=True, capture_output=True,
    )
    return venv_dir


def test_introspect_empty_workspace():
    with tempfile.TemporaryDirectory() as d:
        report = PI.introspect(d)
        assert report["workspace"]
        assert report["venvs"] == []
        assert "summary_one_liner" in report
        assert "no venv" in report["summary_one_liner"]


def test_introspect_invalid_path():
    report = PI.introspect("/no/such/dir/12345")
    assert "error" in report


def test_introspect_detects_existing_venv():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _make_fake_venv(ws, ".venv")
        report = PI.introspect(ws)
        assert len(report["venvs"]) >= 1
        v = report["venvs"][0]
        assert v["path"] == ".venv"
        assert v["python"]
        assert "venv" in report["summary_one_liner"]


def test_introspect_detects_named_venv():
    """venvs with non-standard names like .venv-myproj should still be found."""
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _make_fake_venv(ws, ".venv-myproj")
        report = PI.introspect(ws)
        venv_paths = [v["path"] for v in report["venvs"]]
        assert ".venv-myproj" in venv_paths


def test_introspect_detects_pyproject_with_deps():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_pyproject(ws, deps=["numpy>=1.20", "pytest"])
        report = PI.introspect(ws)
        py = report["manifests"]["pyproject.toml"]
        assert py["present"] is True
        assert py["has_pep621"] is True
        assert any("numpy" in dep for dep in py["deps"])
        assert any("pytest" in dep for dep in py["deps"])


def test_introspect_detects_pytest_via_deps():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_pyproject(ws, deps=["pytest>=7.0"])
        report = PI.introspect(ws)
        assert report["test_framework"] == "pytest"


def test_introspect_detects_requirements_txt():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "requirements.txt").write_text(
            "numpy\npandas\npytest\n# comment\n",
        )
        report = PI.introspect(ws)
        req = report["manifests"]["requirements.txt"]
        assert req["present"] is True
        assert req["lines"] == 3
        assert "numpy" in req["packages"]
        # pytest in requirements counts as test framework
        assert report["test_framework"] == "pytest"


def test_introspect_detects_src_layout():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "src").mkdir()
        (ws / "src" / "__init__.py").touch()
        (ws / "tests").mkdir()
        report = PI.introspect(ws)
        assert report["source_layout"] == "src"


def test_introspect_detects_flat_layout():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "main.py").write_text("print(1)")
        (ws / "helpers.py").write_text("x = 1")
        report = PI.introspect(ws)
        assert report["source_layout"] == "flat"


def test_introspect_detects_git_repo():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        subprocess.run(["git", "init", "--quiet", "--initial-branch=main"],
                        cwd=str(ws), check=True)
        subprocess.run(["git", "config", "user.email", "t@t"], cwd=str(ws), check=True)
        subprocess.run(["git", "config", "user.name", "t"], cwd=str(ws), check=True)
        subprocess.run(["git", "config", "commit.gpgsign", "false"],
                        cwd=str(ws), check=True)
        (ws / "x.txt").write_text("x")
        subprocess.run(["git", "add", "."], cwd=str(ws), check=True)
        subprocess.run(["git", "commit", "--quiet", "-m", "init"],
                        cwd=str(ws), check=True)
        report = PI.introspect(ws)
        assert report["git"]["is_repo"] is True
        assert report["git"]["branch"] == "main"
        assert report["git"]["dirty"] is False


def test_introspect_summary_one_liner_combines():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _make_fake_venv(ws, ".venv")
        _write_pyproject(ws, deps=["pytest"])
        (ws / "src").mkdir()
        (ws / "src" / "__init__.py").touch()
        report = PI.introspect(ws)
        s = report["summary_one_liner"]
        assert "venv" in s
        assert "pytest" in s
        assert "src" in s


def test_introspect_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_pyproject(ws, deps=["numpy"])
        perms = KitToolPermissions(workspace=ws, mode="default")
        out = _DocToolExecutor().execute(
            "project_introspect", {}, permissions=perms,
        )
        payload = json.loads(out)
        assert payload["workspace"] == str(ws.resolve())
        assert payload["manifests"]["pyproject.toml"]["present"] is True


def test_introspect_no_perms_errors():
    out = _DocToolExecutor().execute(
        "project_introspect", {}, permissions=None,
    )
    payload = json.loads(out)
    assert "error" in payload
