"""A CI file's test and lint commands survive the install filter.

Input: a ``.github/workflows/*.yml`` whose ``run:`` lines mix installs
and checks. Output: ``verify_recipe.discover(ws).steps``.

``_INSTALL_HEADS`` holds ``python``, ``npm``, ``yarn`` and ``pnpm``, and
the filter ran on the head word BEFORE classification, so ``npm test``,
``pnpm run lint`` and ``python -m pytest`` were dropped as installs --
although ``_classify_command`` was written to accept npm/yarn/pnpm.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import delfin.agent.verify_recipe as vr


def _ci(ws: Path, run_lines: list[str]) -> None:
    body = "".join(f"      - run: {line}\n" for line in run_lines)
    path = ws / ".github" / "workflows" / "ci.yml"
    path.parent.mkdir(parents=True)
    path.write_text("on: [push]\njobs:\n  check:\n    runs-on: x\n"
                    "    steps:\n" + body, encoding="utf-8")


def _steps(ws: Path) -> list[tuple[str, str]]:
    return [(s.kind, s.command) for s in vr.discover(ws).steps]


@pytest.mark.parametrize("command,kind", [
    ("npm test", "test"),
    ("npm run test", "test"),
    ("yarn test", "test"),
    ("pnpm run lint", "lint"),
    ("python -m pytest -q", "test"),
    ("python3 -m pytest tests/", "test"),
    ("python -m mypy src", "typecheck"),
    ("cargo test --all", "test"),
    ("cargo clippy -- -D warnings", "lint"),
    ("go test ./...", "test"),
    ("go vet ./...", "lint"),
])
def test_a_check_command_survives(tmp_path, command, kind):
    _ci(tmp_path, ["pip install -e .", "npm ci", command])
    assert (kind, command) in _steps(tmp_path)


@pytest.mark.parametrize("command", [
    "pip install -e .", "npm ci", "npm install", "pnpm install",
    "python -m pip install -U pip", "yarn install --frozen-lockfile",
    "python setup.py build", "echo done",
])
def test_an_install_is_still_dropped(tmp_path, command):
    _ci(tmp_path, [command, "pytest -q"])
    assert _steps(tmp_path) == [("test", "pytest -q")]
