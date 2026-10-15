#!/usr/bin/env python3
"""Setup for generalist_conventions: a repository with house rules.

``notes-app/`` -- AGENTS.md and CONTRIBUTING.md state four rules the
prompt does not repeat: unit tests in tests/unit/test_<module>.py; a
``- <module>: <change>`` line under ``## [Unreleased]`` in CHANGELOG.md
for every user-visible change; no print() in library code; a docstring
on every public function. tools/lint.py enforces the last two.
Preconditions: lint and suite green, ``word_count`` absent.
argv[1] is the workspace.
"""
from __future__ import annotations

import subprocess
import sys

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

AGENTS = """# Instructions for coding agents

These rules apply to every change in this repository.

1. Unit tests live in `tests/unit/`, one file per module:
   `notes/<module>.py` is tested in `tests/unit/test_<module>.py`.
   Do not create test files anywhere else.
2. Every user-visible change adds one line under `## [Unreleased]` in
   `CHANGELOG.md`, in the form `- <module>: <what changed>`.
3. Library code (`notes/`) never calls `print()`; use the module logger.
4. Every public function has a docstring.
5. Before you finish, run `python tools/lint.py` and `python -m pytest`;
   both must pass.
"""

CONTRIBUTING = """# Contributing

Read AGENTS.md first; the same rules apply to people. In short: tests in
`tests/unit/test_<module>.py`, a changelog line under `## [Unreleased]`,
no `print()` in `notes/`, docstrings on public functions, and
`python tools/lint.py` plus `python -m pytest` green.
"""

CHANGELOG = """# Changelog

## [Unreleased]

## [0.4.0] - 2026-09-30

- storage: notes are saved atomically.
- stats: add `char_count`.
"""

LINT = '''"""House lint: no print() in notes/, docstrings on public functions.

Exit 0 when clean, 1 with one line per finding.
"""
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def findings():
    for path in sorted((ROOT / "notes").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        rel = path.relative_to(ROOT)
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "print"):
                yield f"{rel}:{node.lineno}: print() in library code"
            if (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and not node.name.startswith("_")
                    and ast.get_docstring(node) is None):
                yield f"{rel}:{node.lineno}: {node.name} has no docstring"


if __name__ == "__main__":
    problems = list(findings())
    for p in problems:
        print(p)
    sys.exit(1 if problems else 0)
'''

STATS = '''"""Statistics about a note."""

import logging

log = logging.getLogger(__name__)


def char_count(text: str) -> int:
    """Number of characters in ``text``, whitespace excluded."""
    return sum(1 for c in text if not c.isspace())
'''

STORAGE = '''"""Saving notes."""

import logging
import os
from pathlib import Path

log = logging.getLogger(__name__)


def save(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` atomically."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
    log.debug("saved %s", path)
'''

FILES = {
    ".gitignore": GITIGNORE_PY,
    "AGENTS.md": AGENTS,
    "CONTRIBUTING.md": CONTRIBUTING,
    "CHANGELOG.md": CHANGELOG,
    "README.md": "# notes-app\n\nA small note keeper.\n",
    "pyproject.toml": ("[project]\nname = \"notes-app\"\nversion = \"0.4.0\"\n\n"
                       "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                       "testpaths = [\"tests\"]\n"),
    "tools/lint.py": LINT,
    "notes/__init__.py": "\"\"\"A small note keeper.\"\"\"\n",
    "notes/stats.py": STATS,
    "notes/storage.py": STORAGE,
    "tests/unit/test_stats.py": (
        "from notes.stats import char_count\n\n\n"
        "def test_char_count_skips_whitespace():\n"
        "    assert char_count(\"a b\\nc\") == 3\n"),
    "tests/unit/test_storage.py": (
        "from notes.storage import save\n\n\n"
        "def test_save_writes_the_text(tmp_path):\n"
        "    p = tmp_path / \"n.md\"\n    save(p, \"hi\")\n"
        "    assert p.read_text(encoding=\"utf-8\") == \"hi\"\n"),
}


def build(workspace: str):
    repo = new_repo(workspace, "notes-app")
    write(repo, FILES)
    mark_base(repo, commit(repo, "notes-app 0.4.0"))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    lint = subprocess.run([sys.executable, "tools/lint.py"], cwd=str(repo),
                          capture_output=True, text=True, timeout=60)
    expect("lint is clean", lint.returncode == 0, lint.stdout)
    return repo


if __name__ == "__main__":
    main(build)
