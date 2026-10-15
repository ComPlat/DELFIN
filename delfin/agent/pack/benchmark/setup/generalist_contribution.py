#!/usr/bin/env python3
"""Setup for generalist_contribution: fix a bug and prepare the PR.

``durfmt/`` -- ``format_duration(seconds)`` never emits hours:
3600 -> "60m", 3661 -> "61m 1s". CONTRIBUTING.md sets the workflow: a
branch ``fix/<topic>`` (never commit to main), Conventional-Commit
subjects of at most 72 characters without a trailing period, no
unrelated changes in a fix, and a PR description with the sections
Summary, Changes, Testing. There is no remote; the prompt asks for the
description in PR_DESCRIPTION.md at the repository root, not committed.
Precondition: the given suite is green (it does not test hours).
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

CONTRIBUTING = """# Contributing to durfmt

## Workflow

- Never commit to `main`. Work on a branch named `fix/<short-topic>` for
  a bug fix or `feat/<short-topic>` for a feature (lower case, dashes).
- Commit subjects follow Conventional Commits:
  `<type>(<optional scope>): <summary>` with type one of `fix`, `feat`,
  `test`, `docs`, `refactor`, `chore`. At most 72 characters, imperative
  mood, no trailing period.
- One logical change per commit. A fix and its regression test may share
  a `fix:` commit or be split into `test:` + `fix:`; never mix unrelated
  edits (reformatting, renames) into a fix.
- Every bug fix comes with a regression test in `tests/`.

## Pull requests

The PR description has exactly these sections:

    ## Summary
    ## Changes
    ## Testing

`Testing` says which command you ran and its result.
"""

DURFMT = '''"""Human-readable durations."""


def format_duration(seconds: int) -> str:
    """``seconds`` as e.g. "1h 2m 3s"; zero parts are left out, 0 -> "0s"."""
    if seconds < 0:
        raise ValueError("negative duration")
    minutes, secs = divmod(seconds, 60)
    parts = []
    if minutes:
        parts.append(f"{minutes}m")
    if secs or not parts:
        parts.append(f"{secs}s")
    return " ".join(parts)
'''

TEST = '''import pytest

from durfmt.core import format_duration


def test_seconds_only():
    assert format_duration(42) == "42s"


def test_minutes_and_seconds():
    assert format_duration(125) == "2m 5s"


def test_zero():
    assert format_duration(0) == "0s"


def test_negative_is_refused():
    with pytest.raises(ValueError):
        format_duration(-1)
'''

FILES = {
    ".gitignore": GITIGNORE_PY,
    "CONTRIBUTING.md": CONTRIBUTING,
    "README.md": ("# durfmt\n\n`format_duration(3725)` -> `\"1h 2m 5s\"`.\n"),
    "pyproject.toml": ("[project]\nname = \"durfmt\"\nversion = \"1.0.1\"\n\n"
                       "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                       "testpaths = [\"tests\"]\n"),
    "durfmt/__init__.py": "\"\"\"Human-readable durations.\"\"\"\n",
    "durfmt/core.py": DURFMT,
    "tests/test_core.py": TEST,
}


def build(workspace: str):
    repo = new_repo(workspace, "durfmt")
    write(repo, FILES)
    mark_base(repo, commit(repo, "chore: release 1.0.1"))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    return repo


if __name__ == "__main__":
    main(build)
