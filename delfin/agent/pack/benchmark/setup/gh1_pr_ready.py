#!/usr/bin/env python3
"""Task setup entry: bring a feature to a PR-ready branch (gh1_pr_ready).

Builds a small Python package on ``main`` with its own tests, plus a
second branch ``docs/notes`` that holds only files the task says the
branch must NOT carry (a probe test file and a scratch note). The
agent must create its feature branch, commit the fix, run the tests,
merge nothing foreign, and leave ``main`` untouched.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gh1_repo import feature_branch, main, mark_base, new_repo, write  # noqa: E402
from generalist_repo import commit  # noqa: E402

PKG = {
    "conftest.py": (
        'import sys\nfrom pathlib import Path\n\n'
        'sys.path.insert(0, str(Path(__file__).parent / "src"))\n'),
    "src/widgetkit/__init__.py": 'from .core import normalise\n\n'
                                 '__all__ = ["normalise"]\n',
    "src/widgetkit/core.py": (
        '"""Normalise widget labels."""\n'
        'import re\n\n'
        '_SEP = re.compile(r"[-_\\s]+")\n\n'
        'def normalise(label: str) -> str:\n'
        '    """A label becomes lower-case, separators single spaces."""\n'
        '    return _SEP.sub(" ", label.strip().lower())\n'),
    "tests/test_core.py": (
        'from widgetkit import normalise\n\n\n'
        'def test_basic():\n'
        '    assert normalise(" Hello-World ") == "hello world"\n'),
}

PROBE = {
    # The probe/scratch files that must NOT ride along on a PR branch.
    "tests/test_zz_probe.py": (
        "def test_probe():\n"
        "    assert True\n"),
    "notes/scratch.txt": "throwaway notes, do not commit\n",
}


def build(workspace: str) -> Path:
    repo = new_repo(workspace, "widgetkit")
    write(repo, PKG)
    base = commit(repo, "initial widgetkit package")
    mark_base(repo, base)
    feature_branch(repo, "docs/notes", PROBE, "scratch notes (probe files)")
    git_switch_main(repo)
    return repo


def git_switch_main(repo: Path) -> None:
    from generalist_repo import git
    git(repo, "checkout", "-q", "main")


if __name__ == "__main__":
    main(build)
