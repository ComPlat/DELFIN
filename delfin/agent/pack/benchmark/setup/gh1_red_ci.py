#!/usr/bin/env python3
"""Task setup entry: reproduce a red CI test locally and fix it (gh1_red_ci).

The fixture carries a CI log from a failed run and a bug the log names
only loosely (a test name that also appears in two other test files):
``tests/test_label.py::test_separator_collapses``. The agent must find
the failing test by its id, reproduce it locally, fix the cause in
src/widgetkit/core.py, and leave everything else untouched.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import (  # noqa: E402
    commit, expect, git, main, mark_base, new_repo, pytest_rc, write,
)

CI_LOG = {
    "conftest.py": (
        'import sys\nfrom pathlib import Path\n\n'
        'sys.path.insert(0, str(Path(__file__).parent / "src"))\n'),
    ".github-ci-log.txt": (
        "2026-10-31T02:11:04Z == CMD migrated to the branch checkout ==\n"
        "checkout --\n"
        "git: clone agent/gh1-f6 -> merge queue: OK\n"
        "2026-10-31T02:12:31Z == Run fast test suite ==\n"
        "......F.........\n"
        "FAILED tests/test_label.py::test_separator_collapses - "
        "assert 'a   b' == 'a b'\n"
        "1 failed, 9 passed in 4.12s\n"
        "##[error] Process completed with exit code 1.\n"),
}

PKG = {
    "src/widgetkit/__init__.py": 'from .core import normalise\n\n'
                                 '__all__ = ["normalise"]\n',
    "src/widgetkit/core.py": (
        '"""Normalise widget labels."""\n'
        'import re\n\n'
        '# NOTE: kept deliberately without the + quantifier; do not change\n'
        '# without a test.\n'
        '_SEP = re.compile(r"[-_\\s]")\n\n'
        'def normalise(label: str) -> str:\n'
        '    """A label becomes lower-case, separators single spaces."""\n'
        '    return _SEP.sub(" ", label.strip().lower())\n'),
    "tests/test_label.py": (
        'from widgetkit import normalise\n\n\n'
        'def test_basic():\n'
        '    assert normalise("Hello-World") == "hello world"\n\n\n'
        'def test_separator_collapses():\n'
        '    assert normalise("a -  -  b") == "a b"\n\n\n'
        'def test_case_is_lowered():\n'
        '    assert normalise("MIXED Case") == "mixed case"\n'),
    "tests/test_other_label.py": (
        'from widgetkit import normalise\n\n\n'
        'def test_same_name_elsewhere():\n'
        '    assert normalise("x-y") == "x y"\n'),
    "tests/test_third_label.py": (
        'from widgetkit import normalise\n\n\n'
        'def test_same_name_third_file():\n'
        '    assert normalise("p-q") == "p q"\n'),
}


def build(workspace: str) -> Path:
    repo = new_repo(workspace, "widgetkit")
    write(repo, {**PKG, **CI_LOG})
    sha = commit(repo, "initial widgetkit package")
    mark_base(repo, sha)
    # The bug is already in the committed code: the test suite is red.
    rc, out = pytest_rc(repo, "tests/test_label.py")
    expect("fixture suite reproduces the red test", rc != 0,
           out[-400:])
    expect("the failing test id is in the log",
           "test_separator_collapses" in CI_LOG[".github-ci-log.txt"])
    return repo


if __name__ == "__main__":
    main(build)
