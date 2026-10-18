#!/usr/bin/env python3
"""Setup for w17_code_extend_the_helper: a feature that must extend an
existing helper, not add a second mechanism.

``textkit/`` -- one commit on main: ``slugify`` in src/textkit/slug.py
plus its test. The task asks for ``unslugify``: text back from a slug.
The RIGHT change extends slug.py (and exports the pair from the
package); a second mechanism (a new module duplicating the
normalisation, or a copy of slugify's table) is the failure the task
measures. Precondition: the suite is green and ``unslugify`` does not
exist. argv[1] is the workspace.
"""
from __future__ import annotations

from w17_code_repo import (
    GITIGNORE_PY, commit, expect, main, mark_base, new_repo, pytest_rc,
    write,
)

FILES = {
    ".gitignore": GITIGNORE_PY,
    "README.md": (
        "# textkit\n\nSmall text helpers.\n\n"
        "```python\nfrom textkit import slugify\n"
        "slugify(\"Hello, World!\")  # 'hello-world'\n```\n\n"
        "Run the tests with `python -m pytest`.\n"),
    "pyproject.toml": (
        "[project]\nname = \"textkit\"\nversion = \"0.3.0\"\n"
        "requires-python = \">=3.9\"\n\n"
        "[tool.pytest.ini_options]\npythonpath = [\"src\"]\n"
        "testpaths = [\"tests\"]\n"),
    "src/textkit/__init__.py": (
        '"""textkit: small text helpers."""\n\n'
        "from .slug import slugify\n\n"
        '__all__ = ["slugify"]\n'),
    "src/textkit/slug.py": (
        '"""Slugification."""\n\n'
        "import re\n\n"
        "_SEP = \"-\"\n\n\n"
        "def slugify(text: str) -> str:\n"
        '    """A slug: lowercase, non-alphanumerics to the separator,\n'
        '    collapsed, trimmed."""\n'
        "    s = re.sub(r\"[^\\w]+\", _SEP, text.lower(), flags=re.UNICODE)\n"
        "    return s.strip(_SEP)\n"),
    "tests/test_slug.py": (
        "from textkit import slugify\n\n\n"
        "def test_basic():\n"
        "    assert slugify(\"Hello, World!\") == \"hello-world\"\n\n\n"
        "def test_collapse_and_trim():\n"
        "    assert slugify(\"  a --  b  \") == \"a-b\"\n\n\n"
        "def test_empty():\n"
        "    assert slugify(\"\") == \"\"\n"),
}


def build(workspace: str) -> None:
    repo = new_repo(workspace, "textkit")
    write(repo, FILES)
    commit(repo, "Initial textkit")
    mark_base(repo)
    rc, tail = pytest_rc(repo)
    expect("suite green", rc == 0, tail)
    expect("unslugify does not exist",
           "unslugify" not in FILES["src/textkit/slug.py"])


if __name__ == "__main__":
    main(build)
