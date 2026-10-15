#!/usr/bin/env python3
"""Setup for generalist_py_feature: a src-layout Python package.

``textkit/`` -- one commit on main: ``slugify`` in src/textkit/slug.py,
its test, pyproject with ``[tool.pytest.ini_options] pythonpath``.
Precondition: the suite is green and ``truncate_words`` does not exist.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

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
        "\"\"\"Small text helpers.\"\"\"\n\n"
        "from .slug import slugify\n\n__all__ = [\"slugify\"]\n"),
    "src/textkit/slug.py": (
        "\"\"\"URL slugs.\"\"\"\n\nimport re\nimport unicodedata\n\n"
        "_NON_WORD = re.compile(r\"[^a-z0-9]+\")\n\n\n"
        "def slugify(text: str) -> str:\n"
        "    \"\"\"Lower-case ASCII words of ``text`` joined by dashes.\"\"\"\n"
        "    ascii_text = (unicodedata.normalize(\"NFKD\", text)\n"
        "                  .encode(\"ascii\", \"ignore\").decode(\"ascii\"))\n"
        "    return _NON_WORD.sub(\"-\", ascii_text.lower()).strip(\"-\")\n"),
    "tests/test_slug.py": (
        "from textkit import slugify\n\n\n"
        "def test_slugify_joins_words_with_dashes():\n"
        "    assert slugify(\"Hello, World!\") == \"hello-world\"\n\n\n"
        "def test_slugify_drops_accents():\n"
        "    assert slugify(\"Crème brûlée\") == \"creme-brulee\"\n"),
}


def build(workspace: str):
    repo = new_repo(workspace, "textkit")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add slugify"))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    return repo


if __name__ == "__main__":
    main(build)
