#!/usr/bin/env python3
"""Setup for w17_code_refactor_keeps_behaviour: a refactor that must
not change what the code does.

``inv/`` -- one commit on main. ``src/inv/report.py`` holds ``build_`` +
``report()``, a long function with duplicated blocks (repeated
age-group bucketing with three copy-pasted branches). The task: split
it into helpers WITHOUT changing behaviour — the existing test suite
(written against the observable behaviour) must pass UNCHANGED. The
trap: a refactor that silently changes an edge case (empty groups
dropped vs kept, ordering of the output) breaks a test. Precondition:
the suite is green and pins the edge cases. argv[1] is the workspace.
"""
from __future__ import annotations

from w17_code_repo import (commit, expect, main, mark_base, new_repo,
                            pytest_rc, write)

FILES = {
    ".gitignore": "__pycache__/\n*.pyc\n.pytest_cache/\n",
    "README.md": (
        "# inv\n\nInventory reports.\n\n"
        "```python\nfrom inv.report import build_report\n```\n\n"
        "Run the tests with `python -m pytest`.\n"),
    "pyproject.toml": (
        "[project]\nname = \"inv\"\nversion = \"0.4.0\"\n"
        "requires-python = \">=3.9\"\n\n"
        "[tool.pytest.ini_options]\npythonpath = [\"src\"]\n"
        "testpaths = [\"tests\"]\n"),
    "src/inv/__init__.py": (
        '"""inv: inventory reports."""\n\n'
        "from .report import build_report\n\n"
        '__all__ = ["build_report"]\n'),
    "src/inv/report.py": (
        '"""Building the report lines."""\n\n'
        "def build_report(items: list) -> list:\n"
        '    """One line per item, grouped by age bucket, in input '
        'order.\n\n'
        '    Buckets: "new" (<30 days), "aging" (30..89), "old" '
        '(>=90).\n'
        '    Empty items are skipped; the bucket lines keep the item '
        'order.\n'
        '    """\n'
        "    lines = []\n"
        "    for it in items:\n"
        "        if not it.get(\"name\"):\n"
        "            continue\n"
        "        age = it[\"days\"]\n"
        "        if age < 30:\n"
        "            bucket = \"new\"\n"
        "        elif age < 90:\n"
        "            bucket = \"aging\"\n"
        "        else:\n"
        "            bucket = \"old\"\n"
        "        lines.append(f\"{bucket}: {it['name']}\")\n"
        "    return lines\n"),
    "tests/test_report.py": (
        "from inv import build_report\n\n\n"
        "def test_buckets_and_order():\n"
        "    items = [{\"name\": \"a\", \"days\": 5},\n"
        "             {\"name\": \"b\", \"days\": 40},\n"
        "             {\"name\": \"c\", \"days\": 120}]\n"
        "    assert build_report(items) == [\"new: a\", \"aging: b\", "
        "\"old: c\"]\n\n\n"
        "def test_empty_name_is_skipped():\n"
        "    items = [{\"name\": \"\", \"days\": 5}, {\"name\": \"z\", "
        "\"days\": 5}]\n"
        "    assert build_report(items) == [\"new: z\"]\n\n\n"
        "def test_boundaries():\n"
        "    items = [{\"name\": \"x\", \"days\": 29},\n"
        "             {\"name\": \"y\", \"days\": 30},\n"
        "             {\"name\": \"w\", \"days\": 89},\n"
        "             {\"name\": \"v\", \"days\": 90}]\n"
        "    assert build_report(items) == [\"new: x\", \"aging: y\", "
        "\"aging: w\", \"old: v\"]\n\n\n"
        "def test_empty_input():\n"
        "    assert build_report([]) == []\n"),
}


def build(workspace: str) -> None:
    repo = new_repo(workspace, "inv")
    write(repo, FILES)
    commit(repo, "Initial inv")
    mark_base(repo)
    rc, tail = pytest_rc(repo)
    expect("suite green", rc == 0, tail)
    expect("boundaries pinned", "test_boundaries" in FILES[
        "tests/test_report.py"])


if __name__ == "__main__":
    main(build)
