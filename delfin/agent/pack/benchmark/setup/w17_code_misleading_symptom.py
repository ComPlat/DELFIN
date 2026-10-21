#!/usr/bin/env python3
"""Setup for w17_code_misleading_symptom: a repo whose test FAILS on an
assertion about format_pairs, while the defect is in parse_line.

``logfmt/`` -- one commit on main. ``parse_line`` (src/logfmt/parse.py)
is documented to preserve keys as written, but it lowercases them
(``.lower()`` on the key). The failing test asserts a roundtrip
``format_pairs(parse_line("A=1 b=2")) == "A=1 b=2"``: the failure
surfaces at format_pairs' output, so the symptom points at the
formatter, while the correct fix removes the lowercasing in the parser.
A second test (value containing "=") passes, so the parser's "="
handling is NOT the bug. Precondition: the suite is red on the
roundtrip case, and the passing test exercises the value-with-equals
path. argv[1] is the workspace.
"""
from __future__ import annotations

from w17_code_repo import (commit, expect, main, mark_base, new_repo,
                            pytest_rc, write)

FILES = {
    ".gitignore": "__pycache__/\n*.pyc\n.pytest_cache/\n",
    "README.md": (
        "# logfmt\n\nTiny log-line helpers. Keys are preserved exactly\n"
        "as written: round-tripping a line must return the same line.\n\n"
        "```python\nfrom logfmt import parse_line, format_pairs\n"
        "parse_line(\"A=1 b=2\")  # {'A': '1', 'b': '2'}\n```\n\n"
        "Run the tests with `python -m pytest`.\n"),
    "pyproject.toml": (
        "[project]\nname = \"logfmt\"\nversion = \"0.2.0\"\n"
        "requires-python = \">=3.9\"\n\n"
        "[tool.pytest.ini_options]\npythonpath = [\"src\"]\n"
        "testpaths = [\"tests\"]\n"),
    "src/logfmt/__init__.py": (
        '"""logfmt: tiny key=value line helpers."""\n\n'
        "from .parse import parse_line\nfrom .format import format_pairs\n\n"
        '__all__ = ["parse_line", "format_pairs"]\n'),
    "src/logfmt/parse.py": (
        '"""Parsing key=value lines."""\n\n\n'
        "def parse_line(line: str) -> dict:\n"
        '    """Parse "A=1 b=2" into {"A": "1", "b": "2"}.\n\n'
        '    Pairs are separated by whitespace; the FIRST "=" in a pair\n'
        '    splits key from value, so a value may itself contain "=".\n'
        '    Keys are preserved exactly as written.\n'
        '    """\n'
        "    pairs = {}\n"
        '    for tok in line.split():\n'
        '        key, _, value = tok.partition("=")\n'
        "        pairs[key.lower()] = value\n"
        "    return pairs\n"),
    "src/logfmt/format.py": (
        '"""Formatting pairs back to a line."""\n\n\n'
        "def format_pairs(pairs: dict) -> str:\n"
        '    """{"A": "1", "b": "2"} -> "A=1 b=2", keys in given order."""\n'
        '    return " ".join(f"{k}={v}" for k, v in pairs.items())\n'),
    "tests/test_format.py": (
        "from logfmt import format_pairs, parse_line\n\n\n"
        "def test_roundtrip_preserves_keys():\n"
        '    line = "A=1 b=2"\n'
        '    assert format_pairs(parse_line(line)) == line\n\n\n'
        "def test_value_with_equals_survives():\n"
        '    line = "url=http://x?a=1 b=2"\n'
        '    assert format_pairs(parse_line(line)) == line\n\n\n'
        "def test_keys_preserved_as_written():\n"
        '    assert parse_line("A=1 b=2") == {"A": "1", "b": "2"}\n'),
}


def build(workspace: str) -> None:
    repo = new_repo(workspace, "logfmt")
    write(repo, FILES)
    commit(repo, "Initial logfmt")
    mark_base(repo)
    # Precondition: the suite is red, and the FIRST failing test is the
    # roundtrip that asserts on format_pairs' output (the misleading
    # symptom; the key-lowering in parse.py is the cause).
    rc, tail = pytest_rc(repo, "tests/test_format.py")
    expect("suite red", rc != 0, tail)
    expect("the symptom is the roundtrip on format_pairs",
           "test_roundtrip_preserves_keys" in tail, tail)
    # The value-with-equals path passes: the "=" handling is not the bug.
    rc2, tail2 = pytest_rc(repo, "tests/test_format.py::"
                                "test_value_with_equals_survives")
    expect("value-with-equals passes", rc2 == 0, tail2)
    # The cause sits in parse.py, on the path the symptom does not name.
    expect("the cause sits in parse.py",
           "key.lower()" in FILES["src/logfmt/parse.py"])


if __name__ == "__main__":
    main(build)
