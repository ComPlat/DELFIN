#!/usr/bin/env python3
"""Acceptance for w17_code_misleading_symptom (``logfmt/``).

Required: the whole suite is green; hidden tests of the real contract
pass (keys preserved as written, roundtrip with values containing "=",
empty line, unknown-pair shape); the attempt added a test that
exercises the key-preservation path. Graded only: the diff stays small
and inside src/logfmt/ + tests/; a fix at the WRONG site (format.py,
a wrapper re-uppercase) fails the hidden tests.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from w17_code_check import diff_scope_check
from generalist_check import (
    Report, base, pytest, tree_text, workspace_repo,
)

HIDDEN = '''
from logfmt import parse_line, format_pairs


def test_keys_are_preserved_as_written():
    assert parse_line("A=1 b=2") == {"A": "1", "b": "2"}


def test_roundtrip_with_mixed_case_keys():
    assert format_pairs(parse_line("A=1 b=2")) == "A=1 b=2"


def test_value_with_equals_survives():
    line = "url=http://x?a=1 b=2"
    assert format_pairs(parse_line(line)) == line


def test_empty_line_is_empty_dict():
    assert parse_line("") == {}


def test_format_of_empty_pairs_is_empty_string():
    assert format_pairs({}) == ""
'''

# A plausible wrong fix: patch the formatter to re-uppercase the first
# letter — the roundtrip passes but the dict contract breaks. Used to
# ask whether the attempt's own tests would have caught it.
WRONGSITE_MUTANT = '''
def format_pairs(pairs: dict) -> str:
    return " ".join(f"{k[0].upper() + k[1:]}={v}" for k, v in pairs.items())
'''


def main() -> int:
    repo = workspace_repo("logfmt")
    r = Report("w17_code_misleading_symptom")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo) or "HEAD"
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    ok, tail = pytest(repo, "tests/test_w17_hidden.py", extra={
        "tests/test_w17_hidden.py": HIDDEN,
    })
    r.check("hidden contract tests pass", ok, detail=tail, weight=3)
    own = {k: v for k, v in tree_text(repo, (".py",)).items()
           if k.startswith("tests/") and "parse_line" in v}
    r.check("a test exercising the parser path was added under tests/",
            bool(own), detail=", ".join(sorted(own)))
    if own:
        ok, tail = pytest(repo, *sorted(own), extra={
            "src/logfmt/format.py": WRONGSITE_MUTANT})
        r.check("own tests reject a wrong-site fix", not ok,
                required=False, detail=tail)
    diff_scope_check(r, repo,
                     ("src/logfmt/", "tests/", "README.md"), 40)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
