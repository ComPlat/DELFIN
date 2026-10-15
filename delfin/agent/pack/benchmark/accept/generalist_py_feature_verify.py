#!/usr/bin/env python3
"""Acceptance for generalist_py_feature (``textkit/``).

Required: the whole suite is green; hidden tests of the specified
behaviour pass (``from textkit import truncate_words`` and
``textkit.words.truncate_words``); the attempt added a test under tests/
that calls ``truncate_words``; tests/test_slug.py is unchanged.
Graded only: the attempt's own tests fail against a plausible wrong
implementation (they have teeth); nothing outside src/textkit/, tests/,
README.md and CHANGELOG.md changed.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, changed_paths, pytest,
                              read, tree_text, workspace_repo)

HIDDEN = '''
import pytest
from textkit import truncate_words
from textkit.words import truncate_words as direct


def test_both_import_paths_name_one_function():
    assert truncate_words is direct


def test_cuts_to_the_limit_and_appends_the_suffix():
    assert truncate_words("one two three", 2) == "one two…"


def test_a_text_within_the_limit_comes_back_unchanged():
    assert truncate_words("one two", 2) == "one two"
    assert truncate_words("  one   two  ", 5) == "  one   two  "


def test_words_are_split_on_any_whitespace_and_joined_by_one_space():
    assert truncate_words("a  b\\tc\\nd", 3) == "a b c…"


def test_a_custom_suffix():
    assert truncate_words("a b c", 1, suffix="...") == "a..."


def test_limit_zero_keeps_only_the_suffix():
    assert truncate_words("a b", 0) == "…"
    assert truncate_words("", 0) == ""


def test_a_negative_limit_is_refused():
    with pytest.raises(ValueError):
        truncate_words("a b", -1)
'''

# Normalises whitespace even when nothing is cut: a plausible slip the
# hidden tests catch, used here to ask whether the attempt's own tests
# would have caught it.
MUTANT = '''
def truncate_words(text, limit, suffix="…"):
    if limit < 0:
        raise ValueError("limit must be >= 0")
    words = text.split()
    if len(words) <= limit:
        return " ".join(words)
    return " ".join(words[:limit + 1]) + suffix
'''


def main() -> int:
    repo = workspace_repo("textkit")
    r = Report("generalist_py_feature")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    ok, tail = pytest(repo, "tests/test_hidden_truncate.py",
                      extra={"tests/test_hidden_truncate.py": HIDDEN})
    r.check("hidden behaviour tests pass", ok, detail=tail, weight=3)
    own = {k: v for k, v in tree_text(repo, (".py",)).items()
           if k.startswith("tests/") and "truncate_words" in v}
    r.check("a test of truncate_words was added under tests/", bool(own),
            detail=", ".join(sorted(own)))
    r.check("tests/test_slug.py unchanged",
            read(repo, "tests/test_slug.py")
            == (blob_at(repo, since, "tests/test_slug.py") or "\0"))
    if own:
        ok, tail = pytest(repo, *sorted(own),
                          extra={"src/textkit/words.py": MUTANT})
        r.check("own tests reject a wrong implementation", not ok,
                required=False, detail=tail)
    allowed = ("src/textkit/", "tests/", "README.md", "CHANGELOG.md")
    stray = sorted(p for p in changed_paths(repo, since)
                   if not p.startswith(allowed))
    r.check("nothing outside the feature changed", not stray,
            required=False, detail=", ".join(stray))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
