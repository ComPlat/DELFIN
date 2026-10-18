#!/usr/bin/env python3
"""Acceptance for w17_code_extend_the_helper (``textkit/``).

Required: whole suite green; hidden tests of ``unslugify`` pass (from
textkit AND from textkit.slug, slug->text is not unique so the API
returns the separator split back into spaces and can raise on an empty
slug); a test of unslugify was added under tests/; tests/test_slug.py
is unchanged. Graded only: the own tests have teeth (they fail against
a plausible wrong implementation that drops the suffix-underscore
rule); the diff stays small and inside src/textkit/, tests/, README,
CHANGELOG. argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from w17_code_check import diff_scope_check
from generalist_check import (
    Report, base, blob_at, pytest, read, tree_text, workspace_repo,
)

HIDDEN = '''
import pytest
from textkit import unslugify
from textkit.slug import unslugify as direct


def test_both_import_paths_name_one_function():
    assert unslugify is direct


def test_separator_becomes_a_space():
    assert unslugify("hello-world") == "hello world"


def test_a_run_of_separators_is_one_space():
    assert unslugify("a--b") == "a b"


def test_leading_and_trailing_separators_are_trimmed():
    assert unslugify("-a-b-").strip() == "a b"


def test_empty_slug_is_refused():
    with pytest.raises(ValueError):
        unslugify("")
'''

# A plausible wrong implementation: handles the empty slug silently
# instead of raising (the API contract says ValueError).
MUTANT = '''
def unslugify(slug: str) -> str:
    if not slug:
        return ""
    return " ".join(p for p in slug.split("-") if p)
'''


def main() -> int:
    repo = workspace_repo("textkit")
    r = Report("w17_code_extend_the_helper")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo) or "HEAD"
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    ok, tail = pytest(repo, "tests/test_w17_hidden.py",
                      extra={"tests/test_w17_hidden.py": HIDDEN})
    r.check("hidden behaviour tests pass", ok, detail=tail, weight=3)
    own = {k: v for k, v in tree_text(repo, (".py",)).items()
           if k.startswith("tests/") and "unslugify" in v}
    r.check("a test of unslugify was added under tests/", bool(own),
            detail=", ".join(sorted(own)))
    r.check("tests/test_slug.py unchanged",
            read(repo, "tests/test_slug.py")
            == (blob_at(repo, since, "tests/test_slug.py") or "\0"))
    if own:
        ok, tail = pytest(repo, *sorted(own),
                          extra={"src/textkit/slug.py": MUTANT})
        r.check("own tests reject a wrong implementation", not ok,
                required=False, detail=tail)
    diff_scope_check(r, repo,
                     ("src/textkit/", "tests/", "README.md",
                      "CHANGELOG.md"), 40)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
