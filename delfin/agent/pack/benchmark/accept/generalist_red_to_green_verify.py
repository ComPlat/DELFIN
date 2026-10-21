#!/usr/bin/env python3
"""Acceptance for generalist_red_to_green (``datewin/``).

Required: the whole suite is green; tests/test_business.py is
byte-identical to the given one (the test was the specification, not
the defect); hidden tests over other dates pass (a fix that matches the
one given case is not a fix).
Graded only: only datewin/business.py changed outside tests/; the fix
is small (at most 4 changed lines in non-test files).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, changed_paths, out,
                              pytest, read, workspace_repo)

HIDDEN = '''from datetime import date

from datewin.business import business_days


def test_hidden_same_day_is_empty():
    assert business_days(date(2024, 1, 3), date(2024, 1, 3)) == 0


def test_hidden_friday_to_monday_is_one():
    assert business_days(date(2024, 1, 5), date(2024, 1, 8)) == 1


def test_hidden_saturday_to_monday_is_none():
    assert business_days(date(2024, 1, 6), date(2024, 1, 8)) == 0


def test_hidden_two_weeks_are_ten():
    assert business_days(date(2024, 3, 4), date(2024, 3, 18)) == 10


def test_hidden_wednesday_to_tuesday_is_four():
    assert business_days(date(2024, 1, 3), date(2024, 1, 9)) == 4
'''


def _changed_lines(repo, since, paths) -> int:
    n = 0
    for line in out(repo, "diff", "--numstat", since, "--",
                    *paths).splitlines():
        parts = line.split("\t")
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            n += max(int(parts[0]), int(parts[1]))
    return n


def main() -> int:
    repo = workspace_repo("datewin")
    r = Report("generalist_red_to_green")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    r.check("the given test file is unchanged",
            read(repo, "tests/test_business.py")
            == (blob_at(repo, since, "tests/test_business.py") or "\0"),
            weight=2)
    ok, tail = pytest(repo, "tests/test_hidden_business.py",
                      extra={"tests/test_hidden_business.py": HIDDEN})
    r.check("hidden cases pass", ok, detail=tail, weight=2)
    code = sorted(p for p in changed_paths(repo, since)
                  if not p.startswith("tests/"))
    r.check("only datewin/business.py changed outside tests/",
            code == ["datewin/business.py"], required=False,
            detail=", ".join(code))
    n = _changed_lines(repo, since, code) if code else 0
    r.check("the fix is small (<= 4 lines)", 0 < n <= 4, required=False,
            detail=f"{n} lines")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
