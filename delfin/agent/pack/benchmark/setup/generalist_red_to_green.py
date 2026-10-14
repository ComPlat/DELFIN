#!/usr/bin/env python3
"""Setup for generalist_red_to_green: one failing test, one defect.

``datewin/`` -- ``business_days(start, end)`` counts weekdays in the
half-open range [start, end) and its loop runs ``while d <= end``, so
an end date on a weekday is counted. tests/test_business.py has three
tests; exactly one fails (Mon 2024-01-01 .. Mon 2024-01-08 -> 6, not 5).
The other two do not see the defect (an end on a Saturday; a reversed
range).
Preconditions: exactly one test fails, and it is that one.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

BUSINESS = '''"""Working-day arithmetic."""

from datetime import date, timedelta


def business_days(start: date, end: date) -> int:
    """Number of weekdays (Mon-Fri) in the half-open range [start, end).

    Raises ValueError when ``end`` is before ``start``.
    """
    if end < start:
        raise ValueError("end before start")
    days = 0
    d = start
    while d <= end:
        if d.weekday() < 5:
            days += 1
        d += timedelta(days=1)
    return days


def is_weekend(d: date) -> bool:
    """True on Saturday and Sunday."""
    return d.weekday() >= 5
'''

TEST = '''from datetime import date

import pytest

from datewin.business import business_days


def test_a_full_week_from_monday_to_monday_has_five_working_days():
    assert business_days(date(2024, 1, 1), date(2024, 1, 8)) == 5


def test_a_range_ending_on_saturday_counts_monday_to_friday():
    assert business_days(date(2024, 1, 1), date(2024, 1, 6)) == 5


def test_a_reversed_range_is_refused():
    with pytest.raises(ValueError):
        business_days(date(2024, 1, 8), date(2024, 1, 1))
'''

FILES = {
    ".gitignore": GITIGNORE_PY,
    "README.md": ("# datewin\n\nWorking-day arithmetic.\n\n"
                  "Tests: `python -m pytest`.\n"),
    "pyproject.toml": ("[project]\nname = \"datewin\"\nversion = \"0.1.0\"\n\n"
                       "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                       "testpaths = [\"tests\"]\n"),
    "datewin/__init__.py": "\"\"\"Working-day arithmetic.\"\"\"\n",
    "datewin/business.py": BUSINESS,
    "tests/test_business.py": TEST,
}


def build(workspace: str):
    repo = new_repo(workspace, "datewin")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add business_days"))
    rc, out = pytest_rc(repo, "-rf")
    expect("the suite is red", rc != 0, out)
    expect("exactly one test fails", "1 failed, 2 passed" in out, out)
    expect("the failing test is the Monday-to-Monday week",
           "test_a_full_week_from_monday_to_monday" in out, out)
    return repo


if __name__ == "__main__":
    main(build)
