#!/usr/bin/env python3
"""Acceptance for generalist_contribution (``durfmt/``).

Required: main still points at the setup's commit; the current branch
is ``fix/<topic>`` or ``feat/<topic>`` (lower case, dashes) with at least
one commit on top of main and no merge commits; every new commit subject
is a Conventional Commit of at most 72 characters without a trailing
period; all tracked changes are committed; the suite and hidden tests of
the hour formatting pass; PR_DESCRIPTION.md exists at the root with
the sections Summary, Changes and Testing and is in no commit.
Graded only: the commits touch only durfmt/, tests/, CHANGELOG.md and
README.md; a commit adds or changes a test; the Testing section names
a pytest command.
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import Report, base, out, pytest, read, workspace_repo

HIDDEN = '''from durfmt.core import format_duration


def test_hidden_whole_hour():
    assert format_duration(3600) == "1h"


def test_hidden_all_parts():
    assert format_duration(3661) == "1h 1m 1s"
    assert format_duration(3725) == "1h 2m 5s"


def test_hidden_hours_and_seconds_without_minutes():
    assert format_duration(7205) == "2h 5s"


def test_hidden_below_an_hour_is_unchanged():
    assert format_duration(125) == "2m 5s"
    assert format_duration(0) == "0s"
'''

_BRANCH = re.compile(r"^(fix|feat)/[a-z0-9]+(-[a-z0-9]+)*$")
_SUBJECT = re.compile(
    r"^(fix|feat|test|docs|refactor|chore)(\([a-z0-9_-]+\))?!?: \S.*[^.]$")


def main() -> int:
    repo = workspace_repo("durfmt")
    r = Report("generalist_contribution")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    r.check("main is unchanged", out(repo, "rev-parse", "main") == since,
            weight=2)
    branch = out(repo, "branch", "--show-current")
    r.check("on a fix/ or feat/ branch", bool(_BRANCH.match(branch)),
            detail=branch)
    shas = out(repo, "rev-list", "--reverse", f"{since}..HEAD").splitlines()
    r.check("at least one commit on top of main", bool(shas))
    merges = out(repo, "rev-list", "--merges", f"{since}..HEAD")
    r.check("no merge commits", not merges)
    subjects = [out(repo, "log", "-1", "--format=%s", s) for s in shas]
    bad = [s for s in subjects if len(s) > 72 or not _SUBJECT.match(s)]
    r.check("commit subjects follow the convention", bool(subjects)
            and not bad, weight=2, detail=" | ".join(bad or subjects))
    dirty = out(repo, "status", "--porcelain", "--untracked-files=no")
    r.check("all tracked changes committed", not dirty, detail=dirty)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    ok, tail = pytest(repo, "tests/test_hidden_hours.py",
                      extra={"tests/test_hidden_hours.py": HIDDEN})
    r.check("hidden hour tests pass", ok, detail=tail, weight=2)
    pr = read(repo, "PR_DESCRIPTION.md")
    heads = [h for h in ("## Summary", "## Changes", "## Testing")
             if re.search(rf"^{re.escape(h)}\s*$", pr, re.M)]
    r.check("PR_DESCRIPTION.md has Summary, Changes, Testing",
            len(heads) == 3, weight=2, detail=", ".join(heads))
    ever = out(repo, "log", "--all", "--format=", "--name-only",
               "--", "PR_DESCRIPTION.md")
    r.check("PR_DESCRIPTION.md is not committed", not ever)
    touched = set(out(repo, "diff", "--name-only", since,
                      "HEAD").splitlines()) if shas else set()
    stray = sorted(p for p in touched if not p.startswith(
        ("durfmt/", "tests/", "CHANGELOG.md", "README.md")))
    r.check("commits touch only the fix", not stray, required=False,
            detail=", ".join(stray))
    r.check("a commit changes a test",
            any(p.startswith("tests/") for p in touched), required=False)
    testing = pr.split("## Testing", 1)[-1] if "## Testing" in pr else ""
    r.check("Testing names the pytest command", "pytest" in testing,
            required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
