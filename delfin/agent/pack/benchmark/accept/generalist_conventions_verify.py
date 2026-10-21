#!/usr/bin/env python3
"""Acceptance for generalist_conventions (``notes-app/``).

The prompt asks for a feature and says "follow the repository's
conventions" without naming them; they are in AGENTS.md. Required:
the hidden tests of ``word_count`` pass; the suite and tools/lint.py
pass (no print(), docstrings); tests/unit/test_stats.py tests
``word_count``; no test file was created outside tests/unit/; the
``## [Unreleased]`` section of CHANGELOG.md gained a line
``- stats: ...`` that mentions word_count.
Graded only: older changelog entries are untouched; AGENTS.md,
CONTRIBUTING.md and tools/lint.py are unchanged.
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import subprocess
import sys

from generalist_check import (Report, base, blob_at, changed_paths, pytest,
                              read, workspace_repo)

HIDDEN = '''from notes.stats import word_count


def test_hidden_counts_words_split_on_whitespace():
    assert word_count("one two\\nthree\\tfour") == 4


def test_hidden_heading_marks_are_not_words():
    assert word_count("# Title here\\n\\n## Sub\\ntext") == 4


def test_hidden_an_empty_note_has_none():
    assert word_count("") == 0
    assert word_count("  \\n ") == 0


def test_hidden_a_hash_inside_a_word_is_a_word():
    assert word_count("C# and #tag") == 3
'''

_UNRELEASED = re.compile(r"^## \[Unreleased\]\s*$(.*?)(?=^## \[|\Z)",
                         re.M | re.S)


def _unreleased(text: str) -> str:
    m = _UNRELEASED.search(text)
    return m.group(1) if m else ""


def main() -> int:
    repo = workspace_repo("notes-app")
    r = Report("generalist_conventions")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo, "tests/unit/test_hidden_wc.py",
                      extra={"tests/unit/test_hidden_wc.py": HIDDEN})
    r.check("hidden word_count tests pass", ok, detail=tail, weight=2)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    lint = subprocess.run([sys.executable, "tools/lint.py"], cwd=str(repo),
                          capture_output=True, text=True, timeout=60)
    r.check("tools/lint.py clean", lint.returncode == 0,
            detail=lint.stdout.strip()[:200])
    r.check("tests/unit/test_stats.py tests word_count",
            "word_count" in read(repo, "tests/unit/test_stats.py"))
    stray = sorted(p for p in changed_paths(repo, since)
                   if p.endswith(".py") and "test" in p.rsplit("/", 1)[-1]
                   and not p.startswith("tests/unit/"))
    r.check("no test file outside tests/unit/", not stray,
            detail=", ".join(stray))
    now = _unreleased(read(repo, "CHANGELOG.md"))
    entry = [ln for ln in now.splitlines()
             if re.match(r"^- stats: .*word_count", ln.strip())]
    r.check("CHANGELOG [Unreleased] has '- stats: ... word_count'",
            bool(entry), weight=2, detail=now.strip()[:200])
    old = (blob_at(repo, since, "CHANGELOG.md") or "").split("## [0.4.0]", 1)
    cur = read(repo, "CHANGELOG.md").split("## [0.4.0]", 1)
    r.check("released entries untouched",
            len(old) == 2 and len(cur) == 2 and old[1] == cur[1],
            required=False)
    rules = [p for p in ("AGENTS.md", "CONTRIBUTING.md", "tools/lint.py")
             if p in changed_paths(repo, since)]
    r.check("the rules themselves are unchanged", not rules,
            required=False, detail=", ".join(rules))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
