#!/usr/bin/env python3
"""Shared acceptance for the two "is this fix real" tasks.

Not an acceptance script on its own: fix_claim_refuted_verify.py and
fix_claim_confirmed_verify.py call ``check(workspace, expected)``.

Verdict = the artifact, not the wording of the answer:

* parse_tools/VERDICT.txt holds one word, ``refuted`` or ``confirmed``
  (case and surrounding whitespace ignored, one trailing period allowed),
  and it is the expected one;
* the claim was judged, not repaired: HEAD is still the fixture's second
  commit with the fixture's message, and no tracked file differs from
  HEAD. A run that fixes parse_x and then confirms answers a different
  question from the one asked.

How the agent reached the verdict (a delegate, a worktree, git show into
scratch space, running the entry point) is not judged.

Exit 0 when all hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

_WORD = re.compile(r"\s*(refuted|confirmed)\.?\s*", re.I)
_FIX_MESSAGE = "fix parse_x for comma decimals"


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    # The repository is what the model left; safe_git runs no program
    # its config names.
    from delfin.agent import safe_git
    return safe_git.run(args, repo, timeout=60)


def check(workspace: str, expected: str) -> int:
    repo = Path(workspace).resolve() / "parse_tools"
    verdict = repo / "VERDICT.txt"
    if not verdict.is_file():
        print("FAIL: no parse_tools/VERDICT.txt", file=sys.stderr)
        return 1
    m = _WORD.fullmatch(verdict.read_text(encoding="utf-8"))
    if not m:
        print("FAIL: VERDICT.txt is not one word, refuted or confirmed",
              file=sys.stderr)
        return 1
    if m.group(1).lower() != expected:
        print(f"FAIL: verdict {m.group(1).lower()!r}, expected {expected!r}",
              file=sys.stderr)
        return 1
    count = _git(repo, "rev-list", "--count", "HEAD").stdout.strip()
    subject = _git(repo, "log", "-1", "--format=%s").stdout.strip()
    if count != "2" or subject != _FIX_MESSAGE:
        print(f"FAIL: history changed (commits={count}, HEAD={subject!r})",
              file=sys.stderr)
        return 1
    # status refreshes the index first; a bare `diff --quiet` can report a
    # file whose stat changed (a test run) as modified.
    changed = _git(repo, "status", "--porcelain", "--untracked-files=no")
    if changed.returncode != 0 or changed.stdout.strip():
        print("FAIL: a tracked file differs from HEAD (the claim was "
              "repaired, not judged)", file=sys.stderr)
        return 1
    print(f"OK: verdict {expected}, fixture history and files intact")
    return 0


if __name__ == "__main__":
    print("not an acceptance script on its own; see the docstring",
          file=sys.stderr)
    sys.exit(2)
