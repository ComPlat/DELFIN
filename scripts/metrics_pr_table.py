#!/usr/bin/env python3
"""Print the before/after metrics table for a pull request, on the reachability half.

Input: a git repository whose ``origin/main`` is the base, and the checked-out
PR tree. Output: the table on stdout, and either a note that says the table
could not be produced or the reason it is empty; exit 0 when the table was
printed, 2 when the measurement could not run.

Why only the reachability half here
-----------------------------------
The row has three parts. The offered-surface and prompt-size parts are
measured without a model and are printed here. The quality part is NOT,
because it needs model calls: measuring it on every pull request would spend
money per push, and this repository's CI holds no key (``ci.yml`` declares
``permissions: contents: read`` and no secrets). A job that pretended to
report quality would report either a zero or a stale number, and both read
as a measurement. So the table says "not measured" there, in the row's own
words, rather than inventing a value.

What it does print is the part that a pull request can change silently: the
number of capabilities the product no longer reaches, and whether that rose
against the base. That is a defect count, not a sample, so it is compared
without a spread and it is the one number with a never-worse rule.

Why a separate script and not the CLI only
------------------------------------------
The CLI records what it measures into one history. CI has no history between
runs (a fresh runner every time) and no write access to a checkout that is
then graded, so this reads both trees in one process and prints the
comparison. It calls the same module the CLI calls -- one measurement, two
front ends.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def _worktree_of(ref: str, into: Path) -> Path | None:
    """A detached worktree at *ref*, or None when it cannot be made.

    A worktree, not a checkout: the tree being graded must not move, and the
    base commit is a different tree from the one CI just built.
    """
    result = subprocess.run(
        ["git", "worktree", "add", "--detach", str(into), ref],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"could not check out {ref}: {result.stderr.strip()}",
              file=sys.stderr)
        return None
    return into


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="metrics_pr_table",
        description="Print the before/after metrics table for a pull request")
    parser.add_argument("--repo", default=".",
                        help="Repository to measure (default: .)")
    parser.add_argument("--base", default="origin/main",
                        help="Base ref to compare against (default: origin/main)")
    parser.add_argument("--json", action="store_true",
                        help="Print the comparison as JSON")
    args = parser.parse_args(argv)

    repo = Path(args.repo).resolve()
    sys.path.insert(0, str(repo))
    from delfin.agent import commit_metrics as cm   # noqa: E402  (path first)

    # Measured under one instrument state: the SAME code reads both trees, so
    # a change in the measurement cannot masquerade as a change in the code.
    with tempfile.TemporaryDirectory() as tmp:
        base_dir = _worktree_of(args.base, Path(tmp) / "base")
        if base_dir is None:
            print("no table: the base ref could not be checked out, so there "
                  "is nothing to compare against", file=sys.stderr)
            return 2
        try:
            old = cm.build_row(base_dir, commit=args.base)
        finally:
            subprocess.run(["git", "worktree", "remove", "--force",
                            str(base_dir)], capture_output=True, text=True)
    new = cm.build_row(repo, commit=cm.git_commit(repo))
    result = cm.compare_rows(old, new)

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(cm.format_comparison(result))
        print()
        print("quality: not measured -- this job holds no model key, so the "
              "benchmark half is not in this table. Record it locally with "
              "`delfin-agent metrics record --bench-runs <file>`.")
        print()
        # The never-worse rule, spelled out: the reachability count is a
        # counted enumeration, so a rise is a defect and not noise.
        violations = result.get("never_worse_violations") or []
        if violations:
            print(f"REGRESSION: {', '.join(violations)}")
        else:
            print("no reachability regression against the base")
    return 0


if __name__ == "__main__":
    sys.exit(main())
