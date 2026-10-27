"""The `delfin prcheck` subcommand body: run pre_pr_check, print, exit.

Kept out of cli.py (protected) so the module is owned and tested here;
cli.py only wires the argparse subcommand to `main` (see .gate/pc1
package patch).  Exit codes: 0 READY, 1 BLOCKED -- so a scripted
hand-over check fails loudly when it is blocked.
"""

from __future__ import annotations

import argparse
import sys

from delfin.agent.pre_pr_check import PrePRReport, run


def main(argv: list[str] | None = None) -> int:
    """Entry point for the `delfin prcheck` subcommand. Returns the
    process exit code (0 = READY, 1 = BLOCKED)."""
    parser = argparse.ArgumentParser(
        prog="delfin prcheck",
        description=(
            "The pre-PR self-check: run the tests the tree selects for "
            "this branch, verify reviewer branches are merged, scan the "
            "diff for credential-shaped strings, and print one verdict."
        ),
    )
    parser.add_argument(
        "--base",
        default="origin/main",
        help="The branch the PR will be based on (default: origin/main)",
    )
    parser.add_argument(
        "--partner",
        action="append",
        default=[],
        metavar="BRANCH",
        help="A partner/reviewer branch that must be merged into this "
        "one before hand-over (repeatable)",
    )
    args = parser.parse_args(argv)
    report = run(worktree=".", base=args.base, partners=args.partner)
    print(report)
    return 0 if report.verdict == "READY" else 1


if __name__ == "__main__":  # pragma: no cover - direct invocation
    sys.exit(main())
