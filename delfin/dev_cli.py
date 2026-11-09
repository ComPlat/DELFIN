"""Developer CLI (``delfin-dev``): maintenance helpers for DELFIN itself.

Currently one subcommand group:

- ``fixtures`` -- the private fixture store for real ORCA output,
  delegated to :mod:`delfin.fixtures_store` (LICENSE_GUARD.md, rule
  family 3).
"""
from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin-dev",
        description="Maintenance helpers for DELFIN developers.",
    )
    sub = parser.add_subparsers(dest="group", required=True)
    sub.add_parser(
        "fixtures",
        help="private fixture store (init/add/verify/where; see "
        "delfin-dev fixtures --help)",
    )
    sub.add_parser(
        "precheck",
        help="run the project-wide pins a PR would trip in CI "
        "(see delfin-dev precheck --help)",
    )
    args, rest = parser.parse_known_args(argv)
    if args.group == "fixtures":
        from delfin.fixtures_store import main as fixtures_main

        return fixtures_main(rest)
    if args.group == "precheck":
        from delfin.dev_precheck import precheck

        return precheck(rest)
    parser.error(f"unknown group {args.group!r}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
