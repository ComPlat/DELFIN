"""CLI front end for ``delfin-agent new-adapter``.

The cli.py hook is minimal (one ``add_parser`` + one ``set_defaults``);
this module owns the argument spec, the handler and the output so the
scaffold logic stays out of the protected-file diff.  If cli.py cannot be
touched, the same handler is importable and testable from here directly.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def build_parser(sub) -> argparse.ArgumentParser:
    """Attach the ``new-adapter`` subcommand to an argparse sub-command set.

    Kept as a standalone function so cli.py needs exactly two lines:
    ``na = cli_new_adapter.build_parser(sub)`` and
    ``na.set_defaults(func=handler)``.
    """
    na = sub.add_parser(
        "new-adapter",
        help="Scaffold a new step adapter (skeleton, contract test, docs)",
    )
    na.add_argument("name", help="Adapter module name (plain identifier)")
    na.add_argument(
        "--kind",
        choices=("program", "parser", "analysis"),
        default="program",
        help="What kind of step to scaffold",
    )
    na.add_argument(
        "--dir",
        default=".",
        help="Target directory (the repo root, or a user adapter tree)",
    )
    return na


def handler(args: argparse.Namespace) -> int:
    """Run the scaffold; print a short report; return a CLI exit code."""
    from delfin.tools.scaffold import write_scaffold

    directory = Path(args.dir).resolve()
    plan = write_scaffold(args.name, args.kind, directory=directory)
    if not plan.ok:
        for err in plan.errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 2
    print(f"Created adapter {plan.adapter_name!r}:")
    for path in sorted(plan.files):
        print(f"  {path}")
    print(
        "Next steps: implement the TODO in execute(); the registry picks it "
        "up on the next import — no hand edit."
    )
    return 0
