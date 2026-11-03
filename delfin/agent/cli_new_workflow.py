"""CLI front end for ``delfin-agent new-workflow``.

Mirrors ``cli_new_adapter``: the cli.py hook is minimal (one ``add_parser`` +
one ``set_defaults``); this module owns the argument spec, the handler and
the output so the scaffold logic stays out of the protected-file diff.

``new-workflow`` scaffolds a NEW workflow in the style of ``new-adapter``: a
module registering a :class:`~delfin.tools.workflows.WorkflowDescriptor`
(which extends the existing Application) plus an acceptance test.  The spec
the workflow runs is validated through the SAME existing checks a real
workflow must pass before anything is written, so the command cannot emit a
workflow that would not run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def build_parser(sub) -> argparse.ArgumentParser:
    """Attach the ``new-workflow`` subcommand to an argparse sub-command set."""
    nw = sub.add_parser(
        "new-workflow",
        help="Scaffold a new workflow (descriptor, steps, acceptance test)",
    )
    nw.add_argument("name", help="Workflow module name (plain identifier)")
    nw.add_argument(
        "--kind",
        choices=("pipeline", "classic-control", "user"),
        default="pipeline",
        help="What kind of workflow to scaffold",
    )
    nw.add_argument(
        "--steps",
        metavar="FILE",
        help="A JSON pipeline spec whose steps seed the new workflow "
             "(default: a smiles_to_xyz -> orca_sp skeleton)",
    )
    nw.add_argument(
        "--dir",
        default=".",
        help="Target directory (the repo root, or a user workflow tree)",
    )
    return nw


def _default_steps() -> list:
    return [
        {"step": "smiles_to_xyz", "label": "structure",
         "kwargs": {"smiles": "CCO"}},
        {"step": "orca_sp", "label": "sp", "kwargs": {"charge": 0, "mult": 1}},
    ]


def handler(args: argparse.Namespace) -> int:
    """Run the scaffold; print a short report; return a CLI exit code."""
    from delfin.tools.scaffold import write_workflow_scaffold

    steps = _default_steps()
    if args.steps:
        try:
            data = json.loads(Path(args.steps).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            print(f"ERROR: cannot read --steps {args.steps}: {exc}",
                  file=sys.stderr)
            return 2
        steps = data.get("steps", data) if isinstance(data, dict) else data

    directory = Path(args.dir).resolve()
    plan = write_workflow_scaffold(
        args.name, kind=args.kind, directory=directory, steps=steps,
    )
    if not plan.ok:
        for err in plan.errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 2
    print(f"Created workflow {plan.name!r} ({plan.kind}):")
    print(f"  {plan.module_path}")
    print(f"  {plan.test_path}")
    print(
        "Next steps: fill in the TODO purpose/entry, then extend the steps "
        "with building blocks from the catalog; the descriptor is registered "
        "on import.  See docs/workflow_integration.md."
    )
    return 0
