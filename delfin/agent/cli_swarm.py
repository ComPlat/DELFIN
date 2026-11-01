"""The `delfin-agent swarm` subcommand body (package SW1).

Kept out of cli.py so the module is owned and tested here; cli.py only
wires the argparse subcommand to ``main``. Every action is a thin call
into delfin/agent/swarm/ -- no orchestration logic lives here, so the
SW2 dashboard and the SW3 lead call the same library functions.

Read-only actions: plan, status, list. Side-effecting: init (writes the
example spec), start (creates worktrees + branches; the user confirms
with --execute after seeing --dry-run).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    """Entry point for the `delfin-agent swarm` subcommand."""
    parser = argparse.ArgumentParser(
        prog="delfin-agent swarm",
        description=(
            "Build and run agent swarms: packages with disjoint write "
            "scopes, role prompts (builder/reviewer/QS), a dry-run start "
            "and a status view. Every member is a full DELFIN session "
            "with its own worktree and the same permission gate as your "
            "own session."
        ),
    )
    sub = parser.add_subparsers(dest="swarm_cmd", required=True)

    p_init = sub.add_parser(
        "init", help="Write a commented example swarm.yaml to edit")
    p_init.add_argument("--force", action="store_true",
                        help="Overwrite an existing swarm.yaml")

    p_plan = sub.add_parser(
        "plan", help="Plan (or propose) a swarm spec from your tree")
    p_plan.add_argument("goal", help="What the swarm should achieve")
    p_plan.add_argument("--spec", default="swarm.yaml",
                        help="Spec file to validate and plan from")

    p_start = sub.add_parser(
        "start", help="Start the swarm: one worktree + branch per member")
    p_start.add_argument("--spec", default="swarm.yaml",
                         help="Spec file to start from")
    p_start.add_argument("--dry-run", action="store_true",
                         help="Show exactly what would start; start nothing")

    p_status = sub.add_parser(
        "status", help="Show swarm status (JSON with --json)")
    p_status.add_argument("--spec", default="swarm.yaml")
    p_status.add_argument("--json", action="store_true", dest="as_json")

    args = parser.parse_args(argv)

    if args.swarm_cmd == "init":
        return _cmd_init(args)
    if args.swarm_cmd == "plan":
        return _cmd_plan(args)
    if args.swarm_cmd == "start":
        return _cmd_start(args)
    if args.swarm_cmd == "status":
        return _cmd_status(args)
    return 2


def _cmd_init(args: argparse.Namespace) -> int:
    from .swarm.launch import LaunchError, write_example_spec

    try:
        path = write_example_spec(Path("swarm.yaml"),
                                  force=bool(args.force))
    except LaunchError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {path}. Edit it, then:")
    print('  delfin-agent swarm plan "<goal>"   # validate + propose')
    print("  delfin-agent swarm start --dry-run # see what would start")
    return 0


def _cmd_plan(args: argparse.Namespace) -> int:
    from .swarm.launch import LaunchError, plan_launch, render_dry_run
    from .swarm.spec import SpecError, load_spec

    try:
        spec = load_spec(args.spec)
        validate_ok = plan_launch(spec)
    except (SpecError, LaunchError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Goal: {args.goal}")
    print(render_dry_run(validate_ok))
    print("Edit swarm.yaml and run `swarm start --dry-run` to confirm.")
    return 0


def _cmd_start(args: argparse.Namespace) -> int:
    from .swarm.launch import LaunchError, plan_launch, render_dry_run
    from .swarm.launch import start_members
    from .swarm.spec import SpecError, load_spec

    try:
        spec = load_spec(args.spec)
        plan = plan_launch(spec)
    except (SpecError, LaunchError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.dry_run:
        print(render_dry_run(plan))
        return 0
    try:
        result = start_members(spec, execute=True)
    except LaunchError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for m in result["started"]:
        print(f"started {m['name']}: branch {m['branch']} at "
              f"{m['worktree']} (permission_mode {m['permission_mode']})")
    print("Pass the prompt to each member session, spaced; "
          "see the role prompts in delfin/agent/swarm/prompts.py.")
    return 0


def _cmd_status(args: argparse.Namespace) -> int:
    from .swarm.spec import SpecError, load_spec
    from .swarm import status as _status
    from .swarm.status import collect_status

    try:
        spec = load_spec(args.spec)
    except SpecError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    report = collect_status(spec)
    if args.as_json:
        print(_status.to_json(report))
    else:
        print(f"swarm: {report['goal']}")
        for pkg in report["packages"]:
            members = ", ".join(
                f"{m['name']}({m['health']})" for m in pkg["members"])
            print(f"  {pkg['name']}: {pkg['phase']} -- {members}")
        if report["waiting_for_user"]:
            print(f"waiting for you: {len(report['waiting_for_user'])} "
                  f"item(s)")
    problems = _status.validate_report(report)
    if problems:
        print("WARNING: status report failed its schema check: "
              + "; ".join(problems), file=sys.stderr)
        return 1
    return 0
