"""The `delfin-agent swarm` subcommand body (package SW1).

Kept out of cli.py so the module is owned and tested here; cli.py only
hands its ``swarm`` subparser to :func:`configure_parser`. Every action
is a thin call into delfin/agent/swarm/ -- no orchestration logic lives
here, so the SW2 dashboard and the SW3 lead call the same library
functions.

Read-only actions: plan, status, jobs. Side-effecting: init (writes the
example spec), start (creates worktrees + branches; on a cluster login
node it submits the swarm as a SLURM job after a check and a
confirmation, see delfin/agent/swarm/job.py), stop (ends a swarm job).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


#: What `delfin-agent swarm --help` says above the subcommand list.
_DESCRIPTION = (
    "Build and run agent swarms: packages with disjoint write "
    "scopes, role prompts (builder/reviewer/QS), a dry-run start "
    "and a status view. Every member is a full DELFIN session "
    "with its own worktree and the same permission gate as your "
    "own session."
)


def main(argv: list[str] | None = None) -> int:
    """Entry point for the `delfin-agent swarm` subcommand, given its own
    argv (the words after ``swarm``)."""
    parser = configure_parser(
        argparse.ArgumentParser(prog="delfin-agent swarm"))
    return run_args(parser.parse_args(argv))


def configure_parser(
        parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Give *parser* the swarm subcommands and route it to :func:`run_args`.

    cli.py calls this on its ``swarm`` subparser, so `delfin-agent swarm
    status --json` and `delfin-agent swarm --help` are parsed by argparse
    itself, as git, gh and kubectl hand a subcommand its own arguments.
    It used to register a bare subparser whose handler was this module's
    ``main``: the handler then received the parsed Namespace instead of an
    argv list (TypeError), and the subparser, which declared no arguments,
    rejected every subcommand ("unrecognized arguments: status").
    """
    parser.description = _DESCRIPTION
    # Not required: a bare `delfin-agent swarm` prints this help, as
    # `gh pr` and `kubectl config` do, instead of a usage error.
    sub = parser.add_subparsers(dest="swarm_cmd", metavar="COMMAND")

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
        "start",
        help="Start the swarm (on a cluster login node: as a SLURM job)",
        description=(
            "Start the swarm a spec describes. On a cluster login node the "
            "swarm is submitted as one SLURM job (at most "
            "8 member sessions per job): it checks the spec, says what "
            "will happen, asks, submits, and tells you how to watch and "
            "stop it. With --here (or hpc.swarm_on_login_node) the members "
            "run on the login node instead, within its session cap, and "
            "send every test and build to SLURM. Inside a job, on a "
            "workstation or a private server it creates one worktree + "
            "branch per member, as before."))
    p_start.add_argument("spec_file", nargs="?", default="",
                         metavar="SPEC",
                         help="Spec file to start from (default: "
                              "swarm.yaml, or --spec)")
    p_start.add_argument("--spec", default="swarm.yaml",
                         help="Spec file to start from")
    p_start.add_argument("--dry-run", action="store_true",
                         help="Show exactly what would start; start nothing")
    p_start.add_argument("--job", action="store_true",
                         help="Submit it as a SLURM job of its own even "
                              "inside a job (a login node does by default)")
    p_start.add_argument("--here", action="store_true",
                         help="On a cluster login node: run the members "
                              "here instead of as a job (or set "
                              "hpc.swarm_on_login_node). Their tests, builds "
                              "and calculations still go to SLURM.")
    p_start.add_argument("-y", "--yes", action="store_true",
                         help="Submit without asking (needed when no "
                              "terminal can answer)")

    sub.add_parser(
        "jobs", help="Swarm jobs: state, node, time left, how to attach")
    p_stop = sub.add_parser(
        "stop", help="Stop a swarm job (its member sessions end)")
    p_stop.add_argument("job", nargs="?", default="",
                        help="The job id (default: the only one running), "
                             "or 'all'")

    p_status = sub.add_parser(
        "status", help="Show swarm status (JSON with --json)")
    p_status.add_argument("--spec", default="swarm.yaml")
    p_status.add_argument("--json", action="store_true", dest="as_json")

    p_triage = sub.add_parser(
        "triage", help="The approval queue for the swarm, sorted")
    p_triage.add_argument("--queue", default="",
                          help="JSON file with the queued items (default: "
                               "read from stdin)")

    p_merge = sub.add_parser(
        "merge", help="Pre-PR check + print the merge command")
    p_merge.add_argument("package", help="The package to merge")

    p_verdict = sub.add_parser(
        "verdict",
        help="Record a real-run verdict (ok|fail). It merges nothing",
        description=("Records a real-run verdict for one member of a wave, "
                     "with the person who recorded it, in the run's event "
                     "log. A 'fail' verdict resumes that member with the "
                     "failure text. Recording a verdict NEVER merges a pull "
                     "request -- that is `swarm merge-run`, which needs "
                     "--confirmed and the head that was checked."))
    p_verdict.add_argument("run", help="The wave run name")
    p_verdict.add_argument("member", help="The member id")
    p_verdict.add_argument("verdict", choices=("ok", "fail"),
                           help="What the real run showed")
    p_verdict.add_argument("--by", required=True,
                           help="Who recorded this verdict (a person)")
    p_verdict.add_argument("--detail", default="",
                           help="The failure text (given to the member)")

    p_merge_run = sub.add_parser(
        "merge-run",
        help="Merge PR <n>, only with --confirmed and the head that was checked",
        description=("The one merge path in this CLI. It refuses without "
                     "--confirmed, refuses when the pull request's head has "
                     "moved away from --head (a check answered for a "
                     "different commit does not count), and refuses "
                     "--admin-override unless you pass it in this command. "
                     "The command it prints is what runs; this verb runs "
                     "nothing."))
    p_merge_run.add_argument("pr", help="The pull request number")
    p_merge_run.add_argument("--head", required=True,
                             help="The head commit the checks were read for")
    p_merge_run.add_argument("--confirmed", action="store_true",
                             help="A person's yes for this merge")
    p_merge_run.add_argument("--checks-ok", action="store_true",
                             help="The checks on --head passed")
    p_merge_run.add_argument("--admin-override", action="store_true",
                             help="Use administrator rights when the checks "
                                  "are red: the printed command carries "
                                  "--admin")

    p_wave = sub.add_parser(
        "wave", help="Waves: a task list of agents, run, watched and steered",
        description=("A wave is a task list; each task is one caged agent session in its "
                     "own worktree, started when what it waits for has ended "
                     "(delfin/agent/swarm/wave.py). The tick keeps it moving; `mcp` "
                     "serves every verb to an orchestrator over MCP."))
    wsub = p_wave.add_subparsers(dest="wave_cmd", metavar="VERB")
    w = wsub.add_parser("start", help="Start a run from a tasks JSON file")
    w.add_argument("run")
    w.add_argument("tasks", help="JSON file with the task list")
    w.add_argument("--repo", default=".", help="The git repository (default: here)")
    w.add_argument("--parallel", type=int, default=0)
    w.add_argument("--no-stages", action="store_true",
                   help="No real runs and resumes after a member ends")
    w = wsub.add_parser("compose", help="Write a wave's task list from a packages file")
    w.add_argument("run")
    w.add_argument("spec", help="JSON: {packages, arms, integrator, lander, reviewer, fixer, base?}")
    w.add_argument("--out", default="", help="Write the tasks here (default: print them)")
    w = wsub.add_parser("add", help="Append tasks from a JSON file to a running wave")
    w.add_argument("run")
    w.add_argument("tasks")
    w = wsub.add_parser("status", help="One row per task of a run")
    w.add_argument("run")
    wsub.add_parser("runs", help="Every run with its members by state")
    w = wsub.add_parser("watch", help="Wait for the next events of the runs")
    w.add_argument("runs", nargs="*")
    w.add_argument("--timeout", type=float, default=1800)
    w = wsub.add_parser("review", help="A member's write-up, commits and real runs")
    w.add_argument("run")
    w.add_argument("member")
    w = wsub.add_parser("resume", help="Continue an ended member with a note")
    w.add_argument("run")
    w.add_argument("member")
    w.add_argument("--note", default="")
    w.add_argument("--fresh", action="store_true")
    w.add_argument("--after-stop", action="store_true",
                   help="Start a member a person stopped: lifting a stop is a "
                        "person's decision, and this is how it is made")
    w = wsub.add_parser("stop", help="End a run's members (or one)")
    w.add_argument("run")
    w.add_argument("member", nargs="?", default="")
    w = wsub.add_parser("tick", help="One pass: start, advance, record, guard")
    w.add_argument("--every", type=float, default=0, metavar="SECONDS",
                   help="Tick forever at this interval (default: one pass)")
    w = wsub.add_parser("emergency-stop", help="Stop every member of every run and hold")
    w.add_argument("why")
    wsub.add_parser("clear-stop", help="Lift an emergency stop (a person's decision)")
    w = wsub.add_parser("friction", help="What the members ran into, counted")
    w.add_argument("runs", nargs="*")
    w = wsub.add_parser("guard", help="Install or remove the tick (scheduler or systemd)")
    w.add_argument("action", choices=("install", "remove"))
    w.add_argument("--every", type=int, default=1, help="Minutes between ticks")
    w.add_argument("--systemd", action="store_true",
                   help="As a systemd user service that survives a reboot")
    w = wsub.add_parser("land", help="Check, review and (with --push) open members' PRs",
                        description="delfin/agent/swarm/land.py: --calibrate [MODEL ...] and "
                                    "--calibrate-pair [A B] measure the reviewers first.")
    w.add_argument("runs", nargs="*", help="Runs to land (default: every run)")
    w.add_argument("--push", action="store_true",
                   help="Push the ready branches and open their PRs (default: a dry run)")
    w.add_argument("--calibrate", nargs="*", metavar="MODEL", default=None,
                   help="Measure reviewers on the corpus (default: the KIT roster)")
    w.add_argument("--calibrate-pair", nargs="*", metavar="MODEL", default=None,
                   help="Score reviewer pairs from stored rows (default: every pair)")
    wsub.add_parser("mcp", help="Serve the wave verbs over MCP on stdio")

    parser.set_defaults(func=run_args, swarm_help=parser.format_help)
    return parser


def run_args(args: argparse.Namespace) -> int:
    """Run the swarm subcommand *args* names (parsed by a parser that
    :func:`configure_parser` set up)."""
    if not getattr(args, "swarm_cmd", None):
        print(args.swarm_help(), end="")
        return 0
    if args.swarm_cmd == "init":
        return _cmd_init(args)
    if args.swarm_cmd == "plan":
        return _cmd_plan(args)
    if args.swarm_cmd == "start":
        return _cmd_start(args)
    if args.swarm_cmd == "status":
        return _cmd_status(args)
    if args.swarm_cmd == "jobs":
        return _cmd_jobs(args)
    if args.swarm_cmd == "stop":
        return _cmd_stop(args)
    if args.swarm_cmd == "triage":
        return _cmd_triage(args)
    if args.swarm_cmd == "merge":
        return _cmd_merge(args)
    if args.swarm_cmd == "verdict":
        return _cmd_verdict(args)
    if args.swarm_cmd == "merge-run":
        return _cmd_merge_run(args)
    if args.swarm_cmd == "wave":
        return _cmd_wave(args)
    return 2


def _cmd_wave(args: argparse.Namespace) -> int:
    """`delfin-agent swarm wave VERB`: each verb is the MCP handler of the
    same name, so the terminal and an orchestrator get the same answer."""
    from .swarm import mcp_server as M
    verb = getattr(args, "wave_cmd", None)
    if verb == "mcp":
        M.run_server([])
        return 0
    if verb == "land":
        from .swarm import land
        if args.calibrate_pair is not None:
            return land.main(["--calibrate-pair", *args.calibrate_pair])
        if args.calibrate is not None:
            return land.main(["--calibrate", *args.calibrate])
        return land.main([*args.runs, *(["--push"] if args.push else [])])
    if verb == "start":
        try:
            tasks = json.loads(Path(args.tasks).read_text())
        except (OSError, ValueError) as exc:
            print(f"cannot read the task list: {exc}")
            return 2
        out = M.h_start(args.run, tasks, str(Path(args.repo).resolve()),
                        not args.no_stages, args.parallel)
    elif verb == "compose":
        try:
            spec = json.loads(Path(args.spec).read_text())
        except (OSError, ValueError) as exc:
            print(f"cannot read the wave spec: {exc}")
            return 2
        out = M.h_compose(args.run, spec.get("packages"), spec.get("arms") or [],
                          spec.get("integrator") or {}, spec.get("lander") or {},
                          spec.get("reviewer") or {}, spec.get("fixer") or {},
                          spec.get("base") or "origin/main", spec.get("stages"),
                          spec.get("arm_refs"))
        if args.out and '"error"' not in out[:40]:
            Path(args.out).write_text(json.dumps(json.loads(out)["tasks"], indent=1))
            out = json.dumps({"written": args.out, "tasks": len(json.loads(out)["tasks"])})
    elif verb == "add":
        try:
            more = json.loads(Path(args.tasks).read_text())
        except (OSError, ValueError) as exc:
            print(f"cannot read the task list: {exc}")
            return 2
        out = M.h_add_tasks(args.run, more)
    elif verb == "status":
        out = M.h_status(args.run)
    elif verb == "runs":
        out = M.h_runs()
    elif verb == "watch":
        out = M.h_watch(args.runs, args.timeout)
    elif verb == "review":
        out = M.h_review(args.run, args.member)
    elif verb == "resume":
        out = M.h_resume(args.run, args.member, args.note, args.fresh, args.after_stop)
    elif verb == "stop":
        out = M.h_stop(args.run, args.member)
    elif verb == "tick":
        if args.every:
            from delfin.agent.swarm import guard
            guard.loop(max(10.0, args.every))
            return 0
        out = M.h_tick()
    elif verb == "emergency-stop":
        out = M.h_emergency_stop(args.why)
    elif verb == "clear-stop":
        out = M.h_clear_stop()
    elif verb == "friction":
        out = M.h_friction(args.runs or None)
    elif verb == "guard":
        out = (M.h_guard_install(args.every, args.systemd) if args.action == "install"
               else M.h_guard_remove())
    else:
        print("delfin-agent swarm wave VERB -- see --help")
        return 2
    print(out)
    return 1 if '"error"' in out[:40] else 0


def _no_spec(path: str) -> int:
    """A missing spec is said, with the command that writes one -- not a
    FileNotFoundError traceback."""
    print(f"ERROR: no swarm spec at {path}. `delfin-agent swarm init` "
          "writes an example to edit (or pass --spec FILE).", file=sys.stderr)
    return 1


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
    except FileNotFoundError:
        return _no_spec(args.spec)
    except (SpecError, LaunchError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Goal: {args.goal}")
    print(render_dry_run(validate_ok))
    print("Edit swarm.yaml and run `swarm start --dry-run` to confirm.")
    return 0


def _cmd_start(args: argparse.Namespace) -> int:
    from .swarm import job as _job
    from .swarm.launch import LaunchError, plan_launch, render_dry_run
    from .swarm.launch import start_members
    from .swarm.spec import SpecError, load_spec

    spec_path = Path(args.spec_file or args.spec)
    try:
        spec = load_spec(spec_path)
        plan = plan_launch(spec)
    except FileNotFoundError:
        print(f"ERROR: no spec at {spec_path}. Write an example to edit "
              "with `delfin-agent swarm init`.", file=sys.stderr)
        return 1
    except (SpecError, LaunchError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if args.here and args.job:
        print("ERROR: --here and --job contradict each other; give one.",
              file=sys.stderr)
        return 2
    where = _job.route(prefer_job=bool(args.job), here=bool(args.here))
    if where == "job":
        return _start_as_job(args, spec_path)
    if _job.host_kind() == "login_node":
        why = _job.login_node_room(len(plan["members"]))
        if why:
            print(f"ERROR: {why}", file=sys.stderr)
            return 1
        print("Members run on this login node. Every test, build and "
              "calculation they start goes to SLURM (enforced on a login "
              "node).")
    if args.dry_run:
        print(render_dry_run(plan))
        print("Run the same command without --dry-run to start.")
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


def _start_as_job(args: argparse.Namespace, spec_path: Path) -> int:
    """Check, explain, ask, submit: a swarm as a SLURM job."""
    from .swarm import job as _job

    text = spec_path.read_text(encoding="utf-8")
    pre = _job.preflight(text, spec_path, route_="job",
                         running=_job.jobs())
    print(_job.confirmation_text(pre))
    if not pre["ok"]:
        return 1
    if args.dry_run:
        print("Nothing submitted (dry run). Run the same command without "
              "--dry-run to submit.")
        return 0
    if not args.yes:
        if not sys.stdin.isatty():
            print("Not submitted: no terminal to confirm. Pass --yes to "
                  "submit without asking.", file=sys.stderr)
            return 2
        try:
            answer = input("Submit this swarm job? [y/N] ").strip().lower()
        except EOFError:
            answer = ""
        if answer not in ("y", "yes"):
            print("Not submitted.")
            return 1
    try:
        rec = _job.submit(text, pre)
    except _job.SwarmJobError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    res = rec["resources"]
    print(f"Submitted swarm job {rec['job']} ({res['partition']}, "
          f"{res['cpus']} cores, {res['mem']}, {res['time']}).")
    rows = {str(r["job"]): r for r in _job.jobs()}
    now = rows.get(rec["job"])
    if now is not None and now.get("state") not in (None, "UNKNOWN"):
        print("  " + _job.status_line(now))
    else:
        print("  It waits for a free compute node; `delfin-agent swarm "
              "jobs` shows the node once it runs.")
    print(f"  Job log: {rec['log']}")
    for line in _job.watch_lines(rec["job"]):
        print("  " + line)
    return 0


def _cmd_jobs(_args: argparse.Namespace) -> int:
    from .swarm import job as _job

    print(_job.format_jobs(_job.jobs()))
    return 0


def _cmd_stop(args: argparse.Namespace) -> int:
    from .swarm import job as _job

    rows = _job.jobs()
    if args.job == "all":
        ids = [str(r["job"]) for r in rows]
    elif args.job:
        ids = [args.job]
    elif len(rows) == 1:
        ids = [str(rows[0]["job"])]
    else:
        print(_job.format_jobs(rows))
        if rows:
            print("Several swarm jobs run -- name one: delfin-agent swarm "
                  "stop JOB  (or: stop all)")
        return 1
    if not ids:
        print("No swarm job runs.")
        return 1
    print(_job.cancel(ids))
    return 0


def _cmd_status(args: argparse.Namespace) -> int:
    from .swarm.spec import SpecError, load_spec
    from .swarm import status as _status
    from .swarm.status import collect_status

    try:
        spec = load_spec(args.spec)
    except FileNotFoundError:
        return _no_spec(args.spec)
    except SpecError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    from .swarm import job as _job
    # The swarm jobs started from here, asked once: a RUNNING job of this
    # spec is where its members' state comes from (its progress file and
    # their session records), and the list is printed below.
    rows = _job.jobs() if _job.host_kind() != "workstation" else None
    report = collect_status(spec, jobs=rows)
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
                  f"item(s): " + "; ".join(
                      f"{w['member']} waits on {w['waits_on']}"
                      for w in report["waiting_for_user"]))
            print("  " + _job.approvals_line())
        if rows:
            print(_job.format_jobs(rows))
    problems = _status.validate_report(report)
    if problems:
        print("WARNING: status report failed its schema check: "
              + "; ".join(problems), file=sys.stderr)
        return 1
    return 0


def _cmd_triage(args: argparse.Namespace) -> int:
    from .swarm.triage import TriageError, classify, triage_list

    if args.queue:
        queue = json.loads(Path(args.queue).read_text(encoding="utf-8"))
    else:
        queue = json.load(sys.stdin)
    try:
        listed = triage_list(queue)
    except TriageError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for item in listed:
        print(f"[{item['action']}] {item.get('id')} {item.get('member')} "
              f"({item.get('kind')}): {item.get('reason')}")
    n_auto = sum(1 for i in listed if i["action"] == "auto_answer")
    print(f"{n_auto} auto-answerable, {len(listed) - n_auto} for you")
    return 0


def _cmd_merge(args: argparse.Namespace) -> int:
    from .swarm.merge import merge_plan
    from .swarm.triage import TriageError

    try:
        plan = merge_plan(args.package)
    except TriageError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"package {plan['package']}: pre-PR check:")
    print(f"  {plan['pre_pr_check']}")
    print("then, if it is READY (this is YOUR action, not the agent's):")
    print(f"  {plan['checkout_command']}")
    print(f"  {plan['merge_command']}")
    return 0


def _cmd_verdict(args: argparse.Namespace) -> int:
    """`delfin-agent swarm verdict RUN MEMBER ok|fail --by WHO`.

    Records the verdict and, for a fail, resumes the member. Merges
    nothing: the printed line says so, because the operator script this
    replaces merged on exactly such a line.
    """
    from .swarm import lead as _lead
    from .swarm.triage import TriageError
    from .swarm.wave import WaveError

    try:
        got = _lead.swarm_verdict(args.run, args.member, args.verdict,
                                  by=args.by, detail=args.detail)
    except (TriageError, WaveError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"recorded: {got['verdict']} for {got['run']}/{got['id']} "
          f"by {got['by']}")
    if got["verdict"] == "fail":
        state = "resumed" if got.get("resumed") else (
            f"not resumed ({got.get('why') or 'unknown'})")
        print(f"  member {state}")
    print("  merges nothing: a merge needs `swarm merge-run --confirmed` "
          "with the head that was checked")
    return 0


def _cmd_merge_run(args: argparse.Namespace) -> int:
    """`delfin-agent swarm merge-run PR --head SHA [--confirmed]`.

    Prints the merge command only when a person confirmed it and the head
    is the one whose checks were read. It runs nothing itself.
    """
    from .swarm import lead as _lead
    from .swarm.triage import TriageError

    try:
        decided = _lead.swarm_merge(
            args.pr, args.head, confirmed=args.confirmed,
            checks_ok=args.checks_ok,
            admin_override=args.admin_override,
            where=str(Path.cwd()))
    except TriageError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"PR {decided['pr']} at head {decided['head']}: {decided['state']}")
    if decided.get("reason"):
        print(f"  {decided['reason']}")
    if decided.get("ask_user"):
        print(f"  {decided['ask_user']}")
    if decided["allow"]:
        print("run this yourself:")
        print(f"  {decided['merge_command']}")
        if decided["admin_override"]:
            print("  (administrator rights, as you asked in this command)")
    return 0 if decided["allow"] else 1
