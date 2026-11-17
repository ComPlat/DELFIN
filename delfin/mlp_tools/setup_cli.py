"""``delfin mlp_setup``: the guided MLIP set-up from a terminal.

Subcommands and what they change:

* ``probe``  -- nothing (read-only inventory).
* ``plan``   -- nothing (proposal, licence notes, commands, settings diff).
* ``result`` -- nothing (reads a verify job directory).
* ``install`` -- runs each proposed command after a y/N answer at the
  terminal, one command per question. Without a terminal to ask (stdin not
  a TTY, e.g. an agent's shell) it refuses; there is no ``--yes``.
* ``configure`` -- prints the settings diff; writes it through
  ``save_settings`` with ``--write`` after a y/N answer (or ``--yes``,
  which the agent's shell gate asks the user about as an acting command).
* ``verify`` -- writes the verify job directory (dry run); ``--submit``
  submits it after a y/N answer (or ``--yes``).
* ``remove`` -- shows the environment removal command and runs it after a
  y/N answer at the terminal; no ``--yes``.

``ask`` is injectable for tests; the default asks on the terminal and
returns None (cannot ask) when stdin is not a TTY.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Callable, Optional

from delfin.mlp_tools import guided_setup as gs


def _tty_ask(question: str) -> Optional[bool]:
    if not (sys.stdin and sys.stdin.isatty()):
        return None
    try:
        answer = input(f"{question} [y/N] ")
    except EOFError:
        return None
    return answer.strip().lower() in ("y", "yes")


def _plan_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--backend", default="", choices=("",) + gs.ORCA_MLIPS_BACKENDS,
                   help="orca-mlips plugin (default: first ungated, commercially usable)")
    p.add_argument("--partition", default="", help="GPU partition(s), comma-separated")
    p.add_argument("--env-name", default=gs.DEFAULT_ENV_NAME)
    p.add_argument("--cuda", default="", help="torch CUDA build, e.g. cu126 or 12.6")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin mlp_setup",
        description="Guided MLIP environment set-up for ORCA ExtOpt via orca-mlips.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("probe", help="read-only inventory of this host")
    p.add_argument("--json", action="store_true")
    p = sub.add_parser("plan", help="proposal with commands and licence notes")
    _plan_args(p)
    p.add_argument("--json", action="store_true")
    p = sub.add_parser("install", help="run the proposed commands, each asked")
    _plan_args(p)
    p = sub.add_parser("configure", help="settings diff; --write to save it")
    _plan_args(p)
    p.add_argument("--write", action="store_true")
    p.add_argument("--yes", action="store_true")
    p = sub.add_parser("verify", help="verify job: dry run, --submit to submit")
    p.add_argument("--backend", default="", choices=("",) + gs.ORCA_MLIPS_BACKENDS)
    p.add_argument("--job-dir", default="mlp_verify")
    p.add_argument("--orca", default="orca")
    p.add_argument("--submit", action="store_true")
    p.add_argument("--yes", action="store_true")
    p = sub.add_parser("result", help="read a verify job's result")
    p.add_argument("--job-dir", default="mlp_verify")
    p.add_argument("--json", action="store_true")
    p = sub.add_parser("remove", help="delete the MLIP environment (asked)")
    p.add_argument("--env-name", default=gs.DEFAULT_ENV_NAME)
    return parser


def _plan(args, probe_fn) -> tuple:
    found = probe_fn()
    plan = gs.propose(found, backend=args.backend, partition=args.partition,
                      env_name=args.env_name, cuda=args.cuda)
    return found, plan


def _print_plan(plan: dict, out) -> None:
    out(f"backend: {plan['backend']} ({plan['backend_reason']})")
    for line in plan["licence"]["lines"]:
        out(f"  {line}")
    out(f"partition: {plan['partition'] or '(none)'} ({plan['partition_reason']})")
    out(f"torch CUDA build: {plan['cuda_build']} ({plan['cuda_reason']})")
    out(f"environment: {plan['env_name']} -> {plan['env_prefix'] or '(prefix known after create)'}")
    out("commands (each asked before it runs):")
    for i, c in enumerate(plan["commands"], 1):
        tag = " [you run this yourself]" if c["kind"] == "manual" else ""
        out(f"  {i}. {c['command']}{tag}")
        out(f"     {c['purpose']}")
    out(f"download: torch {plan['download_estimate']['torch']}; "
        f"weights {plan['download_estimate']['weights']}")
    for b in plan["blockers"]:
        out(f"blocker: {b}")
    if plan.get("remove_command"):
        out(f"remove later with: {plan['remove_command']['command']}")


def main(argv=None, *, ask: Callable[[str], Optional[bool]] = _tty_ask,
         out: Callable[[str], None] = print, probe_fn: Callable[[], dict] = None,
         runner: Optional[Callable] = None, settings_path=None) -> int:
    args = build_parser().parse_args(argv)
    probe_fn = probe_fn or gs.probe

    if args.cmd == "probe":
        found = probe_fn()
        out(json.dumps(found, indent=1, default=str) if args.json
            else "\n".join(found["summary"] + [f"note: {n}" for n in found["notes"]]))
        return 0

    if args.cmd == "plan":
        _found, plan = _plan(args, probe_fn)
        plan["settings_diff"] = gs.configure(plan["settings"], settings_path=settings_path)["diff"]
        if args.json:
            out(json.dumps(plan, indent=1, default=str))
        else:
            _print_plan(plan, out)
            out("settings it would write:")
            out(gs.format_diff(plan["settings_diff"]))
        return 0

    if args.cmd == "install":
        _found, plan = _plan(args, probe_fn)
        _print_plan(plan, out)
        if not plan["commands"]:
            out("nothing to run: " + "; ".join(plan["blockers"]))
            return 1

        def approve(command, purpose):
            return ask(f"run: {command}\n  ({purpose})\n?")

        res = gs.run_install(plan["commands"], approve=approve, runner=runner)
        for s in res["steps"]:
            out(f"  {s['status']}: {s['command']}" + (f" -- {s['note']}" if s.get("note") else ""))
        if res["stopped"] == "no approval channel":
            out("refused: no terminal to ask on. Run `delfin mlp_setup install` in your "
                "own terminal, or run each command above as its own shell call (each is "
                "asked).")
        elif res["ok"]:
            out("done; re-probe with `delfin mlp_setup probe`"
                + (" after `hf auth login`" if plan["licence"]["gated"] else ""))
        return 0 if res["ok"] else 1

    if args.cmd == "configure":
        _found, plan = _plan(args, probe_fn)
        preview = gs.configure(plan["settings"], settings_path=settings_path)
        out(f"settings file: {preview['settings_path']}")
        out(gs.format_diff(preview["diff"]))
        if not args.write or not preview["diff"]:
            return 0
        if not plan["partition"]:
            out("not written: no GPU partition (" + plan["partition_reason"]
                + "); pass --partition")
            return 1
        ok = True if args.yes else ask("write these settings?")
        if not ok:
            out("not written" + (" (no terminal to ask on; pass --yes)" if ok is None else ""))
            return 1
        res = gs.configure(plan["settings"], write=True, settings_path=settings_path)
        out("written" + (f"; not kept: {', '.join(res['not_kept'])}" if res["not_kept"] else ""))
        return 0 if not res["not_kept"] else 1

    if args.cmd == "verify":
        cfg = gs.cluster_config()
        backend = args.backend or (gs.load_mlp_backend() or gs.default_backend()[0])
        prepared = gs.prepare_verify_job(cfg, backend=backend, job_dir=Path(args.job_dir),
                                         orca=args.orca)
        out(f"verify job written (dry run): {prepared['script']}")
        out(f"submit command: {prepared['submit_command']}")
        if not args.submit:
            return 0
        ok = True if args.yes else ask(f"submit to {cfg.gpu_partitions} ({cfg.gres})?")
        if not ok:
            out("not submitted" + (" (no terminal to ask on; pass --yes)" if ok is None else ""))
            return 1
        res = gs.submit_verify_job(prepared, runner=runner)
        if res.get("submitted"):
            out(f"submitted job {res['job_id']}; read it with "
                f"`delfin mlp_setup result --job-dir {prepared['job_dir']}`")
            return 0
        out(res.get("error", "not submitted"))
        for c in res.get("likely_causes") or []:
            out(f"  likely cause [{c['cause']}]: {c['fix']}")
        return 1

    if args.cmd == "result":
        res = gs.read_verify_result(Path(args.job_dir))
        out(json.dumps(res, indent=1, default=str) if args.json else gs.verify_summary(res))
        return 0 if res["status"] == "pass" else 1

    if args.cmd == "remove":
        found = probe_fn()
        cmd = gs.remove_command(found.get("conda") or {}, args.env_name)
        if cmd is None:
            out("no micromamba/mamba/conda found; nothing to remove with")
            return 1
        out(f"{cmd.display}\n  ({cmd.purpose})")
        ok = ask("run it?")
        if not ok:
            out("not run" + (" (no terminal to ask on)" if ok is None else ""))
            return 1
        res = gs.run_install([cmd.to_dict()], approve=lambda c, p: True, runner=runner)
        return 0 if res["ok"] else 1
    return 2
