"""``delfin-agent lab ...``: the person's side of lab control.

Reading (allowed everywhere): ``devices``, ``read DEV``, ``pending``,
``status``, ``manifest check``.

The person's actions (refused in a process an agent session started, and
refused as command text by ``approval_guard``): ``approve ID``, ``clear
DEV``, ``estop --clear``, ``manifest pin``, ``node keygen NODE``.

``estop`` itself is allowed for everyone, an agent included: stopping is
always permitted, starting again is not.

``optimise`` runs the closed loop on the simulator (``--device``, whose
reading depends on the conditions), writing its record under
``<lab_dir>/loops/``: the measurement decides the next conditions, and any
kernel refusal stops the loop.

Daemons: ``daemon --role estop|deadman|both`` (the systemd units run this;
on a workstation run it in a terminal or with ``--detach``). ``smoke`` runs
the simulator plan end to end in a scratch directory and prints
``smoke lab ok``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys


def _print(obj) -> None:
    print(json.dumps(obj, indent=1, default=str, sort_keys=True))


def _person_only(what: str) -> None:
    from delfin.agent.approval_guard import inside_agent_session
    if inside_agent_session():
        print(f"blocked: `lab {what}` is a person's action; this process was "
              "started by an agent session", file=sys.stderr)
        raise SystemExit(3)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="delfin-agent lab",
                                description="Lab device control (safety "
                                            "kernel, e-stop, approvals)")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("devices", help="Devices, capabilities and state")
    r = sub.add_parser("read", help="One device's readings")
    r.add_argument("device")
    sub.add_parser("pending", help="Approval requests waiting for a person")
    sub.add_parser("status", help="E-stop, latches, locks")
    a = sub.add_parser("approve", help="Approve one hazardous action")
    a.add_argument("request_id")
    c = sub.add_parser("clear", help="Clear a latched device (stays safe)")
    c.add_argument("device")
    e = sub.add_parser("estop", help="Emergency stop (global or --module)")
    e.add_argument("--module", default=None)
    e.add_argument("--reason", default="")
    e.add_argument("--clear", action="store_true",
                   help="Lift the stop (a person only); restarts nothing")
    m = sub.add_parser("manifest", help="check | pin | show")
    m.add_argument("action", choices=("check", "pin", "show"))
    n = sub.add_parser("node", help="keygen NODE (enrolment)")
    n.add_argument("action", choices=("keygen",))
    n.add_argument("node")
    d = sub.add_parser("daemon", help="Dead-man / e-stop enforcement loop")
    d.add_argument("--role", choices=("estop", "deadman", "both"),
                   default="both")
    d.add_argument("--interval", type=float, default=1.0)
    d.add_argument("--once", action="store_true")
    d.add_argument("--detach", action="store_true",
                   help="Workstation fallback: run in the background with a "
                        "pid file in the lab directory")
    o = sub.add_parser(
        "optimise",
        help="Closed loop on the simulator: measure, choose the next "
             "conditions, repeat. Stops on any kernel refusal")
    o.add_argument("--device", default="reactor",
                   help="the device whose reading depends on the conditions")
    o.add_argument("--rounds", type=int, default=15,
                   help="how many runs at most (default 15)")
    o.add_argument("--optimiser", choices=("gp", "grid"), default="gp",
                   help="gp: expected improvement on a Gaussian process "
                        "(default); grid: a plain lattice walk, the control")
    o.add_argument("--step", type=float, default=10.0,
                   help="lattice step for --optimiser grid, in each channel's "
                        "unit (default 10)")
    o.add_argument("--setpoint-channel", default="setpoint")
    o.add_argument("--duration-channel", default="duration")
    o.add_argument("--loop-id", default="", help="record name under loops/")
    o.add_argument("--prepare-heating", action="store_true",
                   help="switch the stirrer and fume hood on first (what a "
                        "reactor interlock above 40 degC requires)")
    sub.add_parser("smoke", help="Simulator plan end to end; prints ok")
    w = sub.add_parser(
        "run", help="Judge a workflow document, and run it (--dry-run judges "
                    "only)")
    w.add_argument("workflow", help="workflow document (.yaml/.json)")
    w.add_argument("--dry-run", action="store_true",
                   help="judge it and write nothing (the default for a plan "
                        "nobody has reviewed is to read it first)")
    w.add_argument("--run-id", default="",
                   help="names the run record (default: the document's id)")
    w.add_argument("--resume", action="store_true",
                   help="continue a run that stopped for an approval, from "
                        "the step it stopped at")
    w.add_argument("--parallel", type=int, default=1, metavar="N",
                   help="run up to N steps at once (default 1: document "
                        "order, one at a time); steps that share a device "
                        "never run together, and `after` is always honoured")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from . import store
    from .safety import Kernel
    if args.cmd == "devices":
        k = Kernel.load()
        if k.manifest is None:
            print(k.load_error, file=sys.stderr)
            return 1
        _print(k.devices())
        return 0
    if args.cmd == "read":
        _print(Kernel.load().read(args.device))
        return 0
    if args.cmd == "pending":
        from . import approvals
        from delfin.agent import irreversible
        # What the person reads before `lab approve ID`: the value is set
        # on real hardware, and setting it back does not undo it.
        _print([{**rec, "cannot_be_undone": irreversible.statement(
                    irreversible.lab_consequence(
                        rec.get("device", "?"), rec.get("channel", "?"),
                        rec.get("value", "?"), rec.get("unit")))}
                for rec in approvals.pending()])
        return 0
    if args.cmd == "status":
        from . import estop
        snap = store.snapshot()
        _print({**estop.status(), "locks": snap["locks"],
                "module_locks": snap["module_locks"],
                "active": snap["active"]})
        return 0
    if args.cmd == "approve":
        _person_only("approve")
        from . import approvals
        try:
            rec = approvals.grant(args.request_id)
        except approvals.ApprovalError as exc:
            print(f"refused: {exc}", file=sys.stderr)
            return 1
        print(f"approved {rec['id']}: {rec['device']}.{rec['channel']} = "
              f"{rec['value']} {rec['unit'] or ''} (single use, session "
              f"{rec['session']})")
        return 0
    if args.cmd == "clear":
        _person_only("clear")
        from . import estop
        _print(estop.clear_latch(args.device))
        return 0
    if args.cmd == "estop":
        from . import estop
        if args.clear:
            _person_only("estop --clear")
            _print(estop.clear(module=args.module))
            return 0
        out = estop.trigger(module=args.module, reason=args.reason or
                            "e-stop from the command line")
        _print(out)
        return 0 if out["ok"] else 2
    if args.cmd == "manifest":
        from . import manifest as _m
        from .safety import pin, pinned
        try:
            man = _m.load(store.manifest_path())
        except _m.ManifestError as exc:
            for prob in exc.problems:
                print(f"problem: {prob}", file=sys.stderr)
            return 1
        if args.action == "check":
            p = pinned()
            print(f"ok: {man.installation} v{man.version}, "
                  f"{len(man.devices)} devices, checksum {man.checksum[:16]}, "
                  f"{'pinned' if p.get('checksum') == man.checksum else 'NOT pinned'}")
            return 0
        if args.action == "show":
            _print(man.raw)
            return 0
        _person_only("manifest pin")
        _print(pin(man))
        return 0
    if args.cmd == "node":
        _person_only("node keygen")
        from . import edge
        edge.provision_key(args.node)
        print(f"key for node {args.node}: {edge.key_path(args.node)} (0600); "
              "copy it to the coordinator over a trusted channel")
        return 0
    if args.cmd == "daemon":
        from . import watchdog
        if args.detach:
            return _detach(args)
        return watchdog.serve(args.role, interval=args.interval,
                              once=args.once)
    if args.cmd == "optimise":
        return _optimise(args)
    if args.cmd == "smoke":
        from . import smoke
        return smoke.main([])
    if args.cmd == "run":
        return _run_command(args)
    return 2



def _load_workflow(path: str) -> dict:
    """The document the CLI takes: JSON or YAML, with the ids a resume needs.

    The unique-id check is here rather than in :func:`delfin.lab.plan.validate`
    because a duplicate id makes a RESUME ambiguous, and the resume reads the
    checkpoint this command writes -- the complaint belongs at the door.
    """
    import yaml
    text = open(path, encoding="utf-8").read()
    wf = json.loads(text) if path.endswith(".json") else yaml.safe_load(text)
    if not isinstance(wf, dict) or not wf.get("steps"):
        raise ValueError(f"{path}: a workflow needs a `steps` list")
    seen = set()
    for st in wf["steps"]:
        sid = st.get("id")
        if not sid or sid in seen:
            raise ValueError(f"{path}: every step needs a unique id "
                             f"({sid!r} is missing or repeated)")
        seen.add(sid)
    return wf


def _run_command(args) -> int:
    """``delfin-agent lab run WF``: judge a workflow, then run it.

    One verb, two halves. The document is judged first -- structurally by
    :mod:`delfin.lab.plan`, and every write by the kernel's own dry run -- and
    a document that does not pass is not run: an advisory check on the path to
    a hotplate is a comment. ``--dry-run`` stops after the judgement.

    The runner is the ONE coordinator. ``--parallel N`` overlaps steps that
    share no device, with ``after`` always honoured; the default is document
    order, one at a time. A hazardous write stops the run at
    ``waiting_approval``: a person runs ``lab approve <id>`` and then the same
    command with ``--resume``.
    """
    from . import plan as _plan
    from . import store
    from .safety import Kernel
    try:
        doc = _load_workflow(args.workflow)
    except (OSError, ValueError) as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 1
    k = Kernel.load()
    if k.manifest is None:
        print(k.load_error, file=sys.stderr)
        return 1
    proposal = _plan.validate(k.manifest, doc, k)
    out = proposal.as_dict()
    out["ran"] = False
    if not proposal.ok:
        _print(out)
        print(f"refused: {len(proposal.problems)} problem(s); nothing was run",
              file=sys.stderr)
        return 1
    if args.dry_run and not args.resume:
        _print(out)
        print("judged ok; nothing was written (drop --dry-run to run it)")
        return 0
    run_id = (args.run_id or str(doc.get("id") or "")).strip() or "run"
    try:
        cp = _plan.run(k.manifest, proposal, run_id, person=store.person(),
                       parallel=max(1, int(args.parallel)))
    except ValueError as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 1
    out.update({"ran": True, "run_id": run_id, "status": cp.get("status"),
                "done": cp.get("done"), "pending": cp.get("pending"),
                "waiting": cp.get("waiting") or {},
                "deviations": cp.get("deviations"),
                "measurements": cp.get("measurements"),
                "parallel": max(1, int(args.parallel))})
    _print(out)
    if cp.get("status") == "waiting_approval":
        # A parallel run can park MORE than one step: `pending` names the one
        # this pass reported last, `waiting` every step still owed a grant.
        # Printing only `pending` would send a person to approve one of two
        # and leave the other waiting on an id nobody was told about.
        waiting = cp.get("waiting") or {}
        reqs = " ".join(sorted({str(a) for a in waiting.values() if a})) \
            or str((cp.get("pending") or {}).get("approval") or "")
        print(f"stopped: step {(cp.get('pending') or {}).get('step')!r} needs a "
              "person. It has NOT run.", file=sys.stderr)
        for req in reqs.split():
            print(f"  delfin-agent lab approve {req}", file=sys.stderr)
        print(f"  then: delfin-agent lab run {args.workflow} --run-id {run_id}"
              " --resume", file=sys.stderr)
        return 0
    return 0 if cp.get("status") == "completed" else 1


def _optimise(args) -> int:
    """``delfin-agent lab optimise``: the closed loop, recorded under loops/.

    This is the entry point a person or a scheduled job reaches; the loop
    itself is :mod:`delfin.lab.loop`, and every write in it goes through the
    safety kernel. The record is written whether the loop finished or was
    stopped, because a refusal is data.
    """
    import time as _time

    from . import surface as _surface
    from .loop import GridSearch, Loop, best_measured, rounds_to
    from .safety import Actor, Kernel

    kernel = Kernel.load()
    if kernel.manifest is None:
        print(kernel.load_error, file=sys.stderr)
        return 1
    if args.loop_id:
        session = f"loop-{args.loop_id}"
    else:
        from . import store as _store
        session = f"loop-{_store.now():.0f}"
    actor = Actor.here(session)
    if args.prepare_heating:
        for dev in ("stirrer", "fume_hood"):
            out = kernel.set(actor, dev, "power", "on")
            if not out.get("ok"):
                _print(out)
                return 2
    loop = Loop(kernel, actor, device=args.device,
                setpoint_channel=args.setpoint_channel,
                duration_channel=args.duration_channel, rounds=args.rounds,
                sleep=_time.sleep)
    lo, hi = loop.envelope()
    if args.optimiser == "grid":
        loop.optimiser = GridSearch(
            _surface.lattice(lo, hi, (args.step, args.step)))
    record = loop.run()
    lattice = _surface.lattice(lo, hi, (args.step, args.step))
    best_lattice, best_point = _surface.best_on(lattice)
    record["surface"] = {"lattice_best": best_lattice,
                         "lattice_best_at": list(best_point),
                         "lattice_points": int(len(lattice)),
                         "step": args.step}
    record["rounds_to_99pct"] = rounds_to(record, 0.99 * best_lattice)
    record["best"] = best_measured(record)
    path = _next_record_path(args.loop_id or session)
    path.write_text(json.dumps(record, indent=1, sort_keys=True, default=str),
                    encoding="utf-8")
    record["record"] = str(path)
    _print(record)
    # A stop is reported as a non-zero exit so a scheduled job notices: a
    # search halted at an interlock must not be read as a finished one.
    return 0 if record.get("stopped") == "budget" else 2


def _next_record_path(name: str):
    """A free path under ``<lab_dir>/loops/``: an earlier record is never
    overwritten, because it is a measurement somebody may still be reading."""
    from . import store as _store
    d = _store.ensure_lab_dir() / "loops"
    d.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = d / f"{name}.json"
    n = 1
    while path.exists():
        n += 1
        path = d / f"{name}-{n}.json"
    return path



def _detach(args) -> int:
    """Workstation fallback without systemd: a double-forked daemon."""
    from . import store, watchdog
    pidfile = store.ensure_lab_dir() / f"daemon-{args.role}.pid"
    try:
        old = int(pidfile.read_text().strip())
        if store.pid_alive(old):
            print(f"already running (pid {old})")
            return 0
    except (OSError, ValueError):
        pass
    if os.fork():
        return 0
    os.setsid()
    if os.fork():
        os._exit(0)
    pidfile.write_text(str(os.getpid()))
    log = open(store.lab_dir() / f"daemon-{args.role}.log", "a")
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    try:
        watchdog.serve(args.role, interval=args.interval, out=log)
    finally:
        os._exit(0)
