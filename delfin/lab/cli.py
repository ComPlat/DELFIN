"""``delfin-agent lab ...``: the person's side of lab control.

Reading (allowed everywhere): ``devices``, ``read DEV``, ``pending``,
``status``, ``manifest check``.

The person's actions (refused in a process an agent session started, and
refused as command text by ``approval_guard``): ``approve ID``, ``clear
DEV``, ``estop --clear``, ``manifest pin``, ``node keygen NODE``.

``estop`` itself is allowed for everyone, an agent included: stopping is
always permitted, starting again is not.

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
    sub.add_parser("smoke", help="Simulator plan end to end; prints ok")
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
    if args.cmd == "smoke":
        from . import smoke
        return smoke.main([])
    return 2


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
