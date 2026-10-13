"""``delfin-agent tui`` and ``delfin-agent service status|stop|logs``.

Registered by ``cli.build_parser``; the work happens here, so the CLI module
only gains its parser lines.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

from . import client, launcher, rundir


def _start_settings(args: argparse.Namespace) -> dict:
    """What the TUI's first session starts with: chat's flags, as given."""
    cwd = os.path.abspath(os.path.expanduser(getattr(args, "cwd", "") or os.getcwd()))
    settings = {"cwd": cwd}
    for key, name in (("name", "session_name"), ("model", "model"),
                      ("provider", "provider"), ("effort", "effort"),
                      ("permissionMode", "permission_mode"),
                      ("resume", "session")):
        value = getattr(args, name, "") or ""
        if value:
            settings[key] = str(value)
    for key, name in (("continue", "continue_session"), ("isolate", "isolate"),
                      ("noIsolate", "no_isolate"), ("fake", "fake")):
        if getattr(args, name, False):
            settings[key] = True
    for key, name in (("addDirs", "add_dirs"), ("readDirs", "read_dirs")):
        values = [os.path.abspath(os.path.expanduser(v))
                  for v in (getattr(args, name, None) or [])]
        if values:
            settings[key] = values
    return settings


def _security_preamble() -> None:
    """What `delfin-agent` chat does before a session, said before the TUI
    takes the terminal: the principles check, security notices, and an
    exported provider key put away (cli.main's own steps)."""
    from .. import cli as _cli
    from .. import process_guard as _process_guard

    _cli.guard_principles_or_exit()
    try:
        _cli._show_security_notices()
    except Exception:
        pass
    try:
        exported = _process_guard.exported_provider_keys()
        done = _process_guard.put_exported_keys_away(exported)
        if done:
            print(done, file=sys.stderr)
        elif exported:
            print("Warning: " + _process_guard.exported_key_advice(exported),
                  file=sys.stderr)
    except Exception:
        pass


def cmd_tui(args: argparse.Namespace) -> int:
    _security_preamble()
    try:
        hello = launcher.ensure_service()
    except (launcher.ServiceDidNotStart, rundir.UnsafeRuntimeDir) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    command = launcher.tui_command()
    if command is None:
        print("The terminal UI is not built in this checkout (tui/src/main.tsx "
              "with Bun). The agent service runs: `delfin-agent service "
              f"status` (pid {hello.get('pid')}).", file=sys.stderr)
        return 2
    env = dict(os.environ)
    env["DELFIN_TUI_SOCKET"] = str(rundir.socket_path())
    # Settings first: the workspace is where `delfin-agent tui` was typed.
    env["DELFIN_TUI_START"] = json.dumps(_start_settings(args))
    os.chdir(launcher.tui_root())
    sys.stdout.flush()
    sys.stderr.flush()
    # The TUI takes this process (and its terminal) over: no Python stays
    # behind holding the window.
    os.execvpe(command[0], command, env)
    return 0                                     # pragma: no cover


def _plain(value: object, width: int = 0) -> str:
    """Text from a run file or a session, safe to print: one line, no
    control characters (a forged file must not move the cursor)."""
    text = "".join(ch if ch.isprintable() else "?" for ch in str(value or ""))
    return text[:width] if width else text


def _age(seconds: float) -> str:
    seconds = max(0, int(seconds))
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m"
    return f"{seconds // 3600}h{(seconds % 3600) // 60:02d}m"


def cmd_service(args: argparse.Namespace) -> int:
    action = getattr(args, "action", "status") or "status"
    if action == "logs":
        return _logs(int(getattr(args, "lines", 40) or 40))
    hello = client.reachable()
    if hello is None:
        print(f"no agent service runs on {rundir.hostname()}")
        for other in rundir.other_hosts():
            print(f"  one runs on {_plain(other.get('host'), 60)} "
                  f"(pid {_plain(other.get('pid'), 12)}): "
                  "log in there to reach its sessions")
        return 0 if action == "status" else 1
    if action == "stop":
        return _stop()
    with client.ServiceClient(client="delfin-agent service") as svc:
        status = svc.request("service.status")
        listing = svc.request("sessions.list")
    now = time.time()
    sessions = status.get("sessions") or []
    print(f"agent service on {status.get('host')}: pid {status.get('pid')}, "
          f"up {_age(now - float(status.get('startedAt') or now))}, "
          f"{len(sessions)} session{'s' * (len(sessions) != 1)}, "
          f"{status.get('windows', 0)} window{'s' * (status.get('windows') != 1)}")
    for s in sessions:
        title = _plain(s.get("title") or "(untitled)", 28)
        pending = f", {s['pending']} waiting" if s.get("pending") else ""
        print(f"  {_plain(s.get('sessionId'), 8)}  {title:28}  "
              f"{_plain(s.get('state', '?'), 11):11} {_plain(s.get('workspace'))}{pending}")
    elsewhere = listing.get("elsewhere") or []
    if elsewhere:
        print("elsewhere (dashboards, terminals):")
        for e in elsewhere:
            print(f"  {_plain(e.get('key'), 12):12}  {_plain(e.get('title'), 28):28}  "
                  f"{_plain(e.get('workspace'))}")
    return 0


def _stop() -> int:
    with client.ServiceClient(client="delfin-agent service") as svc:
        reply = svc.request("service.shutdown")
    n = int(reply.get("sessions", 0) or 0)
    print(f"stopping the agent service ({n} session{'s' * (n != 1)} end and "
          "are saved)...", flush=True)
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        if client.reachable(timeout=1.0) is None:
            print("stopped")
            return 0
        time.sleep(0.5)
    print("the service still answers after 60 s; see `delfin-agent service logs`",
          file=sys.stderr)
    return 1


def _logs(lines: int) -> int:
    log = launcher.service_log()
    try:
        text = log.read_text(errors="replace")
    except OSError:
        print(f"no service log yet ({log})")
        return 0
    print(f"== {log}")
    print("\n".join(text.splitlines()[-lines:]))
    workers = sorted(log.parent.glob("worker-*.log"),
                     key=lambda p: p.stat().st_mtime)[-5:]
    if workers:
        print("== the newest session logs:")
        for path in workers:
            print(f"  {path}")
    return 0
