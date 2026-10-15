"""``delfin-agent tui`` and ``delfin-agent service status|end|stop|restart|logs``.

Registered by ``cli.build_parser``; the work happens here, so the CLI module
only gains its parser lines.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from . import client, launcher, rundir

#: How long a service that does not answer gets to end after SIGTERM
#: before SIGKILL: its sessions save within the service's own grace.
UNANSWERING_TERM_WAIT_S = 45.0


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
                      ("noIsolate", "no_isolate"), ("fake", "fake"),
                      ("verbose", "verbose")):
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
    takes the terminal: security notices, and an exported provider key put
    away (cli.main's own steps).

    Not the principles check: it costs every launch 80 ms (it reads the
    digest's copy in api_client), and it is made where it protects -- the
    agent service checks at its start, each session's worker at its own,
    and every agent where it is born (AgentEngine.__init__). A window over
    tampered principles opens, and its session says why it cannot start.
    """
    from .. import cli as _cli
    from .. import process_guard as _process_guard

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


def _set_up_window() -> tuple[list[str] | None, str]:
    """The window's command after a first-use setup (tui_setup), or None
    and why not: refused here (an agent session, root), switched off
    (DELFIN_TUI_AUTO_INSTALL=0, agent.tui_auto_install=false), or failed.
    Runs in this process -- the person's own -- before the window starts."""
    from . import tui_setup

    why = tui_setup.refusal()
    if why:
        return None, tui_setup.not_set_up(why)
    if not tui_setup.auto_install_enabled():
        return None, tui_setup.not_set_up(
            f"automatic setup is off ({tui_setup.ENV_OPT_OUT}=0 or agent.tui_auto_install=false)")
    try:
        tui_setup.set_up(out=sys.stderr)
    except tui_setup.SetupFailed as exc:
        return None, f"The terminal window could not be set up: {exc}"
    command = launcher.tui_command()
    if command is None:
        return None, ("The terminal window was set up but is not found where it was put "
                      f"({launcher.tui_root()}); `delfin-agent doctor` reports what is there.")
    return command, ""


def _setup_only() -> int:
    """`delfin-agent tui --setup`: set the window up and exit. Asked for
    explicitly, so the opt-out of the automatic setup does not apply;
    the refusals do."""
    from . import tui_setup

    command = launcher.tui_command()
    if command is not None:
        print(f"The terminal window is ready: {' '.join(command)}")
        return 0
    try:
        tui_setup.set_up(out=sys.stderr, explicit=True)
    except tui_setup.SetupFailed as exc:
        print(f"The terminal window could not be set up: {exc}", file=sys.stderr)
        return 2
    command = launcher.tui_command()
    if command is None:
        print(f"The terminal window was set up but is not found in {launcher.tui_root()}.",
              file=sys.stderr)
        return 2
    print(f"The terminal window is ready: {' '.join(command)}")
    note = tui_setup.native_note(launcher.tui_root())
    if note:
        print(note)
    return 0


def cmd_tui(args: argparse.Namespace) -> int:
    if getattr(args, "setup", False):
        return _setup_only()
    _security_preamble()
    command = launcher.tui_command()
    missing = ""
    if command is None:
        # The setup first, in this process; the service is started
        # either way, as before, and named when there is no window.
        command, missing = _set_up_window()
    try:
        hello = launcher.ensure_service()
    except (launcher.ServiceDidNotStart, rundir.UnsafeRuntimeDir) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    if command is None:
        print(f"{missing}\nThe agent service runs: `delfin-agent service status` "
              f"(pid {hello.get('pid')}).", file=sys.stderr)
        return 2
    env = dict(os.environ)
    env["DELFIN_TUI_SOCKET"] = str(rundir.socket_path())
    # What starts the next service, should this one vanish under the window.
    env["DELFIN_TUI_SERVICE"] = json.dumps(launcher.service_command())
    # Said in the window: where else sessions of this user are.
    env["DELFIN_TUI_NOTES"] = json.dumps(launcher.other_host_notes())
    # Settings first: the workspace is where `delfin-agent tui` was typed.
    env["DELFIN_TUI_START"] = json.dumps(_start_settings(args))
    # Inline instead of the whole terminal (tui.json "screen" says it for good).
    if getattr(args, "inline", False):
        env["DELFIN_TUI_SCREEN"] = "inline"
    if command[-1].endswith("main.tsx"):
        # From source: Bun reads its configuration in the window's folder.
        os.chdir(str(Path(command[-1]).parent.parent))
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


def state_words(summary: dict) -> str:
    """What a session does, in the words `service status` prints: working,
    waiting for you (a question or an approval is open), idle, sleeping
    (saved, no worker; it wakes when opened or sent a message), ended."""
    state = str(summary.get("state") or "")
    if state in ("sleeping", "ended"):
        return state
    if int(summary.get("pending") or 0) > 0 or state == "needs_input":
        return "waiting for you"
    if state in ("running", "starting"):
        return "working"
    return "idle"


def short_workspace(path: str, width: int = 32) -> str:
    """``path`` with the home as ``~``, cut in the middle to ``width``."""
    text = _plain(path)
    home = os.path.expanduser("~").rstrip("/")
    if home and (text == home or text.startswith(home + "/")):
        text = "~" + text[len(home):]
    if width and len(text) > width:
        head = (width - 1) // 3
        text = text[:head] + "…" + text[len(text) - (width - 1 - head):]
    return text


def matching(rows: list[dict], prefix: str) -> list[dict]:
    """The rows whose id is ``prefix`` or starts with it (an exact id wins)."""
    if not prefix:
        return []
    exact = [r for r in rows if r["id"] == prefix]
    return exact or [r for r in rows if str(r["id"]).startswith(prefix)]


def _session_rows(svc: "client.ServiceClient") -> tuple[list[dict], dict]:
    """Every session of the service, one row each -- the running ones, then
    the sleeping ones -- and the service's status (with ``elsewhere``)."""
    status = svc.request("service.status")
    listing = svc.request("sessions.list")
    now = time.time()
    live = list(status.get("sessions") or [])
    ids = {s.get("sessionId") for s in live}
    asleep = [s for s in listing.get("sessions") or []
              if s.get("state") == "sleeping" and s.get("sessionId") not in ids]
    rows = []
    for s in [*live, *asleep]:
        created = float(s.get("createdAt") or now)
        active = float(s.get("lastActivity") or created)
        rows.append({
            "id": _plain(s.get("sessionId")),
            "title": _plain(str(s.get("title") or s.get("firstPrompt") or "(untitled)").split("\n")[0]),
            "workspace": _plain(s.get("workspace")),
            "state": state_words(s),
            "rawState": _plain(s.get("state")),
            "pending": int(s.get("pending") or 0),
            "windows": int(s.get("windows") or 0),
            "model": _plain(s.get("model")),
            "createdAt": created,
            "lastActivity": active,
            "ageS": round(max(0.0, now - created), 1),
            "lastActivityS": round(max(0.0, now - active), 1),
        })
    return rows, {**status, "elsewhere": listing.get("elsewhere") or []}


def cmd_service(args: argparse.Namespace) -> int:
    action = getattr(args, "action", "status") or "status"
    if action == "logs":
        return _logs(int(getattr(args, "lines", 40) or 40))
    try:
        if action == "end":
            return _end(str(getattr(args, "target", "") or "").strip(),
                        all_idle=bool(getattr(args, "all_idle", False)),
                        yes=bool(getattr(args, "yes", False)))
        return _service(action, as_json=bool(getattr(args, "json", False)))
    except client.ServiceError as exc:
        # Refused (this shell is an agent session's own: watch only), or
        # the service went away in between: said, not a traceback.
        print(f"ERROR: the agent service refused: {exc.message}", file=sys.stderr)
        return 1


def _refused_inside_an_agent(action: str) -> bool:
    """True (and said) when this process was started by an agent session:
    over the socket the service answers it watch-only, and stop's signal
    fallback would end the service past that rule. Only the user stops
    the service or ends a session."""
    from ..approval_guard import inside_agent_session
    if not inside_agent_session():
        return False
    what = "ends sessions" if action == "end" else "stops the agent service"
    print(f"ERROR: `delfin-agent service {action}` is refused: this "
          "process was started by an agent session: it may watch, not "
          f"act. Only the user {what}; ask them to run it.",
          file=sys.stderr)
    return True


def _service(action: str, *, as_json: bool = False) -> int:
    if action in ("stop", "restart") and _refused_inside_an_agent(action):
        return 1
    hello, why = client.probe()
    if hello is None:
        pid = rundir.service_pid()
        if as_json and action == "status":
            print(json.dumps({"host": rundir.hostname(), "service": None, "sessions": [],
                              "unanswering": {"pid": pid, "why": launcher.why_in_words(why)} if pid else None,
                              "otherHosts": [{"host": _plain(o.get("host"), 60), "pid": o.get("pid")}
                                             for o in rundir.other_hosts()]}, indent=2))
            return 1 if pid else 0
        if pid:
            print(f"an agent service runs on {rundir.hostname()} (pid {pid}), "
                  f"but {launcher.why_in_words(why)}")
            if action in ("stop", "restart"):
                code = _stop_unanswering(pid)
                return _start_again() if action == "restart" and code == 0 else code
            print("  `delfin-agent service stop` ends it (its sessions save)")
            return 1
        if action == "restart":
            return _start_again()
        print(f"no agent service runs on {rundir.hostname()}")
        for other in rundir.other_hosts():
            print(f"  one runs on {_plain(other.get('host'), 60)} "
                  f"(pid {_plain(other.get('pid'), 12)}): "
                  "log in there to reach its sessions")
        return 0 if action == "status" else 1
    if action in ("stop", "restart"):
        code = _stop(restart=action == "restart")
        return _start_again() if action == "restart" and code == 0 else code
    with client.ServiceClient(client="delfin-agent service") as svc:
        rows, status = _session_rows(svc)
    elsewhere = status.get("elsewhere") or []
    if as_json:
        print(json.dumps({
            "host": _plain(status.get("host")),
            "service": {"pid": status.get("pid"), "host": _plain(status.get("host")),
                        "startedAt": status.get("startedAt"), "windows": status.get("windows", 0),
                        "socket": _plain(status.get("socket"))},
            "sessions": rows,
            "elsewhere": [{k: _plain(e.get(k)) for k in ("key", "sessionId", "title", "workspace", "host")}
                          for e in elsewhere]}, indent=2))
        return 0
    now = time.time()
    running = [r for r in rows if r["state"] != "sleeping"]
    asleep = len(rows) - len(running)
    windows = int(status.get("windows") or 0)
    print(f"agent service on {_plain(status.get('host'))}: pid {status.get('pid')}, "
          f"up {_age(now - float(status.get('startedAt') or now))}, "
          f"{len(running)} session{'s' * (len(running) != 1)}, "
          f"{windows} window{'s' * (windows != 1)}"
          + (f", {asleep} asleep" if asleep else ""))
    if rows:
        print(f"  {'ID':8}  {'TITLE':28}  {'WORKSPACE':32}  {'STATE':15}  "
              f"{'AGE':6}  {'ACTIVE':8}  WINDOWS")
        for r in rows:
            print(f"  {r['id'][:8]:8}  {r['title'][:28]:28}  {short_workspace(r['workspace']):32}  "
                  f"{r['state']:15}  {_age(r['ageS']):6}  {_age(r['lastActivityS']) + ' ago':8}  "
                  f"{r['windows']}")
        if asleep:
            print("  sleeping: saved, without a worker; it wakes when opened or sent a message")
        print("  `delfin-agent service end <id>` ends one (it is saved and can be resumed)")
    if elsewhere:
        print("elsewhere (dashboards, terminals):")
        for e in elsewhere:
            print(f"  {_plain(e.get('key'), 12):12}  {_plain(e.get('title'), 28):28}  "
                  f"{_plain(e.get('workspace'))}")
    return 0


#: How long `service end` waits for an ended session to have saved.
END_WAIT_S = 60.0


def _end(target: str, *, all_idle: bool, yes: bool) -> int:
    """`delfin-agent service end <id>|--all-idle`: each session ends as the
    agent view's ctrl+x twice ends it (``session.shutdown``): a running
    turn is stopped, open questions are refused, the session is saved and
    can be resumed. The service and the other sessions go on. Never a
    signal: a service that does not answer ends nothing here."""
    if _refused_inside_an_agent("end"):
        return 1
    if all_idle == bool(target):
        print("ERROR: name one session (`delfin-agent service end <id>`, a unique "
              "prefix is enough; `service status` lists them) or --all-idle",
              file=sys.stderr)
        return 2
    hello, why = client.probe()
    if hello is None:
        print(f"ERROR: no agent service answers on {rundir.hostname()} "
              f"({launcher.why_in_words(why)}): no session to end", file=sys.stderr)
        return 1
    with client.ServiceClient(client="delfin-agent service") as svc:
        rows, _status = _session_rows(svc)
        if all_idle:
            chosen = [r for r in rows if r["state"] == "idle"]
            if not chosen:
                print("no idle session to end")
                return 0
        else:
            chosen = matching(rows, target)
            if not chosen:
                print(f"ERROR: no session {_plain(target, 40)!r} in the agent service "
                      "(`delfin-agent service status` lists them)", file=sys.stderr)
                return 2
            if len(chosen) > 1:
                print(f"ERROR: {_plain(target, 40)!r} names {len(chosen)} sessions; "
                      "give more of the id:", file=sys.stderr)
                for r in chosen:
                    print(f"  {r['id']}  {r['title'][:40]}", file=sys.stderr)
                return 2
        if not yes:
            if not (sys.stdin.isatty() and sys.stdout.isatty()):
                print("ERROR: not asked: this is no terminal; `--yes` ends without asking",
                      file=sys.stderr)
                return 2
            if not _confirm(chosen):
                return 1
        mark = len(svc.events)
        ended: list[dict] = []
        kept: list[tuple[dict, str]] = []
        for r in chosen:
            params = {"sessionId": r["id"]}
            if all_idle:
                # Idle when listed; one that started to work since is left.
                params["ifIdle"] = True
            try:
                reply = svc.request("session.shutdown", params, timeout=END_WAIT_S)
            except client.ServiceError as exc:
                kept.append((r, exc.message))
                continue
            if isinstance(reply, dict) and reply.get("closing") is False:
                kept.append((r, str(reply.get("busy") or "it is busy")))
                continue
            ended.append(r)
        saving = _wait_removed(svc, {r["id"] for r in ended}, mark, END_WAIT_S)
    for r in ended:
        still = " (still saving)" if r["id"] in saving else ""
        print(f"ended {r['id'][:8]}  {r['title'][:40]}: saved{still}; "
              f"`delfin-agent tui --resume {r['id']}` resumes it")
    for r, why in kept:
        print(f"not ended {r['id'][:8]}  {r['title'][:40]}: {_plain(why, 200)}", file=sys.stderr)
    return 0 if ended or not kept else 1


def _confirm(chosen: list[dict]) -> bool:
    """y/N on the terminal: what is ended, and what that stops."""
    for r in chosen:
        note = {"working": " -- it works now: its turn is stopped",
                "waiting for you": " -- it waits for you: its open question is refused"}.get(r["state"], "")
        print(f"  {r['id'][:8]}  {r['title'][:40]}  ({r['state']}){note}")
    n = len(chosen)
    try:
        answer = input(f"End {'this session' if n == 1 else f'these {n} sessions'}? "
                       "It is saved and can be resumed. [y/N] ")
    except EOFError:
        answer = ""
    if answer.strip().lower() in ("y", "yes"):
        return True
    print("nothing ended")
    return False


def _wait_removed(svc: "client.ServiceClient", ids: set[str], after: int,
                  seconds: float) -> set[str]:
    """Wait until the service says each of ``ids`` is gone (its worker
    ended and saved); the ones it did not say within ``seconds``."""
    left = set(ids)
    deadline = time.monotonic() + seconds
    while left:
        try:
            event, at = svc.wait_for(lambda e: e.get("type") == "session.removed"
                                     and e.get("session") in left, after=after,
                                     timeout=max(0.0, deadline - time.monotonic()))
        except TimeoutError:
            break
        left.discard(event.get("session"))
    return left


def _start_again() -> int:
    try:
        hello = launcher.ensure_service()
    except (launcher.ServiceDidNotStart, rundir.UnsafeRuntimeDir) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(f"the agent service runs again (pid {hello.get('pid')}): windows "
          "reconnect, and its sessions wake when used")
    return 0


def _stop(*, restart: bool = False) -> int:
    pid = rundir.service_pid()
    with client.ServiceClient(client="delfin-agent service") as svc:
        reply = svc.request("service.shutdown", {"restart": True} if restart else {})
    n = int(reply.get("sessions", 0) or 0)
    print(f"stopping the agent service ({n} session{'s' * (n != 1)} end and "
          "are saved; the next service offers them again)...", flush=True)
    # Stopped when the process is gone, not when it stops listening: it
    # stops listening first, and its sessions take a while to save.
    if _gone(pid, 60):
        print("stopped")
        return 0
    print("the service still runs after 60 s; see `delfin-agent service logs`",
          file=sys.stderr)
    return 1


def _stop_unanswering(pid: int) -> int:
    """A service that does not say hello (another version, or stuck) ends
    by signal: SIGTERM first, which ends and saves its sessions when its
    loop still runs; SIGKILL when it does not react. A worker whose
    service is gone ends on its own and saves (a hangup)."""
    import signal

    print(f"ending it: SIGTERM to {pid}...", flush=True)
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        print("stopped")
        return 0
    if _gone(pid, UNANSWERING_TERM_WAIT_S):
        print("stopped")
        return 0
    if rundir.service_pid() != pid:
        # Gone meanwhile, or the pid is no longer this host's service.
        print("stopped")
        return 0
    print("it did not react: SIGKILL (its sessions end and save on their own)",
          flush=True)
    try:
        os.kill(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    if _gone(pid, 10):
        print("stopped")
        return 0
    print(f"pid {pid} still runs; see `delfin-agent service logs`", file=sys.stderr)
    return 1


def _gone(pid: int, seconds: float) -> bool:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if not pid or rundir.service_pid() != pid:
            return True
        time.sleep(0.25)
    return False


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
