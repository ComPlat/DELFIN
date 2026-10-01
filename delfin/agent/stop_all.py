"""The emergency stop: every agent of this user, on every login node.

The terminal running delfin-voila is the lifeline of what one dashboard
starts (``lifeline.py``). It reaches no further than its own machine, and
on 2026-09-16 that was not far enough: three agent sessions kept on one
login node could not be stopped from the next login, which had landed on
another one.

What every login node shares is the home directory, and every long-lived
DELFIN process already polls its lifeline every few seconds. A stop is
therefore a file, ``~/.delfin/stop_all.json``. Each process that started
before the stop's time ends itself -- dashboard kernels with the agents,
shells and servers in them, and the daemons -- whichever machine it runs
on. On the machine the stop is given from, the kernels and daemons are also
ended directly, so one that is hung, or runs code from before this module,
does not outlive it there.

After a stop nothing starts on its own. Schedules are disabled, and a
session is not woken -- by a finished job, a scheduled wake-up or another
session's message -- until somebody sends it something by hand after the
stop. Cluster jobs keep running: they are calculations, not agents, and
end with their walltime.
"""

from __future__ import annotations

import json
import os
import signal
import socket
import time
from pathlib import Path
from typing import Optional

_PATH = Path.home() / ".delfin" / "stop_all.json"

# (path, mtime_ns, size) -> record. The stop is read by every watcher every
# few seconds; a stat is what that costs, not a read.
_cache: tuple[tuple[str, int, int], dict] | None = None

# What the direct end on this machine looks for.
#
# The daemons and the command-line agent name themselves on the command
# line. Everything a dashboard started -- its kernels, and the shells, MCP
# servers and sandboxes those started -- carries the dashboard's mark in its
# environment: the lifeline (since 2026-09-16) or the port its server was
# started on (since 2026-03), which also finds what a DELFIN from before
# this module left running. The dashboard server itself carries the mark
# too and is left alone: it serves pages, and ends itself once no kernel
# and no window is left.
_COMMAND_MARKERS = (
    "delfin.agent.scheduler_daemon",
    "delfin.agent.job_monitor",
    "delfin.agent.bug_watcher",
    "bin/delfin-agent",
    "delfin.agent.cli",
)
_DASHBOARD_ENV = ("DELFIN_LIFELINE_PID=", "DELFIN_VOILA_PORT=")
_SERVER_MARKERS = ("jupyter-server", "jupyter_server", "voila", "delfin-voila")

# The same exception, for the register. A protected process is found by the
# kind it registered under, not by its command line, and the two paths have
# to agree: the /proc path spares what serves pages, so the register path
# must spare it too.
#
#   "dashboard server"    serves the pages. Ending it turns "reload the page
#                         to start again" into a refused connection.
#   "dashboard launcher"  the server is its child, started with PDEATHSIG,
#                         so ending the launcher ends the server as surely
#                         as signalling the server itself.
#
# Everything else on the register is an agent and ends: the kernels, the
# daemons, the command-line agent.
_SPARED_KINDS = frozenset({"dashboard server", "dashboard launcher"})

# How long the processes that watch for the stop get to end themselves --
# stopping what they started on the way -- before they are ended directly.
_SETTLE_S = 6.0


def _hostname() -> str:
    try:
        return (socket.gethostname() or "").split(".")[0]
    except Exception:
        return ""


def last_stop() -> dict:
    """The most recent stop, or ``{}`` when there has been none. Never raises."""
    global _cache
    try:
        st = os.stat(_PATH)
    except OSError:
        return {}
    key = (str(_PATH), st.st_mtime_ns, st.st_size)
    if _cache is not None and _cache[0] == key:
        return dict(_cache[1])
    try:
        record = json.loads(_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(record, dict):
        return {}
    _cache = (key, record)
    return dict(record)


def stop_time() -> float:
    """When the last stop was given (epoch seconds), 0.0 for never."""
    try:
        return float(last_stop().get("at") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _boot_time() -> Optional[float]:
    try:
        for line in Path("/proc/stat").read_text().splitlines():
            if line.startswith("btime "):
                return float(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return None


def process_started_at(pid: Optional[int] = None) -> Optional[float]:
    """When a process started, in epoch seconds; None without /proc.

    Accurate to about a second (the boot time is kept in whole seconds), so
    a process started within a second of a stop may be taken as older.
    """
    from .lifeline import _start_ticks

    ticks = _start_ticks(os.getpid() if pid is None else pid)
    boot = _boot_time()
    if ticks is None or boot is None:
        return None
    try:
        hz = os.sysconf("SC_CLK_TCK")
    except (ValueError, OSError):
        hz = 100
    return boot + ticks / float(hz or 100)


_STARTED_AT: Optional[float] = None


def stopped_since_start() -> bool:
    """True when a stop was given after this process started. Never raises."""
    global _STARTED_AT
    at = stop_time()
    if not at:
        return False
    if _STARTED_AT is None:
        _STARTED_AT = process_started_at() or time.time()
    return at > _STARTED_AT


def wakes_allowed(armed_at: float) -> bool:
    """May a session start a turn on its own?

    ``armed_at`` is when somebody last sent that session something by hand
    (0.0 for never). After a stop, a session stays quiet until that happens.
    """
    at = stop_time()
    return not at or float(armed_at or 0.0) > at


def held_note() -> str:
    """What a session says instead of waking itself after a stop."""
    stop = last_stop()
    try:
        when = time.strftime("%H:%M", time.localtime(float(stop.get("at") or 0)))
    except (TypeError, ValueError):
        when = "?"
    where = stop.get("host") or "?"
    return (f"⏸ Not started on its own: an emergency stop at {when} "
            f"(on {where}) ended every agent. Nothing here runs again "
            "until you send something yourself.")


def _disable_schedules(reason: str) -> int:
    """Disable every active schedule. Returns how many were disabled."""
    from .scheduler import Scheduler

    sched = Scheduler()
    count = 0
    with sched._lock:
        sched._absorb_external_changes()
        for entry in sched.list_entries():
            if not entry.disabled:
                sched._disable(entry, reason)
                count += 1
        if count:
            sched._save()
    return count


def _ancestors(pid: int) -> set[int]:
    """``pid`` and every process above it: the shell or wrapper a stop was
    typed into mentions the command, and is not what the stop is for."""
    chain: set[int] = set()
    while pid > 1 and pid not in chain:
        chain.add(pid)
        try:
            stat = Path(f"/proc/{pid}/stat").read_text()
            pid = int(stat[stat.rindex(")") + 2:].split()[1])
        except (OSError, ValueError, IndexError):
            break
    return chain


def _own_agent_processes(before: float) -> list[int]:
    """This user's agent processes on this machine that started before
    ``before``: kernels, daemons, command-line agents and whatever a
    dashboard started. Never this process or one above it."""
    uid = os.getuid()
    mine = _ancestors(os.getpid())
    found: list[int] = []
    try:
        entries = os.listdir("/proc")
    except OSError:
        return found
    for entry in entries:
        if not entry.isdigit():
            continue
        pid = int(entry)
        if pid in mine:
            continue
        try:
            if os.stat(f"/proc/{pid}").st_uid != uid:
                continue
            cmdline = Path(f"/proc/{pid}/cmdline").read_bytes()
            if not any(m.encode() in cmdline for m in _COMMAND_MARKERS):
                if any(m.encode() in cmdline for m in _SERVER_MARKERS):
                    continue
                environ = Path(f"/proc/{pid}/environ").read_bytes()
                if not any(b"\0" + m.encode() in b"\0" + environ
                           for m in _DASHBOARD_ENV):
                    continue
        except OSError:
            continue
        started = process_started_at(pid)
        if started is not None and started < before:
            found.append(pid)
    # A protected agent process (process_guard) cannot be read through
    # /proc; it put itself on a register instead.
    try:
        from . import process_guard as _pg
        for rec in _pg.registered_here():
            pid = int(rec.get("pid") or 0)
            if pid in mine or pid in found or _pg.uid_of(pid) != uid:
                continue
            if str(rec.get("kind") or "").strip() in _SPARED_KINDS:
                continue
            started = process_started_at(pid)
            if started is not None and started < before:
                found.append(pid)
    except Exception:
        pass
    return found


def _end_here(pids: list[int], grace_s: float,
              settle_s: float = _SETTLE_S) -> list[int]:
    """End ``pids``: first give them ``settle_s`` to end on their own, then
    SIGTERM, then SIGKILL after ``grace_s``. Returns the pids signalled.

    Each pid is held to the start time it had when it was found, so a pid
    the system hands to a new process meanwhile is left alone.
    """
    from .lifeline import _start_ticks

    targets = [(pid, _start_ticks(pid)) for pid in pids]
    targets = [(pid, ticks) for pid, ticks in targets if ticks is not None]
    deadline = time.monotonic() + settle_s
    while time.monotonic() < deadline and any(_alive(p, t) for p, t in targets):
        time.sleep(0.2)
    ended: list[tuple[int, int]] = []
    for pid, ticks in targets:
        if not _alive(pid, ticks):
            continue
        try:
            os.kill(pid, signal.SIGTERM)
            ended.append((pid, ticks))
        except OSError:
            continue
    deadline = time.monotonic() + grace_s
    while time.monotonic() < deadline and any(_alive(p, t) for p, t in ended):
        time.sleep(0.1)
    for pid, ticks in ended:
        if _alive(pid, ticks):
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
    return [pid for pid, _ticks in ended]


def _alive(pid: int, ticks: Optional[int]) -> bool:
    """The process that had this pid and start time is still running. An
    exited process its parent has not collected yet is not."""
    from . import lifeline

    if not lifeline._alive(pid, ticks):
        return False
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        return stat[stat.rindex(")") + 2:].split()[0] != "Z"
    except (OSError, ValueError, IndexError):
        return True


def request(reason: str = "", *, end_here: bool = True,
            grace_s: float = 5.0, settle_s: float = _SETTLE_S) -> dict:
    """Give the stop. Returns what was done.

    The record is written first: from that moment every watcher on every
    machine ends its process, whatever happens to the rest of this call.
    """
    at = time.time()
    prior = last_stop()
    record = {
        "at": at,
        "host": _hostname(),
        "pid": os.getpid(),
        "reason": str(reason or "")[:300],
        "count": int(prior.get("count") or 0) + 1,
    }
    _PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = _PATH.with_name(f".{_PATH.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(record, indent=1), encoding="utf-8")
    try:
        os.chmod(tmp, 0o600)
    except OSError:
        pass
    os.replace(tmp, _PATH)

    summary = dict(record)
    try:
        summary["schedules_disabled"] = _disable_schedules(
            f"emergency stop at {time.strftime('%Y-%m-%d %H:%M', time.localtime(at))}"
            f" on {record['host']}")
    except Exception as exc:
        summary["schedules_disabled"] = 0
        summary["schedules_error"] = f"{type(exc).__name__}: {exc}"
    found = _own_agent_processes(at) if end_here else []
    summary["found_here"] = found
    summary["ended_here"] = _end_here(found, grace_s, settle_s) if found else []
    return summary


def _environ(pid: int) -> dict[str, str]:
    raw = Path(f"/proc/{pid}/environ").read_bytes().decode(errors="replace")
    return dict(kv.split("=", 1) for kv in raw.split("\0") if "=" in kv)


def _outlived_their_start() -> list[dict]:
    """This user's processes on this machine whose dashboard or terminal is
    gone: the ones that would be working on their own."""
    from .lifeline import ENV_PID, ENV_TICKS, _alive as _lifeline_alive

    uid = os.getuid()
    out: list[dict] = []
    try:
        entries = os.listdir("/proc")
    except OSError:
        return out
    mine = _ancestors(os.getpid())
    try:
        from . import process_guard as _pg
        registered = {int(r.get("pid") or 0): r for r in _pg.registered_here()
                      if _pg.uid_of(int(r.get("pid") or 0)) == uid}
    except Exception:
        registered = {}
    for entry in entries:
        if not entry.isdigit() or int(entry) in mine:
            continue
        pid = int(entry)
        try:
            if pid in registered:
                # Protected: its environment is not readable; the register
                # carries what the judgement needs.
                rec = registered[pid]
                env = {k: str(v) for k, v in (
                    (ENV_PID, rec.get("lifeline_pid") or ""),
                    (ENV_TICKS, rec.get("lifeline_ticks") or ""),
                    ("DELFIN_VOILA_PORT", rec.get("voila_port") or ""),
                    ("JPY_PARENT_PID", rec.get("jpy_parent_pid") or ""),
                ) if v}
            else:
                if os.stat(f"/proc/{pid}").st_uid != uid:
                    continue
                env = _environ(pid)
            cmdline = (Path(f"/proc/{pid}/cmdline").read_bytes()
                       .replace(b"\0", b" ").decode(errors="replace").strip())
        except OSError:
            continue
        why = ""
        root = env.get(ENV_PID, "")
        if root.isdigit() and int(root) != pid:
            raw_ticks = env.get(ENV_TICKS, "")
            ticks = int(raw_ticks) if raw_ticks.isdigit() else None
            if not _lifeline_alive(int(root), ticks):
                why = f"its lifeline {root} is gone"
        server = env.get("JPY_PARENT_PID", "")
        if (not why and "DELFIN_VOILA_PORT" in env and server.isdigit()
                and int(server) != pid and not os.path.exists(f"/proc/{server}")):
            why = f"its dashboard server {server} is gone"
        if why:
            out.append({"pid": pid, "why": why, "command": cmdline[:160]})
    return out


# A dashboard session sends a heartbeat every minute; one silent for longer
# than this has ended (the presence records only count a session as closed
# after fifteen minutes, which is right for the session list and wrong for
# "is anything still running").
HEARTBEAT_FRESH_S = 180


def status() -> dict:
    """What could be running on its own, without changing anything.

    Sessions, kept sessions, schedules and daemons are read from the home
    directory and so cover every login node; processes only this machine.
    """
    from . import session_presence
    from .scheduler import Scheduler

    now = time.time()
    delfin_dir = _PATH.parent
    result: dict = {"host": _hostname(), "at": now, "last_stop": last_stop()}
    result["open_sessions"] = [
        {"host": str(r.get("host") or "?").split(".")[0],
         "title": str(r.get("title") or "")[:60],
         "seconds_since_heartbeat": int(now - float(r.get("updated_at") or 0))}
        for r in session_presence.open_sessions()]
    here = _hostname()
    kept = []
    for path in sorted((delfin_dir / "kept_sessions").glob("*.json")):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(record, dict):
            host = str(record.get("host") or "").split(".")[0]
            pid = record.get("pid")
            kept.append({"session": record.get("session_name"),
                         "host": host,
                         "since": record.get("started_at"),
                         "here": host == here,
                         "alive_here": host == here and str(pid).isdigit()
                         and os.path.exists(f"/proc/{pid}")})
    result["kept_sessions"] = kept
    try:
        result["active_schedules"] = [
            {"id": e.id, "next_fire_at": e.next_fire_at,
             "prompt": (e.reason or e.prompt)[:80]}
            for e in Scheduler().list_entries() if not e.disabled]
    except Exception as exc:
        result["active_schedules"] = []
        result["schedules_error"] = f"{type(exc).__name__}: {exc}"
    daemons = {}
    for name in ("scheduler_daemon", "job_monitor", "bug_watcher"):
        pid_file = delfin_dir / f"{name}.pid"
        try:
            pid = int(pid_file.read_text().strip() or 0)
        except (OSError, ValueError):
            continue
        try:
            running = f"delfin.agent.{name}".encode() in Path(
                f"/proc/{pid}/cmdline").read_bytes()
        except OSError:
            running = False
        daemons[name] = {"pid": pid, "alive_here": running}
    result["daemon_pid_files"] = daemons
    processes = []
    for pid in _own_agent_processes(now + 1.0):
        try:
            cmdline = (Path(f"/proc/{pid}/cmdline").read_bytes()
                       .replace(b"\0", b" ").decode(errors="replace").strip())
        except OSError:
            continue
        processes.append({"pid": pid, "command": cmdline[:160]})
    result["agent_processes_here"] = processes
    result["outlived_their_start_here"] = _outlived_their_start()
    return result
