"""Keep session, for an agent session that runs in a dashboard's kernel.

A session in the agent tab runs in the dashboard's kernel by default, and
the kernel ends with the dashboard (``resume_server``: the window grace,
or the server stopping). With Keep session switched on, the session does
not end there: when its kernel ends, it is resumed in the person's agent
service -- the process ``delfin-agent tui`` and the phone already use,
which lives under no lifeline -- and runs on there, kept. A dashboard
opened later lists it as running; the TUI opens it too.

Input: the switch (``arm`` / ``disarm``) with the session's id, its
working folder and its start settings. Output: one record per kept
session, ``<dir>/<session id>.json`` (owner-only), naming the kernel
process that holds it (pid and start time). ``hand_over`` resumes, in the
agent service, every record whose kernel process is gone, and removes it.

Who calls ``hand_over`` (each claim is an atomic rename, so two callers
never resume one session twice):

* the kernel itself, on its way out (``spawn_hand_over`` from an atexit
  hook and from its lifeline handler): a detached process that waits for
  the kernel's pid to be gone -- the writer lock is released with it --
  and then hands over;
* the dashboard server, after it ended a kernel (``resume_server``): the
  kernel may have been killed after the shutdown wait, with no atexit;
* a dashboard that starts (``agent_sessions.create_tab``): what a killed
  server left behind on this host.

Semantics kept on purpose:

* An emergency stop (``stop_all``) wins: a record armed before the last
  stop is dropped, not resumed -- nothing starts on its own after a stop
  until the person sends something.
* Only this host's records are judged: a pid means nothing on another
  machine sharing the home directory.
* A failed hand-over leaves the record with the reason, for the next
  caller; nothing is lost, the session stays saved in the store.

Rejected: moving the session into the service when the switch is turned
on. The service's view of a session is not the agent tab's full chat
(the take-over in ``agent_sessions.open_from_list`` exists for that
reason), and a switch that changes the chat under the person's hands
surprises. Rejected: keeping the whole kernel for one session -- that is
what "Keep dashboard" is for, and it holds everything else as well.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Optional

#: Overrides the record directory (tests, a second installation).
DIR_ENV = "DELFIN_KEPT_AGENT_SESSIONS_DIR"

#: How long a hand-over waits for the kernel process it was started from.
WAIT_PID_S = 120.0

#: The settings a record may carry into session.create.
_SETTINGS = ("model", "provider", "effort", "mode", "permissionMode")


def record_dir() -> Path:
    raw = os.environ.get(DIR_ENV, "").strip()
    return Path(raw) if raw else Path.home() / ".delfin" / "kept_agent_sessions"


def _host() -> str:
    try:
        return (socket.gethostname() or "").split(".")[0]
    except Exception:
        return ""


def _start_ticks(pid: int) -> Optional[int]:
    try:
        from delfin.agent import proc_identity
        return proc_identity.start_ticks(pid)
    except Exception:
        return None


def _path(session_id: str) -> Path:
    safe = "".join(ch for ch in str(session_id) if ch.isalnum() or ch in "-_")
    if not safe:
        raise ValueError("a kept session needs its id")
    return record_dir() / f"{safe}.json"


def _write(path: Path, payload: dict) -> None:
    directory = path.parent
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        os.chmod(directory, 0o700)
    except OSError:
        pass
    tmp = directory / f".{path.name}.{os.getpid()}.tmp"
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)
    os.replace(tmp, path)


def arm(session_id: str, workspace: str, settings: Optional[dict] = None, *,
        pid: Optional[int] = None, title: str = "") -> Path:
    """Switch Keep session on for ``session_id``, held by this process
    (``pid``: the kernel; default this one). Returns the record's path."""
    pid = int(pid or os.getpid())
    try:
        from delfin.dashboard.session import kernel_id as _kernel_id
        kid = _kernel_id() if pid == os.getpid() else ""
    except Exception:
        kid = ""
    payload = {
        "session_id": str(session_id),
        # The kernel as the dashboard server knows it: the server holds a
        # closed window's kernel for a running turn only when it keeps a
        # session (resume_server).
        "kernel_id": kid,
        "workspace": str(workspace or ""),
        # The name the dashboard showed: the service lists it under that.
        "title": " ".join(str(title or "").split())[:80],
        "settings": {k: v for k, v in (settings or {}).items()
                     if k in _SETTINGS and v not in (None, "")},
        "pid": pid,
        "start_ticks": _start_ticks(pid),
        "host": _host(),
        "armed_at": time.time(),
    }
    path = _path(session_id)
    _write(path, payload)
    return path


def disarm(session_id: str) -> bool:
    """Switch it off: the session ends with its dashboard. True when a
    record was removed."""
    try:
        _path(session_id).unlink()
        return True
    except (OSError, ValueError):
        return False


def is_armed(session_id: str) -> bool:
    try:
        return _path(session_id).is_file()
    except ValueError:
        return False


def records() -> list[dict]:
    """Every record that can be read (unreadable ones are skipped)."""
    out = []
    try:
        entries = sorted(record_dir().glob("*.json"))
    except OSError:
        return out
    for path in entries:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(data, dict) and data.get("session_id"):
            data["_path"] = str(path)
            out.append(data)
    return out


def _alive(pid: int, ticks: Optional[int]) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        pass
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            if fh.read().split(b") ", 1)[1][:1] == b"Z":
                return False              # ended; its parent has not reaped it
    except (OSError, IndexError):
        pass
    if ticks is None:
        return True
    now = _start_ticks(pid)
    return now is None or now == ticks


def holder_gone(record: dict) -> bool:
    """Whether the kernel that holds *record* has ended, as this host can
    tell. A record from another host is never judged here."""
    if str(record.get("host") or "") != _host():
        return False
    try:
        pid = int(record.get("pid") or 0)
    except (TypeError, ValueError):
        return False
    return not _alive(pid, record.get("start_ticks"))


def _stopped_after(armed_at: float) -> bool:
    """An emergency stop was given after the switch was turned on."""
    try:
        from delfin.agent import stop_all
        return not stop_all.wakes_allowed(float(armed_at or 0.0))
    except Exception:
        return False


def _connect() -> Any:
    from delfin.dashboard import live_sessions
    return live_sessions.connect(start=True)


def hand_over(*, wait_pid: int = 0, wait_s: float = WAIT_PID_S,
              connect: Optional[Callable[[], Any]] = None,
              say: Callable[[str], None] = lambda _line: None,
              only: str = "") -> list[dict]:
    """Resume in the agent service every kept session whose kernel ended.

    ``wait_pid``: first wait (up to ``wait_s``) for that process to end --
    the kernel this was started from. ``only``: just that session, now,
    whether or not its kernel still runs -- closed with x and "Keep
    running" its tab has let go of it. Returns one row per record looked
    at: ``session_id`` and ``result`` (``resumed``, ``dropped: <why>``,
    ``failed: <why>``, ``held``: its kernel still runs).
    """
    if wait_pid:
        deadline = time.monotonic() + wait_s
        while time.monotonic() < deadline and _alive(wait_pid, None):
            time.sleep(0.25)
    rows: list[dict] = []
    client = None
    for record in records():
        sid = str(record["session_id"])
        if only and sid != only:
            continue
        if not only and not holder_gone(record):
            rows.append({"session_id": sid, "result": "held"})
            continue
        path = Path(record["_path"])
        claimed = path.with_name(f".{path.name}.claimed-{os.getpid()}")
        try:
            os.rename(path, claimed)
        except OSError:
            continue                      # another caller has it
        if _stopped_after(record.get("armed_at") or 0.0):
            claimed.unlink(missing_ok=True)
            rows.append({"session_id": sid, "result": "dropped: an emergency stop "
                         "was given after it was kept"})
            say(f"[delfin] kept session {sid[:8]} not resumed: an emergency stop "
                "came after it was kept.")
            continue
        try:
            if client is None:
                client = (connect or _connect)()
            running = {str(r.get("sessionId") or "")
                       for r in (client.sessions() or {}).get("sessions") or []}
            if sid not in running:
                settings = dict(record.get("settings") or {})
                if record.get("title"):
                    settings["title"] = str(record["title"])
                client.create(str(record.get("workspace") or ""), resume=sid,
                              keep=True, **settings)
        except Exception as exc:          # noqa: BLE001 -- kept for the next caller
            why = f"{type(exc).__name__}: {getattr(exc, 'message', exc)}"[:300]
            record = {k: v for k, v in record.items() if k != "_path"}
            record["last_error"] = why
            try:
                _write(path, record)
            finally:
                claimed.unlink(missing_ok=True)
            rows.append({"session_id": sid, "result": f"failed: {why}"})
            say(f"[delfin] kept session {sid[:8]} could not be resumed in the "
                f"agent service ({why}); it stays saved and is tried again.")
            continue
        claimed.unlink(missing_ok=True)
        rows.append({"session_id": sid, "result": "resumed"})
        say(f"[delfin] kept session {sid[:8]} runs on in the agent service.")
    if client is not None:
        try:
            client.close()
        except Exception:
            pass
    return rows


def kept_kernel_ids() -> set[str]:
    """The dashboard kernels that hold a kept session and still run."""
    return {str(r.get("kernel_id") or "") for r in records()
            if r.get("kernel_id") and not holder_gone(r)}


def pending_here() -> bool:
    """Whether a record of this host waits for a hand-over (its kernel
    ended)."""
    return any(holder_gone(r) for r in records())


def _child_env() -> dict:
    env = dict(os.environ)
    for name in ("DELFIN_LIFELINE_PID", "DELFIN_LIFELINE_TICKS"):
        env.pop(name, None)
    return env


def spawn_hand_over(*, wait_pid: int = 0, only: str = "") -> Optional[subprocess.Popen]:
    """Run :func:`hand_over` in a process of its own, detached: the caller
    (a kernel on its way out, a server about to stop) does not wait, and
    the hand-over is not ended with it. None when nothing is armed."""
    if not records():
        return None
    argv = [sys.executable, "-m", "delfin.dashboard.kept_agent_sessions",
            "hand-over"]
    if wait_pid:
        argv += ["--wait-pid", str(int(wait_pid))]
    if only:
        argv += ["--only", str(only)]
    log = record_dir() / "hand_over.log"
    try:
        out = open(log, "ab")
        os.chmod(log, 0o600)
    except OSError:
        out = subprocess.DEVNULL
    # Run from the folder this process's own package came from: ``-m``
    # puts the working folder first on the import path, and the kernel's
    # folder is wherever the dashboard was started.
    import delfin
    home = str(Path(delfin.__file__).resolve().parents[1])
    try:
        return subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=out,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                env=_child_env(), close_fds=True, cwd=home)
    except OSError:
        return None
    finally:
        if out is not subprocess.DEVNULL:
            out.close()


def on_kernel_end(*, emergency: bool = False) -> Optional[subprocess.Popen]:
    """This kernel ends: hand its kept sessions over once it is gone. An
    emergency stop hands nothing over (the records are dropped by the
    hand-over that the next caller runs)."""
    if emergency:
        return None
    mine = [r for r in records() if int(r.get("pid") or 0) == os.getpid()]
    if not mine:
        return None
    return spawn_hand_over(wait_pid=os.getpid())


def main(argv: Optional[list[str]] = None) -> int:
    import argparse

    p = argparse.ArgumentParser(prog="python -m delfin.dashboard.kept_agent_sessions")
    sub = p.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("hand-over", help="resume kept sessions whose kernel ended")
    h.add_argument("--wait-pid", type=int, default=0)
    h.add_argument("--only", default="", help="just this session, now")
    args = p.parse_args(argv)
    if args.cmd == "hand-over":
        rows = hand_over(wait_pid=args.wait_pid, only=args.only,
                         say=lambda line: print(line, flush=True))
        print(json.dumps({"at": time.time(), "rows": rows}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
