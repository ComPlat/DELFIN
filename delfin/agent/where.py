"""Where is the dashboard running, and how do I get back to it?

A dashboard on a cluster runs on ONE login node, inside a terminal
multiplexer, behind a port the user forwards. Every one of those three
is easy to lose: the next login lands on a different node, `tmux ls`
says "no server running", and the session that is still working is
simply somewhere else. Reconnecting then means guessing, and the guess
usually starts with "which node was it again?".

Nothing about that needed to be guessed: the home directory is shared
across the login nodes, so a note left in it is readable from all of
them. The sessions already leave one (``session_presence`` records the
host, the pid and the branch of every open session). The dashboard
itself left none — this adds it, and reads both back:

    $ delfin-agent where
    dashboard  uc3n991  port 8867  tmux 'delfin'  started 20:46  (running)
      come back with:  ssh uc3n991 && tmux attach -t delfin
    session    uc3n991  'uc3n991-79b8'  branch session/33c1925d  seen 12s ago

Read-only except for the one record the launcher writes about itself.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

#: Where the running dashboard leaves its note, when nothing redirects it.
#: One file: a machine runs one dashboard per user, and a second one would
#: overwrite the first — which is the truth worth showing, not two stale
#: halves.
DEFAULT_RECORD_NAME = ".delfin/dashboard_here.json"

#: A record older than this is reported as "last seen", never as running:
#: on another node there is no pid to ask.
STALE_S = 24 * 3600


def record_path() -> Path:
    """The note's path, decided now rather than at import.

    A module-level constant freezes ``$HOME`` the moment the module is
    first imported, which is earlier than anything that might redirect
    it. That is not a hypothetical: seven tests that drive the launcher
    wrote into the real note, and the launcher's ``atexit`` hook then
    deleted it on the way out — so an ordinary suite run left a running
    dashboard with no way back, which is the one thing this file exists
    to prevent. Being resolved per call puts it in reach of the state
    table that redirects every other user-wide sink
    (``state_paths.USER_STATE_RESOLVERS``).
    """
    return Path.home() / DEFAULT_RECORD_NAME


def _process_start(pid: int) -> str:
    """A fingerprint that tells one life of a pid from the next, or "".

    Linux hands it over in /proc; everywhere else ``ps`` answers the same
    question, and where neither does the caller falls back to the plain
    pid check -- degraded, never wrong in the other direction.
    """
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return ""
    if pid <= 0:
        return ""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        # The command name may contain spaces and brackets; everything
        # after the closing bracket is fixed-width fields, and the 20th
        # of those is the start time in clock ticks since boot.
        fields = stat.rsplit(")", 1)[1].split()
        return fields[19]
    except (OSError, IndexError, ValueError):
        pass
    try:
        out = subprocess.run(["ps", "-o", "lstart=", "-p", str(pid)],
                             capture_output=True, text=True, timeout=2)
    except Exception:
        return ""
    return out.stdout.strip() if out.returncode == 0 else ""


def _tmux_session() -> str:
    """The tmux session this process runs in, or "".

    Asked of tmux rather than parsed out of $TMUX: the variable holds a
    socket path and a pane id, not the name a person types after
    ``attach -t``.
    """
    if not os.environ.get("TMUX"):
        return ""
    try:
        out = subprocess.run(["tmux", "display-message", "-p", "#S"],
                             capture_output=True, text=True, timeout=2)
    except Exception:
        return ""
    return out.stdout.strip() if out.returncode == 0 else ""


def announce_dashboard(*, port: int = 0, url: str = "", token: str = "",
                       resume_path: str = "") -> str:
    """Leave the note that says where this dashboard runs. Never raises.

    The token is part of it, so the note alone is enough to walk back in
    -- which is the whole point, and also why the file is written 0600
    inside a 0700 directory, exactly like the token file the launcher
    already writes under /run/user. A note that made you hunt for the
    token would send you back to the terminal you were trying to find.

    Returns the path written, or "".
    """
    path = record_path()
    record = {
        "host": socket.gethostname(),
        "pid": os.getpid(),
        "proc_start": _process_start(os.getpid()),
        "port": int(port or 0),
        "url": str(url or ""),
        "token": str(token or ""),
        "resume_path": str(resume_path or ""),
        "tmux": _tmux_session(),
        "cwd": os.getcwd(),
        "started_at": time.time(),
    }
    try:
        path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
        try:
            os.chmod(path.parent, 0o700)
        except OSError:
            pass
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(record, indent=1), encoding="utf-8")
        try:
            os.chmod(tmp, 0o600)
        except OSError:
            pass
        tmp.replace(path)
    except OSError:
        return ""
    return str(path)


def _is_mine(record: dict) -> bool:
    """Did THIS process write that note?

    Pid alone does not answer it: pids are reused, and the number in a
    note left yesterday may belong to something else today. The start
    time of the process pins it to one life of one pid.
    """
    try:
        if int(record.get("pid") or 0) != os.getpid():
            return False
    except (TypeError, ValueError):
        return False
    if str(record.get("host") or "") != socket.gethostname():
        return False
    written = str(record.get("proc_start") or "")
    mine = _process_start(os.getpid())
    # An older note carries no fingerprint; the pid and host still match,
    # so honour it rather than leaving a note nobody will ever remove.
    return not written or not mine or written == mine


def withdraw_dashboard() -> None:
    """Remove the note when the dashboard stops. Never raises.

    Only ever removes THIS process's own note. Anything else deletes the
    way back to a dashboard that is still serving -- which is what an
    ordinary test run used to do on its way out.
    """
    path = record_path()
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if not isinstance(record, dict) or not _is_mine(record):
        return
    try:
        path.unlink()
    except OSError:
        pass


def _still_running(record: dict) -> Optional[bool]:
    """True/False on this machine, None when it cannot be asked."""
    if str(record.get("host") or "") != socket.gethostname():
        return None
    try:
        pid = int(record.get("pid") or 0)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return None
    # The number is in use -- by the same process, or by whatever the
    # system handed it to next. Saying "running" for a stranger sends
    # somebody to a port nothing listens on.
    written = str(record.get("proc_start") or "")
    now = _process_start(pid)
    if written and now and written != now:
        return False
    return True


def dashboard() -> dict:
    """The dashboard note, with ``running`` added, or {}."""
    try:
        record = json.loads(record_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(record, dict):
        return {}
    record["running"] = _still_running(record)
    return record


def reconnect_command(record: dict) -> str:
    """What to type to get back to it, from anywhere."""
    host = str(record.get("host") or "")
    if not host:
        return ""
    session = str(record.get("tmux") or "")
    if session:
        return f"ssh {host} && tmux attach -t {session}"
    return f"ssh {host}"


def forward_command(record: dict) -> str:
    """What to type to reach its port from a laptop."""
    host = str(record.get("host") or "")
    port = int(record.get("port") or 0)
    if not host or not port:
        return ""
    return f"ssh -L {port}:localhost:{port} {host}"


def resume_links(record: dict, rows: list[dict]) -> list[tuple]:
    """(session name, the address that walks back into it).

    A session that is kept has a record of its own; the dashboard note
    has the port and the token. Neither half is enough alone, which is
    why getting back in meant reading two files and pasting them
    together by hand.
    """
    port = int(record.get("port") or 0)
    token = str(record.get("token") or "")
    path = str(record.get("resume_path") or "")
    if not (port and path):
        return []
    out: list[tuple] = []
    for name in kept_session_names():
        query = f"?session={name}"
        if token:
            query = f"?token={token}&session={name}"
        out.append((name, f"http://localhost:{port}{path}{query}"))
    return out


def kept_session_names() -> list[str]:
    """The sessions armed to survive a closed window."""
    try:
        from delfin.dashboard import session as _session

        return sorted(str(r.get("session_name") or "")
                      for r in _session.list_records()
                      if r.get("session_name"))
    except Exception:
        return []


def sessions() -> list[dict]:
    """The agent sessions that are open, as presence recorded them."""
    try:
        from . import session_presence

        return list(session_presence.open_sessions())
    except Exception:
        return []


def _ago(ts: Any) -> str:
    try:
        seconds = time.time() - float(ts)
    except (TypeError, ValueError):
        return "-"
    if seconds < 0:
        return "-"
    if seconds < 90:
        return f"{int(seconds)}s ago"
    if seconds < 5400:
        return f"{int(seconds // 60)}m ago"
    return f"{int(seconds // 3600)}h ago"


def format_text(record: dict, rows: list[dict]) -> str:
    """The answer to "where is it", in the fewest lines that say it."""
    out: list[str] = []
    if not record:
        out.append("No dashboard has left a note in "
                   f"{record_path().parent} — none has run since this was "
                   "added, or it ran as another user.")
    else:
        running = record.get("running")
        state = ("running" if running is True
                 else "gone" if running is False
                 else "on another node — cannot tell from here")
        line = (f"dashboard  {record.get('host', '?')}"
                f"  port {record.get('port') or '?'}")
        if record.get("tmux"):
            line += f"  tmux {record['tmux']!r}"
        line += f"  started {_ago(record.get('started_at'))}  ({state})"
        out.append(line)
        back = reconnect_command(record)
        if back:
            out.append(f"  come back with:  {back}")
        forward = forward_command(record)
        if forward:
            out.append(f"  reach the port:  {forward}")
        links = resume_links(record, rows)
        if links:
            out.append("  walk back into a kept session:")
            for name, url in links:
                out.append(f"    {name}  {url}")
    if rows:
        out.append("")
        for row in rows:
            out.append(
                f"session    {row.get('host', '?')}  {row.get('key', '?')!r}"
                f"  branch {row.get('branch') or '-'}"
                f"  seen {_ago(row.get('updated_at'))}")
    elif record:
        out.append("")
        out.append("No agent session is currently announced.")
    return "\n".join(out)


def main(argv: Optional[list] = None) -> int:
    print(format_text(dashboard(), sessions()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
