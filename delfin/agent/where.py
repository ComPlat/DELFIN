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
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional

from delfin.agent import proc_identity

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


#: A fingerprint that tells one life of a pid from the next, or "".
#: The reading of it lives in ``proc_identity``: the scratch directories
#: ask the same question, and a second reader of /proc drifts from this
#: one. Kept under its old name -- the note format and every caller in
#: this file are written in terms of it.
_process_start = proc_identity.process_start


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
    """True/False on this machine, None when it cannot be asked.

    A note with no host is not answerable either: the pid in it was
    written somewhere, and this is not known to be that machine. The
    process question itself -- including the recycled pid, which would
    otherwise send somebody to a port nothing listens on -- is answered
    in ``proc_identity``, where the scratch sweep asks it too.
    """
    if str(record.get("host") or "") != socket.gethostname():
        return None
    return proc_identity.alive(record.get("pid"),
                               str(record.get("proc_start") or ""))


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
    """Dispatch on the flags; with none of them, the plain print."""
    if argv is None:
        argv = sys.argv[1:]
    import argparse

    p = argparse.ArgumentParser(
        prog="delfin-agent where",
        description="Where the dashboard runs and how to get back into it")
    p.add_argument("--attach", action="store_true",
                   help="Walk into the tmux session the note names")
    p.add_argument("--ssh-config", action="store_true",
                   help="Print a ready ssh config entry for it")
    p.add_argument("--write", action="store_true",
                   help="With --ssh-config: append the entry instead of "
                        "printing it (refuses to change an existing one)")
    args = p.parse_args(argv)

    if args.attach:
        return run_attach(dashboard())
    if args.ssh_config:
        if args.write:
            return write_ssh_config(dashboard())
        print(ssh_config_entry(dashboard()))
        return 0
    print(format_text(dashboard(), sessions()))
    return 0


# -- attach: walk into the session the note names -----------------------------

def run_attach(record: dict) -> int:
    """Replace this process with the way into the note's tmux session.

    Inside that session already: say so, change nothing. On this node:
    tmux directly. On another node: ssh -t with the port forward the note
    filled in, so the URL comes back alive in the same step.

    Nothing here starts anything unattended: the only commands run are
    the ones the user just asked to be placed inside. A host without the
    tool it needs says so in one line and exits 1 -- never a guess.
    """
    if not record:
        print("No dashboard has left a note — there is nothing to attach to.")
        return 1
    session = str(record.get("tmux") or "")
    if not session:
        print("The note names no tmux session; nothing to attach to.")
        return 1
    if session and session == _tmux_session():
        print(f"Already inside tmux session {session!r} — this is it.")
        return 0
    host = str(record.get("host") or "")
    here = socket.gethostname()
    if host and host != here:
        if not shutil.which("ssh"):
            print("ssh is not installed here; cannot reach "
                  f"{host} from this machine.")
            return 1
        port = int(record.get("port") or 0)
        argv = ["ssh", "-t"]
        if port > 0:
            argv += ["-L", f"{port}:localhost:{port}"]
        argv += [host, "tmux", "attach", "-t", session]
    else:
        if not shutil.which("tmux"):
            print("tmux is not installed here; cannot attach to "
                  f"{session!r}.")
            return 1
        argv = ["tmux", "attach", "-t", session]
    # The token stays in the note: this command line is what a terminal
    # and a shell history both record, and the forward alone is enough --
    # the URL on the far side already carries it.
    os.execvp(argv[0], argv)
    return 1          # execvp only returns when it could not run


# -- ssh-config: the entry that makes `ssh delfin` the whole way back ---------

def ssh_config_entry(record: dict) -> str:
    """One ssh config entry for the node the dashboard runs on.

    ``RemoteCommand tmux new -A -s`` attaches when the session exists and
    creates it when it does not -- the exact semantics a reconnect wants.
    The token is NOT part of it: a config file is copied between machines,
    and the note already carries everything private.
    """
    host = str(record.get("host") or "")
    port = int(record.get("port") or 0)
    session = str(record.get("tmux") or "delfin")
    lines = [
        "Host delfin",
        f"    HostName {host or '?'}",
        f"    RemoteCommand tmux new -A -s {session}",
        "    RequestTTY yes",
    ]
    if port > 0:
        lines.append(f"    LocalForward {port} localhost:{port}")
    return "\n".join(lines) + "\n"


def write_ssh_config(record: dict) -> int:
    """Append the entry to ~/.ssh/config, or refuse without changing it.

    An existing ``Host delfin`` block is the user's own work; overwriting
    it silently is the one destructive thing this command could do. The
    refusal says what would change -- old entry and new -- and stops.
    """
    host = str(record.get("host") or "")
    if not host:
        print("No dashboard has left a note — nothing to write.")
        return 1
    config_dir = Path.home() / ".ssh"
    config = config_dir / "config"
    existing = ""
    if config.exists():
        existing = config.read_text(encoding="utf-8")
    for block in existing.split("\n\n"):
        if block.strip().startswith("Host delfin"):
            print("An entry for 'Host delfin' already exists — not "
                  "overwriting it. It says:\n\n"
                  f"{block.strip()}\n\nThe new entry would say:\n\n"
                  f"{ssh_config_entry(record).strip()}\n")
            return 1
    entry = ssh_config_entry(record)
    try:
        config_dir.mkdir(parents=True, mode=0o700, exist_ok=True)
        if config.exists():
            if not existing.endswith("\n"):
                existing += "\n"
            config.write_text(existing + "\n" + entry, encoding="utf-8")
        else:
            config.write_text(entry, encoding="utf-8")
            try:
                os.chmod(config, 0o600)
            except OSError:
                pass
    except OSError as exc:
        print(f"Could not write {config}: {exc}")
        return 1
    print(f"Written to {config} — `ssh delfin` is now the whole way back.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
