"""The one dashboard server of a user on a host: where it is and whether it runs.

``delfin-voila`` runs one server per user and host, the way the agent
service does (``agent/service/rundir.py``). Its launcher -- the process
that started the server and stops it, detached from the terminal or in
the foreground -- writes a record when the server answers:

    ~/.delfin/dashboard_server/run/<host>.json    0600, in a 0700 directory
        pid, proc_start   the launcher, as proc_identity names a process
        server_pid        the jupyter-server it started
        port, ip, scheme  where it listens
        url               the dashboard's render path
        resume_path       the one-cell notebook a kept dashboard returns by
        token             the login token (why the file is the owner's)
        log, started_at, detached, version
    ~/.delfin/dashboard_server/run/<host>.lock    flock'd for the launcher's life
    ~/.delfin/dashboard_server/logs/<host>.log    a detached launcher's output

A record is believed only while its pid is alive AND that pid's start
time is the one written (``proc_identity.alive``): a recycled pid is a
stale record, which the next launch replaces. Per host, because on a
cluster the home directory is shared by the login nodes and a pid means
nothing on another one.

The token is in the record for the same reason it is in the dashboard
note (``agent/where.py``): the record alone is enough to walk back in.
Nobody but the owner can read either; the agent's commands do not see
``~/.delfin`` where their cage runs.

Rejected: the dashboard note alone. It is one file for every host that
shares the home directory (a later dashboard on another login node
overwrites it), so it cannot say whether THIS host runs a server.
Rejected: a lock in /tmp beside the agent service's socket. The suite
starts dashboard servers under private homes on the same host; a
host-wide lock would refuse them or a real one.
"""

from __future__ import annotations

import json
import os
import socket
import time
from pathlib import Path
from typing import Optional

#: Minutes a server with no kernel and no request waits before it ends
#: itself. A kept dashboard and an open window each hold a kernel, so
#: they always hold the server. ``0`` keeps the server up.
IDLE_EXIT_ENV = "DELFIN_DASHBOARD_IDLE_EXIT_MINUTES"

#: The default idle grace, in seconds. Ten minutes: the time the server
#: already waited for a first window (resume_server.NEVER_CONNECTED_SECONDS),
#: so a URL printed in one terminal and opened through an ssh forward set
#: up afterwards still reaches it, and a person who closed the last window
#: and comes back within minutes finds the same address and login. Longer
#: holds a port and a few hundred MB on a shared login node for nobody;
#: shorter costs only the few seconds of the next ``delfin-voila``.
IDLE_EXIT_SECONDS = 600.0

#: Set in a detached launcher and its server: the server then waits out
#: the idle grace instead of stopping when its last kernel ends.
DETACHED_ENV = "DELFIN_VOILA_DETACHED"


def _state_dir() -> Path:
    """Resolved per call (state_paths.USER_STATE_RESOLVERS redirects it)."""
    return Path.home() / ".delfin" / "dashboard_server"


def hostname() -> str:
    return (socket.gethostname() or "localhost").split(".")[0] or "localhost"


def record_file(host: str = "") -> Path:
    return _state_dir() / "run" / f"{host or hostname()}.json"


def lock_file(host: str = "") -> Path:
    return _state_dir() / "run" / f"{host or hostname()}.lock"


def log_file(host: str = "") -> Path:
    return _state_dir() / "logs" / f"{host or hostname()}.log"


def is_detached(env: Optional[dict] = None) -> bool:
    e = os.environ if env is None else env
    return str(e.get(DETACHED_ENV, "")).strip() == "1"


def idle_exit_seconds(env: Optional[dict] = None) -> float:
    """The idle grace in seconds; 0 means the server stays up."""
    e = os.environ if env is None else env
    raw = str(e.get(IDLE_EXIT_ENV, "")).strip()
    if not raw:
        return IDLE_EXIT_SECONDS
    try:
        minutes = float(raw)
    except ValueError:
        return IDLE_EXIT_SECONDS
    if minutes != minutes or minutes < 0:          # NaN, negative
        return IDLE_EXIT_SECONDS
    return minutes * 60.0


# -- the record -----------------------------------------------------------------

def write(record: dict) -> Path:
    """Write this host's record, owner-only. Returns its path."""
    from delfin.agent import state_paths

    path = record_file(str(record.get("host") or ""))
    state_paths.ensure_dir(path.parent)
    state_paths.write_text_atomic(path, json.dumps(record, indent=1))
    state_paths.secure_file(path)
    return path


def read(host: str = "") -> Optional[dict]:
    try:
        data = json.loads(record_file(host).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def alive(record: Optional[dict]) -> bool:
    """The record's launcher runs: same host, same pid, same start time."""
    if not record:
        return False
    from delfin.agent import proc_identity

    if str(record.get("host") or "") != hostname():
        return False
    return proc_identity.alive(record.get("pid"),
                               str(record.get("proc_start") or "")) is True


def running() -> Optional[dict]:
    """This host's record while its launcher runs, else None.

    A record whose launcher is gone (or whose pid now names another
    process) is removed: it is the stale record the next launch replaces.
    """
    record = read()
    if record is None:
        return None
    if alive(record):
        return record
    try:
        record_file().unlink()
    except OSError:
        pass
    return None


def own_record(pid: int = 0) -> dict:
    """The identity part of a record for this process (or *pid*)."""
    from delfin.agent import proc_identity

    pid = int(pid or os.getpid())
    return {"host": hostname(), "pid": pid,
            "proc_start": proc_identity.process_start(pid)}


def withdraw(path: Optional[Path] = None) -> None:
    """Remove the record if this process wrote it. Never raises.

    *path* is fixed when the hook is registered, so an ``atexit`` run
    after the home was redirected back still names the file written.
    """
    target = Path(path) if path is not None else record_file()
    try:
        record = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    mine = own_record()
    if not isinstance(record, dict) or any(
            str(record.get(k) or "") != str(mine[k] or "")
            for k in ("host", "pid", "proc_start")):
        return
    try:
        target.unlink()
    except OSError:
        pass


class Lock:
    """The host's single-server lock: flock on ``<host>.lock``, held by the
    launcher for its life (released by the kernel when it ends)."""

    def __init__(self) -> None:
        self._fd: Optional[int] = None

    def acquire(self) -> bool:
        import fcntl

        from delfin.agent import state_paths

        path = lock_file()
        try:
            state_paths.ensure_dir(path.parent)
            fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o600)
        except OSError:
            # No lock to be had (a read-only home): the record still
            # finds a running server; two racing starts are the rare cost.
            return True
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            os.close(fd)
            return False
        self._fd = fd
        return True

    def release(self) -> None:
        if self._fd is not None:
            try:
                os.close(self._fd)
            except OSError:
                pass
            self._fd = None


def wait_for_record(*, pid: int = 0, timeout: float = 60.0,
                    poll=None, sleep=time.sleep) -> Optional[dict]:
    """The record of a running server (of launcher *pid*, if given), once
    its port answers; None after *timeout* or when *poll* says the
    launcher ended (``poll() is not None``)."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        record = read()
        if record and alive(record) and (not pid or int(record.get("pid") or 0) == pid):
            if port_answers(record):
                return record
        if poll is not None and poll() is not None:
            # One more look: the record may have landed with the exit.
            record = read()
            if record and alive(record) and port_answers(record):
                return record
            return None
        sleep(0.25)
    return None


# -- asking the server -----------------------------------------------------------

def _connect_host(record: dict) -> str:
    ip = str(record.get("ip") or "127.0.0.1")
    return "127.0.0.1" if ip in ("", "0.0.0.0", "localhost", "::") else ip


def port_answers(record: dict, timeout: float = 0.5) -> bool:
    try:
        port = int(record.get("port") or 0)
    except (TypeError, ValueError):
        return False
    if port <= 0:
        return False
    try:
        with socket.create_connection((_connect_host(record), port), timeout=timeout):
            return True
    except OSError:
        return False


def _request(record: dict, method: str, path: str, timeout: float = 5.0):
    """One authenticated request to the server's REST API.

    ``no_track_activity`` keeps the question from counting as activity:
    asking whether a server is idle must not keep it from going idle.
    The token goes in the header, never in a URL.
    """
    import http.client
    import ssl

    host = _connect_host(record)
    port = int(record.get("port") or 0)
    sep = "&" if "?" in path else "?"
    target = f"{path}{sep}no_track_activity=1"
    if str(record.get("scheme") or "http") == "https":
        # Its own certificate on its own loopback port: the token is the
        # authentication here, TLS only keeps it off the wire.
        context = ssl.create_default_context()
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE
        conn = http.client.HTTPSConnection(host, port, timeout=timeout,
                                           context=context)
    else:
        conn = http.client.HTTPConnection(host, port, timeout=timeout)
    try:
        conn.request(method, target, headers={
            "Authorization": f"token {record.get('token') or ''}"})
        response = conn.getresponse()
        return response.status, response.read()
    finally:
        conn.close()


def kernels(record: dict, timeout: float = 5.0) -> Optional[list]:
    """The server's kernels (id, last_activity, connections, ...), or None
    when it cannot be asked."""
    try:
        status, body = _request(record, "GET", "/api/kernels", timeout)
    except (OSError, ValueError):
        return None
    if status != 200:
        return None
    try:
        data = json.loads(body.decode("utf-8"))
    except ValueError:
        return None
    return data if isinstance(data, list) else None


def end_kernel(record: dict, kernel_id: str, timeout: float = 30.0) -> bool:
    """Ask the server to shut one kernel down. True when it did (204)."""
    safe = "".join(ch for ch in str(kernel_id) if ch.isalnum() or ch == "-")
    if not safe:
        return False
    try:
        status, _body = _request(record, "DELETE", f"/api/kernels/{safe}", timeout)
    except (OSError, ValueError):
        return False
    return status in (200, 204)


def from_kernel_environment(env: Optional[dict] = None) -> Optional[dict]:
    """The server a dashboard kernel runs under, as a record-shaped dict,
    from the variables the launcher gives every kernel; None outside one."""
    e = os.environ if env is None else env
    token = str(e.get("JUPYTER_TOKEN") or "")
    try:
        port = int(e.get("DELFIN_VOILA_PORT") or 0)
    except ValueError:
        port = 0
    if not token or port <= 0:
        return None
    return {"port": port, "token": token, "ip": "127.0.0.1", "scheme": "http"}


# -- what it holds ----------------------------------------------------------------

def _seconds_since(iso: str) -> Optional[float]:
    from datetime import datetime, timezone

    try:
        when = datetime.fromisoformat(str(iso).replace("Z", "+00:00"))
    except ValueError:
        return None
    if when.tzinfo is None:
        when = when.replace(tzinfo=timezone.utc)
    return max(0.0, (datetime.now(timezone.utc) - when).total_seconds())


def ago(seconds: Optional[float]) -> str:
    if seconds is None:
        return "-"
    if seconds < 90:
        return f"{int(seconds)} s ago"
    if seconds < 5400:
        return f"{int(seconds // 60)} min ago"
    if seconds < 48 * 3600:
        return f"{seconds / 3600:.0f} h ago"
    return f"{seconds / 86400:.0f} d ago"


def span(seconds: float) -> str:
    seconds = max(0, int(seconds))
    if seconds < 120:
        return f"{seconds} s"
    if seconds < 3600:
        return f"{seconds // 60} min"
    if seconds < 172800:
        return f"{seconds // 3600} h {seconds % 3600 // 60} min"
    return f"{seconds // 86400} d {seconds % 86400 // 3600} h"


def resume_link(record: dict, name: str) -> str:
    """The address back into the kept dashboard *name*, without the token
    (session.resume_url: a token in a URL stays in browser history)."""
    from urllib.parse import quote

    port = int(record.get("port") or 0)
    path = str(record.get("resume_path") or "").split("?")[0]
    if not (port and path and name):
        return ""
    scheme = str(record.get("scheme") or "http")
    return f"{scheme}://localhost:{port}{path}?session={quote(name)}"


def kept_dashboards(record: Optional[dict] = None,
                    kernel_list: Optional[list] = None) -> list[dict]:
    """The kept dashboards (Keep dashboard on) of this host, newest first.

    From the kept-session records (``session.live_records``) -- the
    registry a resume follows -- joined with the server's kernels for
    the windows and the last activity when the server could be asked.
    Each row: name, kernel_id, started_at, url, windows, last_activity
    (seconds ago, or None), state ('running' / 'gone').
    """
    from delfin.dashboard import session as _session

    by_id = {str(k.get("id") or ""): k for k in (kernel_list or [])
             if isinstance(k, dict)}
    rows = []
    for rec in _session.live_records():
        if str(rec.get("host") or "") != hostname():
            continue
        name = str(rec.get("session_name") or "")
        kid = str(rec.get("kernel_id") or "")
        kernel = by_id.get(kid)
        row = {
            "name": name,
            "kernel_id": kid,
            "started_at": rec.get("started_at"),
            "url": resume_link(record or {}, name) if record else "",
            "windows": None,
            "last_activity": None,
            "state": "running",
        }
        if kernel_list is not None:
            if kernel is None:
                row["state"] = "gone"
            else:
                try:
                    row["windows"] = int(kernel.get("connections") or 0)
                except (TypeError, ValueError):
                    row["windows"] = None
                row["last_activity"] = _seconds_since(
                    str(kernel.get("last_activity") or ""))
        rows.append(row)
    return rows


def kept_agent_sessions() -> list[dict]:
    """The agent sessions with Keep session on that this host holds."""
    try:
        from delfin.dashboard import kept_agent_sessions as _kept

        return [r for r in _kept.records()
                if str(r.get("host") or "") == _kept._host()]
    except Exception:
        return []
