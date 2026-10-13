"""Where the agent service runs on this host, and what keeps it single.

The service is one per user and host, like a tmux server:

    /tmp/delfin-<uid>/                   0700, the user's own, never a symlink
        agent.sock                       0600, the service's socket
        lock                             flock'd for the service's lifetime
    ~/.delfin/agent_service/run/<host>.json   0600: pid, socket, version

Why /tmp and not ~/.delfin: on HPC clusters $HOME is NFS shared by every
login node, and a Unix socket can only be connected to on the host that
created it. The run file stays in $HOME on purpose: a window on another
login node finds it and can say where the sessions are running. Why not
$XDG_RUNTIME_DIR: it is removed at the user's last logout, which would
orphan a service whose whole point is to outlive the window.

The lock does two jobs: a second service on the same host sees it held and
leaves, and systemd-tmpfiles skips a locked directory when it ages /tmp.

There is no secret in any of these files: every one of them can be read by
every process of the user, an agent's commands included, so a token here
would protect nothing. Who may connect is the peer's uid; what a connection
may do is decided by where it comes from (supervisor._from_an_agent); and an
agent's commands do not see this directory at all where their cage runs
(api_client._session_doors).
"""

from __future__ import annotations

import errno
import json
import os
import socket
import stat
import struct
import sys
import time
from pathlib import Path

from .eventlog import service_dir

RUNTIME_DIR_ENV = "DELFIN_RUNTIME_DIR"


class UnsafeRuntimeDir(RuntimeError):
    """The runtime directory exists but is not safely the user's own."""


class ServiceBusy(RuntimeError):
    """Another service on this host holds the lock."""


def hostname() -> str:
    return socket.gethostname().split(".")[0] or "localhost"


def runtime_dir() -> Path:
    base = os.environ.get(RUNTIME_DIR_ENV, "").strip() or "/tmp"
    return Path(base) / f"delfin-{os.getuid()}"


def socket_path() -> Path:
    return runtime_dir() / "agent.sock"


def run_file(host: str | None = None) -> Path:
    return service_dir() / "run" / f"{host or hostname()}.json"


def ensure_runtime_dir() -> Path:
    """Create the directory, or accept it only if it is safely ours.

    tmux's rule: a symlink, another owner or a non-directory is refused
    outright. Loose permissions on our own directory are tightened.
    """
    path = runtime_dir()
    try:
        os.mkdir(path, 0o700)
    except FileExistsError:
        pass
    st = os.lstat(path)
    if stat.S_ISLNK(st.st_mode) or not stat.S_ISDIR(st.st_mode):
        raise UnsafeRuntimeDir(f"{path} is not a plain directory")
    if st.st_uid != os.getuid():
        raise UnsafeRuntimeDir(f"{path} belongs to uid {st.st_uid}, not to you")
    if stat.S_IMODE(st.st_mode) & 0o077:
        os.chmod(path, 0o700)
    return path


class ServiceLock:
    """The single-instance lock: flock on a file in the runtime directory
    and (best effort) on the directory itself."""

    def __init__(self) -> None:
        self._fds: list[int] = []

    def acquire(self) -> None:
        import fcntl

        directory = ensure_runtime_dir()
        fd = os.open(directory / "lock", os.O_RDWR | os.O_CREAT, 0o600)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            os.close(fd)
            if exc.errno in (errno.EWOULDBLOCK, errno.EAGAIN, errno.EACCES):
                raise ServiceBusy("another agent service runs on this host") from None
            raise
        self._fds.append(fd)
        try:
            dfd = os.open(directory, os.O_RDONLY)
            fcntl.flock(dfd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self._fds.append(dfd)
        except OSError:
            pass

    def release(self) -> None:
        for fd in self._fds:
            try:
                os.close(fd)
            except OSError:
                pass
        self._fds.clear()


def write_run_file(info: dict) -> Path:
    from .. import state_paths

    path = run_file()
    state_paths.ensure_dir(path.parent)
    state_paths.write_text_atomic(path, json.dumps(info, indent=2))
    state_paths.secure_file(path)
    return path


def read_run_file(host: str | None = None) -> dict | None:
    try:
        data = json.loads(run_file(host).read_text())
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def remove_run_file(pid: int) -> None:
    """Remove this host's run file -- only if it still names ``pid``."""
    current = read_run_file()
    if current and int(current.get("pid", 0) or 0) == pid:
        try:
            run_file().unlink()
        except OSError:
            pass


def other_hosts() -> list[dict]:
    """Run files of services on other hosts (other login nodes)."""
    here = hostname()
    found = []
    for path in sorted((service_dir() / "run").glob("*.json")):
        if path.stem == here:
            continue
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        if isinstance(data, dict):
            found.append(data)
    return found


def peer_uid(sock: socket.socket) -> int | None:
    """The uid of the process at the other end of a Unix socket, if the
    platform says; None where it cannot be asked."""
    if sys.platform.startswith("linux") and hasattr(socket, "SO_PEERCRED"):
        creds = sock.getsockopt(socket.SOL_SOCKET, socket.SO_PEERCRED,
                                struct.calcsize("3i"))
        _pid, uid, _gid = struct.unpack("3i", creds)
        return uid
    if hasattr(os, "getpeereid"):                 # pragma: no cover - BSD/macOS
        try:
            return os.getpeereid(sock.fileno())[0]  # type: ignore[attr-defined]
        except OSError:
            return None
    return None


def started_at() -> float:
    return round(time.time(), 3)
