"""An agent process cannot be read by the commands it runs.

The model's key lives in the agent's memory, and often in its environment.
Every command the agent runs is a process of the same user, and on Linux
the same user may read ``/proc/<pid>/environ`` and, where ptrace is allowed
between siblings (``kernel.yama.ptrace_scope = 0``, the default on RHEL,
Rocky and many HPC systems), attach a debugger and read memory. The
bubblewrap cage hides the agent behind a PID namespace; on a host without
it nothing did.

``protect()`` marks the process non-dumpable (``PR_SET_DUMPABLE = 0``):
its ``/proc`` files become root's, and no process of the user may attach
to it. It is not inherited across ``exec``, so the commands keep their
normal behaviour. Linux only; elsewhere a no-op.

A non-dumpable process's environment cannot be read by the emergency stop
either, which found dashboard kernels that way. So a protected process
also REGISTERS itself -- host, pid, start time and its lifeline -- in a
directory the stop reads (``registered_here``).

Set ``DELFIN_PROCESS_GUARD=off`` to leave a process debuggable.
"""
from __future__ import annotations

import atexit
import json
import os
import socket
import sys
import time
from pathlib import Path
from typing import Optional

_DIR = Path.home() / ".delfin" / "agent_processes"
_PR_SET_DUMPABLE = 4
_PR_GET_DUMPABLE = 3
_DONE: dict = {}


def _libc():
    import ctypes
    return ctypes.CDLL(None, use_errno=True)


def is_protected() -> bool:
    if not sys.platform.startswith("linux"):
        return False
    try:
        return _libc().prctl(_PR_GET_DUMPABLE, 0, 0, 0, 0) == 0
    except Exception:
        return False


def _host() -> str:
    try:
        return (socket.gethostname() or "").strip()
    except Exception:
        return ""


def _record_path(pid: int) -> Path:
    return _DIR / f"{_host() or 'host'}-{pid}.json"


def protect(kind: str) -> bool:
    """Make this process unreadable by the user's other processes and put
    it on the register the emergency stop reads. Returns whether the
    protection holds. Idempotent, never raises."""
    if os.environ.get("DELFIN_PROCESS_GUARD", "").strip().lower() == "off":
        return False
    pid = os.getpid()
    if _DONE.get(pid):
        return is_protected()
    _DONE[pid] = True
    ok = False
    if sys.platform.startswith("linux"):
        try:
            ok = _libc().prctl(_PR_SET_DUMPABLE, 0, 0, 0, 0) == 0
        except Exception:
            ok = False
    try:
        from .lifeline import ENV_PID, ENV_TICKS, _start_ticks
        record = {
            "pid": pid, "ticks": _start_ticks(pid), "host": _host(),
            "uid": os.getuid() if hasattr(os, "getuid") else -1,
            "kind": str(kind or ""), "registered_at": time.time(),
            "lifeline_pid": os.environ.get(ENV_PID, ""),
            "lifeline_ticks": os.environ.get(ENV_TICKS, ""),
            "voila_port": os.environ.get("DELFIN_VOILA_PORT", ""),
            "jpy_parent_pid": os.environ.get("JPY_PARENT_PID", ""),
        }
        from .state_paths import ensure_dir, write_text_atomic
        ensure_dir(_DIR)
        path = _record_path(pid)
        write_text_atomic(path, json.dumps(record))
        atexit.register(_unregister, path, pid)
    except Exception:
        pass
    return ok


def _unregister(path: Path, pid: int) -> None:
    if os.getpid() != pid:
        return                      # a forked child must not remove it
    try:
        path.unlink()
    except OSError:
        pass


def registered_here() -> list[dict]:
    """The live registered agent processes on this machine. A record whose
    process is gone, or whose pid now belongs to another start, is dropped."""
    out: list[dict] = []
    try:
        entries = list(_DIR.glob(f"{_host() or 'host'}-*.json"))
    except OSError:
        return out
    try:
        from .lifeline import _start_ticks
    except Exception:
        return out
    for path in entries:
        try:
            rec = json.loads(path.read_text(encoding="utf-8"))
            pid = int(rec.get("pid") or 0)
        except (OSError, ValueError, TypeError):
            continue
        ticks = _start_ticks(pid) if pid > 0 else None
        if pid <= 0 or ticks is None or (rec.get("ticks") is not None
                                          and int(rec["ticks"]) != ticks):
            try:
                path.unlink()
            except OSError:
                pass
            continue
        out.append(rec)
    return out


def uid_of(pid: int) -> Optional[int]:
    """The real uid of ``pid`` from /proc/<pid>/status, which stays readable
    for a non-dumpable process (its /proc directory belongs to root)."""
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("Uid:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        return None
    return None
