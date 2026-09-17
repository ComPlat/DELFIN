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
    elif sys.platform == "darwin":
        # No process of the user may attach a debugger (PT_DENY_ATTACH).
        try:
            ok = _libc().ptrace(31, 0, 0, 0) == 0
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


#: The model providers' keys DELFIN itself runs on.
PROVIDER_KEYS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY")


def exported_provider_keys() -> list[str]:
    """Names (never values) of provider keys present in this environment."""
    return [k for k in PROVIDER_KEYS if os.environ.get(k)]


def exported_key_advice(names) -> str:
    return (f"{', '.join(names)} is exported in your environment. DELFIN's own "
            "processes are protected, but the shell you exported it in, and "
            "anything else started from it, can be read by every process of "
            "yours through /proc. Store it with `delfin-agent credentials set "
            f"{names[0]}` and remove the export.")


def put_exported_keys_away(names=()) -> str:
    """Do it instead of asking for it. Returns one line, or "".

    The warning above repeated at every start and left four manual steps
    -- store the key, find the export, delete it, open a new shell -- and
    a key sat in a group-readable shell file in the meantime (seen on a
    cluster, 2026-09-17). The key is taken into the 0600 store and the
    export is commented out with a copy of the file kept beside it.

    Silent when there is nothing to do. A stored value that DIFFERS from
    the exported one is left alone and said out loud: which of the two is
    the right key is not this function's to decide.
    """
    try:
        from . import credentials as _cred
        rows = _cred.secure_exported_keys(names or exported_provider_keys())
    except Exception:
        return ""
    stored = [r["name"] for r in rows if r["action"] == "stored"]
    cleaned = [(r["name"], c) for r in rows for c in r.get("cleaned") or []]
    differs = [r["name"] for r in rows if r["action"] == "differs"]
    parts: list[str] = []
    if stored:
        parts.append(f"{', '.join(stored)} taken into "
                     f"{_cred.credentials_path()} (0600)")
    for name, where in cleaned:
        parts.append(f"the line exporting {name} in {where['file']}:"
                     f"{where['line']} is commented out "
                     f"(copy: {where['backup']})")
    for row in rows:
        if row.get("systemd"):
            parts.append(f"{row['name']} removed from the systemd user "
                         "environment (put it back with `systemctl --user "
                         f"set-environment {row['name']}=...` if you meant "
                         "it to be there)")
    unplaced = [r["name"] for r in rows
                if r["action"] in ("stored", "already")
                and not r.get("exports") and not r.get("systemd")]
    if unplaced:
        parts.append(f"{', '.join(unplaced)} is still exported by something "
                     "outside the home directory (a site profile, a job "
                     "script, the parent shell) — the store holds it now, so "
                     "that export can go")
    if differs:
        parts.append(f"{', '.join(differs)} is exported AND stored with a "
                     "different value — neither was changed; decide which is "
                     "current and store that one")
    if not parts:
        return ""
    return ("Key hygiene: " + "; ".join(parts)
            + ". Open a new shell so the export is gone from it too.")

