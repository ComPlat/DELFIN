"""Is that process still the one we meant?

A pid is a number the system hands out again. Anything that writes a pid
down and reads it back later -- a dashboard note, a scratch directory, a
job record -- is asking about a process that may have ended while some
stranger took its number. Answering with ``os.kill(pid, 0)`` alone sends
somebody to a port nothing listens on, or keeps a dead session's scratch
for ever.

The start time settles it: a pid plus the moment that process began is a
name no later process shares. Linux hands it over in /proc; everywhere
else ``ps`` answers the same question; where neither does, the answer
degrades to the plain pid check -- never wrong in the other direction.

This module exists so there is exactly ONE reader of that fact. It was
extracted from ``where.py``, which asked it first about the dashboard
note, when the scratch directories needed the same answer: two copies of
one answer drift, and a test pins the pair.
"""

from __future__ import annotations

import os
import socket
import subprocess
from pathlib import Path
from typing import Optional


def start_ticks(pid) -> Optional[int]:
    """*pid*'s start time in clock ticks since boot, or None off Linux.

    Field 22 of ``/proc/<pid>/stat``, and it must be counted from the
    LAST closing bracket. Field 2 is the command name in brackets and may
    contain spaces, so splitting the whole line and taking index 21 reads
    a different field for those processes. Measured on a binary named
    ``a b``: the whole-line split answered 0 -- the boot moment -- where
    the real value was 26055418, which is a job reported as having run
    for the machine's entire uptime.
    """
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return None
    if pid <= 0:
        return None
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        return int(stat[stat.rindex(")") + 1:].split()[19])
    except (OSError, IndexError, ValueError):
        return None


def process_start(pid: int) -> str:
    """A fingerprint telling one life of *pid* from the next, or "".

    The value is opaque -- compare it, do not read it. It is clock ticks
    since boot on Linux and a date from ``ps`` elsewhere, and neither is
    comparable with the other.
    """
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return ""
    if pid <= 0:
        return ""
    ticks = start_ticks(pid)
    if ticks is not None:
        return str(ticks)
    try:
        out = subprocess.run(["ps", "-o", "lstart=", "-p", str(pid)],
                             capture_output=True, text=True, timeout=2)
    except Exception:
        return ""
    return out.stdout.strip() if out.returncode == 0 else ""


def alive(pid, proc_start: str = "", host: str = "") -> Optional[bool]:
    """True/False on this machine, None when it cannot be asked.

    *host* is the machine the pid was written on. A pid from another one
    means nothing here, and guessing costs somebody their running work.
    """
    if host and host != socket.gethostname():
        return None
    try:
        pid = int(pid or 0)
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
    # system handed it to next.
    written = str(proc_start or "")
    now = process_start(pid)
    if written and now and written != now:
        return False
    return True
