"""A worker keeps the processes its session started as descendants.

Input: none. Output: ``become_subreaper()`` True when the calling process
is now a child subreaper (Linux ``prctl(PR_SET_CHILD_SUBREAPER)``), False
on any other platform or when the call fails; ``Reaper`` collects the
orphans the process adopts that way.

Semantics. Without a subreaper, a process an agent command detaches
(``setsid -f``, a double fork, ``nohup ... &`` from a shell that exits) is
re-parented to pid 1, and the supervisor's ancestry walk
(``supervisor._from_an_agent``) no longer finds the worker above it. With
the worker as subreaper the orphan is re-parented to the worker instead,
so the walk still ends at a worker. What the walk cannot see -- a process
whose worker has itself ended, a process started outside the worker's
tree -- is unchanged.

A subreaper must also reap what it adopts, or each adopted orphan stays a
zombie until the worker exits. ``waitpid(-1)`` cannot be used: it would
take the exit status of children the engine's ``subprocess.Popen`` objects
wait for (Popen then reports 0 on ECHILD, so a failed command would read
as a success). ``Reaper`` therefore reaps a direct child only when an
earlier scan saw it below another of the worker's children (it was
adopted), with the same start time (the pid was not reused), and it is a
zombie. An orphan whose parent ended between two scans is not recognised
and stays a zombie until the worker ends; init reaps it then.

Cost: one scan reads ``/proc/<pid>/task/*/children`` and ``/proc/<pid>/stat``
for the worker's own subtree only, every ``interval`` seconds (5 s).
Rejected: a SIGCHLD handler (the same Popen problem, and Python runs
handlers on the main thread, which the session loop occupies); reaping in
the supervisor (orphans would then be descendants of the service, not of a
worker, and windows started from a test process would look agent-born).
"""

from __future__ import annotations

import os
import sys
import threading

PR_SET_CHILD_SUBREAPER = 36


def become_subreaper() -> bool:
    """Make this process the reaper of its orphaned descendants (Linux)."""
    if not sys.platform.startswith("linux"):
        return False
    try:
        import ctypes
        libc = ctypes.CDLL(None, use_errno=True)
        return int(libc.prctl(PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0)) == 0
    except Exception:
        return False


def _children(pid: int) -> list[int]:
    out: list[int] = []
    try:
        tasks = os.listdir(f"/proc/{pid}/task")
    except OSError:
        return out
    for tid in tasks:
        try:
            with open(f"/proc/{pid}/task/{tid}/children") as fh:
                out.extend(int(x) for x in fh.read().split())
        except (OSError, ValueError):
            continue
    return out


def _stat(pid: int) -> tuple[str, str] | None:
    """(state, start time) from /proc/<pid>/stat, or None."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            text = fh.read()
        rest = text[text.rindex(")") + 2:].split()
        return rest[0], rest[19]
    except (OSError, ValueError, IndexError):
        return None


class Reaper:
    """Reaps the orphans a subreaper adopted, and nothing else."""

    def __init__(self, interval: float = 5.0):
        self.interval = interval
        #: pid -> start time of every process seen below a direct child.
        self.adoptable: dict[int, str] = {}
        self.reaped: list[int] = []

    def reap_once(self) -> None:
        me = os.getpid()
        direct = _children(me)
        seen: dict[int, str] = {}
        for pid in direct:
            known = self.adoptable.get(pid)
            if known is None:
                continue
            st = _stat(pid)
            if st is None or st[1] != known:
                continue
            if st[0] == "Z":
                try:
                    done, _ = os.waitpid(pid, os.WNOHANG)
                except ChildProcessError:
                    done = 0
                if done == pid:
                    self.reaped.append(pid)
                continue
            seen[pid] = known            # adopted and still running
        stack = [c for pid in direct for c in _children(pid)]
        hops = 0
        while stack and hops < 100_000:
            pid = stack.pop()
            hops += 1
            st = _stat(pid)
            if st is None:
                continue
            seen[pid] = st[1]
            stack.extend(_children(pid))
        self.adoptable = seen

    def _loop(self) -> None:
        import time
        while True:
            try:
                self.reap_once()
            except Exception:
                pass
            time.sleep(self.interval)

    def start(self) -> "Reaper":
        threading.Thread(target=self._loop, name="subreaper", daemon=True).start()
        return self


def keep_descendants() -> bool:
    """Become a subreaper and reap what is adopted; False where the
    platform has no subreaper (then nothing is started)."""
    if not become_subreaper():
        return False
    Reaper().start()
    return True
