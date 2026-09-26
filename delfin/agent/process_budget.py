"""Per-tool-call process and memory budget.

On 2026-09-25 a single tool call spawned ~75 nested ``bash -lc
"module -t avail"`` processes; nothing in DELFIN bounded it. This module
is that bound: limits as data, a /proc-only census of a call's
descendants, a breach verdict, and a guard that ends the call's own
process group when the verdict fires. It is deliberately free of any
import from api_client or sandbox so the security core stays untouched.
"""
from __future__ import annotations

import errno
import os
import signal
import threading
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BudgetLimits:
    max_processes: int = 64
    max_rss_mb: float = 8192.0
    max_cpu_seconds: float = 1200.0


@dataclass(frozen=True)
class Breach:
    limit: str
    value: float
    ceiling: float
    message: str


def _default_limits() -> BudgetLimits:
    """Defaults for a shared login node, overridable via settings.

    Precedence mirrors the other ``agent.*`` keys (see
    api_client._resolve_max_tool_rounds): an explicit
    ``agent.process_budget`` mapping in settings.json wins, individual
    subkeys merge over the built-in defaults, 0 or negative disables the
    respective limit.

    Rationale for the built-in numbers — a login node is shared by dozens
    of interactive users, so the budget must cover a legitimate tool call
    without letting one swallow the node:
      * 64 processes: an install or a small build fans out (make -j,
        pip, conda) but stays comfortably under ~64 concurrent live
        descendants; the observed runaway needed ~78 and would have been
        stopped at 64.
      * 8192 MB RSS summed over descendants: a modest compile or import
        fan-out fits; one runaway process cannot eat the node's memory.
      * 1200 CPU-seconds total (not per process): the existing wall-clock
        timeout bounds a call that *waits*; this bounds one that *burns*
        CPU, e.g. a spin loop, on a node where CPU is the shared good.
    """
    limits = BudgetLimits()
    try:
        from delfin import user_settings
        raw = ((user_settings.load_settings() or {}).get("agent") or {}
               ).get("process_budget")
    except Exception:
        raw = None
    if not isinstance(raw, dict):
        return limits
    def _num(key: str, current: float) -> float:
        val = raw.get(key, current)
        try:
            val = float(val)
        except (TypeError, ValueError):
            return current
        return 100_000_000.0 if val <= 0 else val
    return BudgetLimits(
        max_processes=int(_num("max_processes", limits.max_processes)),
        max_rss_mb=_num("max_rss_mb", limits.max_rss_mb),
        max_cpu_seconds=_num("max_cpu_seconds", limits.max_cpu_seconds),
    )


# --- Census: /proc only, no ps, no filesystem tree walk -------------------

_PAGE = os.sysconf("SC_PAGE_SIZE") if hasattr(os, "sysconf") else 4096
_CLK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
_HZ_MB = _PAGE / (1024.0 * 1024.0)


def _stat_fields(pid: int) -> list[str] | None:
    """Fields of /proc/<pid>/stat, or None when the pid is gone.

    The comm field can contain spaces and parentheses, so the parse
    anchors on the LAST ')' in the line, as proc(5) prescribes.
    """
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            line = fh.read().decode("utf-8", "replace")
    except (OSError, ValueError):
        return None
    close = line.rfind(")")
    if close < 0:
        return None
    fields = line[close + 2:].split()
    return [line[:close].split("(", 1)[0]] + fields


def _is_zombie(fields: list[str]) -> bool:
    # The list is [pid] + stat fields from 3 on, so state is fields[1].
    return len(fields) > 1 and fields[1] == "Z"


def descendants_of(root_pid: int) -> set[int]:
    """All live descendants of *root_pid*, zombies excluded.

    Enumerates /proc once and builds a ppid map (fields[4] of stat),
    then walks it in memory. ``/proc/<pid>/task/*/children`` exists but
    needs one syscall per thread of every process; the single scan of
    ``/proc/*/stat`` covers the same ground and cannot miss a process
    that changed parent while we were reading a children file.
    """
    ppid: dict[int, int] = {}
    proc_root = Path("/proc")
    try:
        entries = list(proc_root.iterdir())
    except OSError:
        return set()
    for entry in entries:
        if not entry.name.isdigit():
            continue
        fields = _stat_fields(int(entry.name))
        if fields is None or _is_zombie(fields) or len(fields) < 5:
            continue
        try:
            ppid[int(entry.name)] = int(fields[2])  # ppid is stat field 4
        except (IndexError, ValueError):
            continue
    # A zombie root or one that already exited contributes nothing.
    if ppid.get(root_pid) is None and root_pid not in ppid:
        return set()
    # Also drop the root itself from the map when we only want children.
    found: set[int] = set()
    frontier = [root_pid]
    seen = {root_pid}
    while frontier:
        current = frontier.pop()
        for pid, parent in ppid.items():
            if parent == current and pid not in seen:
                seen.add(pid)
                found.add(pid)
                frontier.append(pid)
    return found


def count_descendants(root_pid: int) -> int:
    """Number of live (non-zombie) descendants, root excluded."""
    return len(descendants_of(root_pid))


def _rss_mb(pid: int) -> float:
    # statm field 2 (0-based index 1) is resident pages, in pages.
    try:
        with open(f"/proc/{pid}/statm", "rb") as fh:
            parts = fh.read().split()
        return int(parts[1]) * _HZ_MB
    except (OSError, IndexError, ValueError):
        return 0.0


def _cpu_seconds(fields: list[str]) -> float:
    # utime is stat field 14, stime field 15; in this list (pid + fields
    # from 3 on) that is indices 12 and 13.
    try:
        return (int(fields[12]) + int(fields[13])) / float(_CLK)
    except (IndexError, ValueError):
        return 0.0


def sample_descendants(root_pid: int) -> dict:
    """One snapshot: live descendant count, summed RSS, summed CPU time.

    The root process itself is NOT counted as a process (the limit is on
    the fan-out), but its memory and CPU are included in the sums — a
    single bloated process must breach the memory budget even without
    children.
    """
    procs = descendants_of(root_pid)
    rss = _rss_mb(root_pid)
    cpu = 0.0
    root_fields = _stat_fields(root_pid)
    if root_fields is not None and not _is_zombie(root_fields):
        cpu += _cpu_seconds(root_fields)
    for pid in procs:
        rss += _rss_mb(pid)
        fields = _stat_fields(pid)
        if fields is not None:
            cpu += _cpu_seconds(fields)
    return {"processes": len(procs), "rss_mb": rss, "cpu_s": cpu}


def verdict(sample: dict, limits: BudgetLimits) -> Breach | None:
    """Which limit the sample breaks, or None when it is within budget."""
    checks = [
        ("max_processes", sample.get("processes", 0), limits.max_processes,
         "processes"),
        ("max_rss_mb", sample.get("rss_mb", 0.0), limits.max_rss_mb,
         "MB resident memory summed over the call's processes"),
        ("max_cpu_seconds", sample.get("cpu_s", 0.0), limits.max_cpu_seconds,
         "CPU-seconds summed over the call's processes"),
    ]
    for name, value, ceiling, human in checks:
        if value > ceiling:
            return Breach(
                limit=name, value=value, ceiling=ceiling,
                message=(
                    f"Process budget exceeded: {name} reached {value:g} "
                    f"(limit {ceiling:g}) — {human}. The call's process "
                    f"group was terminated; re-run with narrower work, "
                    f"fewer parallel children, or raise "
                    f"agent.process_budget.{name} in settings if the "
                    f"workload legitimately needs more."))
    return None


def child_rlimits(limits: BudgetLimits):
    """A preexec_fn-style mapping of rlimits for the call's own child.

    Second line of defence: where the kernel enforces a limit directly,
    no polling race exists. RLIMIT_CPU bounds total CPU of the direct
    child AND its descendants (the limit is inherited), SIGKILL/SIGXCPU
    on excess; RLIMIT_AS bounds per-process virtual memory, which
    approximates — conservatively, since virtual >= resident — the RSS
    budget.

    Deliberately NOT RLIMIT_NPROC: it counts every process of the whole
    USER on the machine (the kernel check is against the user's total
    process count in the user namespace, not the children of this call),
    so on a shared login node with dozens of the user's sessions it
    would either be set so high it never fires, or kill unrelated work
    by failing innocent fork() calls elsewhere. The process count is
    therefore watched by the census in this module, not by a rlimit.

    Not applied here — sandbox.py (the caller) passes this into its
    preexec_fn; this module must not import or touch the sandbox core.
    """
    import resource
    return {
        resource.RLIMIT_CPU: (int(limits.max_cpu_seconds),
                              int(limits.max_cpu_seconds)),
        resource.RLIMIT_AS: (int(limits.max_rss_mb * 1024 * 1024),
                             int(limits.max_rss_mb * 1024 * 1024)),
    }


class BudgetGuard:
    """Watches one call's process group and ends it on a breach.

    The guard signals ONLY the process group of *root_pid* — never any
    other process on the node. TERM first, then KILL after a short
    grace so a well-behaved tree can clean up. Daemon thread so an
    abandoned guard never keeps a session alive.
    """

    def __init__(self, root_pid: int, limits: BudgetLimits | None = None,
                 poll_s: float = 0.25, term_grace_s: float = 5.0,
                 sample_fn=sample_descendants):
        self.root_pid = root_pid
        self.limits = limits if limits is not None else _default_limits()
        self.poll_s = poll_s
        self.term_grace_s = term_grace_s
        self._sample_fn = sample_fn
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.breach: Breach | None = None

    def start(self) -> "BudgetGuard":
        self._thread = threading.Thread(
            target=self._watch, name=f"process-budget-{self.root_pid}",
            daemon=True)
        self._thread.start()
        return self

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)

    # Exposed for callers that prefer polling to a thread.
    def poll_once(self) -> Breach | None:
        sample = self._sample_fn(self.root_pid)
        breach = verdict(sample, self.limits)
        if breach is not None:
            self.breach = breach
            self._kill_group()
        return breach

    def _kill_group(self) -> None:
        try:
            pgid = os.getpgid(self.root_pid)
        except (ProcessLookupError, PermissionError, OSError):
            return
        if pgid <= 1:
            # Never signal a group we did not create (init's group).
            return
        try:
            os.killpg(pgid, signal.SIGTERM)
        except OSError:
            pass
        deadline = time.monotonic() + self.term_grace_s
        while time.monotonic() < deadline:
            try:
                os.killpg(pgid, 0)
            except ProcessLookupError:
                return
            except PermissionError:
                break
            except OSError as exc:
                if exc.errno == errno.ESRCH:
                    return
                break
            time.sleep(min(0.1, self.poll_s))
        try:
            os.killpg(pgid, signal.SIGKILL)
        except OSError:
            pass

    def _watch(self) -> None:
        while not self._stop.wait(self.poll_s):
            if self.poll_once() is not None:
                return
