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


def _children(pid: int) -> list[int] | None:
    """Direct children of *pid* from /proc/<pid>/task/*/children, or None
    when the kernel does not offer the file (CONFIG_PROC_CHILDREN off)."""
    kids: list[int] = []
    try:
        tasks = os.listdir(f"/proc/{pid}/task")
    except OSError:
        return []
    for tid in tasks:
        try:
            with open(f"/proc/{pid}/task/{tid}/children", "rb") as fh:
                kids.extend(int(x) for x in fh.read().split())
        except FileNotFoundError:
            return None
        except (OSError, ValueError):
            continue
    return kids


def _descendants_by_scan(root_pid: int) -> set[int]:
    """Fallback: one pass over /proc/*/stat into a parent -> children map."""
    children: dict[int, list[int]] = {}
    try:
        entries = os.listdir("/proc")
    except OSError:
        return set()
    for name in entries:
        if not name.isdigit():
            continue
        fields = _stat_fields(int(name))
        if fields is None or len(fields) < 3:
            continue
        try:
            children.setdefault(int(fields[2]), []).append(int(name))
        except ValueError:
            continue
    found: set[int] = set()
    frontier = [root_pid]
    while frontier:
        for kid in children.get(frontier.pop(), ()):
            if kid not in found and kid != root_pid:
                found.add(kid)
                frontier.append(kid)
    return found


def descendants_of(root_pid: int) -> set[int]:
    """All live descendants of *root_pid*, zombies excluded.

    Walks /proc/<pid>/task/*/children from the root, so the cost is the
    size of the call's own tree, not of every process on a shared login
    node (thousands, sampled every poll). Only where the kernel lacks the
    children file does it fall back to one scan of /proc/*/stat.
    """
    found: set[int] = set()
    frontier = [root_pid]
    while frontier:
        kids = _children(frontier.pop())
        if kids is None:
            found = _descendants_by_scan(root_pid)
            break
        for kid in kids:
            if kid not in found and kid != root_pid:
                found.add(kid)
                frontier.append(kid)
    live = set()
    for pid in found:
        fields = _stat_fields(pid)
        if fields is not None and not _is_zombie(fields):
            live.add(pid)
    return live


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


# No rlimits on the child. RLIMIT_NPROC counts every process of the
# whole USER on the machine, not this call's children, so on a shared
# login node it either never fires or fails innocent fork() calls
# elsewhere. RLIMIT_AS limits virtual address space per process, and
# tools that reserve large address ranges without touching them (MPI,
# OpenMP runtimes, JVMs, glibc arenas per thread) break far below their
# real memory use. RLIMIT_CPU is per process, not per call. The census
# above measures what the budget means; it is the one line.


class BudgetGuard:
    """Watches one call's process group and ends it on a breach.

    The guard signals ONLY the process group of *root_pid* — never any
    other process on the node. TERM first, then KILL after a short
    grace so a well-behaved tree can clean up. Daemon thread so an
    abandoned guard never keeps a session alive.
    """

    def __init__(self, root_pid: int, limits: BudgetLimits | None = None,
                 poll_s: float = 1.0, term_grace_s: float = 5.0,
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
        if pgid <= 1 or pgid == os.getpgrp():
            # Never init's group, and never our own: a call that somehow
            # shares the agent's group must not take the agent with it.
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
