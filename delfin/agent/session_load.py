"""The load an open agent session puts on this host — reported, never acted on.

A shared login node is for everyone. An agent that starts too many
processes, eats memory or parks cores hurts every other user, and until
now only the operator's own script could see it. This module is DELFIN's
answer: read ``/proc`` (own uid only, no filesystem tree walks), add up
what a session's process tree costs, and say it in plain English.

It never kills, never renices, never writes anything: ``delfin-agent
sessions --load`` shows the numbers and the alarms, and that is all.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

# Defaults for `alarms`. Overridable per settings file as
# ``agent.load_limits`` (same keys as DEFAULT_LIMITS). A login node is
# shared: these mark the point where one session stops being invisible
# to its neighbours.
DEFAULT_LIMITS = {
    "max_procs": 60,       # processes in the session's tree
    "max_rss_gb": 12.0,    # resident memory of the whole tree
    "max_cores": 4.0,      # average cores over the measurement window
    "max_d_state": 2,      # processes stuck in uninterruptible sleep
}

_TICKS_PER_SECOND = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
_MY_UID = os.getuid()


def _parse_stat(text: str) -> dict | None:
    """pid, ppid, state, cpu ticks from a /proc/<pid>/stat line.

    ``comm`` sits in parentheses and may itself contain spaces and
    parentheses, so the fields are what follows the LAST ``)`` — the
    same rule every reader of /proc uses. Unparseable → None; one bad
    line never takes down the whole picture.
    """
    close = text.rfind(")")
    if close < 0:
        return None
    try:
        fields = text[close + 2:].split()
        pid = int(text[: text.find("(")].strip())
        return {
            "pid": pid,
            "state": fields[0],
            "ppid": int(fields[1]),
            "ticks": int(fields[11]) + int(fields[12]),  # utime + stime
        }
    except (IndexError, ValueError):
        return None


def _read_stat(proc_root: Path, pid: int) -> dict | None:
    try:
        return _parse_stat(
            (proc_root / str(pid) / "stat").read_text(encoding="utf-8"))
    except OSError:
        return None


def _uid_of(proc_root: Path, pid: int) -> int | None:
    """The real uid from /proc/<pid>/status, or None when unreadable."""
    try:
        for line in (proc_root / str(pid) / "status").read_text(
                encoding="utf-8").splitlines():
            if line.startswith("Uid:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return None


def process_tree(pid: int, *, proc_root: str | Path = "/proc") -> list[int]:
    """``pid`` and every descendant of it that runs under OUR uid.

    One pass over the proc root builds pid → ppid, then the tree is
    collected from the root down. Other users' processes are never
    entered, so a foreign pid cannot pull its children into our count.
    """
    root = Path(proc_root)
    children: dict[int, list[int]] = {}
    try:
        entries = list(root.iterdir())
    except OSError:
        return []
    pids: set[int] = set()
    for entry in entries:
        if not entry.name.isdigit():
            continue
        other = _read_stat(root, int(entry.name))
        if other is None:
            continue
        pids.add(other["pid"])
        children.setdefault(other["ppid"], []).append(other["pid"])
    out: list[int] = []
    stack = [int(pid)]
    seen = {int(pid)}
    while stack:
        current = stack.pop()
        if current in pids and _uid_of(root, current) == _MY_UID:
            out.append(current)
        for child in children.get(current, ()):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    out.sort()
    return out


def _rss_bytes(proc_root: Path, pid: int) -> int:
    """VmRSS of one process in bytes, 0 when unreadable."""
    try:
        for line in (proc_root / str(pid) / "status").read_text(
                encoding="utf-8").splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024  # kB → bytes
    except (OSError, ValueError, IndexError):
        pass
    return 0


def load_of(pid: int, previous: dict | None = None, *,
            proc_root: str | Path = "/proc") -> dict:
    """What the process tree of ``pid`` costs right now.

    ``previous`` is an earlier ``load_of`` result: CPU cores are the
    growth of the tree's CPU time divided by the elapsed wall time, so
    one core fully busy reads 1.0. Without a previous measurement there
    is no window and ``cpu_cores`` is None — no invented number.

    Returns ``{"procs", "rss_bytes", "cpu_cores", "d_state",
    "cpu_seconds", "monotonic"}``; the last two are for the NEXT call's
    ``previous``.
    """
    root = Path(proc_root)
    tree = process_tree(pid, proc_root=root)
    rss = 0
    d_state = 0
    ticks = 0
    for one in tree:
        info = _read_stat(root, one)
        if info is None:
            continue
        rss += _rss_bytes(root, one)
        ticks += info["ticks"]
        if info["state"] == "D":
            d_state += 1
    now = time.monotonic()
    cpu_seconds = ticks / _TICKS_PER_SECOND
    cores = None
    if previous:
        elapsed = now - float(previous.get("monotonic") or now)
        before = float(previous.get("cpu_seconds") or cpu_seconds)
        if elapsed > 0:
            cores = (cpu_seconds - before) / elapsed
    return {
        "procs": len(tree),
        "rss_bytes": rss,
        "cpu_cores": cores,
        "d_state": d_state,
        "cpu_seconds": cpu_seconds,
        "monotonic": now,
    }


def alarms(load: dict, limits: dict | None = None) -> list[str]:
    """English findings for a load that outgrew its limits. Empty = quiet.

    Pure reporting: the caller decides what a finding means. ``limits``
    overrides the defaults (``DEFAULT_LIMITS`` keys); unknown keys are
    ignored rather than fatal, so a settings file may grow new knobs.
    """
    caps = dict(DEFAULT_LIMITS)
    for key, value in (limits or {}).items():
        if key in caps:
            caps[key] = value
    found: list[str] = []
    procs = int(load.get("procs") or 0)
    if procs > caps["max_procs"]:
        found.append(f"{procs} processes (limit {caps['max_procs']})")
    rss_gb = (load.get("rss_bytes") or 0) / (1 << 30)
    if rss_gb > caps["max_rss_gb"]:
        found.append(f"{rss_gb:.1f} GB resident (limit "
                     f"{caps['max_rss_gb']:.0f} GB)")
    cores = load.get("cpu_cores")
    if cores is not None and cores > caps["max_cores"]:
        found.append(f"{cores:.1f} cores on average "
                     f"(limit {caps['max_cores']:.0f})")
    d_state = int(load.get("d_state") or 0)
    if d_state > caps["max_d_state"]:
        found.append(f"{d_state} processes in D state "
                     f"(limit {caps['max_d_state']})")
    return found

