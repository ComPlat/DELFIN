"""Is any process using a worktree? Asked of the processes, not of a registry.

Input: directory paths. Output of :func:`in_use`: ``{path: why}`` for every
path that a process of this user uses -- its working directory, its
executable, an open file or a mapped file lies under the path. ``why``
names the process (pid and command) and the kind of use.

Sources, in this order:

1. ``/proc`` (Linux): ``cwd``, ``exe`` and ``fd/*`` links and the ``maps``
   file of every process owned by this user. One pass over ``/proc`` for
   all paths together.
2. ``lsof`` (hosts without ``/proc``, e.g. macOS): ``lsof -n -w -Fpc +D
   <path>`` per path, bounded by a timeout; any process listed uses it.
3. Neither: a file modified within ``recent_s`` (one hour) counts as use,
   from a walk bounded to ``max_files`` entries.

Semantics: an unreadable process (another user's, or one that ended during
the scan) is skipped, so a worktree is reported as used only on evidence;
the callers (Clean up, Delete branch, Close preview) keep a worktree that
is used and ask again at the moment they would remove it.

Rejected: DELFIN's own job registry as the source (a benchmark or a shell
started from a terminal registers nowhere; measured on the development
machine, a running benchmark's worktree was offered for removal); ``fuser``
(not on every host, and it does not report the working directory on all
of them).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable, Optional

RECENT_S = 3600.0
MAX_FILES = 20000


def _norm(path) -> str:
    try:
        return str(Path(path).resolve())
    except OSError:
        return os.path.abspath(str(path))


def _owner(path: str, roots: list[str]) -> Optional[str]:
    """The root *path* lies under (or is), or None."""
    if not path:
        return None
    path = path.split(" (deleted)")[0]
    for root in roots:
        if path == root or path.startswith(root.rstrip("/") + "/"):
            return root
    return None


def _readlink(p: str) -> str:
    try:
        return os.readlink(p)
    except OSError:
        return ""


def _comm(proc_dir: str) -> str:
    try:
        return Path(proc_dir, "comm").read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return "?"


def _from_proc(roots: list[str], proc_root: str,
               pids: Optional[Iterable[str]] = None) -> dict[str, str]:
    found: dict[str, str] = {}
    uid = os.getuid() if hasattr(os, "getuid") else None
    try:
        names = list(pids) if pids is not None else os.listdir(proc_root)
    except OSError:
        return found
    for pid in sorted((n for n in names if n.isdigit()), key=int):
        d = os.path.join(proc_root, pid)
        try:
            if uid is not None and os.stat(d).st_uid != uid:
                continue
        except OSError:
            continue
        hits: list[tuple[str, str]] = []
        for kind, link in (("working directory", "cwd"), ("executable", "exe")):
            root = _owner(_readlink(os.path.join(d, link)), roots)
            if root:
                hits.append((root, kind))
        try:
            fds = os.listdir(os.path.join(d, "fd"))
        except OSError:
            fds = []
        for fd in fds:
            root = _owner(_readlink(os.path.join(d, "fd", fd)), roots)
            if root:
                hits.append((root, "open file"))
                break
        try:
            with open(os.path.join(d, "maps"), encoding="utf-8",
                      errors="replace") as fh:
                for line in fh:
                    parts = line.split(None, 5)
                    if len(parts) == 6:
                        root = _owner(parts[5].strip(), roots)
                        if root:
                            hits.append((root, "mapped file"))
                            break
        except OSError:
            pass
        for root, kind in hits:
            found.setdefault(root, f"in use by process {pid} ({_comm(d)}): {kind}")
    return found


def _from_lsof(roots: list[str], lsof: str, timeout: float = 30.0) -> dict[str, str]:
    found: dict[str, str] = {}
    for root in roots:
        try:
            r = subprocess.run([lsof, "-n", "-w", "-Fpc", "+D", root],
                               capture_output=True, text=True, timeout=timeout,
                               check=False)
        except subprocess.TimeoutExpired:
            found[root] = "in use: lsof did not finish, so it is treated as used"
            continue
        except OSError:
            continue
        pid, cmd = "", ""
        for line in (r.stdout or "").splitlines():
            if line.startswith("p") and not pid:
                pid = line[1:]
            elif line.startswith("c") and not cmd:
                cmd = line[1:]
        if pid:
            found[root] = f"in use by process {pid} ({cmd or '?'}): lsof"
    return found


def _from_mtime(roots: list[str], recent_s: float, now: float,
                max_files: int) -> dict[str, str]:
    found: dict[str, str] = {}
    for root in roots:
        seen = 0
        for dirpath, dirnames, filenames in os.walk(root):
            if ".git" in dirnames:
                dirnames.remove(".git")
            for name in filenames:
                seen += 1
                try:
                    m = os.lstat(os.path.join(dirpath, name)).st_mtime
                except OSError:
                    continue
                if now - m < recent_s:
                    found[root] = (f"modified {int((now - m) // 60)} min ago "
                                   "(no /proc and no lsof to ask the processes)")
                    break
            if root in found:
                break
            if seen > max_files:
                found[root] = ("too large to check without /proc or lsof, so "
                               "it is treated as used")
                break
    return found


def in_use(paths: Iterable, *, proc_root: str = "/proc", lsof: Optional[str] = "auto",
           recent_s: float = RECENT_S, now: Optional[float] = None,
           max_files: int = MAX_FILES,
           pids: Optional[Iterable[str]] = None) -> dict[str, str]:
    """``{path: why}`` for the paths some process uses (module docstring).

    Keys are the paths as given."""
    given = [str(p) for p in paths]
    norm = {_norm(p): p for p in given}
    roots = list(norm)
    if not roots:
        return {}
    if os.path.isdir(proc_root):
        hits = _from_proc(roots, proc_root, pids)
    else:
        tool = shutil.which("lsof") if lsof == "auto" else lsof
        if tool:
            hits = _from_lsof(roots, tool)
        else:
            hits = _from_mtime(roots, recent_s,
                               time.time() if now is None else now, max_files)
    return {norm[r]: why for r, why in hits.items()}
