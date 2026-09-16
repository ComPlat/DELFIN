"""Which agent sessions are open, and where they work.

Several sessions can run at once: side by side in one dashboard, or in a
dashboard on another login node that shares the home directory. Each open
session keeps a small record here -- its title, working directory, and the
git repository and branch that directory is in. Sessions in one repository
read it to stay out of each other's way, and to find each other.

A record is keyed by the session's place in its session list, which exists
from the moment the session is opened; the conversation id comes later.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "session_presence"
# A record not refreshed for this long belongs to a session that is gone.
_STALE_S = 15 * 60.0
# An unchanged record is rewritten at most this often (a heartbeat), so a
# session list refreshing every few seconds does not write every few seconds.
_HEARTBEAT_S = 60.0
_GIT_TTL_S = 60.0

_git_cache: dict[str, tuple[float, dict]] = {}
_last_written: dict[str, tuple[str, float]] = {}


def _git(workspace: str, *args: str) -> str:
    try:
        out = subprocess.run(["git", "-C", workspace, *args],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def repository_of(workspace: str) -> dict:
    """``{"root", "common_dir", "branch"}`` for the git checkout ``workspace``
    is in, empty strings outside one. ``common_dir`` is shared by a
    repository and all of its worktrees; ``root`` is the checkout itself."""
    ws = str(workspace or "")
    now = time.monotonic()
    hit = _git_cache.get(ws)
    if hit and now - hit[0] < _GIT_TTL_S:
        return dict(hit[1])
    info = {"root": "", "common_dir": "", "branch": ""}
    if ws and Path(ws).is_dir():
        root = _git(ws, "rev-parse", "--show-toplevel")
        if root:
            common = _git(ws, "rev-parse", "--path-format=absolute",
                          "--git-common-dir")
            info = {
                "root": str(Path(root).resolve()),
                "common_dir": str(Path(common).resolve()) if common else "",
                "branch": _git(ws, "rev-parse", "--abbrev-ref", "HEAD"),
            }
    _git_cache[ws] = (now, info)
    return dict(info)


def _safe(key: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in key)[:80]


def _path(key: str) -> Path:
    return _DIR / f"{_safe(key)}.json"


def announce(key: str, *, session_id: str = "", title: str = "",
             workspace: str = "") -> None:
    """Say that session ``key`` is open and where it works. Never raises."""
    key = str(key or "").strip()
    if not key:
        return
    repo = repository_of(workspace)
    record = {
        "key": key,
        "session_id": str(session_id or ""),
        "title": str(title or "")[:80],
        "workspace": str(workspace or ""),
        "repo_root": repo["root"],
        "common_dir": repo["common_dir"],
        "branch": repo["branch"],
        "host": socket.gethostname(),
        "pid": os.getpid(),
    }
    body = json.dumps(record, sort_keys=True)
    last = _last_written.get(key)
    now = time.time()
    if last and last[0] == body and now - last[1] < _HEARTBEAT_S:
        return
    try:
        from .state_paths import ensure_dir, write_text
        ensure_dir(_DIR)
        write_text(_path(key), json.dumps({**record, "updated_at": now}))
        _last_written[key] = (body, now)
    except Exception:
        pass


def withdraw(key: str) -> None:
    """Say that session ``key`` has closed. Never raises."""
    _last_written.pop(str(key or ""), None)
    try:
        _path(str(key or "")).unlink()
    except Exception:
        pass


def _alive(record: dict, now: float) -> bool:
    if now - float(record.get("updated_at") or 0) > _STALE_S:
        return False
    if record.get("host") == socket.gethostname():
        pid = int(record.get("pid") or 0)
        if pid <= 0:
            return False
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except Exception:
            pass          # someone else's process: alive
    return True


# A record this much older than stale is reaped, not merely skipped. Crashed
# kernels never withdraw their record; the directory grew by one file per
# crash and was read whole on every refresh (review 2026-09-16).
_REAP_AFTER_S = 4 * _STALE_S


def _reap(now: float) -> int:
    """Remove records of sessions that are long gone. Never raises."""
    removed = 0
    try:
        files = list(_DIR.glob("*.json"))
    except Exception:
        return 0
    for f in files:
        try:
            record = json.loads(f.read_text(encoding="utf-8"))
            updated = float((record or {}).get("updated_at") or 0)
            dead_here = (isinstance(record, dict)
                         and record.get("host") == socket.gethostname()
                         and not _alive(record, now))
            if now - updated > _REAP_AFTER_S or dead_here:
                f.unlink()
                removed += 1
        except FileNotFoundError:
            continue
        except Exception:
            continue
    return removed


def open_sessions(*, exclude_key: str = "") -> list[dict]:
    """The records of the sessions that are open, other than ``exclude_key``."""
    now = time.time()
    out: list[dict] = []
    _reap(now)
    try:
        files = sorted(_DIR.glob("*.json"))
    except Exception:
        return out
    for f in files:
        try:
            record = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        if (isinstance(record, dict) and record.get("key") != exclude_key
                and _alive(record, now)):
            out.append(record)
    return out


def in_same_repository(workspace: str, *, exclude_key: str = "") -> list[dict]:
    """Open sessions working in the repository ``workspace`` is in -- in the
    same checkout or in another worktree of it."""
    common = repository_of(workspace)["common_dir"]
    if not common:
        return []
    return [r for r in open_sessions(exclude_key=exclude_key)
            if r.get("common_dir") == common]
