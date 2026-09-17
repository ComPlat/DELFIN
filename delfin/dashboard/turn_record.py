"""The record that says a turn is running in this kernel.

The server culls a dashboard kernel that has had no window for the
grace. That rule reads windows on purpose: a kernel is idle to Jupyter
while its agent works in a thread, and busy with widget traffic while
nobody watches, so the server's own activity signal says nothing here.

Windows alone are not enough either. A browser tab that is closed -- or
a WebSocket that a ping timeout drops on a slow link -- takes the window
away while the agent is mid-turn, and the grace then ended the kernel
with the work still in flight (seen on 2026-09-17: a kernel ended 99
seconds after its window went, in the middle of a run).

So the kernel leaves a record here while a turn runs, and the culler
reads it. One small JSON file per kernel, refreshed by a thread while
the turn lasts and removed when it ends, so a kernel that died mid-turn
cannot hold itself alive: the record is honoured only while it is fresh
and its process is still there.

Nothing here keeps a kernel alive past its turn. Once the turn ends the
record goes and the ordinary grace applies, so an unwatched session
still stops on its own.
"""

from __future__ import annotations

import json
import os
import socket
import threading
import time
from typing import Optional

#: Where a running turn announces itself.
RECORD_DIR = os.environ.get("DELFIN_TURN_RECORD_DIR") or os.path.join(
    os.path.expanduser("~"), ".delfin", "running_turns")

#: How often the record is refreshed while the turn lasts.
TOUCH_SECONDS = 30.0

#: How old a record may be and still count. Four touches of headroom:
#: a kernel busy in a tight loop can be late, and the cost of counting a
#: stale record is one more poll of the grace, not a kernel kept forever.
FRESH_SECONDS = 120.0

_lock = threading.Lock()
#: kernel id -> (stop event, thread), for the turns this process runs.
_touchers: dict[str, tuple] = {}


def _hostname() -> str:
    try:
        return socket.gethostname()
    except OSError:
        return ""


def record_path(kid: str, *, root: str = "") -> str:
    return os.path.join(root or RECORD_DIR, f"{kid}.json")


def _write(kid: str, root: str) -> str:
    directory = root or RECORD_DIR
    payload = {
        "kernel_id": kid,
        "pid": os.getpid(),
        "host": _hostname(),
        "updated_at": time.time(),
    }
    try:
        os.makedirs(directory, exist_ok=True)
        path = record_path(kid, root=directory)
        tmp = f"{path}.{os.getpid()}.tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
        os.replace(tmp, path)
    except OSError:
        # A record that cannot be written costs the turn its protection,
        # never the turn itself.
        return ""
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass
    return path


def _remove(kid: str, root: str) -> None:
    try:
        os.unlink(record_path(kid, root=root))
    except OSError:
        pass


def _keep_fresh(kid: str, root: str, stop: threading.Event) -> None:
    while not stop.wait(TOUCH_SECONDS):
        if not _write(kid, root):
            return


def mark(on: bool, *, kid: str = "", root: str = "") -> str:
    """Say whether a turn is running in this kernel.

    Returns the record path while a turn runs, "" otherwise -- including
    outside a kernel, where there is no id to announce and the culler
    this speaks to does not exist.
    """
    ident = kid or _kernel_id()
    if not ident:
        return ""
    with _lock:
        entry = _touchers.pop(ident, None)
        if entry is not None:
            entry[0].set()
        if not on:
            _remove(ident, root)
            return ""
        path = _write(ident, root)
        if not path:
            return ""
        stop = threading.Event()
        thread = threading.Thread(
            target=_keep_fresh, args=(ident, root, stop),
            name=f"delfin-turn-{ident[:8]}", daemon=True)
        _touchers[ident] = (stop, thread)
        thread.start()
        return path


def _kernel_id() -> str:
    from delfin.dashboard import session as _session

    try:
        return _session.kernel_id()
    except Exception:
        return ""


def _still_there(record: dict) -> bool:
    """Whether the process that left this record is still running.

    Only answerable on the machine that wrote it; elsewhere the record's
    own age is the whole answer.
    """
    if str(record.get("host") or "") != _hostname():
        return True
    try:
        pid = int(record.get("pid") or 0)
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
        return True
    return True


def running_kernel_ids(*, root: str = "", now: Optional[float] = None) -> set[str]:
    """Every kernel whose record says a turn is running right now."""
    directory = root or RECORD_DIR
    moment = time.time() if now is None else now
    found: set[str] = set()
    try:
        names = os.listdir(directory)
    except OSError:
        return found
    for name in names:
        if not name.endswith(".json"):
            continue
        try:
            with open(os.path.join(directory, name), encoding="utf-8") as handle:
                record = json.load(handle)
        except (OSError, ValueError):
            continue
        if not isinstance(record, dict):
            continue
        kid = str(record.get("kernel_id") or "")
        if not kid:
            continue
        try:
            updated = float(record.get("updated_at") or 0.0)
        except (TypeError, ValueError):
            continue
        if moment - updated > FRESH_SECONDS:
            continue
        if not _still_there(record):
            continue
        found.add(kid)
    return found


def _reset_for_tests() -> None:
    with _lock:
        for stop, _thread in _touchers.values():
            stop.set()
        _touchers.clear()
