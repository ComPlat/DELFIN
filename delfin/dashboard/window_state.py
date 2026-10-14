"""Which dashboard kernels have no window, as the server sees it.

Input: the server's connection count per kernel (``notify_connect`` /
``notify_disconnect`` in resume_server's manager). Output: one small JSON
file per kernel that currently has no window, removed when a window
connects or the kernel ends.

Why a file: the kernels are separate processes and cannot see the server's
count. A reload starts a NEW kernel while the server keeps the old one for
the grace; the old one still holds the conversation's writer lock, and the
new one has to know that the holder is a kernel nobody is looking at
before it may take the conversation over (``may_take_over``). The server
is the only party that knows, so it writes it down.

Rejected: a heartbeat from the page into the kernel (the page cannot beat
once it is gone, and a missed beat is indistinguishable from a slow one);
asking the server over HTTP from the kernel (needs the token and the base
url in every kernel, for one bit).
"""

from __future__ import annotations

import json
import os
import socket
import time

#: Where the server writes the windowless kernels.
RECORD_DIR = os.environ.get("DELFIN_WINDOWLESS_DIR") or os.path.join(
    os.path.expanduser("~"), ".delfin", "windowless_kernels")


def _hostname() -> str:
    try:
        return socket.gethostname()
    except OSError:
        return ""


def _safe(kid: str) -> str:
    return "".join(c for c in str(kid) if c.isalnum() or c in "-_")[:80]


def _path(kid: str, root: str = "") -> str:
    return os.path.join(root or RECORD_DIR, f"{_safe(kid)}.json")


def mark(kid: str, windowless: bool, *, root: str = "") -> None:
    """Record that kernel ``kid`` has (``windowless``) no window, or drop
    the record. Never raises: a lost record costs a reload its takeover,
    which then refuses as before."""
    if not _safe(kid):
        return
    path = _path(kid, root)
    if not windowless:
        try:
            os.unlink(path)
        except OSError:
            pass
        return
    try:
        directory = root or RECORD_DIR
        os.makedirs(directory, mode=0o700, exist_ok=True)
        tmp = f"{path}.{os.getpid()}.tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump({"kernel_id": str(kid), "host": _hostname(),
                       "since": time.time()}, handle)
        os.replace(tmp, path)
    except OSError:
        pass


def is_windowless(kid: str, *, root: str = "") -> bool:
    """Whether the server on this host recorded ``kid`` as having no window.

    A record from another host, or one owned by another user, is not
    evidence: it says nothing about a kernel this host could hand over."""
    if not _safe(kid):
        return False
    try:
        with open(_path(kid, root), encoding="utf-8") as handle:
            if os.fstat(handle.fileno()).st_uid != os.getuid():
                return False
            record = json.load(handle)
    except (OSError, ValueError):
        return False
    return (isinstance(record, dict)
            and str(record.get("kernel_id") or "") == str(kid)
            and str(record.get("host") or "") == _hostname())


def may_take_over(holder: dict | None) -> bool:
    """Whether a new page may take a conversation from its lock ``holder``.

    Input: ``session_store.session_lock_holder``'s answer. Output: True only
    when the holder is a dashboard kernel on this host (it named its kernel
    id in the lock), the server recorded that kernel as windowless, and no
    turn is running in it. Everything else -- a CLI, another login node, a
    kernel with a window, a kernel mid-turn -- keeps the conversation.
    """
    if not holder or holder.get("elsewhere"):
        return False
    kid = str(holder.get("kernel") or "")
    if not kid or not is_windowless(kid):
        return False
    try:
        from delfin.dashboard import turn_record
        if kid in turn_record.running_kernel_ids():
            return False
    except Exception:
        return False
    return True
