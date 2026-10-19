"""A read of a file the model already has, unchanged and still in view.

Input: each read_file call (its arguments and the workspace), the file's
size and modification time, the result row the loop appended, and at the
next read of the same range the request list about to be sent. Output:
:meth:`ReadLedger.repeat_note` -- a short answer instead of the file, or
None when the file must be read.

The short answer is given only when all three hold:

* the same file and the same lines (``offset``/``limit``) were read;
* the file's size and nanosecond mtime are what they were then -- any
  write, by a tool, the shell or the user, changes them;
* the earlier result row is still in the request list, byte for byte. A
  compaction, a context edit or a new session that dropped or rewrote it
  makes the model read the file again, so the note can never point at
  text the model no longer has.

Measured 2026-10-19 (session 78c36f8c): a 441-line HTML file read whole
again and again in one session, the requests growing past 400k input
tokens; with the rows already there, every repeat was paid twice.

Pure apart from ``os.stat``; never raises.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any, Optional


def _key(args: Any, workspace: Any) -> Optional[tuple]:
    if not isinstance(args, dict):
        return None
    raw = str(args.get("path") or args.get("file_path") or "")
    if not raw:
        return None
    try:
        p = Path(raw).expanduser()
        if not p.is_absolute():
            p = Path(str(workspace or ".")) / p
        resolved = str(p.resolve())
    except Exception:
        return None
    return (resolved, str(args.get("offset") or 0), str(args.get("limit") or ""))


def _stamp(path: str) -> Optional[tuple]:
    try:
        st = os.stat(path)
        return (st.st_size, st.st_mtime_ns)
    except OSError:
        return None


def _digest(content: Any) -> str:
    return hashlib.sha256(str(content).encode("utf-8", "replace")).hexdigest()


class ReadLedger:
    """The reads whose results the conversation may still hold."""

    def __init__(self) -> None:
        self._seen: dict[tuple, tuple] = {}

    def record(self, args: Any, workspace: Any, call_id: str,
               content: Any) -> None:
        """A read_file result row was appended under *call_id*."""
        try:
            key = _key(args, workspace)
            stamp = _stamp(key[0]) if key else None
            if key and stamp and call_id:
                self._seen[key] = (call_id, stamp, _digest(content))
        except Exception:
            pass

    def repeat_note(self, args: Any, workspace: Any,
                    messages: Any) -> Optional[str]:
        """The short answer for a read the model already has, else None."""
        try:
            key = _key(args, workspace)
            entry = self._seen.get(key) if key else None
            if entry is None:
                return None
            call_id, stamp, digest = entry
            if _stamp(key[0]) != stamp:
                return None
            for m in reversed(list(messages or [])):
                if (isinstance(m, dict) and m.get("role") == "tool"
                        and m.get("tool_call_id") == call_id):
                    if _digest(m.get("content")) != digest:
                        return None
                    shown = str((args or {}).get("path")
                                or (args or {}).get("file_path") or key[0])
                    return (
                        f"{shown} is unchanged since you read these lines "
                        "earlier in this conversation; that result is still "
                        "above -- work from it. For other lines pass "
                        "offset/limit; to find one, use grep_file.")
            return None
        except Exception:
            return None


__all__ = ["ReadLedger"]
