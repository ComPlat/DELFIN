"""Pause/resume: a per-session flag the engine checks between tool calls.

Wave-12 finding 2: there is no pause — SIGSTOP on the agent process has no
effect; the only stop was Esc + /exit, and even then a partner message woke
the session again. The fix is a DELFIN-level pause: a flag file for a session,
checked by the engine between tool calls. While the flag is present, no tool
call may start and no wake may fire.

This module owns ONLY the state: the flag-file path, the atomic create/remove,
and the two decision functions the engine/caller ask. The engine wiring (check
between tool calls) and the ``delfin-agent pause|resume`` CLI commands live in
api_client.py / cli.py, which are protected — they go to the operator as a
.gate patch. Everything here is a pure file operation, testable with a temp
dir.

Keys are the session address used elsewhere (repl._presence_key): the ``-n``
name, else the first 8 chars of the session id. ``pause`` is idempotent and
atomic (tmp + os.replace, under the same cross-process lock the message inbox
uses); ``resume`` only ever removes this session's own flag.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "session_pause"


def _path(key: str) -> Path:
    """The flag file for ``key``, sanitised so it can never escape _DIR.

    Same rule as the session message inbox: only alphanumerics and ``-_.``
    survive, everything else becomes ``_``. A key like ``../../escape``
    therefore becomes one single filename under _DIR, never a path that
    climbs out of it.
    """
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(key))[:80]
    return _DIR / f"{safe}.pause.json"


def is_paused(key: str) -> bool:
    """Whether the session ``key`` is paused right now."""
    try:
        return _path(key).exists()
    except OSError:
        return False


def pause_gate(key: str) -> bool:
    """Whether a tool call may START for ``key`` (False while paused)."""
    return not is_paused(key)


def wake_blocked(key: str) -> bool:
    """Whether the wake look must NOT fire for ``key`` (True while paused)."""
    return is_paused(key)


def pause(key: str, reason: str = "") -> dict:
    """Mark session ``key`` paused (idempotent, atomic).

    Writes the flag under the cross-process lock so two pauses cannot tear
    each other. Returns the flag record. Reason is for the operator's
    benefit when it lists pauses; it is not read by the engine.
    """
    path = _path(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {"key": str(key), "reason": str(reason or ""),
              "paused_at": time.time()}
    blob = (json.dumps(record) + "\n").encode("utf-8")
    from .bash_jobs import cross_process_lock
    with cross_process_lock(path):
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "wb") as fh:
            fh.write(blob)
        os.replace(tmp, path)
    return record


def resume(key: str) -> bool:
    """Un-pause session ``key``. Returns True if a flag was removed."""
    path = _path(key)
    try:
        from .bash_jobs import cross_process_lock
        with cross_process_lock(path):
            try:
                os.remove(path)
                return True
            except FileNotFoundError:
                return False
    except OSError:
        return False


def paused_keys() -> list[str]:
    """Every currently paused session key, sorted."""
    try:
        keys = []
        for child in _DIR.iterdir():
            if child.is_file() and child.name.endswith(".pause.json"):
                key = child.name[:-len(".pause.json")]
                keys.append(key)
        return sorted(keys)
    except FileNotFoundError:
        return []
