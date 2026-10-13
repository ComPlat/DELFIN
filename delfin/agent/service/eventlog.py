"""The durable record of one session, for the windows that attach later.

Every event the worker emits gets the next sequence number; the durable
ones (protocol.DURABLE_EVENTS) are also appended to a JSONL file. A window
that attaches to a session already in progress asks for the last few turns
and gets exactly what a window watching from the start would have kept:
finished messages, tool calls with their results, notices, approvals.

Written by one process only (the session's worker); others may read at any
time. A line is written with one call and readers drop a last line that is
not finished yet, so a reader never sees half an event.
"""

from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from typing import Any

from .protocol import DURABLE_EVENTS

SERVICE_DIR_ENV = "DELFIN_AGENT_SERVICE_DIR"


def service_dir() -> Path:
    """Where the service keeps its state: run files, logs, session records."""
    override = os.environ.get(SERVICE_DIR_ENV, "").strip()
    if override:
        return Path(override).expanduser()
    return Path.home() / ".delfin" / "agent_service"


def log_path(session_id: str) -> Path:
    safe = "".join(c for c in str(session_id) if c.isalnum() or c in "-_")
    if not safe:
        raise ValueError(f"not a usable session id: {session_id!r}")
    return service_dir() / "sessions" / f"{safe}.jsonl"


def read_events(path: Path) -> list[dict]:
    """Every complete event in ``path``; a torn last line is left out."""
    try:
        raw = path.read_bytes()
    except FileNotFoundError:
        return []
    events: list[dict] = []
    for line in raw.split(b"\n"):
        if not line.strip():
            continue
        try:
            event = json.loads(line)
        except ValueError:
            continue           # a line still being written, or damage
        if isinstance(event, dict):
            events.append(event)
    return events


def tail_turns(events: list[dict], turns: int) -> list[dict]:
    """The events from the start of the ``turns``-th last turn onwards.

    The latest ``session.info`` always comes first, so a replay can draw
    the session's header even when its first turn is long gone.
    ``turns <= 0`` returns everything.
    """
    info = next((e for e in reversed(events)
                 if e.get("type") == "session.info"), None)
    if turns <= 0:
        start = 0
    else:
        starts = [i for i, e in enumerate(events)
                  if e.get("type") == "user.message"]
        start = starts[-turns] if len(starts) >= turns else 0
    picked = events[start:]
    if info is not None and info not in picked:
        picked = [info, *picked]
    return picked


class EventLog:
    """Sequence numbers for every event; a file line for the durable ones."""

    def __init__(self, path: Path):
        from .. import state_paths

        self.path = path
        state_paths.ensure_dir(path.parent)
        self._lock = threading.Lock()
        existing = read_events(path)
        self._seq = max((int(e.get("seq", 0) or 0) for e in existing),
                        default=0)

    def stamp(self, event: dict) -> dict:
        """Give ``event`` the next sequence number and keep it if durable."""
        from .. import state_paths

        with self._lock:
            self._seq += 1
            event["seq"] = self._seq
            if event.get("type") in DURABLE_EVENTS:
                line = json.dumps(event, ensure_ascii=False,
                                  separators=(",", ":")) + "\n"
                with state_paths.open_append(self.path) as fh:
                    fh.write(line)
        return event

    def replay(self, turns: int = 0) -> list[dict[str, Any]]:
        return tail_turns(read_events(self.path), turns)
