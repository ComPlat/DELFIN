"""The durable record of one session, for the windows that attach later.

The durable events (protocol.DURABLE_EVENTS) get the next sequence number
and are appended to a JSONL file. A window that attaches to a session
already in progress asks for the last few turns and gets exactly what a
window watching from the start would have kept: finished messages, tool
calls with their results, notices, approvals.

Sequence numbers belong to durable events only, so they are never reused:
a worker that resumes a session continues after the highest number in the
file, and every number below it is in the file too. Live-only events
(deltas, usage samples, state changes) carry none -- nothing replays them,
so nothing has to tell them apart.

Written by one process only (the session's worker); others may read at any
time. A line is written with one call and readers drop a line that is not
complete JSON, so a reader never sees half an event.
"""

from __future__ import annotations

import json
import os
import threading
from collections import deque
from pathlib import Path
from typing import Any

from .protocol import DURABLE_EVENTS, encode

SERVICE_DIR_ENV = "DELFIN_AGENT_SERVICE_DIR"

#: The most event bytes a replay carries (an attach, a create, a
#: session.replay): the newest whole turns that fit. Every reader bounds
#: its lines (protocol.MAX_LINE_BYTES), so an answer must stay well below.
REPLAY_BUDGET_BYTES = 16 * 1024 * 1024


def service_dir() -> Path:
    """Where the service keeps its state: run files, logs, session records.

    Follows ``DELFIN_SCRATCH_STATE`` the way every other ``~/.delfin`` sink
    does (state_paths), so a probe or a test run keeps its records out of
    the user's real home.
    """
    override = os.environ.get(SERVICE_DIR_ENV, "").strip()
    if override:
        return Path(override).expanduser()
    scratch = os.environ.get("DELFIN_SCRATCH_STATE", "").strip()
    home = Path(scratch).expanduser() if scratch else Path.home()
    return home / ".delfin" / "agent_service"


def log_path(session_id: str) -> Path:
    safe = "".join(c for c in str(session_id) if c.isalnum() or c in "-_")
    if not safe:
        raise ValueError(f"not a usable session id: {session_id!r}")
    return service_dir() / "sessions" / f"{safe}.jsonl"


def read_events(path: Path) -> list[dict]:
    """Every complete event in ``path``; a torn or damaged line is left out."""
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


def tail_turns(events: list[dict], turns: int,
               header: dict | None = None, max_bytes: int = 0) -> list[dict]:
    """The events from the start of the ``turns``-th last turn onwards.

    A turn starts at its ``turn.started`` -- together with the messages
    sent right before it, which are what started it. A message sent while
    an earlier turn ran (queued, steered) stays inside that earlier turn,
    where it happened, and never counts as a turn of its own.

    A session header comes first, so a replay can draw it even when the
    first turn is long gone: ``header`` when given (the session as it is
    now), else the latest ``session.info`` in the log. ``turns <= 0``
    returns everything.

    ``max_bytes`` bounds the events (the header aside): the newest whole
    turns that fit, with a notice in front when older ones were left out.
    Every reader bounds its lines, so an answer that carries a replay must
    stay well below that.
    """
    events = [e for e in events if e.get("type") != "session.info"] \
        if header is not None else events
    info = header if header is not None else next(
        (e for e in reversed(events) if e.get("type") == "session.info"), None)
    start = 0
    if turns > 0:
        starts = [i for i, e in enumerate(events)
                  if e.get("type") == "turn.started"]
        if len(starts) >= turns:
            start = starts[-turns]
            while start > 0 and events[start - 1].get("type") == "user.message":
                start -= 1
    picked = events[start:]
    if max_bytes > 0:
        window = RecentEvents(max_bytes)
        for event in picked:
            if event is not info:
                window.add(event)
        if window.dropped:
            picked = [OLDER_LEFT_OUT, *window.events()]
    if info is not None:
        picked = [info, *[e for e in picked if e is not info]]
    return picked


#: Where a replay leaves out turns that did not fit its budget.
OLDER_LEFT_OUT = {"type": "notice", "ts": 0,
                  "text": "(older turns of this session are not shown here; "
                          "`delfin-agent report` has them all)"}


class RecentEvents:
    """The newest durable events of one session within a byte budget,
    beginning at the start of a turn where one fits: what an attach
    replays, kept in memory so no attach reads the log in the event loop."""

    def __init__(self, budget: int):
        self.budget = int(budget)
        self._events: deque[tuple[dict, int]] = deque()
        self._bytes = 0
        #: Older events were let go to stay within the budget.
        self.dropped = False

    def add(self, event: dict) -> None:
        size = len(encode(event))
        self._events.append((event, size))
        self._bytes += size
        self._trim()

    def merge_older(self, older: list[dict]) -> None:
        """Put the events that came before the first one held in front
        (a resumed session's log), newest kept when they do not all fit."""
        held = [e for e, _ in self._events]
        cutoff = int(held[0].get("seq", 0) or 0) if held else 0

        def earlier(event: dict) -> bool:
            seq = int(event.get("seq", 0) or 0)
            return seq > 0 and (not held or seq < cutoff)

        before = sorted(filter(earlier, older), key=lambda e: int(e.get("seq", 0) or 0))
        if not before:
            return
        sized = [(e, len(encode(e))) for e in before]
        self._events.extendleft(reversed(sized))
        self._bytes += sum(size for _, size in sized)
        self._trim()

    def _trim(self) -> None:
        if self._bytes <= self.budget:
            return
        events = self._events
        while self._bytes > self.budget and len(events) > 1:
            _, size = events.popleft()
            self._bytes -= size
            self.dropped = True
        # Start at a turn (its turn.started, with the messages that sent
        # it): a replay that opens mid-turn shows a tool row without the
        # question that led to it. A single turn larger than the budget is
        # kept from wherever it fits.
        start = self._first_turn_start()
        for _ in range(start):
            _, size = events.popleft()
            self._bytes -= size

    def _first_turn_start(self) -> int:
        messages = 0
        for i, (event, _) in enumerate(self._events):
            kind = event.get("type")
            if kind == "turn.started":
                return i - messages
            messages = messages + 1 if kind == "user.message" else 0
        return 0

    def events(self) -> list[dict]:
        return [e for e, _ in self._events]

    def turns(self) -> int:
        return sum(1 for e, _ in self._events if e.get("type") == "turn.started")


class EventLog:
    """Sequence numbers for the durable events, and a file line for each."""

    def __init__(self, path: Path, on_error=None):
        from .. import state_paths

        self.path = path
        self._lock = threading.Lock()
        self._on_error = on_error
        self._fd: int | None = None
        self.writable = True
        try:
            # Every level owner-only, not just the last one.
            for parent in reversed(path.parents[:3]):
                state_paths.ensure_dir(parent)
            self._fd = os.open(str(path),
                               os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
            state_paths.secure_file(path)
            self._end_torn_line()
        except OSError as exc:
            self._fail(exc)
        existing = read_events(path)
        self._seq = max((int(e.get("seq", 0) or 0) for e in existing),
                        default=0)
        #: Turns the log already holds: a resumed session numbers on.
        self.turns = sum(1 for e in existing if e.get("type") == "turn.started")

    def _end_torn_line(self) -> None:
        """A worker killed mid-write leaves a last line without its newline;
        the next event appended to it would be lost with it. End it first."""
        size = os.fstat(self._fd).st_size
        if size == 0:
            return
        with open(self.path, "rb") as fh:
            fh.seek(size - 1)
            if fh.read(1) != b"\n":
                self._write(b"\n")

    def _write(self, data: bytes) -> None:
        view = memoryview(data)
        while view:
            n = os.write(self._fd, view)
            view = view[n:]

    def _fail(self, exc: Exception) -> None:
        """A log that cannot be written costs the replay, never the session:
        a full disk or an HPC home over quota must not end the turn."""
        if self.writable and self._on_error is not None:
            try:
                self._on_error(exc)
            except Exception:
                pass
        self.writable = False

    def stamp(self, event: dict) -> dict:
        """Number ``event`` and keep it, when it is durable. Never raises."""
        if event.get("type") not in DURABLE_EVENTS:
            return event
        with self._lock:
            self._seq += 1
            event["seq"] = self._seq
            if self.writable and self._fd is not None:
                try:
                    self._write(encode(event))
                except (OSError, ValueError) as exc:
                    self._fail(exc)
        return event

    def replay(self, turns: int = 0, header: dict | None = None,
               max_bytes: int = 0) -> list[dict[str, Any]]:
        return tail_turns(read_events(self.path), turns, header, max_bytes)

    def close(self) -> None:
        with self._lock:
            if self._fd is not None:
                try:
                    os.close(self._fd)
                except OSError:
                    pass
                self._fd = None
