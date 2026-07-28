"""Schedule and cron-style recurring agent invocations.

Two kinds of entries:

  - ``once`` — fires after ``delay_seconds`` and self-deletes.
    Persistent scheduling primitive — ScheduleWakeup.
  - ``interval`` — fires every ``every_seconds`` until removed.
    The minimal cron substitute (full cron expression parsing is
    overkill for the in-dashboard use cases).

Persistence lives at ``~/.delfin/cron.json`` so entries survive
restart. The Scheduler runs a background thread that polls every
30 s and fires due entries via a user-supplied ``fire_callback``;
the callback is set by the dashboard once it knows how to dispatch
into the agent loop. Without a callback, due entries are still
recorded as overdue but not executed — they fire as soon as a
callback is bound. Headless execution without any dashboard is
provided by :mod:`delfin.agent.scheduler_daemon`, which passes its
own fire callback straight into :meth:`Scheduler.tick`.

Failures are isolated: a buggy callback never crashes the
scheduler thread. Consecutive callback failures are counted per
entry; after ``_MAX_CONSECUTIVE_FAILURES`` the entry is disabled
(persisted, with a reason) instead of retrying forever.
"""

from __future__ import annotations

import json
import os
import threading
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Optional


_DEFAULT_PATH = Path.home() / ".delfin" / "cron.json"
_POLL_S = 30.0
_MAX_CONSECUTIVE_FAILURES = 3


class DisableEntry(Exception):
    """Raised by a fire callback to disable the entry immediately.

    ``tick`` persists the ``disabled`` flag + reason instead of counting
    the raise as an ordinary (retryable) failure. Used by the headless
    daemon, e.g. when an entry's workspace no longer exists.
    """

    def __init__(self, reason: str = ""):
        super().__init__(reason)
        self.reason = reason


@dataclass
class ScheduleEntry:
    id: str
    kind: str                       # "once" | "interval"
    prompt: str
    reason: str = ""
    delay_seconds: int = 0          # used by "once"
    every_seconds: int = 0          # used by "interval"
    created_at: float = field(default_factory=time.time)
    next_fire_at: float = 0.0
    last_fired_at: float = 0.0
    fire_count: int = 0
    # Headless-execution metadata (scheduler_daemon):
    workspace: str = ""             # directory the prompt runs in headlessly
    budget_usd: float = 0.0         # optional per-run cost cap (0 = default)
    fail_count: int = 0             # consecutive fire failures
    disabled: bool = False          # never fired while set (persisted)
    disabled_reason: str = ""

    def is_due(self, now: float | None = None) -> bool:
        return (now or time.time()) >= self.next_fire_at

    def reschedule(self, fired_at: float | None = None) -> bool:
        """Update next_fire_at after a fire. Return True if still active."""
        fired_at = fired_at or time.time()
        self.last_fired_at = fired_at
        self.fire_count += 1
        if self.kind == "once":
            return False
        if self.kind == "interval" and self.every_seconds > 0:
            self.next_fire_at = fired_at + self.every_seconds
            return True
        return False


class Scheduler:
    def __init__(self, path: Path | None = None):
        self.path = path or _DEFAULT_PATH
        self._entries: dict[str, ScheduleEntry] = {}
        self._lock = threading.RLock()
        self._fire_callback: Optional[Callable[[ScheduleEntry], None]] = None
        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._load()

    # --- persistence ------------------------------------------------------

    def _load(self) -> None:
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return
        for raw in data.get("entries", []):
            try:
                ent = ScheduleEntry(**raw)
                self._entries[ent.id] = ent
            except (TypeError, ValueError):
                continue

    def _save(self) -> None:
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            payload = {"entries": [asdict(e) for e in self._entries.values()]}
            self.path.write_text(
                json.dumps(payload, indent=2), encoding="utf-8",
            )
        except OSError:
            pass

    # --- API --------------------------------------------------------------

    def schedule_once(
        self, *, delay_seconds: int, prompt: str, reason: str = "",
        workspace: str = "", budget_usd: float = 0.0,
    ) -> ScheduleEntry:
        if delay_seconds < 1:
            raise ValueError("delay_seconds must be >= 1")
        ent = ScheduleEntry(
            id=uuid.uuid4().hex[:10],
            kind="once",
            prompt=prompt,
            reason=reason,
            delay_seconds=delay_seconds,
            next_fire_at=time.time() + delay_seconds,
            # Recorded so the headless daemon knows where to run the
            # prompt; the creating process's cwd is the best default.
            workspace=str(workspace or os.getcwd()),
            budget_usd=max(0.0, float(budget_usd or 0.0)),
        )
        with self._lock:
            self._entries[ent.id] = ent
            self._save()
        return ent

    def schedule_interval(
        self, *, every_seconds: int, prompt: str, reason: str = "",
        fire_immediately: bool = False,
        workspace: str = "", budget_usd: float = 0.0,
    ) -> ScheduleEntry:
        if every_seconds < 60:
            raise ValueError("every_seconds must be >= 60 (be sensible)")
        first_fire = time.time() + (0 if fire_immediately else every_seconds)
        ent = ScheduleEntry(
            id=uuid.uuid4().hex[:10],
            kind="interval",
            prompt=prompt,
            reason=reason,
            every_seconds=every_seconds,
            next_fire_at=first_fire,
            workspace=str(workspace or os.getcwd()),
            budget_usd=max(0.0, float(budget_usd or 0.0)),
        )
        with self._lock:
            self._entries[ent.id] = ent
            self._save()
        return ent

    def list_entries(self) -> list[ScheduleEntry]:
        with self._lock:
            return list(self._entries.values())

    def delete(self, entry_id: str) -> bool:
        with self._lock:
            removed = self._entries.pop(entry_id, None) is not None
            if removed:
                self._save()
            return removed

    def set_fire_callback(self, cb: Callable[[ScheduleEntry], None] | None) -> None:
        self._fire_callback = cb

    # --- background thread -----------------------------------------------

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop.clear()
        t = threading.Thread(target=self._run, name="delfin-scheduler", daemon=True)
        t.start()
        self._thread = t

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2)

    def tick(
        self,
        fire_callback: Callable[[ScheduleEntry], None] | None = None,
    ) -> int:
        """Run a single polling pass. Returns number of successful fires.

        ``fire_callback`` overrides the bound callback for this pass — the
        headless daemon passes its executor here without touching the
        dashboard binding. Disabled entries are never fired. A callback
        exception counts as a consecutive failure (entry retried next
        pass; disabled with a reason after ``_MAX_CONSECUTIVE_FAILURES``);
        raising :class:`DisableEntry` disables the entry immediately.

        Exposed separately so tests don't have to wait for the thread.
        """
        fired = 0
        changed = False
        now = time.time()
        with self._lock:
            for ent in list(self._entries.values()):
                if ent.disabled or not ent.is_due(now):
                    continue
                cb = fire_callback or self._fire_callback
                if cb is None:
                    continue
                try:
                    cb(ent)
                except DisableEntry as exc:
                    ent.disabled = True
                    ent.disabled_reason = (
                        exc.reason or "disabled by fire callback")
                    changed = True
                    continue
                except Exception as exc:
                    ent.fail_count += 1
                    if ent.fail_count >= _MAX_CONSECUTIVE_FAILURES:
                        ent.disabled = True
                        ent.disabled_reason = (
                            f"disabled after {ent.fail_count} consecutive "
                            f"failures (last: {str(exc)[:160]})")
                    changed = True
                    continue
                ent.fail_count = 0
                fired += 1
                still_active = ent.reschedule(now)
                if not still_active:
                    self._entries.pop(ent.id, None)
            if fired or changed:
                self._save()
        return fired

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self.tick()
            except Exception:
                pass
            self._stop.wait(_POLL_S)


def parse_interval_seconds(token: str) -> int | None:
    """Parse a ``/loop`` interval token (``5m`` / ``2h`` / ``1d`` / ``90s``)
    into seconds, clamped to the scheduler minimum (60s). None when the token
    is not a valid ``<N><s|m|h|d>`` interval."""
    import re as _re
    m = _re.fullmatch(r"(\d+)\s*([smhd])", (token or "").strip().lower())
    if not m:
        return None
    n = int(m.group(1))
    mult = {"s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]
    return max(60, n * mult)


_GLOBAL_LOCK = threading.Lock()
_GLOBAL: Scheduler | None = None


def get_scheduler(path: Path | None = None) -> Scheduler:
    global _GLOBAL
    with _GLOBAL_LOCK:
        if _GLOBAL is None:
            _GLOBAL = Scheduler(path=path)
            _GLOBAL.start()
    return _GLOBAL


def reset_scheduler() -> None:
    """Stop the global scheduler. Used by tests."""
    global _GLOBAL
    with _GLOBAL_LOCK:
        if _GLOBAL is not None:
            _GLOBAL.stop()
        _GLOBAL = None


__all__ = [
    "DisableEntry",
    "ScheduleEntry",
    "Scheduler",
    "get_scheduler",
    "reset_scheduler",
    "parse_interval_seconds",
]
