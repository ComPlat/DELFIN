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

Every route to a disabled entry goes through ``Scheduler._disable``,
which emits an attention event as well as persisting the reason: a
schedule that stops without saying so cannot be told apart from one
that is simply not due yet. The save is atomic and merges with what
is on disk, so a second writer's entries are neither torn nor lost.
"""

from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Optional


_DEFAULT_PATH = Path.home() / ".delfin" / "cron.json"
_POLL_S = 30.0
_MAX_CONSECUTIVE_FAILURES = 3
# How late a ONE-SHOT wake-up may still be honoured. A one-shot is a
# statement about a moment ("in two hours, check the calculation"); firing
# it days later, on a machine whose state has moved on, is not a late
# version of that request. Observed live: starting the daemon fired an
# entry that had been due for 85 minutes, against a /tmp workspace left
# by a test, and spent a real agent turn on it. Wide enough that a closed
# laptop or a rebooted node still gets its wake-up. Recurring entries are
# deliberately exempt — firing late IS the schedule.
_STALE_ONCE_GRACE_S = 4 * 3600


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

    def is_stale(self, now: float | None = None) -> bool:
        """Whether this entry is too late to be worth firing.

        Only one-shots can go stale. A recurring entry firing late is the
        schedule catching up, which is what it is for; a one-shot names a
        moment, and the moment is gone.
        """
        if self.kind != "once":
            return False
        return ((now or time.time()) - self.next_fire_at) > _STALE_ONCE_GRACE_S

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
        self._removed: set[str] = set()
        self._lock = threading.RLock()
        self._fire_callback: Optional[Callable[[ScheduleEntry], None]] = None
        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()
        self._load()

    # --- persistence ------------------------------------------------------

    def _read_file(self) -> dict[str, ScheduleEntry]:
        out: dict[str, ScheduleEntry] = {}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return out
        for raw in data.get("entries", []):
            try:
                ent = ScheduleEntry(**raw)
                out[ent.id] = ent
            except (TypeError, ValueError):
                continue
        return out

    def _load(self) -> None:
        self._entries.update(self._read_file())

    def _save(self) -> None:
        """Persist atomically, and merge rather than overwrite.

        Two problems, one write. A plain ``write_text`` is not atomic: this
        was the one state file in the agent that a reader could catch
        half-written or empty, which for a schedule means entries simply
        gone. And an in-process Scheduler loads once and never reloads, so
        its save flattened whatever another process (the headless daemon,
        the CLI, a second dashboard) had created in the meantime —
        scheduling into a file that another writer would silently revert.

        Merging on write keeps entries this instance never knew about,
        while entries it deliberately removed stay removed.
        """
        try:
            merged = self._read_file()
            merged.update(self._entries)
            for entry_id in self._removed:
                merged.pop(entry_id, None)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            payload = {"entries": [asdict(e) for e in merged.values()]}
            fd, tmp_name = tempfile.mkstemp(
                prefix=f".{self.path.name}.", suffix=".tmp",
                dir=str(self.path.parent))
            tmp = Path(tmp_name)
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                    json.dump(payload, fh, indent=2)
                os.replace(tmp, self.path)
            finally:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass
        except OSError:
            pass

    # --- disable notification (the single choke point) --------------------

    def _disable(self, ent: ScheduleEntry, reason: str) -> None:
        """Disable an entry AND say so where a user will see it.

        Every route to a disabled entry goes through here. The workspace
        gate used to raise outside the block that emitted attention, so the
        one case most worth telling about — a moved or deleted workspace
        killing the schedule outright — was also the only silent one, while
        an ordinary run failure did notify. A schedule that stops without
        saying so is indistinguishable from one that is simply not due yet.
        """
        ent.disabled = True
        ent.disabled_reason = reason
        try:
            from .attention import emit_attention
            emit_attention(
                "run_failed",
                title=f"Schedule {ent.id} disabled",
                detail=f"{ent.reason or ent.prompt[:120]} — {reason}"[:400],
                workspace=str(ent.workspace or ""),
            )
        except Exception:
            pass        # notification is best-effort; the disable is not

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
        """Every entry in the schedule, including ones another process
        created since this instance was built."""
        with self._lock:
            self._entries.update(
                {k: v for k, v in self._read_file().items()
                 if k not in self._entries and k not in self._removed})
            return list(self._entries.values())

    def delete(self, entry_id: str) -> bool:
        with self._lock:
            removed = self._entries.pop(entry_id, None) is not None
            if not removed and entry_id in self._read_file():
                removed = True          # created by another process
            if removed:
                self._removed.add(entry_id)
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
            # Reload first: a long-lived in-process scheduler otherwise
            # never sees an entry the CLI, the daemon or a second front end
            # wrote after it was constructed.
            self._entries.update(
                {k: v for k, v in self._read_file().items()
                 if k not in self._entries and k not in self._removed})
            for ent in list(self._entries.values()):
                if ent.disabled or not ent.is_due(now):
                    continue
                # A one-shot that went unseen for hours is not fired late;
                # it is disabled with the reason, so the entry can still be
                # found and understood rather than vanishing.
                if ent.is_stale(now):
                    self._disable(ent, (
                        f"not fired: overdue by "
                        f"{int(now - ent.next_fire_at)}s, past the "
                        f"{int(_STALE_ONCE_GRACE_S)}s grace for a one-shot "
                        "wake-up. Schedule it again if it is still wanted."))
                    changed = True
                    continue
                cb = fire_callback or self._fire_callback
                if cb is None:
                    continue
                try:
                    cb(ent)
                except DisableEntry as exc:
                    self._disable(
                        ent, exc.reason or "disabled by fire callback")
                    changed = True
                    continue
                except Exception as exc:
                    ent.fail_count += 1
                    if ent.fail_count >= _MAX_CONSECUTIVE_FAILURES:
                        self._disable(ent, (
                            f"disabled after {ent.fail_count} consecutive "
                            f"failures (last: {str(exc)[:160]})"))
                    changed = True
                    continue
                ent.fail_count = 0
                fired += 1
                still_active = ent.reschedule(now)
                if not still_active:
                    self._entries.pop(ent.id, None)
                    self._removed.add(ent.id)
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
