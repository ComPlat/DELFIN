"""Runtime layer — async runs, a persistent run store, and an event stream.

Turns the platform from "call ``run_application`` and block" into a managed
service: submit a run, get a handle/id back immediately, poll status, read the
named outputs, follow a structured event stream, list past runs, and cancel.
This is what lets a human or an agent fold a long DELFIN workflow into their own
program without blocking on it.

Run records persist as JSON under ``$DELFIN_RUNS_DIR`` (default
``~/.delfin/runs``), so history survives the process.  Execution is synchronous
*inside a background thread*; this module is the seam — callers see only handles.
"""

from __future__ import annotations

import json
import os
import threading
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


class RunStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    CANCELLED = "cancelled"


_TERMINAL = {RunStatus.SUCCESS.value, RunStatus.FAILED.value, RunStatus.CANCELLED.value}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_runs_dir() -> Path:
    env = os.environ.get("DELFIN_RUNS_DIR")
    return Path(env) if env else (Path.home() / ".delfin" / "runs")


class _Cancelled(Exception):
    """Raised internally to abort a run between steps."""


@dataclass
class RunRecord:
    id: str
    kind: str                       # "application" | "pipeline"
    name: str
    inputs: Dict[str, Any] = field(default_factory=dict)
    status: str = RunStatus.PENDING.value
    created_at: str = ""
    started_at: str = ""
    finished_at: str = ""
    outputs: Dict[str, Any] = field(default_factory=dict)
    error: str = ""
    events: List[Dict[str, Any]] = field(default_factory=list)
    work_dir: str = ""
    metrics: Dict[str, Any] = field(default_factory=dict)

    @property
    def done(self) -> bool:
        return self.status in _TERMINAL


class RunStore:
    """Thread-safe JSON-file persistence for :class:`RunRecord`."""

    def __init__(self, base: Optional[Path | str] = None):
        self.base = Path(base) if base else _default_runs_dir()
        self.base.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _path(self, run_id: str) -> Path:
        return self.base / f"{run_id}.json"

    def save(self, rec: RunRecord) -> None:
        with self._lock:
            self._path(rec.id).write_text(
                json.dumps(asdict(rec), indent=2, default=str), encoding="utf-8"
            )

    def get(self, run_id: str) -> Optional[RunRecord]:
        p = self._path(run_id)
        if not p.is_file():
            return None
        try:
            return RunRecord(**json.loads(p.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError, TypeError):
            return None

    def list(self) -> List[RunRecord]:
        recs: List[RunRecord] = []
        for p in self.base.glob("*.json"):
            try:
                recs.append(RunRecord(**json.loads(p.read_text(encoding="utf-8"))))
            except (OSError, json.JSONDecodeError, TypeError):
                continue
        return sorted(recs, key=lambda r: r.created_at, reverse=True)


@dataclass
class RunHandle:
    """A reference to a submitted run."""

    id: str
    _runtime: "Runtime"

    def record(self) -> Optional[RunRecord]:
        return self._runtime.store.get(self.id)

    def status(self) -> Optional[str]:
        rec = self.record()
        return rec.status if rec else None

    def done(self) -> bool:
        rec = self.record()
        return bool(rec and rec.done)

    def outputs(self) -> Dict[str, Any]:
        rec = self.record()
        return rec.outputs if rec else {}

    def events(self) -> List[Dict[str, Any]]:
        rec = self.record()
        return rec.events if rec else []

    def cancel(self) -> None:
        self._runtime.cancel(self.id)

    def wait(self, *, timeout: Optional[float] = None, poll: float = 0.05) -> Optional[RunRecord]:
        """Block until the run reaches a terminal state (or *timeout*)."""
        import time as _time
        deadline = None if timeout is None else (_time.monotonic() + timeout)
        while True:
            rec = self.record()
            if rec and rec.done:
                return rec
            if deadline is not None and _time.monotonic() >= deadline:
                return rec
            _time.sleep(poll)


class Runtime:
    """Submits runs to background threads and tracks them in a :class:`RunStore`."""

    def __init__(self, store: Optional[RunStore] = None):
        self.store = store or RunStore()
        self._cancels: Dict[str, threading.Event] = {}
        self._threads: Dict[str, threading.Thread] = {}
        self._lock = threading.Lock()

    # -- event helpers -------------------------------------------------

    def _append_event(self, run_id: str, kind: str, data: Dict[str, Any]) -> None:
        rec = self.store.get(run_id)
        if rec is None:
            return
        rec.events.append({"t": _now(), "event": kind, **data})
        self.store.save(rec)

    # -- submission ----------------------------------------------------

    def submit_application(
        self,
        name: str,
        *,
        cores: int = 1,
        geometry: Optional[str | Path] = None,
        work_dir: Optional[str | Path] = None,
        inputs: Optional[Dict[str, Any]] = None,
    ) -> RunHandle:
        """Submit an application run; returns immediately with a handle."""
        run_id = uuid.uuid4().hex[:12]
        rec = RunRecord(
            id=run_id, kind="application", name=name,
            inputs=dict(inputs or {}), created_at=_now(),
            work_dir=str(work_dir) if work_dir else "",
        )
        self.store.save(rec)

        cancel = threading.Event()
        with self._lock:
            self._cancels[run_id] = cancel
        thread = threading.Thread(
            target=self._run_application,
            args=(run_id, name, cores, geometry, work_dir, dict(inputs or {}), cancel),
            daemon=True,
        )
        with self._lock:
            self._threads[run_id] = thread
        thread.start()
        return RunHandle(run_id, self)

    def _run_application(self, run_id, name, cores, geometry, work_dir, inputs, cancel):
        rec = self.store.get(run_id)
        if rec is None:
            return
        rec.status = RunStatus.RUNNING.value
        rec.started_at = _now()
        rec.events.append({"t": _now(), "event": "run_started", "name": name})
        self.store.save(rec)

        try:
            from delfin.tools._application import extract_outputs, get_application

            app = get_application(name)
            if app is None:
                raise ValueError(f"unknown application {name!r}")
            missing = app.missing_inputs(inputs)
            if missing:
                raise ValueError(f"missing required input(s): {', '.join(missing)}")

            pipeline = app.build(**inputs)
            report = pipeline.validate(geometry=bool(geometry))
            if not report.ok:
                raise ValueError("application failed static validation")

            def _on_step(r):
                if cancel.is_set():
                    raise _Cancelled()
                self._append_event(
                    run_id, "step_finished",
                    {"step": r.step_name, "ok": r.ok,
                     "elapsed_s": round(r.elapsed_seconds, 3)},
                )

            pipeline.on_step(_on_step)
            result = pipeline.run(
                cores=cores, geometry=geometry,
                work_dir=Path(work_dir) if work_dir else None,
            )
            outputs = extract_outputs(app, pipeline, result)

            rec = self.store.get(run_id)
            rec.status = RunStatus.SUCCESS.value if result.ok else RunStatus.FAILED.value
            rec.outputs = outputs
            rec.finished_at = _now()
            rec.metrics = {
                "steps": len(result.results),
                "elapsed_s": round(sum(r.elapsed_seconds for r in result.all_results), 3),
            }
            rec.events.append({"t": _now(), "event": "run_finished", "status": rec.status})
            self.store.save(rec)

        except _Cancelled:
            rec = self.store.get(run_id)
            if rec is not None:
                rec.status = RunStatus.CANCELLED.value
                rec.finished_at = _now()
                rec.events.append({"t": _now(), "event": "run_cancelled"})
                self.store.save(rec)
        except Exception as exc:  # noqa: BLE001 — record any failure
            rec = self.store.get(run_id)
            if rec is not None:
                rec.status = RunStatus.FAILED.value
                rec.error = str(exc)
                rec.finished_at = _now()
                rec.events.append({"t": _now(), "event": "run_failed", "error": str(exc)})
                self.store.save(rec)
        finally:
            with self._lock:
                self._cancels.pop(run_id, None)
                self._threads.pop(run_id, None)

    # -- control / query ----------------------------------------------

    def cancel(self, run_id: str) -> bool:
        """Request cooperative cancellation (takes effect before the next step)."""
        with self._lock:
            ev = self._cancels.get(run_id)
        if ev is None:
            return False
        ev.set()
        return True

    def get(self, run_id: str) -> Optional[RunRecord]:
        return self.store.get(run_id)

    def list_runs(self) -> List[RunRecord]:
        return self.store.list()


# --- module-level default runtime -----------------------------------------

_RUNTIME: Optional[Runtime] = None
_RT_LOCK = threading.Lock()


def get_runtime() -> Runtime:
    """The process-wide default runtime (lazily created)."""
    global _RUNTIME
    if _RUNTIME is None:
        with _RT_LOCK:
            if _RUNTIME is None:
                _RUNTIME = Runtime()
    return _RUNTIME


__all__ = [
    "RunStatus",
    "RunRecord",
    "RunStore",
    "RunHandle",
    "Runtime",
    "get_runtime",
]
