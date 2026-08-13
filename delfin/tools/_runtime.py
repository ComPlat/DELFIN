"""Runtime layer — async runs, a persistent run store, and an event stream.

Turns the platform from "call ``run_application`` and block" into a managed
service: submit a run, get a handle/id back immediately, poll status, read the
named outputs, follow a structured event stream, list past runs, and cancel.
This is what lets a human or an agent fold a long DELFIN workflow into their own
program without blocking on it.

Run records persist as JSON under ``$DELFIN_RUNS_DIR`` (default
``~/.delfin/runs``), so history survives the process.  Execution is synchronous
*inside a background thread*; this module is the seam — callers see only handles.

For ``backend="slurm"`` the executing process lives on a compute node, so the
login-node record cannot learn anything by itself.  Every non-terminal
slurm-backed record is therefore reconciled against ``sacct`` when it is read,
and the answer is tri-state on purpose: a state, "the scheduler does not know
this id", or "the scheduler could not be asked".
"""

from __future__ import annotations

import json
import os
import threading
import time
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

# --- SLURM reconciliation -------------------------------------------------
#
# A ``backend="slurm"`` run is executed by ``execute_run`` on a compute node,
# and that compute-node process was the ONLY writer able to move the record
# off PENDING. Every way a batch job can end without running its payload --
# OUT_OF_MEMORY, TIMEOUT, NODE_FAIL, an admin ``scancel``, a held job that
# never starts -- therefore left the record PENDING for ever, indistinguishable
# from "queued" and from "running fine". On a shared node that is not a
# cosmetic bug: the run store is what tells a user whether their allocation is
# still being spent, and a record that can only ever say "pending" says nothing.
#
# The scheduler is the authority on a scheduler job, so a non-terminal
# slurm-backed record is reconciled against ``sacct`` whenever it is read.

# The scheduler answered and the job is in this state.
_SLURM_RUNNING_STATES = frozenset({
    "RUNNING", "COMPLETING", "STAGE_OUT", "SUSPENDED",
})
_SLURM_QUEUED_STATES = frozenset({
    "PENDING", "CONFIGURING", "REQUEUED", "REQUEUE_HOLD", "REQUEUE_FED",
    "RESIZING", "SIGNALING", "RESV_DEL_HOLD",
})
_SLURM_FAILED_STATES = frozenset({
    "FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL",
    "DEADLINE", "PREEMPTED", "REVOKED", "SPECIAL_EXIT",
})
_SLURM_CANCELLED_STATES = frozenset({"CANCELLED"})
_SLURM_OK_STATES = frozenset({"COMPLETED"})

# How long an id may be unknown to ``sacct`` before the record is closed.
# Accounting can lag a fresh submission by seconds; a job that is still
# missing a quarter of an hour later is not going to appear.
_SLURM_ABSENT_GRACE_S = 900.0

# ``RunHandle.wait`` used to default to no timeout at all with a 0.05 s poll:
# a 20 Hz busy-wait on the login node, beside the queue it was waiting on,
# that could not return while a record was stuck non-terminal. It now has a
# default cap and backs off, and a wait that expires says so through
# ``RunHandle.state()`` rather than pretending the run is still going.
_DEFAULT_WAIT_TIMEOUT_S = 3600.0
_WAIT_POLL_START_S = 0.5
_WAIT_POLL_MAX_S = 15.0


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_runs_dir() -> Path:
    env = os.environ.get("DELFIN_RUNS_DIR")
    return Path(env) if env else (Path.home() / ".delfin" / "runs")


def _default_calc_dir() -> Path:
    """DELFIN's standard calculations directory (``$DELFIN_CALC_DIR`` → ``~/calc``).

    Platform runs land here (one ``<app>_<run_id>`` folder each) so their outputs
    appear in the Calculations browser, consistent with the classic DELFIN
    workflow.
    """
    env = os.environ.get("DELFIN_CALC_DIR")
    return Path(env) if env else (Path.home() / "calc")


class _Cancelled(Exception):
    """Raised internally to abort a run between steps."""


def _default_slurm_query(job_id: str) -> Optional[str]:
    """Ask the scheduler for one job's state. Tri-state, deliberately.

    * ``"RUNNING"`` / ``"TIMEOUT"`` / … — the scheduler answered.
    * ``""`` — the scheduler answered and does not know this id.
    * ``None`` — the scheduler could NOT be asked (no ``sacct``, non-zero
      exit, timeout).

    Collapsing the last two into one empty answer is what turns "SLURM says
    nothing changed" into "SLURM could not be reached", and a run that is
    burning an allocation into one that looks merely quiet.
    """
    import shutil
    import subprocess

    sacct = shutil.which("sacct")
    if sacct is None:
        return None
    try:
        proc = subprocess.run(
            [sacct, "-j", str(job_id), "-n", "-X", "-P", "-o", "State"],
            capture_output=True, text=True, timeout=20,
        )
    except Exception:      # noqa: BLE001 — any failure means "not asked"
        return None
    if proc.returncode != 0:
        return None
    for line in (proc.stdout or "").splitlines():
        text = line.strip().upper()
        if text:
            # "CANCELLED by 1234" → CANCELLED; "COMPLETED+" → COMPLETED.
            return text.split()[0].rstrip("+")
    return ""


# Injection seam: tests replace this instead of touching a real queue.
_slurm_query = _default_slurm_query


def _reconcile_slurm_record(rec: "RunRecord", store: "RunStore") -> "RunRecord":
    """Fold the scheduler's view of a slurm-backed run into its record.

    A no-op for local runs, for records without a job id, and for records
    that already reached a terminal state. Persists only on an actual
    change, so reading a run stays cheap.
    """
    if rec is None or rec.done:
        return rec
    job_id = str((rec.metrics or {}).get("slurm_job_id") or "").strip()
    if not job_id:
        return rec
    try:
        state = _slurm_query(job_id)
    except Exception:      # noqa: BLE001 — a broken probe is "not asked"
        state = None

    prev_note = (rec.metrics or {}).get("slurm_state")
    new_status = ""
    error = ""
    note = ""

    if state is None:
        # Visible degradation: the record keeps its status, but it now says
        # the status is unverified rather than implying it was checked.
        note = "unavailable"
    elif state == "":
        note = "absent"
        first_missing = float((rec.metrics or {}).get("slurm_absent_since") or 0.0)
        now = time.time()
        if not first_missing:
            rec.metrics["slurm_absent_since"] = now
        elif (now - first_missing) > _SLURM_ABSENT_GRACE_S:
            new_status = RunStatus.FAILED.value
            error = (f"SLURM job {job_id} is unknown to the scheduler and the "
                     f"run never reported a result — the job is gone")
    else:
        note = state
        rec.metrics.pop("slurm_absent_since", None)
        if state in _SLURM_RUNNING_STATES:
            if rec.status == RunStatus.PENDING.value:
                new_status = RunStatus.RUNNING.value
        elif state in _SLURM_QUEUED_STATES:
            pass
        elif state in _SLURM_CANCELLED_STATES:
            new_status = RunStatus.CANCELLED.value
        elif state in _SLURM_FAILED_STATES:
            new_status = RunStatus.FAILED.value
            error = f"SLURM job {job_id} ended in state {state}"
        elif state in _SLURM_OK_STATES:
            # The batch job is over but the payload never wrote a result:
            # either it never started or it died before saving. Reporting
            # success here would invent outputs that do not exist.
            new_status = RunStatus.FAILED.value
            error = (f"SLURM job {job_id} finished (COMPLETED) without the run "
                     f"recording a result — check the slurm_*.out/.err files "
                     f"in {rec.work_dir or 'the work directory'}")
        else:
            new_status = RunStatus.FAILED.value
            error = f"SLURM job {job_id} ended in state {state}"

    if not new_status and note == prev_note:
        return rec

    # Re-read before writing: the compute node may have finished between the
    # query and now, and its result always wins over our inference.
    fresh = store.get(rec.id)
    if fresh is not None:
        if fresh.done:
            return fresh
        fresh.metrics = {**(fresh.metrics or {}), **(rec.metrics or {})}
        rec = fresh
    rec.metrics["slurm_state"] = note
    if new_status:
        rec.status = new_status
        rec.finished_at = rec.finished_at or _now()
        if error:
            rec.error = error
        rec.events.append({"t": _now(), "event": "slurm_reconciled",
                           "job_id": job_id, "state": note,
                           "status": new_status})
    store.save(rec)
    return rec


def _apply_resources(pipeline, resources: Dict[str, Any]) -> None:
    """Inject run-level resource kwargs (e.g. maxcore) into every step.

    Explicit per-step kwargs always win; non-QM steps ignore the extra kwargs.
    """
    if not resources:
        return

    def _inject(trunk) -> None:
        for spec in trunk:
            if getattr(spec, "kwargs", None) is not None:
                spec.kwargs = {**resources, **spec.kwargs}

    _inject(pipeline._trunk)
    for branch in getattr(pipeline, "_branches", {}).values():
        _inject(branch._trunk)


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
        # Goes through the runtime, not straight to the store: a slurm-backed
        # record is reconciled against the scheduler on every read.
        return self._runtime.get(self.id)

    def status(self) -> Optional[str]:
        rec = self.record()
        return rec.status if rec else None

    def state(self) -> Dict[str, Any]:
        """Tri-state view of the run — including "we could not find out".

        ``known`` is the point. A slurm-backed run whose scheduler could not
        be reached is NOT running-as-far-as-we-know; it is unverified, and a
        caller that treats the two the same will wait for ever on a job that
        was cancelled hours ago.

        Keys: ``status``, ``done``, ``known`` (False when the last
        reconciliation could not reach the scheduler, or the run has no
        record at all), ``slurm_state`` and ``detail``.
        """
        rec = self.record()
        if rec is None:
            return {"status": "", "done": False, "known": False,
                    "slurm_state": "", "detail": f"no record for run {self.id}"}
        slurm_state = str((rec.metrics or {}).get("slurm_state") or "")
        known = rec.done or slurm_state != "unavailable"
        detail = ""
        if not known:
            detail = ("the scheduler could not be asked, so this status is "
                      "the last thing we were told, not the current state")
        elif slurm_state == "absent":
            detail = ("the scheduler does not (yet) know this job id")
        return {"status": rec.status, "done": rec.done, "known": known,
                "slurm_state": slurm_state, "detail": detail}

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

    def wait(
        self,
        *,
        timeout: Optional[float] = None,
        poll: float = _WAIT_POLL_START_S,
    ) -> Optional[RunRecord]:
        """Block until the run reaches a terminal state, or the wait expires.

        ``timeout=None`` means the default cap (:data:`_DEFAULT_WAIT_TIMEOUT_S`),
        not "for ever": an unbounded wait on a run whose state nothing could
        advance is an idle core held on a login node next to the queue it is
        waiting on. The poll interval backs off from *poll* to
        :data:`_WAIT_POLL_MAX_S` for the same reason.

        On expiry the returned record is simply NOT terminal — ask
        :meth:`state` whether that means "still going" or "we could not find
        out". ``None`` means there is no such run.
        """
        limit = _DEFAULT_WAIT_TIMEOUT_S if timeout is None else float(timeout)
        deadline = time.monotonic() + max(0.0, limit)
        interval = max(0.01, float(poll))
        while True:
            rec = self.record()
            if rec and rec.done:
                return rec
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return rec
            time.sleep(min(interval, remaining))
            interval = min(_WAIT_POLL_MAX_S, interval * 2)


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
        maxcore: Optional[int] = None,
        geometry: Optional[str | Path] = None,
        work_dir: Optional[str | Path] = None,
        inputs: Optional[Dict[str, Any]] = None,
        backend: str = "local",
    ) -> RunHandle:
        """Submit an application run; returns immediately with a handle.

        ``cores`` is the parallelism (ORCA PAL); ``maxcore`` the memory per core
        (MB) injected into the QM steps. ``backend="local"`` runs in a background
        thread on this machine; ``backend="slurm"`` submits an sbatch job that
        executes the run on a compute node and writes its result back into the
        (shared) run store.
        """
        run_id = uuid.uuid4().hex[:12]
        # Default each run into DELFIN's standard calculations directory
        # (~/calc/<app>_<id>) so its outputs show up in the Calculations browser,
        # just like the classic DELFIN workflow.
        if work_dir is None:
            work_dir = _default_calc_dir() / f"{name}_{run_id}"
        resources: Dict[str, Any] = {}
        if maxcore is not None:
            resources["maxcore"] = maxcore
        rec = RunRecord(
            id=run_id, kind="application", name=name,
            inputs=dict(inputs or {}), created_at=_now(),
            work_dir=str(work_dir),
            metrics={"cores": cores, "backend": backend, "resources": resources},
        )
        self.store.save(rec)

        if backend == "slurm":
            self._submit_slurm(run_id, name, cores, work_dir)
            return RunHandle(run_id, self)

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

    def _submit_slurm(self, run_id, name, cores, work_dir):
        """Submit the run as a SLURM batch job (executed via execute_run)."""
        import shlex
        import shutil
        import subprocess
        import sys

        rec = self.store.get(run_id)
        sbatch = shutil.which("sbatch")
        if sbatch is None:
            if rec is not None:
                rec.status = RunStatus.FAILED.value
                rec.error = "sbatch not found — SLURM is not available on this host"
                rec.finished_at = _now()
                self.store.save(rec)
            return

        wd = Path(work_dir) if work_dir else (self.store.base / run_id)
        wd.mkdir(parents=True, exist_ok=True)
        store_base = str(self.store.base)
        code = (f"from delfin.tools._runtime import execute_run; "
                f"execute_run({run_id!r}, {store_base!r})")
        script = (
            "#!/bin/bash\n"
            f"#SBATCH --job-name=delfin-app-{name}\n"
            "#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={cores}\n"
            f"#SBATCH --output={wd}/slurm_%j.out\n"
            f"#SBATCH --error={wd}/slurm_%j.err\n\n"
            f"{shlex.quote(sys.executable)} -c {shlex.quote(code)}\n"
        )
        script_path = wd / "submit.sh"
        script_path.write_text(script)

        proc = subprocess.run([sbatch, str(script_path)], capture_output=True, text=True)
        rec = self.store.get(run_id)
        if rec is None:
            return
        if proc.returncode == 0:
            job_id = (proc.stdout or "").strip().split()[-1] if proc.stdout else ""
            rec.metrics["slurm_job_id"] = job_id
            rec.events.append({"t": _now(), "event": "slurm_submitted", "job_id": job_id})
        else:
            rec.status = RunStatus.FAILED.value
            rec.error = f"sbatch failed: {(proc.stderr or '')[:300]}"
            rec.finished_at = _now()
        self.store.save(rec)

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
            _apply_resources(pipeline, (rec.metrics or {}).get("resources") or {})
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
            rec.metrics.update(
                steps=len(result.results),
                elapsed_s=round(sum(r.elapsed_seconds for r in result.all_results), 3),
            )
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
        """Cancel a run: ``scancel`` a SLURM job, else cooperative local stop.

        Returns True only when the cancellation was actually requested. A
        missing ``scancel`` or a failing one used to be swallowed: the record
        was marked CANCELLED and True returned, and because CANCELLED is
        terminal nothing ever looked at that record again — while the job ran
        its full allocation on a node somebody else was queued for. The
        terminal status is now written only after the scheduler confirms it.
        """
        rec = self.store.get(run_id)
        slurm_job = rec.metrics.get("slurm_job_id") if rec else None
        if slurm_job:
            import shutil
            import subprocess
            scancel = shutil.which("scancel")
            if scancel is None:
                return False           # nothing was cancelled; say so
            try:
                proc = subprocess.run(
                    [scancel, str(slurm_job)], capture_output=True, text=True,
                    timeout=20,
                )
            except Exception:          # noqa: BLE001 — could not even run it
                return False
            if proc.returncode != 0:
                return False
            if rec is not None and not rec.done:
                rec.events.append({"t": _now(), "event": "slurm_cancel_requested",
                                   "job_id": str(slurm_job)})
                self.store.save(rec)
            # One state read decides whether the job is really over. scancel
            # is asynchronous, so "not terminal yet" is normal — the record
            # stays non-terminal and the next read reconciles it.
            self.get(run_id)
            return True
        with self._lock:
            ev = self._cancels.get(run_id)
        if ev is None:
            return False
        ev.set()
        return True

    def get(self, run_id: str) -> Optional[RunRecord]:
        return _reconcile_slurm_record(self.store.get(run_id), self.store)

    def list_runs(self) -> List[RunRecord]:
        return [_reconcile_slurm_record(r, self.store) for r in self.store.list()]


# --- compute-node executor (used by SLURM jobs) ---------------------------


def execute_run(run_id: str, store_base: Optional[str] = None) -> None:
    """Execute a previously-submitted application run by id.

    The login node creates the :class:`RunRecord` (with name + inputs); a SLURM
    job calls this on the compute node to run the application locally and write
    its status / outputs back into the shared run store.
    """
    store = RunStore(store_base) if store_base else RunStore()
    rec = store.get(run_id)
    if rec is None:
        raise SystemExit(f"unknown run {run_id!r}")

    rec.status = RunStatus.RUNNING.value
    rec.started_at = _now()
    rec.events.append({"t": _now(), "event": "run_started", "name": rec.name})
    store.save(rec)

    try:
        from delfin.tools._application import extract_outputs, get_application

        app = get_application(rec.name)
        if app is None:
            raise ValueError(f"unknown application {rec.name!r}")
        cores = int(rec.metrics.get("cores", 1) or 1)
        pipeline = app.build(**rec.inputs)
        _apply_resources(pipeline, (rec.metrics or {}).get("resources") or {})
        result = pipeline.run(
            cores=cores,
            work_dir=Path(rec.work_dir) if rec.work_dir else None,
        )
        outputs = extract_outputs(app, pipeline, result)

        rec = store.get(run_id)
        rec.status = RunStatus.SUCCESS.value if result.ok else RunStatus.FAILED.value
        rec.outputs = outputs
        rec.finished_at = _now()
        rec.metrics.update(
            steps=len(result.results),
            elapsed_s=round(sum(r.elapsed_seconds for r in result.all_results), 3),
        )
        rec.events.append({"t": _now(), "event": "run_finished", "status": rec.status})
        store.save(rec)
    except Exception as exc:  # noqa: BLE001
        rec = store.get(run_id)
        if rec is not None:
            rec.status = RunStatus.FAILED.value
            rec.error = str(exc)
            rec.finished_at = _now()
            store.save(rec)
        raise


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
    "execute_run",
]
