"""Parallel execution helpers for DELFIN classic and manually modes."""

from __future__ import annotations

import logging
import os
import re
import time
import statistics
import threading
from pathlib import Path
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, Optional, Set, List

from delfin.common.logging import get_logger
from delfin.dynamic_pool import PoolJob, JobPriority
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.imag import run_IMAG
from delfin.xyz_io import (
    create_s1_optimization_input,
    read_and_modify_file_1,
    read_xyz_and_create_input2,
    read_xyz_and_create_input3,
)

logger = get_logger(__name__)


JOB_DURATION_HISTORY: Dict[str, deque[float]] = defaultdict(lambda: deque(maxlen=8))


@dataclass
class WorkflowJob:
    """Represents a single ORCA task with dependency metadata."""

    job_id: str
    work: Callable[[int], None]
    description: str
    dependencies: Set[str] = field(default_factory=set)
    cores_min: int = 1
    cores_optimal: int = 2
    cores_max: int = 2
    priority: JobPriority = JobPriority.NORMAL
    memory_mb: Optional[int] = None
    estimated_duration: float = 3600.0


@dataclass
class WorkflowRunResult:
    """Summary of a workflow scheduler execution."""

    completed: Set[str] = field(default_factory=set)
    failed: Dict[str, str] = field(default_factory=dict)
    skipped: Dict[str, List[str]] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        return not self.failed and not self.skipped


class _WorkflowManager:
    """Schedules dependent ORCA jobs on the shared global dynamic core pool."""

    def __init__(self, config: Dict[str, Any], label: str, *, max_jobs_override: Optional[int] = None):
        self.config = config
        self.label = label

        global_mgr = get_global_manager()
        if not global_mgr.is_initialized():
            raise RuntimeError(
                f"[{label}] Global job manager not initialized; call get_global_manager().initialize(config) first."
            )

        self.pool = global_mgr.get_pool()
        pool_id = id(self.pool)

        self.total_cores = max(1, global_mgr.total_cores)
        self.maxcore_mb = max(256, _parse_int(config.get('maxcore'), fallback=1000))

        if max_jobs_override is not None and max_jobs_override > 0:
            desired_jobs = max(1, max_jobs_override)
        else:
            desired_jobs = max(1, global_mgr.max_jobs)

        self.max_jobs = max(1, min(desired_jobs, self.pool.max_concurrent_jobs))

        logger.info(
            "[%s] ✓ USING GLOBAL SHARED POOL (pool_id=%d, %d cores)",
            label,
            pool_id,
            self.total_cores,
        )

        self._sync_parallel_flag()

        self._jobs: Dict[str, WorkflowJob] = {}
        self._completed: Set[str] = set()
        self._failed: Dict[str, str] = {}
        self._skipped: Dict[str, List[str]] = {}
        self._lock = threading.RLock()
        self._event = threading.Event()

    def derive_core_bounds(self, preferred_opt: Optional[int] = None, *, hint: Optional[str] = None) -> tuple[int, int, int]:
        cores_min = 1 if self.total_cores == 1 else 2
        cores_max = self.total_cores

        if not self._parallel_enabled:
            default_opt = cores_max
        else:
            default_opt = self._base_share()
            if hint:
                default_opt = self._suggest_optimal_from_hint(
                    hint,
                    default_opt,
                    cores_max,
                    cores_min,
                    None,
                )

        if preferred_opt is not None:
            default_opt = preferred_opt

        preferred = max(cores_min, min(default_opt, cores_max))
        return cores_min, preferred, cores_max

    def add_job(self, job: WorkflowJob) -> None:
        if job.job_id in self._jobs:
            raise ValueError(f"Duplicate workflow job id '{job.job_id}'")

        deps = set(job.dependencies)
        job.dependencies = deps

        job.cores_min = max(1, min(job.cores_min, self.total_cores))
        job.cores_max = max(job.cores_min, min(job.cores_max, self.total_cores))
        job.cores_optimal = max(job.cores_min, min(job.cores_optimal, job.cores_max))

        self._auto_tune_job(job)

        if job.memory_mb is None:
            job.memory_mb = job.cores_optimal * self.maxcore_mb

        self._jobs[job.job_id] = job
        logger.info(
            "[%s] Registered job %s (%s); deps=%s",
            self.label,
            job.job_id,
            job.description,
            ",".join(sorted(job.dependencies)) or "none",
        )

    def has_jobs(self) -> bool:
        return bool(self._jobs)

    def run(self) -> None:
        if not self._jobs:
            logger.info("[%s] No jobs to schedule", self.label)
            return

        pending: Dict[str, WorkflowJob] = dict(self._jobs)
        self._sync_parallel_flag()
        logger.info(
            "[%s] Scheduling %d jobs across %d cores using GLOBAL SHARED pool (pool_id=%d)",
            self.label,
            len(pending),
            self.total_cores,
            id(self.pool),
        )

        while pending:
            ready = [job for job in pending.values() if job.dependencies <= self._completed]

            if not ready:
                blocked = {
                    job.job_id: sorted(job.dependencies - self._completed)
                    for job in pending.values()
                }

                failed_ids = set(self._failed)
                skipped_any = False
                if failed_ids:
                    for job_id, missing in list(blocked.items()):
                        if any(dep in failed_ids for dep in missing):
                            self._mark_skipped(job_id, missing)
                            pending.pop(job_id, None)
                            skipped_any = True
                if skipped_any:
                    continue

                try:
                    status = self.pool.get_status()
                except Exception:
                    status = None

                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        "[%s] Waiting; blocked jobs=%s | completed=%s",
                        self.label,
                        blocked,
                        sorted(self._completed),
                    )

                running = (status and (
                    status.get('running_jobs', 0) > 0 or
                    status.get('queued_jobs', 0) > 0
                ))

                if running:
                    self._event.wait(timeout=0.5)
                    self._event.clear()
                    continue

                for job_id, missing in blocked.items():
                    self._mark_skipped(job_id, missing or ['unresolved dependency'])
                    pending.pop(job_id, None)
                break

            ready.sort(key=self._job_order_key)

            try:
                status = self.pool.get_status()
            except Exception:
                status = None

            logger.debug(
                "[%s] Ready jobs=%s | completed=%s",
                self.label,
                [job.job_id for job in ready],
                sorted(self._completed),
            )

            if (status and len(ready) == 1 and
                    status.get('running_jobs', 0) == 0 and
                    status.get('queued_jobs', 0) == 0):
                self._apply_exclusive_allocation(ready[0])

            for job in ready:
                self._submit(job)
                pending.pop(job.job_id, None)

            # After submitting ready jobs, wait briefly for state changes
            # This allows the loop to react immediately when dependencies are fulfilled
            # instead of only checking every 0.5s when blocked
            if pending:
                self._event.wait(timeout=0.1)
                self._event.clear()

        self.pool.wait_for_completion()

        if self._failed:
            logger.warning(
                "[%s] Workflow encountered failures: %s",
                self.label,
                ", ".join(f"{job_id} ({msg})" for job_id, msg in self._failed.items()),
            )
        if self._skipped:
            logger.warning(
                "[%s] Workflow skipped jobs due to unmet dependencies: %s",
                self.label,
                ", ".join(
                    f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
                    for job_id, deps in self._skipped.items()
                ),
            )

    def shutdown(self) -> None:
        logger.debug("[%s] Global pool in use - shutdown handled by GlobalJobManager", self.label)

    def _sync_parallel_flag(self) -> None:
        self._parallel_enabled = (
            self.pool.max_concurrent_jobs > 1 and self.total_cores > 1
        )

    def _submit(self, job: WorkflowJob) -> None:
        def runner(*_args, **kwargs):
            cores = kwargs.get('cores', job.cores_optimal)
            pool_snapshot = kwargs.get('pool_snapshot')
            start_time = time.time()
            usage_suffix = ""
            if pool_snapshot:
                used, total = pool_snapshot
                try:
                    used_int = int(used)
                    total_int = int(total)
                except (TypeError, ValueError):
                    usage_suffix = ""
                else:
                    usage_suffix = f"; {used_int}/{total_int} cores used"
            logger.info(
                "[%s] Starting %s with %d cores (%s%s)",
                self.label,
                job.job_id,
                cores,
                job.description,
                usage_suffix,
            )
            try:
                job.work(cores)
            except Exception as exc:  # noqa: BLE001
                self._mark_failed(job.job_id, exc)
                raise
            else:
                duration = time.time() - start_time
                self._record_duration(job, duration)
                logger.info("[%s] Job %s completed", self.label, job.job_id)
                self._mark_completed(job.job_id)

        pool_job = PoolJob(
            job_id=job.job_id,
            cores_min=job.cores_min,
            cores_optimal=job.cores_optimal,
            cores_max=job.cores_max,
            memory_mb=job.memory_mb,
            priority=job.priority,
            execute_func=runner,
            args=(),
            kwargs={},
            estimated_duration=job.estimated_duration,
        )
        pool_job.suppress_pool_logs = True

        self.pool.submit_job(pool_job)

    @staticmethod
    def _job_order_key(job: WorkflowJob):
        """Provide a stable sort key so lower-index FoBs dispatch first."""
        digits = re.findall(r"\d+", job.job_id)
        numeric = int(digits[0]) if digits else 0
        priority_value = getattr(job.priority, "value", 0)
        return (priority_value, numeric, job.job_id)

    def _mark_completed(self, job_id: str) -> None:
        with self._lock:
            self._completed.add(job_id)
            self._event.set()

    def _mark_failed(self, job_id: str, exc: Exception) -> None:
        message = f"{exc.__class__.__name__}: {exc}"
        with self._lock:
            self._failed[job_id] = message
            self._event.set()

    def _mark_skipped(self, job_id: str, missing: Iterable[str]) -> None:
        with self._lock:
            self._skipped[job_id] = list(missing)
            self._event.set()

    def _format_failure(self) -> str:
        parts = [f"{job_id}: {message}" for job_id, message in self._failed.items()]
        return f"Workflow failures ({self.label}): " + "; ".join(parts)

    @property
    def completed_jobs(self) -> Set[str]:
        with self._lock:
            return set(self._completed)

    @property
    def failed_jobs(self) -> Dict[str, str]:
        with self._lock:
            return dict(self._failed)

    @property
    def skipped_jobs(self) -> Dict[str, List[str]]:
        with self._lock:
            return {job_id: list(deps) for job_id, deps in self._skipped.items()}

    def _base_share(self) -> int:
        if not self._parallel_enabled:
            return self.total_cores
        share = max(1, self.total_cores // max(1, self.max_jobs))
        if self.total_cores > 2:
            share = max(2, share)
        return min(self.total_cores, share)

    def _apply_exclusive_allocation(self, job: WorkflowJob) -> None:
        target = max(1, self.total_cores)
        if job.cores_optimal >= target and job.cores_max >= target:
            return

        job.cores_max = target
        job.cores_optimal = target
        job.memory_mb = target * self.maxcore_mb

        logger.info(
            "[%s] Exclusive allocation for %s (%s) → %d cores",
            self.label,
            job.job_id,
            job.description,
            target,
        )

    def _auto_tune_job(self, job: WorkflowJob) -> None:
        if not self._parallel_enabled:
            job.cores_min = job.cores_max = job.cores_optimal = self.total_cores
            job.memory_mb = job.cores_optimal * self.maxcore_mb
            return

        base_share = self._base_share()
        hint = f"{job.job_id} {job.description}".lower()
        duration_hint = self._get_duration_hint(job)
        suggestion = self._suggest_optimal_from_hint(
            hint,
            base_share,
            job.cores_max,
            job.cores_min,
            duration_hint,
        )

        min_required = self._minimum_required_cores(hint, base_share, duration_hint)
        job.cores_min = max(job.cores_min, min(min_required, suggestion, job.cores_max))

        if job.cores_optimal >= job.cores_max:
            job.cores_optimal = suggestion
        else:
            job.cores_optimal = max(job.cores_min, min(job.cores_optimal, suggestion))

    def enforce_sequential_allocation(self) -> None:
        """Force all jobs to use the full PAL when running sequentially."""
        self._parallel_enabled = False
        for job in self._jobs.values():
            full = max(1, self.total_cores)
            job.cores_min = full
            job.cores_max = full
            job.cores_optimal = full
            job.memory_mb = job.cores_optimal * self.maxcore_mb

    def _suggest_optimal_from_hint(
        self,
        hint: str,
        base_share: int,
        cores_max: int,
        cores_min: int,
        duration_hint: Optional[float],
    ) -> int:
        hint_lc = hint.lower()

        light_tokens = ("absorption", "emission", "spectrum", "td-dft", "td dft", "tddft")
        heavy_tokens = ("optimization", "freq", "frequency", "geometry", "ox", "red", "initial")

        suggestion = base_share

        if "fob" in hint_lc or "occ_" in hint_lc:
            cap = self._foB_cap()
            suggestion = min(cap, base_share)
        elif any(token in hint_lc for token in light_tokens):
            suggestion = max(1, base_share // 2)
        elif any(token in hint_lc for token in heavy_tokens):
            extra = max(1, self.total_cores // max(2, self.max_jobs))
            suggestion = base_share + extra

        suggestion = self._apply_duration_bias(
            suggestion,
            cores_min,
            cores_max,
            duration_hint,
            hint_lc,
        )

        return max(cores_min, min(cores_max, suggestion))

    def _foB_cap(self) -> int:
        if self.total_cores <= 8:
            return max(2, self.total_cores)

        default_cap = 16
        pal = max(1, _parse_int(self.config.get('PAL'), fallback=self.total_cores))
        if pal >= 48:
            default_cap = 24
        if pal >= 64:
            default_cap = 32

        return min(default_cap, self.total_cores)

    def _apply_duration_bias(
        self,
        suggestion: int,
        cores_min: int,
        cores_max: int,
        duration_hint: Optional[float],
        hint_lc: str,
    ) -> int:
        if duration_hint is None:
            return suggestion

        if duration_hint <= 45:
            adjusted = max(cores_min, suggestion // 2)
            if "fob" in hint_lc or "occ_" in hint_lc:
                adjusted = max(cores_min, min(adjusted, self._foB_cap()))
            return adjusted

        if duration_hint >= 180:
            boost = max(1, suggestion // 2)
            adjusted = min(cores_max, suggestion + boost)
            if "fob" in hint_lc or "occ_" in hint_lc:
                adjusted = min(adjusted, self._foB_cap())
            return adjusted

        return suggestion

    def _minimum_required_cores(
        self,
        hint: str,
        base_share: int,
        duration_hint: Optional[float],
    ) -> int:
        hint_lc = hint.lower()
        light_tokens = ("absorption", "emission", "spectrum", "td-dft", "td dft", "tddft")
        heavy_tokens = ("optimization", "freq", "frequency", "geometry", "ox", "red", "initial")

        min_required = 2

        if "fob" in hint_lc or "occ_" in hint_lc:
            min_required = max(4, min(self._foB_cap(), base_share))
        elif any(token in hint_lc for token in heavy_tokens):
            min_required = max(4, base_share)
        elif any(token in hint_lc for token in light_tokens):
            min_required = 2
        else:
            min_required = max(2, base_share // 2)

        if duration_hint is not None and duration_hint >= 180:
            min_required = max(min_required, min(self.total_cores // 2, base_share))

        return min_required

    def _get_duration_hint(self, job: WorkflowJob) -> Optional[float]:
        history = JOB_DURATION_HISTORY.get(self._duration_key(job))
        if not history:
            return None
        if len(history) == 1:
            return history[0]
        try:
            return statistics.median(history)
        except statistics.StatisticsError:  # pragma: no cover - defensive fallback
            return sum(history) / len(history)

    def _duration_key(self, job: WorkflowJob) -> str:
        return f"{self.label}:{job.job_id}:{job.description}".lower()

    def _record_duration(self, job: WorkflowJob, duration: float) -> None:
        key = self._duration_key(job)
        JOB_DURATION_HISTORY[key].append(duration)
        logger.debug(
            "[%s] Duration recorded for %s: %.1fs (samples=%d)",
            self.label,
            job.job_id,
            duration,
            len(JOB_DURATION_HISTORY[key]),
        )


def normalize_parallel_token(value: Any, default: str = "auto") -> str:
    token = str(value).strip().lower() if value not in (None, "") else default
    if token in ("no", "false", "off", "0"):  # explicit disable
        return "disable"
    if token in ("yes", "true", "on", "1"):  # explicit enable
        return "enable"
    return "auto"


def estimate_parallel_width(jobs: Iterable[WorkflowJob]) -> int:
    job_list = list(jobs)
    if not job_list:
        return 0

    job_ids = {job.job_id for job in job_list}
    dependency_map: Dict[str, Set[str]] = {
        job.job_id: set(dep for dep in job.dependencies if dep in job_ids)
        for job in job_list
    }

    completed: Set[str] = set()
    remaining = set(job_ids)
    max_width = 0
    guard = 0

    while remaining and guard <= len(job_ids) * 2:
        ready = {job_id for job_id in remaining if dependency_map[job_id] <= completed}
        if not ready:
            break
        max_width = max(max_width, len(ready))
        completed.update(ready)
        remaining -= ready
        guard += 1

    if remaining:
        return max_width or 1

    return max(max_width, 1)


def jobs_have_parallel_potential(jobs: Iterable[WorkflowJob]) -> bool:
    return estimate_parallel_width(jobs) > 1


def determine_effective_slots(
    total_cores: int,
    jobs: Iterable[WorkflowJob],
    requested_slots: int,
    width: int,
) -> int:
    width = max(1, width)
    requested = requested_slots if requested_slots > 0 else width
    baseline = max(1, min(width, requested))

    job_list = list(jobs)
    if not job_list or baseline >= width:
        return baseline

    min_opt = min(max(job.cores_optimal, job.cores_min) for job in job_list)
    capacity_limit = max(1, total_cores // max(1, min_opt))

    light_threshold = max(2, total_cores // 8)
    light_jobs = sum(1 for job in job_list if job.cores_optimal <= light_threshold)
    bonus = max(0, light_jobs // 2)

    candidate = min(width, baseline + bonus, capacity_limit)
    return max(1, max(baseline, candidate))


def execute_classic_workflows(config: Dict[str, Any], *, allow_parallel: bool, **kwargs) -> WorkflowRunResult:
    """Run classic oxidation/reduction steps via the shared workflow scheduler."""

    manager = _WorkflowManager(config, label="classic")

    try:
        _populate_classic_jobs(manager, config, kwargs)
        if not manager.has_jobs():
            logger.info("[classic] No oxidation/reduction jobs queued for execution")
            return WorkflowRunResult()

        width = estimate_parallel_width(manager._jobs.values())
        pal_jobs_cap = _parse_int(config.get('pal_jobs'), fallback=0)

        if allow_parallel:
            effective = determine_effective_slots(
                manager.total_cores,
                manager._jobs.values(),
                pal_jobs_cap,
                width,
            )
        else:
            effective = 1
            manager.enforce_sequential_allocation()

        if effective <= 0:
            effective = 1

        if manager.pool.max_concurrent_jobs != effective:
            logger.info(
                "[classic] Adjusting scheduler slots to %d (width=%d, pal_jobs=%s, allow_parallel=%s)",
                effective,
                width,
                config.get('pal_jobs'),
                allow_parallel,
            )
            manager.pool.max_concurrent_jobs = effective
            manager.max_jobs = effective
            manager._sync_parallel_flag()

        manager.run()
        return WorkflowRunResult(
            completed=set(manager.completed_jobs),
            failed=dict(manager.failed_jobs),
            skipped={job_id: list(deps) for job_id, deps in manager.skipped_jobs.items()},
        )

    except Exception as exc:  # noqa: BLE001
        logger.error("Classic workflows failed: %s", exc)
        result = WorkflowRunResult(
            completed=set(getattr(manager, 'completed_jobs', set())),
            failed=dict(getattr(manager, 'failed_jobs', {}) or {}),
            skipped={
                job_id: list(deps)
                for job_id, deps in (getattr(manager, 'skipped_jobs', {}) or {}).items()
            },
        )
        result.failed.setdefault('scheduler_error', f"{exc.__class__.__name__}: {exc}")
        return result

    finally:
        manager.shutdown()


def execute_manually_workflows(config: Dict[str, Any], *, allow_parallel: bool, **kwargs) -> WorkflowRunResult:
    """Run manual oxidation/reduction steps via the shared workflow scheduler."""

    manager = _WorkflowManager(config, label="manually")

    try:
        _populate_manual_jobs(manager, config, kwargs)
        if not manager.has_jobs():
            logger.info("[manually] No oxidation/reduction jobs queued for execution")
            return WorkflowRunResult()

        width = estimate_parallel_width(manager._jobs.values())
        pal_jobs_cap = _parse_int(config.get('pal_jobs'), fallback=0)

        if allow_parallel:
            effective = determine_effective_slots(
                manager.total_cores,
                manager._jobs.values(),
                pal_jobs_cap,
                width,
            )
        else:
            effective = 1
            manager.enforce_sequential_allocation()

        if effective <= 0:
            effective = 1

        if manager.pool.max_concurrent_jobs != effective:
            logger.info(
                "[manually] Adjusting scheduler slots to %d (width=%d, pal_jobs=%s, allow_parallel=%s)",
                effective,
                width,
                config.get('pal_jobs'),
                allow_parallel,
            )
            manager.pool.max_concurrent_jobs = effective
            manager.max_jobs = effective
            manager._sync_parallel_flag()

        manager.run()
        return WorkflowRunResult(
            completed=set(manager.completed_jobs),
            failed=dict(manager.failed_jobs),
            skipped={job_id: list(deps) for job_id, deps in manager.skipped_jobs.items()},
        )

    except Exception as exc:  # noqa: BLE001
        logger.error("Manual workflows failed: %s", exc)
        result = WorkflowRunResult(
            completed=set(getattr(manager, 'completed_jobs', set())),
            failed=dict(getattr(manager, 'failed_jobs', {}) or {}),
            skipped={
                job_id: list(deps)
                for job_id, deps in (getattr(manager, 'skipped_jobs', {}) or {}).items()
            },
        )
        result.failed.setdefault('scheduler_error', f"{exc.__class__.__name__}: {exc}")
        return result

    finally:
        manager.shutdown()


def _populate_classic_jobs(manager: _WorkflowManager, config: Dict[str, Any], kwargs: Dict[str, Any]) -> None:
    solvents = kwargs['solvent']
    metals = kwargs['metals']
    metal_basis = kwargs['metal_basisset']
    main_basis = kwargs['main_basisset']
    additions = kwargs['additions']
    total_electrons_txt = kwargs['total_electrons_txt']
    include_excited = bool(kwargs.get('include_excited_jobs', False))

    base_charge = _parse_int(config.get('charge'))
    base_multiplicity = _parse_int(kwargs.get('ground_multiplicity'), fallback=1)

    initial_job_id: Optional[str] = None

    def _add_job(job_id: str, description: str, work: Callable[[int], None],
                 dependencies: Optional[Set[str]] = None,
                 preferred_opt: Optional[int] = None) -> None:
        cores_min, cores_opt, cores_max = manager.derive_core_bounds(preferred_opt)
        manager.add_job(
            WorkflowJob(
                job_id=job_id,
                work=work,
                description=description,
                dependencies=dependencies or set(),
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )

    if include_excited and str(config.get('calc_initial', 'yes')).strip().lower() == 'yes':
        input_path = kwargs.get('input_file_path')
        output_initial = kwargs.get('output_initial', 'initial.inp')

        def run_initial(cores: int) -> None:
            if not input_path:
                raise RuntimeError("Input file path for classic initial job missing")
            read_and_modify_file_1(
                input_path,
                output_initial,
                base_charge,
                base_multiplicity,
                solvents,
                metals,
                metal_basis,
                main_basis,
                config,
                additions,
            )
            _update_pal_block(output_initial, cores)
            run_orca(output_initial, 'initial.out')
            if not _verify_orca_output('initial.out'):
                raise RuntimeError('ORCA terminated abnormally for initial.out')
            run_IMAG(
                'initial.out',
                'initial',
                base_charge,
                base_multiplicity,
                solvents,
                metals,
                config,
                main_basis,
                metal_basis,
                additions,
            )

        initial_job_id = 'classic_initial'
        _add_job(initial_job_id, 'initial frequency & geometry optimization', run_initial)

    if include_excited and str(config.get('absorption_spec', 'no')).strip().lower() == 'yes':
        additions_tddft = config.get('additions_TDDFT', '')

        def run_absorption(cores: int) -> None:
            read_xyz_and_create_input2(
                'initial.xyz',
                'absorption_td.inp',
                base_charge,
                1,
                solvents,
                metals,
                config,
                main_basis,
                metal_basis,
                additions_tddft,
            )
            _update_pal_block('absorption_td.inp', cores)
            run_orca('absorption_td.inp', 'absorption_spec.out')
            if not _verify_orca_output('absorption_spec.out'):
                raise RuntimeError('ORCA terminated abnormally for absorption_spec.out')

        deps = {initial_job_id} if initial_job_id else set()
        _add_job('classic_absorption', 'TD-DFT absorption spectrum', run_absorption, deps, preferred_opt=manager.total_cores // 2 or None)

    if include_excited and str(config.get('E_00', 'no')).strip().lower() == 'yes':
        excitation = str(config.get('excitation', '')).lower()
        additions_tddft = config.get('additions_TDDFT', '')

        if 't' in excitation:

            def run_t1_state(cores: int) -> None:
                read_xyz_and_create_input3(
                    'initial.xyz',
                    't1_state_opt.inp',
                    base_charge,
                    3,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block('t1_state_opt.inp', cores)
                run_orca('t1_state_opt.inp', 't1_state_opt.out')
                if not _verify_orca_output('t1_state_opt.out'):
                    raise RuntimeError('ORCA terminated abnormally for t1_state_opt.out')

            deps = {initial_job_id} if initial_job_id else set()
            _add_job('classic_t1_state', 'T1 geometry optimization', run_t1_state, deps)

            if str(config.get('emission_spec', 'no')).strip().lower() == 'yes':

                def run_t1_emission(cores: int) -> None:
                    read_xyz_and_create_input2(
                        't1_state_opt.xyz',
                        'emission_t1.inp',
                        base_charge,
                        1,
                        solvents,
                        metals,
                        config,
                        main_basis,
                        metal_basis,
                        additions_tddft,
                    )
                    _update_pal_block('emission_t1.inp', cores)
                    run_orca('emission_t1.inp', 'emission_t1.out')
                    if not _verify_orca_output('emission_t1.out'):
                        raise RuntimeError('ORCA terminated abnormally for emission_t1.out')

                _add_job('classic_t1_emission', 'T1 emission spectrum', run_t1_emission, {'classic_t1_state'}, preferred_opt=manager.total_cores // 2 or None)

        if 's' in excitation:

            def run_s1_state(cores: int) -> None:
                create_s1_optimization_input(
                    'initial.xyz',
                    's1_state_opt.inp',
                    base_charge,
                    1,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block('s1_state_opt.inp', cores)
                run_orca('s1_state_opt.inp', 's1_state_opt.out')
                if not _verify_orca_output('s1_state_opt.out'):
                    raise RuntimeError('ORCA terminated abnormally for s1_state_opt.out')

            deps = {initial_job_id} if initial_job_id else set()
            _add_job('classic_s1_state', 'S1 geometry optimization', run_s1_state, deps)

            if str(config.get('emission_spec', 'no')).strip().lower() == 'yes':

                def run_s1_emission(cores: int) -> None:
                    read_xyz_and_create_input2(
                        's1_state_opt.xyz',
                        'emission_s1.inp',
                        base_charge,
                        1,
                        solvents,
                        metals,
                        config,
                        main_basis,
                        metal_basis,
                        additions_tddft,
                    )
                    _update_pal_block('emission_s1.inp', cores)
                    run_orca('emission_s1.inp', 'emission_s1.out')
                    if not _verify_orca_output('emission_s1.out'):
                        raise RuntimeError('ORCA terminated abnormally for emission_s1.out')

                _add_job('classic_s1_emission', 'S1 emission spectrum', run_s1_emission, {'classic_s1_state'}, preferred_opt=manager.total_cores // 2 or None)

    ox_sources = {1: kwargs['xyz_file'], 2: kwargs['xyz_file4'], 3: kwargs['xyz_file8']}
    ox_inputs = {1: kwargs['output_file5'], 2: kwargs['output_file9'], 3: kwargs['output_file10']}
    ox_outputs = {1: "ox_step_1.out", 2: "ox_step_2.out", 3: "ox_step_3.out"}

    red_sources = {1: kwargs['xyz_file'], 2: kwargs['xyz_file2'], 3: kwargs['xyz_file3']}
    red_inputs = {1: kwargs['output_file6'], 2: kwargs['output_file7'], 3: kwargs['output_file8']}
    red_outputs = {1: "red_step_1.out", 2: "red_step_2.out", 3: "red_step_3.out"}

    for step in (1, 2, 3):
        if not _step_enabled(config.get('oxidation_steps'), step):
            continue

        dependencies = {f"classic_ox{step - 1}"} if step > 1 else set()
        if initial_job_id:
            dependencies.add(initial_job_id)
        cores_min, cores_opt, cores_max = manager.derive_core_bounds()

        def make_work(idx: int) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                charge = base_charge + idx
                total_electrons = total_electrons_txt - charge
                multiplicity = 1 if total_electrons % 2 == 0 else 2

                read_xyz_and_create_input3(
                    ox_sources[idx],
                    ox_inputs[idx],
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block(ox_inputs[idx], cores)
                run_orca(ox_inputs[idx], ox_outputs[idx])
                if not _verify_orca_output(ox_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {ox_outputs[idx]}")

            return _work

        manager.add_job(
            WorkflowJob(
                job_id=f"classic_ox{step}",
                work=make_work(step),
                description=f"oxidation step {step}",
                dependencies=dependencies,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )

    for step in (1, 2, 3):
        if not _step_enabled(config.get('reduction_steps'), step):
            continue

        dependencies = {f"classic_red{step - 1}"} if step > 1 else set()
        if initial_job_id:
            dependencies.add(initial_job_id)
        cores_min, cores_opt, cores_max = manager.derive_core_bounds()

        def make_work(idx: int) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                charge = base_charge - idx
                total_electrons = total_electrons_txt - charge
                multiplicity = 1 if total_electrons % 2 == 0 else 2

                read_xyz_and_create_input3(
                    red_sources[idx],
                    red_inputs[idx],
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block(red_inputs[idx], cores)
                run_orca(red_inputs[idx], red_outputs[idx])
                if not _verify_orca_output(red_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {red_outputs[idx]}")

            return _work

        manager.add_job(
            WorkflowJob(
                job_id=f"classic_red{step}",
                work=make_work(step),
                description=f"reduction step {step}",
                dependencies=dependencies,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )


def _populate_manual_jobs(manager: _WorkflowManager, config: Dict[str, Any], kwargs: Dict[str, Any]) -> None:
    solvents = kwargs['solvent']
    metals = kwargs['metals']
    metal_basis = kwargs['metal_basisset']
    main_basis = kwargs['main_basisset']
    total_electrons_txt = kwargs['total_electrons_txt']
    include_excited = bool(kwargs.get('include_excited_jobs', False))
    base_charge = _parse_int(config.get('charge'))
    base_multiplicity = _parse_int(kwargs.get('ground_multiplicity'), fallback=1)
    ground_additions = kwargs.get('ground_additions', '')
    initial_job_id: Optional[str] = None

    def _add_job(job_id: str, description: str, work: Callable[[int], None],
                 dependencies: Optional[Set[str]] = None,
                 preferred_opt: Optional[int] = None) -> None:
        cores_min, cores_opt, cores_max = manager.derive_core_bounds(preferred_opt)
        manager.add_job(
            WorkflowJob(
                job_id=job_id,
                work=work,
                description=description,
                dependencies=dependencies or set(),
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )

    if include_excited:
        input_path = kwargs.get('input_file_path')
        output_initial = kwargs.get('output_initial', 'initial.inp')

        def run_initial(cores: int) -> None:
            if not input_path:
                raise RuntimeError('Input file path missing for manual initial job')
            read_and_modify_file_1(
                input_path,
                output_initial,
                base_charge,
                base_multiplicity,
                solvents,
                metals,
                metal_basis,
                main_basis,
                config,
                ground_additions,
            )
            _update_pal_block(output_initial, cores)
            run_orca(output_initial, 'initial.out')
            if not _verify_orca_output('initial.out'):
                raise RuntimeError('ORCA terminated abnormally for initial.out')
            run_IMAG(
                'initial.out',
                'initial',
                base_charge,
                base_multiplicity,
                solvents,
                metals,
                config,
                main_basis,
                metal_basis,
                ground_additions,
            )

        initial_job_id = 'manual_initial'
        _add_job(initial_job_id, 'manual initial frequency job', run_initial)

        additions_td = config.get('additions_TDDFT', '')

        def run_absorption(cores: int) -> None:
            read_xyz_and_create_input2(
                'initial.xyz',
                'absorption_td.inp',
                base_charge,
                1,
                solvents,
                metals,
                config,
                main_basis,
                metal_basis,
                additions_td,
            )
            _update_pal_block('absorption_td.inp', cores)
            run_orca('absorption_td.inp', 'absorption_spec.out')
            if not _verify_orca_output('absorption_spec.out'):
                raise RuntimeError('ORCA terminated abnormally for absorption_spec.out')

        deps_abs = {initial_job_id} if initial_job_id else set()
        _add_job('manual_absorption', 'manual TD-DFT absorption', run_absorption, deps_abs, preferred_opt=manager.total_cores // 2 or None)

        if str(config.get('E_00', 'no')).strip().lower() == 'yes':
            excitation = str(config.get('excitation', '')).lower()

            if 't' in excitation:
                add_t1 = _extract_manual_additions(config.get('additions_T1', '')) or ground_additions

                def run_t1_state(cores: int) -> None:
                    read_xyz_and_create_input3(
                        'initial.xyz',
                        't1_state_opt.inp',
                        base_charge,
                        3,
                        solvents,
                        metals,
                        metal_basis,
                        main_basis,
                        config,
                        add_t1,
                    )
                    _update_pal_block('t1_state_opt.inp', cores)
                    run_orca('t1_state_opt.inp', 't1_state_opt.out')
                    if not _verify_orca_output('t1_state_opt.out'):
                        raise RuntimeError('ORCA terminated abnormally for t1_state_opt.out')

                deps_t1 = {initial_job_id} if initial_job_id else set()
                _add_job('manual_t1_state', 'manual T1 geometry optimization', run_t1_state, deps_t1)

                if str(config.get('emission_spec', 'no')).strip().lower() == 'yes':

                    def run_t1_emission(cores: int) -> None:
                        read_xyz_and_create_input2(
                            't1_state_opt.xyz',
                            'emission_t1.inp',
                            base_charge,
                            1,
                            solvents,
                            metals,
                            config,
                            main_basis,
                            metal_basis,
                            additions_td,
                        )
                        _update_pal_block('emission_t1.inp', cores)
                        run_orca('emission_t1.inp', 'emission_t1.out')
                        if not _verify_orca_output('emission_t1.out'):
                            raise RuntimeError('ORCA terminated abnormally for emission_t1.out')

                    _add_job('manual_t1_emission', 'manual T1 emission spectrum', run_t1_emission, {'manual_t1_state'}, preferred_opt=manager.total_cores // 2 or None)

            if 's' in excitation:
                add_s1 = _extract_manual_additions(config.get('additions_S1', '')) or ground_additions

                def run_s1_state(cores: int) -> None:
                    create_s1_optimization_input(
                        'initial.xyz',
                        's1_state_opt.inp',
                        base_charge,
                        1,
                        solvents,
                        metals,
                        metal_basis,
                        main_basis,
                        config,
                        add_s1,
                    )
                    _update_pal_block('s1_state_opt.inp', cores)
                    run_orca('s1_state_opt.inp', 's1_state_opt.out')
                    if not _verify_orca_output('s1_state_opt.out'):
                        raise RuntimeError('ORCA terminated abnormally for s1_state_opt.out')

                deps_s1 = {initial_job_id} if initial_job_id else set()
                _add_job('manual_s1_state', 'manual S1 geometry optimization', run_s1_state, deps_s1)

                if str(config.get('emission_spec', 'no')).strip().lower() == 'yes':

                    def run_s1_emission(cores: int) -> None:
                        read_xyz_and_create_input2(
                            's1_state_opt.xyz',
                            'emission_s1.inp',
                            base_charge,
                            1,
                            solvents,
                            metals,
                            config,
                            main_basis,
                            metal_basis,
                            additions_td,
                        )
                        _update_pal_block('emission_s1.inp', cores)
                        run_orca('emission_s1.inp', 'emission_s1.out')
                        if not _verify_orca_output('emission_s1.out'):
                            raise RuntimeError('ORCA terminated abnormally for emission_s1.out')

                    _add_job('manual_s1_emission', 'manual S1 emission spectrum', run_s1_emission, {'manual_s1_state'}, preferred_opt=manager.total_cores // 2 or None)

    ox_sources = {1: kwargs['xyz_file'], 2: kwargs['xyz_file4'], 3: kwargs['xyz_file8']}
    ox_inputs = {1: kwargs['output_file5'], 2: kwargs['output_file9'], 3: kwargs['output_file10']}
    ox_outputs = {1: "ox_step_1.out", 2: "ox_step_2.out", 3: "ox_step_3.out"}

    red_sources = {1: kwargs['xyz_file'], 2: kwargs['xyz_file2'], 3: kwargs['xyz_file3']}
    red_inputs = {1: kwargs['output_file6'], 2: kwargs['output_file7'], 3: kwargs['output_file8']}
    red_outputs = {1: "red_step_1.out", 2: "red_step_2.out", 3: "red_step_3.out"}

    for step in (1, 2, 3):
        if not _step_enabled(config.get('oxidation_steps'), step):
            continue

        dependencies = {f"manual_ox{step - 1}"} if step > 1 else set()
        if initial_job_id:
            dependencies.add(initial_job_id)
        cores_min, cores_opt, cores_max = manager.derive_core_bounds()

        additions_key = f"additions_ox{step}"
        multiplicity_key = f"multiplicity_ox{step}"

        def make_work(idx: int, add_key: str, mult_key: str) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                charge = base_charge + idx
                multiplicity = _parse_int(config.get(mult_key), fallback=1)
                additions = _extract_manual_additions(config.get(add_key, ""))

                read_xyz_and_create_input3(
                    ox_sources[idx],
                    ox_inputs[idx],
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block(ox_inputs[idx], cores)
                run_orca(ox_inputs[idx], ox_outputs[idx])
                if not _verify_orca_output(ox_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {ox_outputs[idx]}")

            return _work

        manager.add_job(
            WorkflowJob(
                job_id=f"manual_ox{step}",
                work=make_work(step, additions_key, multiplicity_key),
                description=f"manual oxidation step {step}",
                dependencies=dependencies,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )

    for step in (1, 2, 3):
        if not _step_enabled(config.get('reduction_steps'), step):
            continue

        dependencies = {f"manual_red{step - 1}"} if step > 1 else set()
        if initial_job_id:
            dependencies.add(initial_job_id)
        cores_min, cores_opt, cores_max = manager.derive_core_bounds()

        additions_key = f"additions_red{step}"
        multiplicity_key = f"multiplicity_red{step}"

        def make_work(idx: int, add_key: str, mult_key: str) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                charge = base_charge - idx
                multiplicity = _parse_int(config.get(mult_key), fallback=1)
                additions = _extract_manual_additions(config.get(add_key, ""))

                read_xyz_and_create_input3(
                    red_sources[idx],
                    red_inputs[idx],
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    additions,
                )
                _update_pal_block(red_inputs[idx], cores)
                run_orca(red_inputs[idx], red_outputs[idx])
                if not _verify_orca_output(red_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {red_outputs[idx]}")

            return _work

        manager.add_job(
            WorkflowJob(
                job_id=f"manual_red{step}",
                work=make_work(step, additions_key, multiplicity_key),
                description=f"manual reduction step {step}",
                dependencies=dependencies,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )


def _parse_int(value: Any, fallback: int = 0) -> int:
    try:
        return int(str(value).strip())
    except Exception:  # noqa: BLE001
        return fallback


def _normalize_tokens(raw: Any) -> Set[str]:
    if not raw:
        return set()
    if isinstance(raw, str):
        parts = re.split(r"[;,\s]+", raw.strip())
    elif isinstance(raw, Iterable):
        parts = []
        for item in raw:
            if item is None:
                continue
            parts.extend(re.split(r"[;,\s]+", str(item)))
    else:
        parts = [str(raw)]
    return {token for token in (part.strip() for part in parts) if token}


def _step_enabled(step_config: Any, step: int) -> bool:
    tokens = _normalize_tokens(step_config)
    return str(step) in tokens


def _update_pal_block(input_path: str, cores: int) -> None:
    try:
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as stream:
            lines = stream.readlines()
    except FileNotFoundError as exc:
        raise RuntimeError(f"Input file '{input_path}' missing") from exc

    pal_line = f"%pal nprocs {cores} end\n"
    replaced = False

    for idx, line in enumerate(lines):
        if line.strip().startswith('%pal'):
            lines[idx] = pal_line
            replaced = True
            break

    if not replaced:
        insert_idx = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('%') and not stripped.startswith('%pal'):
                insert_idx = idx + 1
            elif stripped and not stripped.startswith('%'):
                break
        lines.insert(insert_idx, pal_line)

    with open(input_path, 'w', encoding='utf-8') as stream:
        stream.writelines(lines)


def _verify_orca_output(path: str) -> bool:
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as stream:
            return "ORCA TERMINATED NORMALLY" in stream.read()
    except FileNotFoundError:
        return False


def _extract_manual_additions(raw: Any) -> str:
    if raw is None:
        return ""
    if isinstance(raw, str):
        value = raw.strip()
        if not value:
            return ""
        if re.fullmatch(r"\d+,\d+", value):
            return f"%scf BrokenSym {value} end"
        return value
    if isinstance(raw, Iterable):
        values = [str(item).strip() for item in raw if str(item).strip()]
        if not values:
            return ""
        return f"%scf BrokenSym {','.join(values)} end"
    return str(raw)
