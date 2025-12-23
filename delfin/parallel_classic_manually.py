"""Parallel execution helpers for DELFIN classic and manually modes."""

from __future__ import annotations

import logging
import math
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
from delfin.esd_input_generator import append_properties_of_interest_jobs
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.imag import run_IMAG
from delfin.xyz_io import read_and_modify_file_1, read_xyz_and_create_input3

logger = get_logger(__name__)

JOB_DURATION_HISTORY: Dict[str, deque[float]] = defaultdict(lambda: deque(maxlen=8))

# Global registry of active workflow managers for nested workflow coordination
_ACTIVE_MANAGERS: Dict[int, '_WorkflowManager'] = {}
_ACTIVE_MANAGERS_LOCK = threading.Lock()


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
    inline: bool = False  # Run inline without reserving pool cores
    preserve_cores_optimal: bool = False  # Skip auto-tuning for explicit core allocations
    working_dir: Optional[Path] = None  # Working directory for job-specific process tracking

    # Cache original core preferences so dynamic scheduling can adjust per run.
    base_cores_min: int = field(init=False, repr=False)
    base_cores_optimal: int = field(init=False, repr=False)
    base_cores_max: int = field(init=False, repr=False)

    def __post_init__(self) -> None:
        # The actual values will be set by _WorkflowManager.add_job once PAL is known.
        self.base_cores_min = self.cores_min
        self.base_cores_optimal = self.cores_optimal
        self.base_cores_max = self.cores_max


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

        # Fallback safety: if user explicitly disabled parallel workflows,
        # force a single-slot scheduler even if pool/max_jobs suggest otherwise.
        if normalize_parallel_token(self.config.get('parallel_workflows', 'auto')) == "disable":
            if self.pool.max_concurrent_jobs != 1 or self.max_jobs != 1:
                logger.info("[%s] parallel_workflows=disable → enforcing sequential scheduling", self.label)
            self.pool.max_concurrent_jobs = 1
            self.max_jobs = 1
            self._sync_parallel_flag()

        self._jobs: Dict[str, WorkflowJob] = {}
        self._completed: Set[str] = set()
        self._failed: Dict[str, str] = {}
        self._skipped: Dict[str, List[str]] = {}
        self._inflight: Set[str] = set()
        self._job_start_times: Dict[str, float] = {}  # Track when jobs started
        self._lock = threading.RLock()
        self._event = threading.Event()
        self._completion_listeners: List[Callable[[str], None]] = []
        self._parent_manager: Optional['_WorkflowManager'] = None

        # Batch logging: accumulate jobs for consolidated output
        self._pending_job_registrations: List[str] = []
        self._last_registration_log_time: float = 0.0

        # Register this manager globally for nested workflow coordination
        # If there's already a manager registered, this is a nested workflow
        with _ACTIVE_MANAGERS_LOCK:
            # Find parent manager (the most recently created one that isn't us)
            if _ACTIVE_MANAGERS:
                # Get the most recent manager as parent (excluding self if somehow already there)
                for mgr_id, mgr in reversed(list(_ACTIVE_MANAGERS.items())):
                    if mgr_id != id(self) and mgr.label != label:
                        self._parent_manager = mgr
                        logger.debug("[%s] Detected parent manager: [%s]", label, mgr.label)
                        break
            _ACTIVE_MANAGERS[id(self)] = self

        callback = self.config.pop('_post_attach_callback', None)
        if callable(callback):
            try:
                callback(self)
            except Exception:  # noqa: BLE001
                logger.debug("[%s] Post-attach callback raised", self.label, exc_info=True)

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

    def _validate_no_cycles(self) -> None:
        """Detect circular dependencies in the job graph using DFS.

        Raises:
            ValueError: If a circular dependency is detected.
        """
        visited: Set[str] = set()
        rec_stack: Set[str] = set()
        path: List[str] = []

        def has_cycle(job_id: str) -> bool:
            """DFS to detect cycles. Returns True if cycle found."""
            visited.add(job_id)
            rec_stack.add(job_id)
            path.append(job_id)

            job = self._jobs.get(job_id)
            if job:
                for dep in job.dependencies:
                    if dep not in visited:
                        if has_cycle(dep):
                            return True
                    elif dep in rec_stack:
                        # Cycle detected! Build cycle path
                        cycle_start = path.index(dep)
                        cycle_path = path[cycle_start:] + [dep]
                        raise ValueError(
                            f"[{self.label}] Circular dependency detected: {' → '.join(cycle_path)}"
                        )

            path.pop()
            rec_stack.remove(job_id)
            return False

        # Check all jobs for cycles
        for job_id in self._jobs:
            if job_id not in visited:
                has_cycle(job_id)

    def add_job(self, job: WorkflowJob) -> None:
        with self._lock:
            if job.job_id in self._jobs:
                raise ValueError(f"Duplicate workflow job id '{job.job_id}'")

            deps = set(job.dependencies)
            job.dependencies = deps

            if job.inline:
                # Zero-cost jobs: no pool cores reserved, keep metadata at 0
                job.cores_min = job.cores_optimal = job.cores_max = 0
                job.base_cores_min = job.base_cores_optimal = job.base_cores_max = 0
                job.memory_mb = 0
            else:
                job.cores_min = max(1, min(job.cores_min, self.total_cores))
                job.cores_max = max(job.cores_min, min(job.cores_max, self.total_cores))
                job.cores_optimal = max(job.cores_min, min(job.cores_optimal, job.cores_max))

                # Remember the original preferences so we can recompute dynamic allocations per dispatch.
                job.base_cores_min = job.cores_min
                job.base_cores_optimal = job.cores_optimal
                job.base_cores_max = job.cores_max

                self._auto_tune_job(job)

                if job.memory_mb is None:
                    job.memory_mb = job.cores_optimal * self.maxcore_mb

            self._jobs[job.job_id] = job

            # Validate no circular dependencies after adding job
            try:
                self._validate_no_cycles()
            except ValueError as e:
                # Remove the job that caused the cycle
                self._jobs.pop(job.job_id, None)
                raise ValueError(f"Cannot add job {job.job_id}: {e}") from e

            # Accumulate job registration for batch logging
            self._pending_job_registrations.append(job.job_id)

            # Detailed per-job logging only in DEBUG mode
            logger.debug(
                "[%s] Registered job %s (%s); deps=%s",
                self.label,
                job.job_id,
                job.description,
                ",".join(sorted(job.dependencies)) or "none",
            )

            # Flush if we haven't logged in a while (500ms threshold)
            # This ensures quick feedback while still batching rapid additions
            time_since_last_log = time.time() - self._last_registration_log_time
            if time_since_last_log > 0.5:
                self._flush_pending_registrations()

            # Wake up scheduler to process new job
            self._event.set()

    def _flush_pending_registrations(self) -> None:
        """Log accumulated job registrations in a consolidated message."""
        if not self._pending_job_registrations:
            return

        job_count = len(self._pending_job_registrations)
        if job_count == 1:
            # Single job - log with details
            job_id = self._pending_job_registrations[0]
            job = self._jobs[job_id]
            logger.info(
                "[%s] Registered 1 job: %s (%s)",
                self.label,
                job_id,
                job.description,
            )
        else:
            # Multiple jobs - show compact summary
            job_names = ", ".join(self._pending_job_registrations[:5])
            if job_count > 5:
                job_names += f", ... (+{job_count - 5} more)"
            logger.info(
                "[%s] Registered %d jobs: %s",
                self.label,
                job_count,
                job_names,
            )

        self._pending_job_registrations.clear()
        self._last_registration_log_time = time.time()

    def register_completion_listener(self, listener: Callable[[str], None]) -> None:
        with self._lock:
            self._completion_listeners.append(listener)

    def unregister_completion_listener(self, listener: Callable[[str], None]) -> None:
        with self._lock:
            try:
                self._completion_listeners.remove(listener)
            except ValueError:
                pass

    def _notify_completion(self, job_id: str) -> None:
        listeners = list(self._completion_listeners)
        for listener in listeners:
            try:
                listener(job_id)
            except Exception:  # noqa: BLE001
                logger.debug("[%s] Completion listener raised", self.label, exc_info=True)

    def reschedule_pending(self) -> None:
        with self._lock:
            self._event.set()

    def has_jobs(self) -> bool:
        return bool(self._jobs)

    def _adjust_priorities_for_bottlenecks(self) -> None:
        """Boost priority for jobs that block many downstream jobs."""
        try:
            from delfin.job_priority import adjust_job_priorities
            boosted = adjust_job_priorities(self._jobs, bottleneck_threshold=3)
            if boosted > 0:
                logger.info(
                    "[%s] Boosted priority for %d bottleneck job(s)",
                    self.label,
                    boosted,
                )
        except Exception as e:
            logger.debug("[%s] Failed to adjust priorities: %s", self.label, e, exc_info=True)

    def _check_exclusive_bottleneck_boost(self, job: WorkflowJob) -> Optional[int]:
        """
        Check if job is an exclusive bottleneck and should get all cores.

        Returns forced_cores if boost applies, None otherwise.
        """
        try:
            from delfin.job_priority import is_exclusive_bottleneck

            # Get pending jobs
            with self._lock:
                finished = self._completed | set(self._failed) | set(self._skipped)
                inflight = set(self._inflight)
                pending_ids = set(self._jobs.keys()) - finished - inflight

            # Check if this job is exclusive bottleneck
            if is_exclusive_bottleneck(job.job_id, self._jobs, pending_ids):
                logger.info(
                    "[%s] Job %s is exclusive bottleneck → allocating max cores (%d)",
                    self.label,
                    job.job_id,
                    job.cores_max,
                )
                return job.cores_max

        except Exception as e:
            logger.debug("[%s] Failed to check exclusive bottleneck: %s", self.label, e, exc_info=True)

        return None

    def _should_wait_for_exclusive_bottleneck(self, ready_jobs: List[WorkflowJob], status: Optional[Dict]) -> bool:
        """
        Check if we should wait for an exclusive bottleneck to run alone.

        Returns True if:
        - There is an exclusive bottleneck in ready jobs
        - Other jobs are currently running
        - We should wait for them to finish so bottleneck gets all cores

        Args:
            ready_jobs: List of jobs ready to run
            status: Pool status dict

        Returns:
            True if we should wait, False otherwise
        """
        if not ready_jobs:
            return False

        # Check if any jobs are running
        jobs_running = False
        if status:
            jobs_running = (
                status.get('running_jobs', 0) > 0 or
                status.get('queued_jobs', 0) > 0
            )
        if not jobs_running and not self._inflight:
            # No jobs running, no need to wait
            return False

        try:
            from delfin.job_priority import is_exclusive_bottleneck

            # Get pending jobs (including ready ones)
            with self._lock:
                finished = self._completed | set(self._failed) | set(self._skipped)
                inflight = set(self._inflight)
                pending_ids = set(self._jobs.keys()) - finished - inflight

            # Check if any ready job is an exclusive bottleneck
            has_exclusive_bottleneck = False
            bottleneck_job = None
            for job in ready_jobs:
                if is_exclusive_bottleneck(job.job_id, self._jobs, pending_ids):
                    has_exclusive_bottleneck = True
                    bottleneck_job = job
                    break

            if not has_exclusive_bottleneck:
                return False

            # Check if waiting makes sense based on job progress
            # For long-running jobs (4-12h), only wait if almost done
            current_time = time.time()
            max_wait_time = 15 * 60  # Maximum 15 minutes willing to wait (jobs are typically 4-12h)

            # Check running job start times and estimated completion
            with self._lock:
                for job_id in list(self._inflight):
                    start_time = self._job_start_times.get(job_id)
                    if start_time is None:
                        continue

                    elapsed = current_time - start_time

                    # Get job estimated duration
                    job = self._jobs.get(job_id)
                    estimated_duration = job.estimated_duration if job else 14400.0  # Default 4h

                    # Calculate estimated remaining time
                    estimated_remaining = estimated_duration - elapsed
                    progress_ratio = elapsed / estimated_duration if estimated_duration > 0 else 0

                    # Only wait if job will finish within max_wait_time (30 min)
                    # AND has completed at least 90% of its time
                    if estimated_remaining > max_wait_time or progress_ratio < 0.90:
                        # Job would take too long to finish
                        estimated_remaining_min = estimated_remaining / 60
                        logger.debug(
                            "[%s] Not waiting for exclusive bottleneck %s: "
                            "job %s is %.0f%% done (est. %.1f min remaining, max wait: %.0f min)",
                            self.label,
                            bottleneck_job.job_id,
                            job_id,
                            progress_ratio * 100,
                            estimated_remaining_min,
                            max_wait_time / 60,
                        )
                        return False

            # Check available cores
            # For long-running jobs (4-12h), only wait if very few cores available
            if status:
                allocated = status.get('allocated_cores', 0)
                available = self.total_cores - allocated
                # If >= 25% cores already available, don't wait (too long for benefit)
                if available >= self.total_cores * 0.25:
                    logger.debug(
                        "[%s] Not waiting for exclusive bottleneck %s: "
                        "%d cores (%.0f%%) already available (jobs typically 4-12h)",
                        self.label,
                        bottleneck_job.job_id,
                        available,
                        (available / self.total_cores) * 100,
                    )
                    return False

            # All checks passed: wait for jobs to finish
            logger.info(
                "[%s] Waiting for running jobs to finish before starting exclusive bottleneck %s "
                "(will maximize core allocation)",
                self.label,
                bottleneck_job.job_id,
            )
            return True

        except Exception as e:
            logger.debug("[%s] Failed to check exclusive bottleneck wait: %s", self.label, e, exc_info=True)

        return False

    def run(self) -> None:
        # Flush any pending job registrations before starting
        self._flush_pending_registrations()

        if not self._jobs:
            logger.info("[%s] No jobs to schedule", self.label)
            return

        # Adjust job priorities based on bottleneck detection
        self._adjust_priorities_for_bottlenecks()

        pending: Dict[str, WorkflowJob] = {}
        self._sync_parallel_flag()
        logger.info(
            "[%s] Scheduling %d jobs across %d cores using GLOBAL SHARED pool (pool_id=%d)",
            self.label,
            len(self._jobs),
            self.total_cores,
            id(self.pool),
        )

        # Deadlock detection: Track last progress to detect infinite blocking
        last_progress_time = time.time()
        last_completed_count = 0
        DEADLOCK_TIMEOUT = 300  # 5 minutes without progress

        while True:
            for job_id, job in list(self._jobs.items()):
                if job_id in pending:
                    continue
                if job_id in self._completed or job_id in self._failed or job_id in self._skipped:
                    continue
                if job_id in self._inflight:
                    continue
                pending[job_id] = job

            if not pending:
                any_running = bool(self._inflight)
                if not any_running:
                    try:
                        status = self.pool.get_status()
                    except Exception:
                        status = None
                    if status:
                        any_running = (
                            status.get('running_jobs', 0) > 0
                            or status.get('queued_jobs', 0) > 0
                        )
                if any_running:
                    self._event.wait(timeout=0.1)
                    self._event.clear()
                    continue
                break

            finished_tokens = self._completed | set(self._failed) | set(self._skipped)
            ready = [
                job for job in pending.values()
                if job.dependencies <= finished_tokens
            ]

            if not ready:
                blocked = {
                    job.job_id: sorted(job.dependencies - finished_tokens)
                    for job in pending.values()
                }

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
                    # Check if we have capacity to start jobs despite other jobs running
                    # This handles the case where nested workflows are using the pool
                    if status:
                        allocated = status.get('allocated_cores', 0)
                        available = self.total_cores - allocated
                        # If we have enough cores available and there are blocked jobs,
                        # check more frequently instead of waiting the full timeout
                        if available >= 2 and blocked:
                            self._event.wait(timeout=0.1)
                            self._event.clear()
                            continue
                    self._event.wait(timeout=0.5)
                    self._event.clear()
                    continue

                if self._inflight:
                    self._event.wait(timeout=0.1)
                    self._event.clear()
                    continue

                # Deadlock detection: Check if we're stuck with pending jobs but no progress
                current_completed_count = len(self._completed) + len(self._failed) + len(self._skipped)
                if current_completed_count > last_completed_count:
                    # Progress made - reset timeout
                    last_progress_time = time.time()
                    last_completed_count = current_completed_count
                elif time.time() - last_progress_time > DEADLOCK_TIMEOUT:
                    # No progress for DEADLOCK_TIMEOUT seconds and jobs are blocked
                    logger.error(
                        "[%s] DEADLOCK DETECTED: %d jobs blocked for %.1f minutes without progress. "
                        "Blocked jobs: %s. Aborting workflow to prevent infinite hang.",
                        self.label,
                        len(blocked),
                        (time.time() - last_progress_time) / 60,
                        blocked,
                    )
                    # Mark all blocked jobs as skipped to cleanly exit
                    for job_id, missing in blocked.items():
                        self._mark_skipped(job_id, missing or ['deadlock - timeout waiting for dependencies'])
                        pending.pop(job_id, None)
                    break

                for job_id, missing in blocked.items():
                    self._mark_skipped(job_id, missing or ['unresolved dependency'])
                    pending.pop(job_id, None)
                continue

            ready.sort(key=self._job_order_key)

            failed_ids = set(self._failed)
            skipped_ids = set(self._skipped)

            # CRITICAL FIX: Get fresh pool status to ensure allocated_cores is up-to-date
            # This prevents using stale snapshots after job completions
            try:
                status = self.pool.get_status()
            except Exception:
                status = None

            # ROBUST FIX: If pool snapshot shows more cores allocated than we have inflight jobs,
            # it's stale. Use a conservative estimate based on our tracking.
            if status and self._inflight:
                snapshot_allocated = status.get('allocated_cores', 0)
                # Assume inflight jobs use at least cores_min each as lower bound
                with self._lock:
                    inflight_count = len(self._inflight)
                # If snapshot shows unreasonably high allocation, cap it
                # This handles race conditions where pool hasn't updated after completions
                if snapshot_allocated > self.total_cores or (
                    inflight_count == 0 and snapshot_allocated > 0
                ):
                    logger.debug(
                        "[%s] Pool snapshot stale (allocated=%d, inflight=%d), refreshing",
                        self.label,
                        snapshot_allocated,
                        inflight_count,
                    )
                    # Force a small delay to let pool update, then re-fetch
                    time.sleep(0.05)
                    try:
                        status = self.pool.get_status()
                    except Exception:
                        pass

            logger.debug(
                "[%s] Ready jobs=%s | completed=%s",
                self.label,
                [job.job_id for job in ready],
                sorted(self._completed),
            )

            inline_ready = [job for job in ready if job.inline]
            if inline_ready:
                for job in inline_ready:
                    self._submit(job)
                    pending.pop(job.job_id, None)
                ready = [job for job in ready if not job.inline]
                if not ready:
                    # Inline jobs finished; loop will pick up newly ready work immediately.
                    continue

            # Check for exclusive bottleneck that should run alone
            if self._should_wait_for_exclusive_bottleneck(ready, status):
                # Wait for running jobs to finish before starting exclusive bottleneck
                self._event.wait(timeout=1.0)
                self._event.clear()
                continue

            allocations, caps = self._plan_core_allocations(ready, status)

            jobs_submitted_this_round = 0
            for job in ready:
                if job.job_id not in allocations:
                    # Not enough cores to start this job in the current packing round
                    continue
                failed_prereqs = sorted(job.dependencies & failed_ids)
                skipped_prereqs = sorted(job.dependencies & skipped_ids)
                if failed_prereqs or skipped_prereqs:
                    logger.warning(
                        "[%s] Job %s starting although prerequisites failed/skipped: failed=%s, skipped=%s",
                        self.label,
                        job.job_id,
                        failed_prereqs or "none",
                        skipped_prereqs or "none",
                    )
                self._submit(
                    job,
                    allocations.get(job.job_id),
                    lock_to_forced=caps.get(job.job_id, False),
                )
                pending.pop(job.job_id, None)
                jobs_submitted_this_round += 1

            # Deadlock detection: Reset progress timer when jobs are submitted
            if jobs_submitted_this_round > 0:
                last_progress_time = time.time()

            # After submitting ready jobs, wait briefly for state changes
            # This allows the loop to react immediately when dependencies are fulfilled
            # instead of only checking every 0.5s when blocked
            if pending:
                self._event.wait(timeout=0.1)
                self._event.clear()

        self.pool.wait_for_completion()

        # Summary logging
        total_jobs = len(self._jobs)
        completed_count = len(self._completed)
        failed_count = len(self._failed)
        skipped_count = len(self._skipped)

        if failed_count == 0 and skipped_count == 0:
            logger.info(
                "[%s] ✓ Workflow completed successfully: %d/%d jobs",
                self.label,
                completed_count,
                total_jobs,
            )
        else:
            logger.info(
                "[%s] Workflow finished: %d completed, %d failed, %d skipped (total: %d)",
                self.label,
                completed_count,
                failed_count,
                skipped_count,
                total_jobs,
            )

        if self._failed:
            logger.warning(
                "[%s] Failed jobs: %s",
                self.label,
                ", ".join(f"{job_id} ({msg})" for job_id, msg in self._failed.items()),
            )
        if self._skipped:
            logger.warning(
                "[%s] Skipped jobs: %s",
                self.label,
                ", ".join(
                    f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
                    for job_id, deps in self._skipped.items()
                ),
            )

    def shutdown(self) -> None:
        logger.debug("[%s] Global pool in use - shutdown handled by GlobalJobManager", self.label)
        # Unregister from global registry
        with _ACTIVE_MANAGERS_LOCK:
            _ACTIVE_MANAGERS.pop(id(self), None)

    def _sync_parallel_flag(self) -> None:
        self._parallel_enabled = (
            self.pool.max_concurrent_jobs > 1 and self.total_cores > 1
        )

    def _submit(
        self,
        job: WorkflowJob,
        forced_cores: Optional[int] = None,
        *,
        lock_to_forced: bool = False,
    ) -> None:
        if job.inline:
            start_time = time.time()
            logger.info(
                "[%s] Running %s inline (zero-cost; %s)",
                self.label,
                job.job_id,
                job.description,
            )
            with self._lock:
                self._job_start_times[job.job_id] = start_time
            try:
                job.work(0)
            except Exception as exc:  # noqa: BLE001
                self._mark_failed(job.job_id, exc)
                raise
            else:
                duration = time.time() - start_time
                self._record_duration(job, duration)
                self._mark_completed(job.job_id)
            return

        # Check if this is an exclusive bottleneck and boost cores
        if forced_cores is None:
            forced_cores = self._check_exclusive_bottleneck_boost(job)

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

            # Only log start for non-inline jobs or in DEBUG mode
            if not job.inline:
                logger.info(
                    "[%s] Starting %s with %d cores (%s%s)",
                    self.label,
                    job.job_id,
                    cores,
                    job.description,
                    usage_suffix,
                )
            else:
                logger.debug(
                    "[%s] Starting %s (inline, zero-cost)",
                    self.label,
                    job.job_id,
                )

            try:
                job.work(cores)
            except Exception as exc:  # noqa: BLE001
                self._mark_failed(job.job_id, exc)
                raise
            else:
                duration = time.time() - start_time
                self._record_duration(job, duration)

                # Only log completion for jobs that took significant time (>5s) or failed/inline
                # This reduces noise from quick recalc jobs
                if duration > 5.0 or job.inline:
                    logger.debug("[%s] Job %s completed (%.1fs)", self.label, job.job_id, duration)
                else:
                    logger.debug("[%s] Job %s completed", self.label, job.job_id)

                self._mark_completed(job.job_id)

        pool_job = PoolJob(
            job_id=job.job_id,
            cores_min=self._resolve_min_cores(job, forced_cores, lock_to_forced),
            cores_optimal=self._resolve_opt_cores(job, forced_cores, lock_to_forced),
            cores_max=self._resolve_max_cores(job, forced_cores, lock_to_forced),
            memory_mb=self._resolve_memory(job, forced_cores, lock_to_forced),
            priority=job.priority,
            execute_func=runner,
            args=(),
            kwargs={},
            estimated_duration=job.estimated_duration,
            working_dir=job.working_dir,
        )
        pool_job.suppress_pool_logs = True

        self.pool.submit_job(pool_job)
        with self._lock:
            self._inflight.add(job.job_id)
            self._job_start_times[job.job_id] = time.time()

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
            self._inflight.discard(job_id)
            self._job_start_times.pop(job_id, None)  # Clean up start time
            self._event.set()
        self._notify_completion(job_id)
        # Notify parent manager if this is a nested workflow
        if self._parent_manager is not None:
            self._parent_manager.reschedule_pending()

    def _mark_failed(self, job_id: str, exc: Exception) -> None:
        message = f"{exc.__class__.__name__}: {exc}"
        with self._lock:
            self._failed[job_id] = message
            self._inflight.discard(job_id)
            self._job_start_times.pop(job_id, None)  # Clean up start time
            self._event.set()
        # Notify parent manager if this is a nested workflow
        if self._parent_manager is not None:
            self._parent_manager.reschedule_pending()

    def _mark_skipped(self, job_id: str, missing: Iterable[str]) -> None:
        with self._lock:
            self._skipped[job_id] = list(missing)
            self._inflight.discard(job_id)
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

    def _auto_tune_job(self, job: WorkflowJob) -> None:
        if job.inline:
            return

        # Skip auto-tuning for jobs with explicit core allocations (e.g., weighted FoB jobs)
        # This must come BEFORE parallel check to preserve asymmetric allocations
        if job.preserve_cores_optimal:
            logger.debug(
                "[%s] Preserving FoB weight-based allocation: %s=%d cores",
                self.label,
                job.job_id,
                job.cores_optimal,
            )
            return

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

    def _plan_core_allocations(
        self,
        ready_jobs: List[WorkflowJob],
        status_snapshot: Optional[Dict[str, Any]],
    ) -> tuple[Dict[str, int], Dict[str, bool]]:
        """
        Determine per-job core targets with a packing step that maximizes the
        number of concurrent ready jobs while respecting available cores.

        Returns:
            allocations: job_id -> planned cores
            caps: job_id -> whether to hard-cap job to the planned cores
        """
        if not ready_jobs:
            return {}, {}

        status_snapshot = status_snapshot or {}
        running_jobs = max(0, int(status_snapshot.get('running_jobs', 0)))
        queued_jobs = max(0, int(status_snapshot.get('queued_jobs', 0)))
        allocated_cores = max(0, int(status_snapshot.get('allocated_cores', 0)))
        total = self.total_cores

        ready_count = len(ready_jobs)

        # ROBUST FIX: Better available cores calculation with stale snapshot detection
        with self._lock:
            actual_inflight = len(self._inflight)

        # If nothing else is running, the ready jobs may use the full PAL.
        if running_jobs == 0 and queued_jobs == 0:
            available = total
        elif actual_inflight == 0 and allocated_cores > 0:
            # CRITICAL FIX: Snapshot is stale - we have no inflight jobs but snapshot shows allocation
            # This happens when jobs just completed. Use full pool.
            logger.debug(
                "[%s] Stale snapshot detected (allocated=%d but inflight=0), using full pool",
                self.label,
                allocated_cores,
            )
            available = total
        else:
            # Use snapshot's allocated_cores, but sanity-check it
            # allocated_cores should never exceed total, and should be reasonable given inflight count
            available = max(1, total - allocated_cores)
            # Guard against pathological snapshots where allocated exceeds total
            available = min(total, available)

            # Additional sanity check: if available is suspiciously low and we have few inflight jobs,
            # the snapshot is likely stale. Be more optimistic about available cores.
            if actual_inflight > 0 and available < total * 0.2:
                # Estimate a more realistic available based on typical job sizes
                # Assume each inflight job uses ~optimal cores on average
                estimated_used = min(total, actual_inflight * 12)  # 12 cores typical for ESD jobs
                estimated_available = max(1, total - estimated_used)
                if estimated_available > available:
                    logger.debug(
                        "[%s] Snapshot shows low availability (%d cores), but inflight=%d suggests more (%d cores). Using optimistic estimate.",
                        self.label,
                        available,
                        actual_inflight,
                        estimated_available,
                    )
                    available = estimated_available

        if not self._parallel_enabled:
            return {job.job_id: total for job in ready_jobs}, {
                job.job_id: True for job in ready_jobs
            }

        if ready_count == 1:
            # Single ready job: default to its optimal share; only grab the full pool
            # when nothing else can run (no other pending work and no running/queued jobs).
            with self._lock:
                unfinished = (
                    set(self._jobs.keys())
                    - self._completed
                    - set(self._failed)
                    - set(self._skipped)
                    - self._inflight
                )
            blocked_other = max(0, len(unfinished) - ready_count)
            job = ready_jobs[0]

            # If this job is an exclusive bottleneck (all pending depend on it), give it max cores.
            forced_exclusive = self._check_exclusive_bottleneck_boost(job)
            if forced_exclusive is not None:
                target = max(job.cores_min, min(available, job.cores_max, forced_exclusive))
                cap = False
                logger.debug(
                    "[%s] Exclusive bottleneck %s (%s) → %d cores (all pending depend on it)",
                    self.label,
                    job.job_id,
                    job.description,
                    target,
                )
                return {job.job_id: target}, {job.job_id: cap}

            base_target = max(job.cores_min, min(job.cores_optimal, available))
            nothing_else_runnable = (
                blocked_other == 0 and running_jobs == 0 and queued_jobs == 0
            )
            if nothing_else_runnable:
                target = max(base_target, min(available, job.cores_max))
                cap = False
                logger.debug(
                    "[%s] Exclusive allocation for %s (%s) → %d cores (full pool, nothing else runnable)",
                    self.label,
                    job.job_id,
                    job.description,
                    target,
                )
            else:
                target = base_target
                cap = False
                logger.debug(
                    "[%s] Single ready job %s (%s) capped to optimal %d cores "
                    "(blocked_other=%d, running=%d, queued=%d, available=%d, max=%d)",
                    self.label,
                    job.job_id,
                    job.description,
                    target,
                    blocked_other,
                    running_jobs,
                    queued_jobs,
                    available,
                    job.cores_max,
                )
            return {job.job_id: target}, {job.job_id: cap}

        # Multiple jobs ready: distribute available capacity as evenly as possible.
        # Fall back to base share if available cores are insufficient.
        if available < ready_count:
            available = max(ready_count, self._base_share() * ready_count)
            available = min(total, available)

        # For jobs with preserve_cores_optimal (e.g., weighted FoB jobs),
        # use their explicit cores_optimal instead of runtime-based weights
        allocations: Dict[str, int] = {}
        allocation_pool = available

        # First pass: assign preserved allocations
        preserved_jobs = [job for job in ready_jobs if job.preserve_cores_optimal]
        dynamic_jobs = [job for job in ready_jobs if not job.preserve_cores_optimal]

        for job in preserved_jobs:
            # Use the explicit core allocation from FoB weights
            share = max(job.cores_min, min(job.cores_optimal, job.cores_max, total, available))
            allocations[job.job_id] = share
            allocation_pool -= share

        if allocation_pool < 0:
            allocation_pool = 0

        # Second pass: distribute remaining cores to dynamic jobs using runtime weights
        if dynamic_jobs:
            weights: Dict[str, float] = {}
            weight_sum = 0.0
            for job in dynamic_jobs:
                duration = self._get_duration_hint(job)
                weight = max(1.0, math.sqrt(duration / 300.0)) if duration else 1.0
                weights[job.job_id] = weight
                weight_sum += weight

            per_job = max(1, allocation_pool // len(dynamic_jobs)) if allocation_pool > 0 else 1
            remainder = max(0, allocation_pool - per_job * len(dynamic_jobs))

            for idx, job in enumerate(dynamic_jobs):
                base = per_job
                if remainder > 0:
                    base += 1
                    remainder -= 1

                share = base
                if weight_sum > 0 and allocation_pool > 0:
                    proportional = int(round(allocation_pool * (weights[job.job_id] / weight_sum)))
                    share = max(base, proportional)

                # Respect job's core limits: min <= share <= max
                share = max(job.cores_min, min(share, job.cores_max, total))
                allocations[job.job_id] = share

        # Distribute remaining cores to jobs that can use more (up to their max)
        if allocation_pool > 0 and allocations:
            for job in ready_jobs:
                if allocation_pool <= 0:
                    break
                current = allocations[job.job_id]
                can_add = min(allocation_pool, job.cores_max - current)
                if can_add > 0:
                    allocations[job.job_id] = current + can_add
                    allocation_pool -= can_add
                    logger.debug(
                        "[%s] Allocated +%d cores to %s (now %d/%d cores)",
                        self.label,
                        can_add,
                        job.job_id,
                        allocations[job.job_id],
                        job.cores_max,
                    )

        # Pack allocations to ensure we start as many ready jobs as possible
        packed_allocations, caps = self._pack_allocations(
            ready_jobs,
            available,
            allocations,
        )

        return packed_allocations, caps

    def _pack_allocations(
        self,
        ready_jobs: List[WorkflowJob],
        available: int,
        proposed: Dict[str, int],
    ) -> tuple[Dict[str, int], Dict[str, bool]]:
        """
        Bin-pack ready jobs into the available cores.

        We keep jobs with smaller cores_min first to maximize concurrency,
        then distribute any remaining cores towards the proposed targets.
        """
        if available <= 0 or not ready_jobs:
            return {}, {}

        jobs_sorted = sorted(
            ready_jobs,
            key=lambda j: (self._job_order_key(j), j.cores_min),
        )

        targets: Dict[str, int] = {}
        for job in jobs_sorted:
            target = proposed.get(job.job_id, job.cores_optimal)
            target = max(job.cores_min, min(target, job.cores_max, self.total_cores))
            targets[job.job_id] = target

        # First, try to fit all jobs with their target allocations (preserves asymmetric weights)
        total_target_cores = sum(targets[job.job_id] for job in jobs_sorted)

        if total_target_cores <= available:
            # We have enough cores to give everyone their target allocation
            # This preserves the asymmetric FoB weight-based allocation
            allocations = {job.job_id: targets[job.job_id] for job in jobs_sorted}

            # CRITICAL FIX: Redistribute unused cores to maximize utilization
            # Even when all targets fit, we should use ALL available cores
            total_allocated = sum(allocations.values())
            remaining_cores = available - total_allocated
            if remaining_cores > 0:
                logger.debug(
                    "[%s] Redistributing %d unused cores to jobs that can use more",
                    self.label,
                    remaining_cores,
                )
                for job in jobs_sorted:
                    if remaining_cores <= 0:
                        break
                    current = allocations[job.job_id]
                    can_add = min(remaining_cores, job.cores_max - current)
                    if can_add > 0:
                        allocations[job.job_id] = current + can_add
                        remaining_cores -= can_add
                        logger.debug(
                            "[%s] Allocated +%d cores to %s (now %d/%d cores)",
                            self.label,
                            can_add,
                            job.job_id,
                            allocations[job.job_id],
                            job.cores_max,
                        )

            # Cap jobs with preserve_cores_optimal to enforce asymmetric allocation
            # Otherwise dynamic pool might expand them
            caps = {
                job.job_id: job.preserve_cores_optimal
                for job in jobs_sorted
            }

            logger.info(
                "[%s] ✓ Using asymmetric allocations: %s (total: %d/%d cores)",
                self.label,
                allocations,
                sum(allocations.values()),
                available,
            )

            return allocations, caps

        # Not enough cores for all targets - distribute evenly (simpler management)
        logger.debug(
            "[%s] Insufficient cores for targets (%d < %d), using even distribution",
            self.label,
            available,
            total_target_cores,
        )

        selected: List[WorkflowJob] = []
        remaining = available
        for job in jobs_sorted:
            if job.cores_min <= remaining:
                selected.append(job)
                remaining -= job.cores_min

        if not selected:
            return {}, {}

        # Distribute available cores evenly across selected jobs
        per_job = max(1, available // len(selected))
        remainder = available - (per_job * len(selected))

        allocations: Dict[str, int] = {}
        for idx, job in enumerate(selected):
            share = per_job
            if idx < remainder:
                share += 1  # Distribute remainder to first N jobs

            # Respect job's core limits
            share = max(job.cores_min, min(share, job.cores_max, self.total_cores))
            allocations[job.job_id] = share

        # After initial distribution, redistribute any unused cores to jobs that can use more
        total_allocated = sum(allocations.values())
        remaining_cores = available - total_allocated
        if remaining_cores > 0:
            for job in selected:
                if remaining_cores <= 0:
                    break
                current = allocations[job.job_id]
                can_add = min(remaining_cores, job.cores_max - current)
                if can_add > 0:
                    allocations[job.job_id] = current + can_add
                    remaining_cores -= can_add

        # Cap when we had to pack tightly (either we deferred jobs or could not reach targets)
        could_not_reach_targets = any(
            allocations[job.job_id] < targets[job.job_id] for job in selected
        )
        cap_needed = (
            could_not_reach_targets
            or len(selected) < len(ready_jobs)
        )
        caps = {job.job_id: cap_needed for job in selected}

        return allocations, caps

    def _resolve_min_cores(self, job: WorkflowJob, forced: Optional[int], lock: bool) -> int:
        # Minimum requirement stays at the job's baseline unless we explicitly lock to a packed value.
        if lock and forced is not None:
            return max(1, min(forced, self.total_cores))
        return max(1, job.base_cores_min)

    def _resolve_opt_cores(self, job: WorkflowJob, forced: Optional[int], lock: bool) -> int:
        if forced is None:
            return job.cores_optimal
        forced = max(1, min(forced, self.total_cores))
        if lock:
            return forced
        return max(job.base_cores_min, forced)

    def _resolve_max_cores(self, job: WorkflowJob, forced: Optional[int], lock: bool) -> int:
        base_max = job.base_cores_max
        if forced is None:
            return job.cores_max
        forced = max(1, min(forced, self.total_cores))
        if lock:
            return forced
        return max(base_max, job.cores_max, forced)

    def _resolve_memory(self, job: WorkflowJob, forced: Optional[int], lock: bool) -> int:
        if forced is None:
            return job.memory_mb
        resolved = max(job.base_cores_min, min(forced, self.total_cores))
        return resolved * self.maxcore_mb


    def enforce_sequential_allocation(self) -> None:
        """Force all jobs to use the full PAL when running sequentially."""
        self._parallel_enabled = False
        for job in self._jobs.values():
            full = max(1, self.total_cores)
            job.cores_min = full
            job.cores_max = full
            job.cores_optimal = full
            job.memory_mb = job.cores_optimal * self.maxcore_mb
            job.base_cores_min = full
            job.base_cores_optimal = full
            job.base_cores_max = full

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
            suggestion = base_share
            cap = self._foB_cap()
            suggestion = min(cap, suggestion)
        elif any(token in hint_lc for token in light_tokens):
            suggestion = max(1, base_share // 2)
        elif any(token in hint_lc for token in heavy_tokens):
            suggestion = base_share

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
            min_required = max(4, base_share)
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
    if token in ("no", "false", "off", "0", "disable", "disabled"):  # explicit disable
        return "disable"
    if token in ("yes", "true", "on", "1", "enable", "enabled"):  # explicit enable
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


def execute_classic_workflows(
    config: Dict[str, Any],
    *,
    allow_parallel: bool,
    scheduler: Optional["GlobalOrcaScheduler"] = None,
    **kwargs,
) -> WorkflowRunResult:
    """Run classic oxidation/reduction steps via the shared workflow scheduler."""

    # If a scheduler is provided, use its manager (which may already have ESD jobs)
    # Otherwise create a new manager
    if scheduler is not None:
        manager = scheduler.manager
    else:
        manager = _WorkflowManager(config, label="classic")

    try:
        _populate_classic_jobs(manager, config, kwargs)

        # If using a shared scheduler, check if scheduler has any jobs
        # (e.g., ESD jobs may have been added before execute_classic_workflows was called)
        if not manager.has_jobs():
            if scheduler is not None and scheduler.manager.has_jobs():
                # No classic jobs, but scheduler has other jobs (e.g., ESD) - run them
                logger.info("[classic] No oxidation/reduction jobs queued for execution")
                return scheduler.run()
            else:
                # No jobs at all
                logger.info("[classic] No oxidation/reduction jobs queued for execution")
                return WorkflowRunResult()

        jobs_snapshot = list(manager._jobs.values())

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

        if scheduler is not None:
            # Apply sequential allocation to scheduler's manager if needed
            if not allow_parallel:
                scheduler.manager.enforce_sequential_allocation()
                if scheduler.manager.pool.max_concurrent_jobs != 1:
                    scheduler.manager.pool.max_concurrent_jobs = 1
                    scheduler.manager.max_jobs = 1
                    scheduler.manager._sync_parallel_flag()
            # Jobs were already added to scheduler.manager (since manager = scheduler.manager)
            # No need to add them again
            return scheduler.run()

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


def execute_manually_workflows(
    config: Dict[str, Any],
    *,
    allow_parallel: bool,
    scheduler: Optional["GlobalOrcaScheduler"] = None,
    **kwargs,
) -> WorkflowRunResult:
    """Run manual oxidation/reduction steps via the shared workflow scheduler."""

    # If a scheduler is provided, use its manager (which may already have ESD jobs)
    # Otherwise create a new manager
    if scheduler is not None:
        manager = scheduler.manager
    else:
        manager = _WorkflowManager(config, label="manually")

    try:
        _populate_manual_jobs(manager, config, kwargs)
        if not manager.has_jobs():
            logger.info("[manually] No oxidation/reduction jobs queued for execution")
            return WorkflowRunResult()

        jobs_snapshot = list(manager._jobs.values())

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

        if scheduler is not None:
            # Apply sequential allocation to scheduler's manager if needed
            if not allow_parallel:
                scheduler.manager.enforce_sequential_allocation()
                if scheduler.manager.pool.max_concurrent_jobs != 1:
                    scheduler.manager.pool.max_concurrent_jobs = 1
                    scheduler.manager.max_jobs = 1
                    scheduler.manager._sync_parallel_flag()
            # Jobs were already added to scheduler.manager (since manager = scheduler.manager)
            # No need to add them again
            return scheduler.run()

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
            if not run_orca(output_initial, 'initial.out'):
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
                source_input='initial.inp',
            )

        initial_job_id = 'classic_initial'
        _add_job(initial_job_id, 'initial frequency & geometry optimization', run_initial)

    # Note: Excited state calculations (E_00, S1, T1) are now handled by the ESD module
    # Set ESD_modul=yes and states=[S0,S1,T1] in CONTROL.txt to calculate excited states

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
        elif "esd_S0" in manager._jobs:
            # If initial is skipped but ESD is enabled, wait for S0 to provide initial.xyz
            dependencies.add("esd_S0")
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

                # Add IP/EA jobs for ox_step_1 if properties_of_interest is set, calc_prop_of_interest=yes, and method != classic
                if idx == 1:
                    method = str(config.get('method', '')).strip().lower()
                    calc_prop = str(config.get('calc_prop_of_interest', 'no')).strip().lower()
                    properties = config.get('properties_of_interest', '')
                    if calc_prop in ('yes', 'true', '1', 'on') and properties and method != 'classic':
                        xyz_file = ox_sources[idx]  # Use the same xyz file as the main job
                        append_properties_of_interest_jobs(
                            inp_file=ox_inputs[idx],
                            xyz_file=xyz_file,
                            base_charge=charge,
                            base_multiplicity=multiplicity,
                            properties=properties,
                            config=config,
                            solvent=solvents,
                            metals=metals,
                            main_basisset=main_basis,
                            metal_basisset=metal_basis,
                        )

                if not run_orca(ox_inputs[idx], ox_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {ox_outputs[idx]}")
                run_IMAG(
                    ox_outputs[idx],
                    f"ox_step_{idx}",
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    additions,
                    step_name=f"ox_step_{idx}",
                    source_input=ox_inputs[idx],
                )

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
        elif "esd_S0" in manager._jobs:
            # If initial is skipped but ESD is enabled, wait for S0 to provide initial.xyz
            dependencies.add("esd_S0")
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

                # Add IP/EA jobs for red_step_1 if properties_of_interest is set, calc_prop_of_interest=yes, and method != classic
                if idx == 1:
                    method = str(config.get('method', '')).strip().lower()
                    calc_prop = str(config.get('calc_prop_of_interest', 'no')).strip().lower()
                    properties = config.get('properties_of_interest', '')
                    if calc_prop in ('yes', 'true', '1', 'on') and properties and method != 'classic':
                        xyz_file = red_sources[idx]  # Use the same xyz file as the main job
                        append_properties_of_interest_jobs(
                            inp_file=red_inputs[idx],
                            xyz_file=xyz_file,
                            base_charge=charge,
                            base_multiplicity=multiplicity,
                            properties=properties,
                            config=config,
                            solvent=solvents,
                            metals=metals,
                            main_basisset=main_basis,
                            metal_basisset=metal_basis,
                        )

                if not run_orca(red_inputs[idx], red_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {red_outputs[idx]}")
                run_IMAG(
                    red_outputs[idx],
                    f"red_step_{idx}",
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    additions,
                    step_name=f"red_step_{idx}",
                    source_input=red_inputs[idx],
                )

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
            if not run_orca(output_initial, 'initial.out'):
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
                source_input=output_initial,
            )

        initial_job_id = 'manual_initial'
        _add_job(initial_job_id, 'manual initial frequency job', run_initial)

        # Note: Excited state calculations (E_00, S1, T1) are now handled by the ESD module
        # Set ESD_modul=yes and states=[S0,S1,T1] in CONTROL.txt to calculate excited states

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

                # Add IP/EA jobs for ox_step_1 if properties_of_interest is set, calc_prop_of_interest=yes, and method != classic
                if idx == 1:
                    method = str(config.get('method', '')).strip().lower()
                    calc_prop = str(config.get('calc_prop_of_interest', 'no')).strip().lower()
                    properties = config.get('properties_of_interest', '')
                    if calc_prop in ('yes', 'true', '1', 'on') and properties and method != 'classic':
                        xyz_file = ox_sources[idx]  # Use the same xyz file as the main job
                        append_properties_of_interest_jobs(
                            inp_file=ox_inputs[idx],
                            xyz_file=xyz_file,
                            base_charge=charge,
                            base_multiplicity=multiplicity,
                            properties=properties,
                            config=config,
                            solvent=solvents,
                            metals=metals,
                            main_basisset=main_basis,
                            metal_basisset=metal_basis,
                        )

                if not run_orca(ox_inputs[idx], ox_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {ox_outputs[idx]}")
                run_IMAG(
                    ox_outputs[idx],
                    f"ox_step_{idx}",
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    additions,
                    step_name=f"ox_step_{idx}",
                    source_input=ox_inputs[idx],
                )

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

                # Add IP/EA jobs for red_step_1 if properties_of_interest is set, calc_prop_of_interest=yes, and method != classic
                if idx == 1:
                    method = str(config.get('method', '')).strip().lower()
                    calc_prop = str(config.get('calc_prop_of_interest', 'no')).strip().lower()
                    properties = config.get('properties_of_interest', '')
                    if calc_prop in ('yes', 'true', '1', 'on') and properties and method != 'classic':
                        xyz_file = red_sources[idx]  # Use the same xyz file as the main job
                        append_properties_of_interest_jobs(
                            inp_file=red_inputs[idx],
                            xyz_file=xyz_file,
                            base_charge=charge,
                            base_multiplicity=multiplicity,
                            properties=properties,
                            config=config,
                            solvent=solvents,
                            metals=metals,
                            main_basisset=main_basis,
                            metal_basisset=metal_basis,
                        )

                if not run_orca(red_inputs[idx], red_outputs[idx]):
                    raise RuntimeError(f"ORCA terminated abnormally for {red_outputs[idx]}")
                run_IMAG(
                    red_outputs[idx],
                    f"red_step_{idx}",
                    charge,
                    multiplicity,
                    solvents,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    additions,
                    step_name=f"red_step_{idx}",
                    source_input=red_inputs[idx],
                )

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


def _add_moinp_block(input_path: str, gbw_path: str) -> None:
    """Add %moinp block and MOREAD keyword to reuse wavefunction from OCCUPIER GBW file."""
    try:
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as stream:
            lines = stream.readlines()
    except FileNotFoundError as exc:
        raise RuntimeError(f"Input file '{input_path}' missing") from exc

    moinp_line = f'%moinp "{gbw_path}"\n'

    # Check if %moinp already exists
    has_moinp = False
    for line in lines:
        if line.strip().startswith('%moinp'):
            has_moinp = True
            break

    if not has_moinp:
        # Insert %moinp before %maxcore (or before first % block if no maxcore)
        insert_idx = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('%maxcore'):
                insert_idx = idx
                break
            elif stripped.startswith('%') and insert_idx == 0:
                insert_idx = idx

        lines.insert(insert_idx, moinp_line)

    # Add MOREAD keyword if not already present (required when using MOINP with .gbw)
    for idx, line in enumerate(lines):
        if line.strip().startswith('!'):
            if 'MOREAD' not in line:
                lines[idx] = line.rstrip() + ' MOREAD\n'
            break

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
