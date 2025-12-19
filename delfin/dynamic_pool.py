"""Dynamic Core Pool Management for optimal resource utilization."""

import threading
import time
import queue
from concurrent.futures import ThreadPoolExecutor, Future
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum
from itertools import count

from delfin.common.logging import get_logger
from delfin.cleanup import _collect_orca_processes, _terminate_process_groups

logger = get_logger(__name__)

# Thread-local storage for tracking current job context
_job_context = threading.local()


def get_current_job_id() -> Optional[str]:
    """Get the ID of the currently executing pool job (if any).

    Returns None if not currently executing within a pool job.
    """
    return getattr(_job_context, 'current_job_id', None)


def _set_current_job_id(job_id: Optional[str]) -> None:
    """Set the current job ID in thread-local context (internal use only)."""
    _job_context.current_job_id = job_id


class JobPriority(Enum):
    HIGH = 1
    NORMAL = 2
    LOW = 3


@dataclass
class PoolJob:
    """Represents a job in the dynamic pool."""
    job_id: str
    cores_min: int          # Minimum cores needed
    cores_optimal: int      # Optimal cores for best performance
    cores_max: int          # Maximum cores that can be used
    memory_mb: int          # Memory requirement
    priority: JobPriority
    execute_func: Callable
    args: tuple
    kwargs: dict
    estimated_duration: float = 3600.0  # Default 1 hour
    max_duration: Optional[float] = None  # Max allowed duration (None = 3× estimated)
    actual_start_time: Optional[float] = None
    allocated_cores: int = 0
    retry_count: int = 0
    next_retry_time: float = 0.0
    suppress_pool_logs: bool = False
    working_dir: Optional[Path] = None  # Working directory for job-specific process tracking

    # Parent-Child Job Tracking for nested jobs
    parent_job_id: Optional[str] = None  # ID of parent job if this is a child
    child_jobs: List[str] = field(default_factory=list)  # IDs of child jobs spawned by this job
    borrowed_cores: int = 0  # Cores borrowed from parent job
    lent_cores: int = 0  # Cores lent to child jobs
    last_progress_time: Optional[float] = None  # Last time job showed progress

    # Starvation Prevention
    queue_time: Optional[float] = None  # When job was queued
    priority_boost: int = 0  # Number of priority levels boosted to prevent starvation


class DynamicCorePool:
    """Dynamic core pool with intelligent resource allocation."""

    def __init__(self, total_cores: int, total_memory_mb: int, max_jobs: int = 4, config: Optional[Dict[str, Any]] = None):
        self.total_cores = total_cores
        self.total_memory_mb = total_memory_mb
        self.max_concurrent_jobs = max_jobs
        self._original_max_jobs = max_jobs  # Store original for adaptive parallelism

        # Configuration
        self.config = config or {}

        # Configurable timeouts (in seconds)
        self.enable_job_timeouts = self._parse_bool(self.config.get('enable_job_timeouts', 'yes'))
        self.job_timeout_default = int(self.config.get('job_timeout_hours', 24)) * 3600
        self.opt_timeout = int(self.config.get('opt_timeout_hours', 12)) * 3600
        self.frequency_timeout = int(self.config.get('frequency_timeout_hours', 36)) * 3600
        self.sp_timeout = int(self.config.get('sp_timeout_hours', 6)) * 3600

        # Feature flags
        self.enable_adaptive_parallelism = self._parse_bool(self.config.get('enable_adaptive_parallelism', 'yes'))
        self.enable_performance_metrics = self._parse_bool(self.config.get('enable_performance_metrics', 'yes'))

        # Core allocation tracking
        self._allocated_cores = 0
        self._allocated_memory = 0
        self._lock = threading.RLock()
        self._condition = threading.Condition(self._lock)

        # Job management
        self._running_jobs: Dict[str, PoolJob] = {}
        self._job_queue = queue.PriorityQueue()
        self._completed_jobs: List[str] = []
        self._failed_jobs: List[str] = []
        self._job_counter = count()

        # Performance metrics tracking
        from collections import defaultdict
        self._job_metrics = defaultdict(list) if self.enable_performance_metrics else None

        # Execution
        self._executor = ThreadPoolExecutor(max_workers=max_jobs)
        self._futures: Dict[str, Future] = {}
        self._shutdown = False
        self._resource_event = threading.Event()

        # Track warnings for stuck jobs: job_id -> {warned_near_shutdown}
        self._stuck_job_warnings: Dict[str, Dict[str, bool]] = {}

        # Track accounting mismatch warnings to avoid spam: job_id -> last_warning_time
        self._accounting_warnings: Dict[str, float] = {}

        # Track termination attempts: job_id -> {soft_kill_time, warned_before_hard_kill}
        self._termination_state: Dict[str, Dict] = {}

        # Start resource monitor
        self._monitor_thread = threading.Thread(target=self._resource_monitor, daemon=True)
        self._monitor_thread.start()

        # Scheduler thread to react immediately to new work
        self._scheduler_thread = threading.Thread(target=self._scheduler_loop, daemon=True)
        self._scheduler_thread.start()

        logger.info(f"Dynamic Core Pool initialized: {total_cores} cores, {total_memory_mb}MB memory")

    @staticmethod
    def _parse_bool(value: Any) -> bool:
        """Parse boolean value from config."""
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() in ('yes', 'true', '1', 'on')
        return bool(value)

    def submit_job(self, job: PoolJob) -> str:
        """Submit a job to the dynamic pool.

        Automatically detects if being called from within another pool job
        and sets up parent-child relationship for core borrowing.
        """
        if self._shutdown:
            logger.warning(f"Pool is shutting down, cannot accept job {job.job_id}")
            raise RuntimeError("Cannot submit job: pool is shutting down")

        # Auto-detect parent job if not explicitly set
        if job.parent_job_id is None:
            current_job_id = get_current_job_id()
            if current_job_id is not None:
                job.parent_job_id = current_job_id
                # Register this job as a child of the parent
                with self._lock:
                    if current_job_id in self._running_jobs:
                        parent_job = self._running_jobs[current_job_id]
                        parent_job.child_jobs.append(job.job_id)
                        logger.debug(
                            f"Job {job.job_id} registered as child of {current_job_id}"
                        )

        job.retry_count = 0
        now = time.time()
        job.next_retry_time = now
        if job.queue_time is None:
            job.queue_time = now  # Track when job was first queued
        # Apply priority boost for starvation prevention
        priority_value = job.priority.value - job.priority_boost
        self._job_queue.put((priority_value, job.next_retry_time, next(self._job_counter), job))
        if not job.suppress_pool_logs:
            parent_info = f", parent={job.parent_job_id}" if job.parent_job_id else ""
            logger.info(
                f"Job {job.job_id} queued (cores: {job.cores_min}-{job.cores_optimal}-{job.cores_max}{parent_info})"
            )
        self._signal_state_change()
        return job.job_id

    def try_borrow_cores_from_parent(self, child_job: PoolJob, requested_cores: int) -> int:
        """Try to borrow cores from parent job for a child job.

        Includes circular borrowing prevention and validation.

        Args:
            child_job: The child job requesting cores
            requested_cores: Number of cores requested

        Returns:
            Number of cores actually borrowed (may be less than requested or 0)
        """
        with self._lock:
            if child_job.parent_job_id is None:
                return 0

            parent_job_id = child_job.parent_job_id
            if parent_job_id not in self._running_jobs:
                logger.debug(
                    f"Cannot borrow cores: parent job {parent_job_id} not found or not running"
                )
                return 0

            parent_job = self._running_jobs[parent_job_id]

            # ROBUSTNESS: Prevent circular borrowing (A borrows from B, B borrows from C, C borrows from A)
            if parent_job.borrowed_cores > 0:
                logger.warning(
                    f"Circular borrowing prevented: Parent {parent_job_id} itself has borrowed cores. "
                    f"Child {child_job.job_id} cannot borrow from a borrower."
                )
                return 0

            # ROBUSTNESS: Validate parent state
            if parent_job.allocated_cores <= 0:
                logger.error(
                    f"Invalid parent state: {parent_job_id} has {parent_job.allocated_cores} allocated cores"
                )
                return 0

            # Calculate how many cores the parent can lend
            # Parent must keep at least cores_min for itself
            available_to_lend = parent_job.allocated_cores - parent_job.lent_cores - parent_job.cores_min
            cores_to_lend = min(requested_cores, max(0, available_to_lend))

            # ROBUSTNESS: Validate lending won't cause negative cores
            if cores_to_lend > 0:
                if parent_job.lent_cores + cores_to_lend > parent_job.allocated_cores - parent_job.cores_min:
                    logger.error(
                        f"Lending validation failed: Would lend too many cores. "
                        f"Parent {parent_job_id} allocated={parent_job.allocated_cores}, "
                        f"lent={parent_job.lent_cores}, trying to lend={cores_to_lend}"
                    )
                    return 0

                parent_job.lent_cores += cores_to_lend
                child_job.borrowed_cores = cores_to_lend
                logger.debug(
                    f"Job {child_job.job_id} borrowed {cores_to_lend} cores from parent {parent_job_id} "
                    f"(parent has {parent_job.allocated_cores - parent_job.lent_cores} cores remaining)"
                )

                # ROBUSTNESS: Validate accounting after lending
                self._validate_resource_accounting()

            return cores_to_lend

    def return_borrowed_cores(self, child_job: PoolJob) -> None:
        """Return borrowed cores from child back to parent.

        Args:
            child_job: The child job returning cores
        """
        with self._lock:
            if child_job.borrowed_cores == 0 or child_job.parent_job_id is None:
                return

            parent_job_id = child_job.parent_job_id
            cores_to_return = child_job.borrowed_cores

            if parent_job_id in self._running_jobs:
                parent_job = self._running_jobs[parent_job_id]

                # ROBUSTNESS: Validate parent has enough lent cores
                if parent_job.lent_cores < cores_to_return:
                    # This can happen transiently if validation thread corrected parent but child not yet synced
                    if parent_job.lent_cores == 0:
                        # Parent already corrected to 0 - child should skip return (no-op)
                        logger.debug(
                            f"Child {child_job.job_id} skipping core return: parent {parent_job_id} already corrected to 0 lent"
                        )
                        child_job.borrowed_cores = 0
                        return
                    logger.warning(
                        f"Accounting mismatch: Child {child_job.job_id} returning {cores_to_return} cores, "
                        f"but parent {parent_job_id} only has {parent_job.lent_cores} lent. Auto-correcting."
                    )
                    cores_to_return = parent_job.lent_cores

                parent_job.lent_cores -= cores_to_return
                logger.debug(
                    f"Job {child_job.job_id} returned {cores_to_return} cores to parent {parent_job_id}"
                )

                # ROBUSTNESS: Validate lent_cores doesn't go negative
                if parent_job.lent_cores < 0:
                    logger.error(f"BUG: Parent {parent_job_id} lent_cores went negative! Forcing to 0.")
                    parent_job.lent_cores = 0
            else:
                # Parent no longer running - child was orphaned
                logger.debug(
                    f"Child {child_job.job_id} returning cores but parent {parent_job_id} no longer running (orphaned)"
                )

            child_job.borrowed_cores = 0

    def _validate_resource_accounting(self):
        """Validate resource accounting consistency (debug/robustness check).

        Checks:
        - allocated_cores matches sum of running job allocations
        - parent lent_cores matches sum of children borrowed_cores
        - no negative values
        - pool capacity not exceeded

        Auto-corrects accounting errors when safe to do so.
        """
        try:
            # Count pool-allocated cores (excluding borrowed cores)
            total_allocated = sum(
                job.allocated_cores
                for job in self._running_jobs.values()
                if job.borrowed_cores == 0
            )

            # CRITICAL: Pool accounting mismatch check
            if total_allocated != self._allocated_cores:
                job_details = [
                    f"{job.job_id}={job.allocated_cores}c"
                    + (f" (borrowed={job.borrowed_cores})" if job.borrowed_cores > 0 else "")
                    + (f" (lent={job.lent_cores})" if job.lent_cores > 0 else "")
                    for job in self._running_jobs.values()
                ]
                logger.error(
                    f"CRITICAL RESOURCE ACCOUNTING MISMATCH:\n"
                    f"  Pool allocated:  {self._allocated_cores} cores\n"
                    f"  Jobs sum:        {total_allocated} cores\n"
                    f"  Leak/Deficit:    {self._allocated_cores - total_allocated} cores\n"
                    f"  Running jobs:    {', '.join(job_details)}\n"
                    f"  Pool capacity:   {self.total_cores} cores"
                )

            # CRITICAL: Pool overflow check
            if self._allocated_cores > self.total_cores:
                logger.error(
                    f"CRITICAL: Pool overflow detected! "
                    f"{self._allocated_cores} cores allocated exceeds capacity of {self.total_cores} cores. "
                    f"This should never happen!"
                )

            # Validate parent-child accounting
            for job_id, job in self._running_jobs.items():
                if job.lent_cores > 0:
                    # Calculate actual borrowed by children
                    actual_borrowed = sum(
                        child.borrowed_cores
                        for child_id in job.child_jobs
                        if child_id in self._running_jobs
                        for child in [self._running_jobs[child_id]]
                    )
                    if actual_borrowed != job.lent_cores:
                        # Only log warning once every 10 seconds to avoid spam
                        now = time.time()
                        last_warned = self._accounting_warnings.get(job_id, 0)
                        if now - last_warned > 10.0:
                            logger.warning(
                                f"Parent-child accounting mismatch: {job_id} shows {job.lent_cores} lent, "
                                f"but children have {actual_borrowed} borrowed. Auto-correcting."
                            )
                            self._accounting_warnings[job_id] = now
                        # ROBUSTNESS: Auto-correct stale lending counters to avoid repeated noise and
                        # overly conservative scheduling when children already returned cores.
                        job.lent_cores = actual_borrowed

                # ROBUSTNESS: Check for negative values
                if job.lent_cores < 0:
                    logger.error(f"BUG: Job {job_id} has negative lent_cores={job.lent_cores}! Forcing to 0.")
                    job.lent_cores = 0
                if job.borrowed_cores < 0:
                    logger.error(f"BUG: Job {job_id} has negative borrowed_cores={job.borrowed_cores}! Forcing to 0.")
                    job.borrowed_cores = 0

        except Exception as e:
            logger.error(f"Resource validation failed (this is a bug): {e}", exc_info=True)

    def _calculate_optimal_allocation(self, job: PoolJob) -> Optional[int]:
        """Calculate optimal core allocation for a job.

        For child jobs, tries to borrow cores from parent if pool resources are insufficient.
        """
        with self._lock:
            available_cores = self.total_cores - self._allocated_cores
            available_memory = self.total_memory_mb - self._allocated_memory

            # If this is a child job and pool doesn't have enough cores,
            # try to borrow from parent
            if job.parent_job_id is not None and available_cores < job.cores_min:
                borrowed = self.try_borrow_cores_from_parent(job, job.cores_optimal)
                if borrowed >= job.cores_min:
                    # Successfully borrowed enough cores from parent
                    # Child job doesn't allocate from pool, uses borrowed cores
                    logger.debug(
                        f"Job {job.job_id} will run with {borrowed} borrowed cores from parent"
                    )
                    return borrowed
                elif borrowed > 0:
                    # Partially borrowed, return them and try normal allocation
                    self.return_borrowed_cores(job)

            # Check if job can run at all from pool resources
            if (available_cores < job.cores_min or
                available_memory < job.memory_mb or
                len(self._running_jobs) >= self.max_concurrent_jobs):
                return None

            # Calculate allocation based on current load and job requirements
            effective_available = max(job.cores_min, available_cores)

            if len(self._running_jobs) == 0:
                # No other jobs running, give optimal cores while keeping reserved capacity if required
                cores = min(job.cores_optimal, effective_available)
            else:
                # Balance between current jobs and new job
                cores = self._calculate_balanced_allocation(job, effective_available)

            # Ensure we stay within job constraints
            cores = max(job.cores_min, min(cores, job.cores_max, effective_available))

            return cores

    def _calculate_balanced_allocation(self, new_job: PoolJob, available_cores: int) -> int:
        """
        Allocate cores intelligently based on available resources and pending jobs.

        If no other jobs are waiting, give more than cores_optimal (up to cores_max).
        If jobs are waiting BUT can't start (max_concurrent_jobs reached), use more cores.
        Otherwise, stick to cores_optimal to leave room for other jobs.
        """
        # Check if other jobs are waiting AND can actually start
        has_pending = not self._job_queue.empty()
        can_start_more = len(self._running_jobs) < self.max_concurrent_jobs

        if has_pending and can_start_more:
            # Jobs waiting AND can start - use cores_optimal to leave resources
            return max(
                new_job.cores_min,
                min(new_job.cores_optimal, new_job.cores_max, available_cores)
            )
        else:
            # No jobs waiting OR max_concurrent_jobs reached - use all available cores
            return max(
                new_job.cores_min,
                min(new_job.cores_max, available_cores)
            )

    def _try_start_next_job(self) -> bool:
        """Try to start the next job from the queue."""
        try:
            # Get next job (blocks briefly)
            _, scheduled_time, _, job = self._job_queue.get(timeout=0.1)

            now = time.time()
            if scheduled_time > now:
                # Not ready yet; put it back and wait for the next cycle
                self._job_queue.put((job.priority.value, scheduled_time, next(self._job_counter), job))
                return False

            cores = self._calculate_optimal_allocation(job)
            if cores is not None:
                job.retry_count = 0
                job.next_retry_time = 0.0
                self._start_job(job, cores)
                return True

            # Unable to start due to resource limits; backoff briefly
            job.retry_count += 1
            # Prevent overflow by limiting retry_count in the exponential calculation
            limited_retry_count = min(job.retry_count - 1, 10)  # Cap at 2^10 = 1024

            # Use shorter delays for the first few retries, then cap at 1s instead of 2s
            if job.retry_count <= 3:
                # Fast retries for first attempts (100ms, 200ms, 400ms)
                delay = min(0.1 * (2 ** limited_retry_count), 0.5)
            else:
                # Slower retries after initial attempts, max 1s
                delay = min(0.2 * (2 ** limited_retry_count), 1.0)

            job.next_retry_time = now + delay
            self._job_queue.put((job.priority.value, job.next_retry_time, next(self._job_counter), job))
            return False

        except queue.Empty:
            return False

    def _start_job(self, job: PoolJob, allocated_cores: int):
        """Start a job with allocated resources (either from pool or borrowed from parent)."""
        with self._lock:
            # Check if pool is shutting down
            if self._shutdown:
                logger.debug(f"Pool shutting down, cannot start job {job.job_id}")
                return

            # Reserve resources only if not using borrowed cores
            if job.borrowed_cores == 0:
                self._allocated_cores += allocated_cores
                self._allocated_memory += job.memory_mb

            job.allocated_cores = allocated_cores
            job.actual_start_time = time.time()
            self._running_jobs[job.job_id] = job

            if not job.suppress_pool_logs:
                source = "borrowed from parent" if job.borrowed_cores > 0 else "from pool"
                logger.info(
                    f"Starting job {job.job_id} with {allocated_cores} cores {source} "
                    f"({self._allocated_cores}/{self.total_cores} pool cores used)"
                )

            # Create modified execution function with allocated cores and job context
            def execute_with_cores():
                # Set job context so child jobs can detect their parent
                _set_current_job_id(job.job_id)
                try:
                    # Add allocated cores and pool usage snapshot to kwargs
                    modified_kwargs = job.kwargs.copy()
                    modified_kwargs['cores'] = allocated_cores
                    modified_kwargs['pool_snapshot'] = (
                        self._allocated_cores,
                        self.total_cores,
                    )

                    return job.execute_func(*job.args, **modified_kwargs)
                finally:
                    # Clear job context when done
                    _set_current_job_id(None)

            # Submit to executor
            future = self._executor.submit(execute_with_cores)
            self._futures[job.job_id] = future

            # Set completion callback
            future.add_done_callback(lambda f: self._job_completed(job.job_id, f))

    def _job_completed(self, job_id: str, future: Future):
        """Handle job completion and resource release.

        Includes orphaned child handling and resource cleanup validation.
        """
        with self._lock:
            if job_id not in self._running_jobs:
                return

            job = self._running_jobs.pop(job_id)
            self._futures.pop(job_id, None)

            # Cleanup warning and termination trackers for this job
            self._stuck_job_warnings.pop(job_id, None)
            self._termination_state.pop(job_id, None)
            self._accounting_warnings.pop(job_id, None)

            # Handle orphaned children if parent is completing
            if job.child_jobs:
                self._handle_orphaned_children(job)

            # Return borrowed cores to parent if applicable
            if job.borrowed_cores > 0:
                self.return_borrowed_cores(job)

            # CRITICAL: Validate that all lent cores have been returned
            if job.lent_cores > 0:
                # Calculate actual borrowed cores by children
                actual_borrowed = sum(
                    child.borrowed_cores
                    for child_id in job.child_jobs
                    if child_id in self._running_jobs
                    for child in [self._running_jobs[child_id]]
                )

                logger.error(
                    f"CRITICAL: Job {job_id} completing with {job.lent_cores} cores still lent "
                    f"(actual borrowed by children: {actual_borrowed}). "
                    f"This indicates a resource accounting leak. Auto-recovering."
                )

                # Auto-recovery: If children still have borrowed cores, they were orphaned
                # This should have been handled by _handle_orphaned_children, but double-check
                if actual_borrowed > 0:
                    logger.error(
                        f"BUG: Orphaned children still have {actual_borrowed} borrowed cores! "
                        f"Force-converting to pool allocation to prevent leak."
                    )
                    for child_id in job.child_jobs:
                        if child_id in self._running_jobs:
                            child = self._running_jobs[child_id]
                            if child.borrowed_cores > 0:
                                self._allocated_cores += child.borrowed_cores
                                logger.error(
                                    f"Auto-recovery: Child {child_id} converted {child.borrowed_cores} "
                                    f"borrowed cores → pool allocation"
                                )
                                child.borrowed_cores = 0
                                child.parent_job_id = None

                # Force cleanup to prevent core leak
                job.lent_cores = 0

            # Release pool resources only if job was using pool cores
            if job.borrowed_cores == 0:
                self._allocated_cores -= job.allocated_cores
                self._allocated_memory -= job.memory_mb

            duration = time.time() - job.actual_start_time if job.actual_start_time else 0

            try:
                result = future.result()
                self._completed_jobs.append(job_id)

                # PERFORMANCE METRICS: Track job completion
                if self._job_metrics is not None:
                    self._record_job_metrics(job_id, duration, job.allocated_cores, success=True)

                if not job.suppress_pool_logs:
                    source = "borrowed" if job.borrowed_cores > 0 else "pool"
                    logger.info(
                        f"Job {job_id} completed in {duration:.1f}s, "
                        f"freed {job.allocated_cores} {source} cores"
                    )
            except Exception as e:
                self._failed_jobs.append(job_id)

                # PERFORMANCE METRICS: Track job failure
                if self._job_metrics is not None:
                    self._record_job_metrics(job_id, duration, job.allocated_cores, success=False)

                logger.error(f"Job {job_id} failed after {duration:.1f}s: {e}")

        # Try to start waiting jobs now that resources are free
        self._try_rebalance_resources()
        self._signal_state_change()

    def _handle_orphaned_children(self, parent_job: PoolJob):
        """Handle children whose parent is completing/dying.

        Converts borrowed cores to pool-allocated cores for orphaned children
        to prevent resource accounting errors.

        This is a critical auto-recovery mechanism that prevents resource leaks
        when parent jobs complete before their children.
        """
        orphaned_count = 0
        total_cores_converted = 0

        for child_id in parent_job.child_jobs:
            if child_id in self._running_jobs:
                child = self._running_jobs[child_id]
                if child.borrowed_cores > 0:
                    # Child is orphaned with borrowed cores
                    # Convert borrowed cores to pool allocation
                    cores_to_convert = child.borrowed_cores

                    # ROBUSTNESS: Validate pool won't overflow
                    new_allocation = self._allocated_cores + cores_to_convert
                    if new_allocation > self.total_cores:
                        logger.error(
                            f"CRITICAL: Converting {cores_to_convert} borrowed cores from orphaned child {child_id} "
                            f"would exceed pool capacity ({new_allocation} > {self.total_cores}). "
                            f"This indicates a serious accounting bug!"
                        )
                        # Continue anyway to prevent worse corruption

                    logger.error(
                        f"AUTO-RECOVERY: Child job {child_id} orphaned by parent {parent_job.job_id}. "
                        f"Converting {cores_to_convert} borrowed cores → pool allocation. "
                        f"Pool: {self._allocated_cores} → {self._allocated_cores + cores_to_convert}/{self.total_cores}"
                    )

                    # Allocate from pool to replace borrowed cores
                    self._allocated_cores += cores_to_convert
                    child.borrowed_cores = 0
                    child.parent_job_id = None  # Orphan the child
                    orphaned_count += 1
                    total_cores_converted += cores_to_convert

        if orphaned_count > 0:
            logger.warning(
                f"Handled {orphaned_count} orphaned child(ren) from parent {parent_job.job_id}. "
                f"Converted {total_cores_converted} cores: borrowed → pool allocation."
            )
            # Validate accounting after orphan handling
            self._validate_resource_accounting()

    def drain_completed_jobs(self) -> List[str]:
        """Return and clear list of recently completed job ids."""
        with self._lock:
            done = list(self._completed_jobs)
            self._completed_jobs.clear()
            return done

    def _try_rebalance_resources(self):
        """Attempt to rebalance resources and start new jobs."""
        if self._shutdown:
            return
        self._schedule_pending_jobs()

    def _resource_monitor(self):
        """Background thread to monitor and optimize resource usage.

        Also detects stuck jobs and handles timeouts.
        """
        last_validation_time = 0
        validation_interval = 60  # Validate accounting every 60 seconds
        last_starvation_check = 0
        starvation_check_interval = 30  # Check for starving jobs every 30 seconds

        while not self._shutdown:
            try:
                # Wake periodically or when state changes
                self._resource_event.wait(timeout=5)
                self._resource_event.clear()

                with self._lock:
                    utilization = (self._allocated_cores / self.total_cores) * 100

                    if len(self._running_jobs) > 0:
                        avg_cores = self._allocated_cores / len(self._running_jobs)
                        logger.debug(f"Pool status: {utilization:.1f}% cores used, "
                                   f"{len(self._running_jobs)} jobs, avg {avg_cores:.1f} cores/job")

                    # ROBUSTNESS: Detect stuck/timeout jobs
                    self._detect_stuck_jobs()

                    # ROBUSTNESS: Periodically validate resource accounting
                    now = time.time()
                    if now - last_validation_time >= validation_interval:
                        self._validate_resource_accounting()
                        last_validation_time = now

                    # ROBUSTNESS: Prevent job starvation
                    if now - last_starvation_check >= starvation_check_interval:
                        self._prevent_starvation()
                        last_starvation_check = now

                    # ADAPTIVE PARALLELISM: Adjust parallelism based on failure rate
                    if self.enable_adaptive_parallelism:
                        self._check_failure_rate_and_adapt()

                # Try to optimize resource allocation
                self._try_rebalance_resources()

            except Exception as e:
                logger.error(f"Resource monitor error: {e}")

    def _get_job_timeout(self, job_id: str) -> Optional[float]:
        """Determine timeout for a job based on its type.

        Args:
            job_id: The job identifier

        Returns:
            Timeout in seconds, or None if timeouts are disabled
        """
        # If timeouts are globally disabled, return None (no timeout)
        if not self.enable_job_timeouts:
            return None

        job_id_lower = job_id.lower()

        # Check job type and return appropriate timeout
        if 'frequency' in job_id_lower or 'freq' in job_id_lower:
            return self.frequency_timeout
        elif 'opt' in job_id_lower or 'optimization' in job_id_lower:
            return self.opt_timeout
        elif 'sp' in job_id_lower or 'single_point' in job_id_lower:
            return self.sp_timeout
        else:
            return self.job_timeout_default

    def _detect_stuck_jobs(self):
        """Detect and handle stuck or timed-out jobs.

        Called by resource monitor (expects lock is held).
        Automatically terminates jobs that are severely stuck.
        Uses configurable timeouts based on job type.

        Warns only twice per job:
        1. When job first exceeds max_duration (3× estimated)
        2. One hour before auto-terminate threshold
        """
        now = time.time()
        jobs_to_kill = []

        for job_id, job in self._running_jobs.items():
            if job.actual_start_time is None:
                continue

            runtime = now - job.actual_start_time

            # Get job-type-specific timeout
            stuck_kill_threshold = self._get_job_timeout(job_id)

            # Skip timeout checks if timeouts are disabled
            if stuck_kill_threshold is None:
                continue

            # Calculate max allowed duration
            max_duration = job.max_duration if job.max_duration is not None else (job.estimated_duration * 3)

            # Initialize warning tracker for this job
            if job_id not in self._stuck_job_warnings:
                self._stuck_job_warnings[job_id] = {'warned_near_shutdown': False}

            warnings = self._stuck_job_warnings[job_id]

            # Warning: One hour before auto-terminate (first time only)
            time_until_shutdown = stuck_kill_threshold - runtime
            if time_until_shutdown <= 3600 and time_until_shutdown > 0 and not warnings['warned_near_shutdown']:
                logger.warning(
                    f"Job {job_id} will be auto-terminated in {time_until_shutdown/60:.0f} minutes "
                    f"(running {runtime/3600:.1f}h / max {stuck_kill_threshold/3600:.1f}h)"
                )
                warnings['warned_near_shutdown'] = True

            # Check if job exceeded timeout → Start staged termination
            if runtime > stuck_kill_threshold:
                # Initialize termination state if not exists
                if job_id not in self._termination_state:
                    self._termination_state[job_id] = {
                        'soft_kill_time': now,
                        'warned_before_hard_kill': False
                    }
                    logger.error(
                        f"Job {job_id} exceeded timeout: running for {runtime/3600:.1f}h, "
                        f"max allowed {stuck_kill_threshold/3600:.1f}h. "
                        f"Sending SIGTERM (soft kill), will hard-kill after 1h grace period."
                    )
                    jobs_to_kill.append((job_id, job, 'soft'))
                else:
                    # Check grace period (1 hour = 3600 seconds)
                    term_state = self._termination_state[job_id]
                    time_since_soft_kill = now - term_state['soft_kill_time']
                    grace_period = 3600  # 1 hour

                    if time_since_soft_kill >= grace_period:
                        # Grace period expired → HARD KILL
                        logger.error(
                            f"Job {job_id} did not terminate after 1h grace period. "
                            f"Sending SIGKILL (hard kill) to force termination."
                        )
                        jobs_to_kill.append((job_id, job, 'hard'))
                    elif time_since_soft_kill >= grace_period - 600 and not term_state['warned_before_hard_kill']:
                        # Warn 10 minutes before hard kill
                        logger.warning(
                            f"Job {job_id} will be force-killed in {(grace_period - time_since_soft_kill)/60:.0f} minutes "
                            f"if it does not terminate gracefully."
                        )
                        term_state['warned_before_hard_kill'] = True

        # Terminate stuck jobs (soft or hard)
        if jobs_to_kill:
            for job_id, job, kill_type in jobs_to_kill:
                self._terminate_stuck_job(job_id, job, force=kill_type == 'hard')

    def _terminate_stuck_job(self, job_id: str, job: PoolJob, force: bool = False):
        """Terminate a stuck job and cleanup its resources.

        Called when job exceeds the stuck kill threshold.

        Args:
            job_id: Job identifier
            job: PoolJob instance
            force: If True, use SIGKILL (hard kill). If False, use SIGTERM (soft kill).

        Soft kill (force=False): Sends SIGTERM, allows graceful shutdown
        Hard kill (force=True): Sends SIGKILL, forceful termination
        """
        # Get the future for this job
        future = self._futures.get(job_id)
        if future is None:
            logger.warning(f"Cannot terminate {job_id}: future not found")
            return

        try:
            # Attempt to cancel the future
            # Note: This may not work if the job is already running
            cancelled = future.cancel()
            if cancelled:
                logger.info(f"Successfully cancelled stuck job {job_id}")
            else:
                # Job is already running and can't be cancelled
                kill_type = "SIGKILL (hard)" if force else "SIGTERM (soft)"
                logger.warning(
                    f"Could not cancel stuck job {job_id} (already running). "
                    f"Attempting {kill_type} on ORCA processes."
                )

        except Exception as e:
            logger.error(f"Error attempting to terminate stuck job {job_id}: {e}")
            return

        # If we could not cancel, terminate ORCA processes
        # Soft kill (force=False): SIGTERM only (graceful, allows 1h grace period)
        # Hard kill (force=True): SIGTERM → wait → SIGKILL (forceful after grace period)
        if not cancelled:
            try:
                roots: List[Path] = []
                if job.working_dir is not None:
                    roots.append(Path(job.working_dir).resolve())
                else:
                    # Fallback to current working directory if no job-specific path is available
                    roots.append(Path.cwd())
                    logger.warning(
                        "Job %s missing working_dir; falling back to current directory for termination scope",
                        job_id,
                    )

                procs = _collect_orca_processes(roots)
                if not procs:
                    logger.debug("No ORCA processes found to terminate for %s", job_id)
                    return

                if force:
                    # Hard kill: Full termination sequence (SIGTERM → SIGKILL)
                    summaries = _terminate_process_groups(procs)
                    logger.info("Force-terminated (SIGKILL) ORCA processes for %s: %s", job_id, summaries)
                else:
                    # Soft kill: Only SIGTERM, no escalation
                    # This allows the job to gracefully shutdown during the 1h grace period
                    import signal
                    import os
                    for proc in procs:
                        try:
                            # ProcessInfo is a NamedTuple, not a Process object
                            # Send SIGTERM directly to the process group
                            os.killpg(proc.pgid, signal.SIGTERM)
                            logger.info(f"Sent SIGTERM to process group {proc.pgid} (PID {proc.pid}) for job {job_id}")
                        except ProcessLookupError:
                            logger.debug(f"Process group {proc.pgid} already exited")
                        except Exception as exc:
                            logger.warning(f"Failed to send SIGTERM to PGID {proc.pgid}: {exc}")

            except Exception as exc:  # noqa: BLE001
                logger.warning("Failed to terminate processes for %s: %s", job_id, exc)

    def _prevent_starvation(self):
        """Boost priority of jobs waiting too long in queue to prevent starvation.

        Called periodically by resource monitor (expects lock is held).
        Jobs waiting > 2 hours get priority boost.
        """
        # Based on real data: typical jobs run 0.9h, max 2.9h
        # Jobs waiting > 2× max observed (5.8h) definitely need a boost
        # But we'll be more aggressive and boost after 2 hours of waiting
        starvation_threshold = 2 * 3600  # 2 hours

        now = time.time()

        # Get all queued jobs
        try:
            pending = [item[3] for item in list(self._job_queue.queue)]
        except Exception as e:
            logger.debug(f"Could not inspect queue for starvation prevention: {e}")
            return

        boosted_count = 0
        for job in pending:
            if job.queue_time is None:
                continue

            wait_time = now - job.queue_time

            # Boost priority if waiting too long and not already boosted
            if wait_time > starvation_threshold and job.priority_boost == 0:
                job.priority_boost = 1  # Boost by 1 priority level
                boosted_count += 1
                logger.warning(
                    f"Job {job.job_id} has been waiting {wait_time/3600:.1f}h "
                    f"(priority={job.priority.name}). "
                    f"Boosting priority to prevent starvation."
                )

        if boosted_count > 0:
            logger.info(
                f"Starvation prevention: Boosted priority for {boosted_count} job(s). "
                f"Re-queuing to apply new priorities."
            )
            # Note: The priority boost will take effect on the next retry
            # when the job is re-queued in _try_start_next_job()

    def _check_failure_rate_and_adapt(self):
        """Adapt parallelism based on failure rate (graceful degradation).

        Called periodically by resource monitor (expects lock is held).
        Reduces parallelism if high failure rate detected (resource exhaustion).
        Gradually increases back when failure rate improves.
        """
        total_jobs = len(self._completed_jobs) + len(self._failed_jobs)

        # Need enough data for statistical significance
        if total_jobs < 10:
            return

        failure_rate = len(self._failed_jobs) / total_jobs

        # Determine action based on failure rate
        if failure_rate > 0.50:  # >50% failures - CRITICAL
            new_max = 1
            reason = f"CRITICAL failure rate ({failure_rate:.0%})"

        elif failure_rate > 0.30:  # >30% failures - HIGH
            new_max = max(1, self.max_concurrent_jobs // 2)
            reason = f"HIGH failure rate ({failure_rate:.0%})"

        elif failure_rate > 0.15:  # >15% failures - MODERATE
            new_max = max(1, int(self.max_concurrent_jobs * 0.75))
            reason = f"MODERATE failure rate ({failure_rate:.0%})"

        else:
            # Success rate OK - try to restore parallelism gradually
            if self.max_concurrent_jobs < self._original_max_jobs:
                new_max = min(
                    self._original_max_jobs,
                    self.max_concurrent_jobs + 1
                )
                reason = f"recovering from previous failures ({failure_rate:.0%})"
            else:
                return  # All good, no action needed

        # Apply changes if different from current
        if new_max != self.max_concurrent_jobs:
            old_max = self.max_concurrent_jobs
            self.max_concurrent_jobs = new_max

            logger.warning(
                f"⚠️  ADAPTIVE PARALLELISM: {reason}\n"
                f"   Adjusting max concurrent jobs: {old_max} → {new_max}\n"
                f"   Stats: {len(self._completed_jobs)} completed, {len(self._failed_jobs)} failed"
            )

    def wait_for_completion(self, timeout: Optional[float] = None) -> bool:
        """Wait for all jobs to complete."""
        deadline = time.time() + timeout if timeout else None

        def _all_done() -> bool:
            return (not self._running_jobs and
                    self._job_queue.empty() and
                    not any(f.running() for f in self._futures.values()))

        with self._condition:
            while not _all_done():
                remaining = None if deadline is None else max(0.0, deadline - time.time())
                if deadline is not None and remaining <= 0:
                    logger.warning("Timeout waiting for job completion")
                    return False
                self._condition.wait(timeout=remaining)

        # Auto-print performance metrics if enabled
        if self.enable_performance_metrics and self._job_metrics:
            self.print_performance_metrics()

        return len(self._failed_jobs) == 0

    def get_status(self) -> Dict[str, Any]:
        """Get current pool status."""
        with self._lock:
            return {
                'total_cores': self.total_cores,
                'allocated_cores': self._allocated_cores,
                'utilization_percent': (self._allocated_cores / self.total_cores) * 100,
                'running_jobs': len(self._running_jobs),
                'queued_jobs': self._job_queue.qsize(),
                'completed_jobs': len(self._completed_jobs),
                'failed_jobs': len(self._failed_jobs),
                'job_details': {
                    job_id: {
                        'cores': job.allocated_cores,
                        'priority': job.priority.name,
                        'duration': time.time() - job.actual_start_time if job.actual_start_time else 0
                    }
                    for job_id, job in self._running_jobs.items()
                }
            }

    def shutdown(self, wait: bool = True, cancel_pending: bool = False):
        """Shutdown the pool gracefully.

        Args:
            wait: If True, wait for running jobs to complete
            cancel_pending: If True, cancel all pending jobs (not running)
        """
        logger.info("Shutting down dynamic core pool...")
        self._shutdown = True

        if cancel_pending:
            # Clear pending queue
            while not self._job_queue.empty():
                try:
                    self._job_queue.get_nowait()
                except queue.Empty:
                    break
            logger.info("Cancelled all pending jobs")

        self._signal_state_change()
        self._resource_event.set()
        if self._scheduler_thread.is_alive():
            self._scheduler_thread.join(timeout=2)
        self._executor.shutdown(wait=wait, cancel_futures=not wait)
        if self._monitor_thread.is_alive():
            self._monitor_thread.join(timeout=2)

        # Print performance metrics before final shutdown
        if self.enable_performance_metrics and self._job_metrics:
            self.print_performance_metrics()

    def _classify_job_type(self, job_id: str) -> str:
        """Classify job by type for metrics tracking.

        Args:
            job_id: Job identifier

        Returns:
            Job type string (e.g., 'initial', 'ox_step_1', 'frequency')
        """
        job_id_lower = job_id.lower()

        if 'initial' in job_id_lower:
            return 'initial'
        elif 'frequency' in job_id_lower or 'freq' in job_id_lower:
            return 'frequency'
        elif 'ox' in job_id_lower:
            # Extract step number if present
            import re
            match = re.search(r'ox.*?(\d+)', job_id_lower)
            if match:
                return f'ox_step_{match.group(1)}'
            return 'oxidation'
        elif 'red' in job_id_lower:
            # Extract step number if present
            import re
            match = re.search(r'red.*?(\d+)', job_id_lower)
            if match:
                return f'red_step_{match.group(1)}'
            return 'reduction'
        elif 'occupier' in job_id_lower:
            return 'occupier'
        elif 'sp' in job_id_lower:
            return 'single_point'
        elif 'opt' in job_id_lower:
            return 'optimization'
        else:
            return 'other'

    def _record_job_metrics(self, job_id: str, duration: float, cores: int, success: bool):
        """Record performance metrics for a completed/failed job.

        Args:
            job_id: Job identifier
            duration: Job runtime in seconds
            cores: Number of cores used
            success: Whether job completed successfully
        """
        job_type = self._classify_job_type(job_id)
        self._job_metrics[job_type].append({
            'job_id': job_id,
            'duration': duration,
            'cores': cores,
            'success': success,
        })

    def print_performance_metrics(self, output_file: Optional[str] = None):
        """Print performance metrics report to logger and optionally to file.

        Args:
            output_file: Optional path to write metrics to. If None, uses
                         'performance_metrics_<timestamp>.txt' in current directory.
        """
        if self._job_metrics is None or not self._job_metrics:
            return

        import statistics
        from datetime import datetime

        # Build metrics report as string
        lines = []
        lines.append("="*70)
        lines.append("📊 JOB PERFORMANCE METRICS")
        lines.append("="*70)
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")

        total_jobs = 0
        total_successes = 0
        total_failures = 0
        total_runtime = 0.0

        for job_type in sorted(self._job_metrics.keys()):
            metrics = self._job_metrics[job_type]
            if not metrics:
                continue

            runtimes = [m['duration'] for m in metrics]
            cores_used = [m['cores'] for m in metrics]
            successes = sum(1 for m in metrics if m['success'])
            failures = len(metrics) - successes

            avg_time = statistics.mean(runtimes)
            std_time = statistics.stdev(runtimes) if len(runtimes) > 1 else 0
            min_time = min(runtimes)
            max_time = max(runtimes)
            avg_cores = statistics.mean(cores_used)
            total_time = sum(runtimes)

            # Accumulate totals
            total_jobs += len(metrics)
            total_successes += successes
            total_failures += failures
            total_runtime += total_time

            lines.append(f"{job_type}:")
            lines.append(f"  Jobs:         {len(metrics)} ({successes} ✓, {failures} ✗)")
            lines.append(f"  Avg time:     {avg_time/60:.1f} min ±{std_time/60:.1f} min")
            lines.append(f"  Range:        {min_time/60:.1f} - {max_time/60:.1f} min")
            lines.append(f"  Avg cores:    {avg_cores:.1f}")
            lines.append(f"  Total time:   {total_time/60:.1f} min")
            lines.append("")

        # Add summary
        lines.append("="*70)
        lines.append("SUMMARY:")
        lines.append(f"  Total jobs:       {total_jobs} ({total_successes} ✓, {total_failures} ✗)")
        lines.append(f"  Success rate:     {(total_successes/total_jobs*100):.1f}%" if total_jobs > 0 else "  Success rate:     N/A")
        lines.append(f"  Total runtime:    {total_runtime/60:.1f} min ({total_runtime/3600:.1f} hours)")
        lines.append("="*70)

        # Write to file (no terminal output)
        if output_file is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_file = f"performance_metrics_{timestamp}.txt"

        try:
            from pathlib import Path
            output_path = Path(output_file)
            output_path.write_text('\n'.join(lines), encoding='utf-8')
            logger.info(f"📁 Performance metrics written to: {output_path.absolute()}")
        except Exception as e:
            logger.error(f"Failed to write metrics to file {output_file}: {e}")

    def _schedule_pending_jobs(self):
        """Start as many queued jobs as resources allow."""
        if self._shutdown:
            return
        # Keep trying to pack jobs until we can't start any more
        while True:
            started_any = False
            # Inner loop: try to start jobs until no more will fit
            while True:
                with self._lock:
                    if self._job_queue.empty() or len(self._running_jobs) >= self.max_concurrent_jobs:
                        break
                if self._try_start_next_job():
                    started_any = True
                else:
                    # Failed to start - either resources exhausted or job needs retry
                    break

            # If we started nothing in this round, we're done
            if not started_any:
                break

            # We started something - try another round to pack more jobs
            # This handles cases where starting a small job frees up space for others

    def _scheduler_loop(self):
        """Reactively schedule jobs when resources or queue state changes."""
        while not self._shutdown:
            with self._condition:
                self._condition.wait_for(lambda: self._shutdown or self._should_attempt_schedule(), timeout=1.0)
                if self._shutdown:
                    return
            self._schedule_pending_jobs()

    def _should_attempt_schedule(self) -> bool:
        return (not self._job_queue.empty() and
                len(self._running_jobs) < self.max_concurrent_jobs and
                self._allocated_cores < self.total_cores)

    def _signal_state_change(self):
        with self._condition:
            self._condition.notify_all()
        self._resource_event.set()
