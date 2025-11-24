"""Dynamic Core Pool Management for optimal resource utilization."""

import threading
import time
import queue
from concurrent.futures import ThreadPoolExecutor, Future
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum
from itertools import count

from delfin.common.logging import get_logger

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
    actual_start_time: Optional[float] = None
    allocated_cores: int = 0
    retry_count: int = 0
    next_retry_time: float = 0.0
    suppress_pool_logs: bool = False

    # Parent-Child Job Tracking for nested jobs
    parent_job_id: Optional[str] = None  # ID of parent job if this is a child
    child_jobs: List[str] = field(default_factory=list)  # IDs of child jobs spawned by this job
    borrowed_cores: int = 0  # Cores borrowed from parent job
    lent_cores: int = 0  # Cores lent to child jobs


class DynamicCorePool:
    """Dynamic core pool with intelligent resource allocation."""

    def __init__(self, total_cores: int, total_memory_mb: int, max_jobs: int = 4):
        self.total_cores = total_cores
        self.total_memory_mb = total_memory_mb
        self.max_concurrent_jobs = max_jobs

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

        # Execution
        self._executor = ThreadPoolExecutor(max_workers=max_jobs)
        self._futures: Dict[str, Future] = {}
        self._shutdown = False
        self._resource_event = threading.Event()

        # Start resource monitor
        self._monitor_thread = threading.Thread(target=self._resource_monitor, daemon=True)
        self._monitor_thread.start()

        # Scheduler thread to react immediately to new work
        self._scheduler_thread = threading.Thread(target=self._scheduler_loop, daemon=True)
        self._scheduler_thread.start()

        logger.info(f"Dynamic Core Pool initialized: {total_cores} cores, {total_memory_mb}MB memory")

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
        job.next_retry_time = time.time()
        priority_value = job.priority.value
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

            # Calculate how many cores the parent can lend
            # Parent must keep at least cores_min for itself
            available_to_lend = parent_job.allocated_cores - parent_job.lent_cores - parent_job.cores_min
            cores_to_lend = min(requested_cores, max(0, available_to_lend))

            if cores_to_lend > 0:
                parent_job.lent_cores += cores_to_lend
                child_job.borrowed_cores = cores_to_lend
                logger.debug(
                    f"Job {child_job.job_id} borrowed {cores_to_lend} cores from parent {parent_job_id} "
                    f"(parent has {parent_job.allocated_cores - parent_job.lent_cores} cores remaining)"
                )

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
            if parent_job_id in self._running_jobs:
                parent_job = self._running_jobs[parent_job_id]
                parent_job.lent_cores -= child_job.borrowed_cores
                logger.debug(
                    f"Job {child_job.job_id} returned {child_job.borrowed_cores} cores to parent {parent_job_id}"
                )

            child_job.borrowed_cores = 0

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
        Otherwise, stick to cores_optimal to leave room for other jobs.
        """
        # Check if other jobs are waiting
        has_pending = not self._job_queue.empty()

        if has_pending:
            # Other jobs waiting - use cores_optimal to leave resources
            return max(
                new_job.cores_min,
                min(new_job.cores_optimal, new_job.cores_max, available_cores)
            )
        else:
            # No other jobs waiting - use all available cores up to cores_max
            return max(
                new_job.cores_min,
                min(new_job.cores_max, available_cores)
            )

    def _min_waiting_cores_locked(self, default: int = 0) -> int:
        """Return minimum cores_min among queued jobs (expects caller holds lock)."""
        if self._job_queue.empty():
            return 0
        try:
            pending = [item[3] for item in list(self._job_queue.queue)]
        except Exception:
            return default
        mins = [job.cores_min for job in pending if hasattr(job, 'cores_min')]
        return min(mins) if mins else default

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
        """Handle job completion and resource release."""
        with self._lock:
            if job_id not in self._running_jobs:
                return

            job = self._running_jobs.pop(job_id)
            self._futures.pop(job_id, None)

            # Return borrowed cores to parent if applicable
            if job.borrowed_cores > 0:
                self.return_borrowed_cores(job)

            # Release pool resources only if job was using pool cores
            if job.borrowed_cores == 0:
                self._allocated_cores -= job.allocated_cores
                self._allocated_memory -= job.memory_mb

            duration = time.time() - job.actual_start_time if job.actual_start_time else 0

            try:
                result = future.result()
                self._completed_jobs.append(job_id)
                if not job.suppress_pool_logs:
                    source = "borrowed" if job.borrowed_cores > 0 else "pool"
                    logger.info(
                        f"Job {job_id} completed in {duration:.1f}s, "
                        f"freed {job.allocated_cores} {source} cores"
                    )
            except Exception as e:
                self._failed_jobs.append(job_id)
                logger.error(f"Job {job_id} failed after {duration:.1f}s: {e}")

        # Try to start waiting jobs now that resources are free
        self._try_rebalance_resources()
        self._signal_state_change()

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
        # NOTE: Core reallocation is disabled for ORCA jobs because ORCA cannot
        # change core count during execution (it's fixed in the .inp file at start).
        # Reallocating cores to running jobs would waste resources instead of
        # allowing new jobs to start.
        # self._consider_reallocation()

    def _consider_reallocation(self):
        """Consider reallocating cores to running jobs for better utilization."""
        with self._lock:
            if not self._running_jobs:
                return

            # Do not reallocate while there are queued jobs waiting for resources.
            if not self._job_queue.empty():
                return

            spare_cores = self.total_cores - self._allocated_cores
            if spare_cores <= 0:
                return

            # Find jobs that could benefit from more cores
            candidates = []
            for job in self._running_jobs.values():
                if job.allocated_cores < job.cores_max:
                    potential_gain = min(job.cores_max - job.allocated_cores, spare_cores)
                    if potential_gain > 0:
                        # Estimate benefit (simple heuristic)
                        benefit_score = potential_gain / (job.allocated_cores + 1)
                        candidates.append((job, potential_gain, benefit_score))

            if not candidates:
                return

            # Sort by benefit score (highest first)
            candidates.sort(key=lambda x: x[2], reverse=True)

            # Allocate spare cores to most beneficial jobs
            for job, potential_gain, _ in candidates:
                if spare_cores <= 0:
                    break

                additional_cores = min(potential_gain, spare_cores)
                job.allocated_cores += additional_cores
                self._allocated_cores += additional_cores
                spare_cores -= additional_cores

                logger.info(f"Reallocated +{additional_cores} cores to job {job.job_id} "
                           f"(now {job.allocated_cores} cores)")

    def _resource_monitor(self):
        """Background thread to monitor and optimize resource usage."""
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

                # Try to optimize resource allocation
                self._try_rebalance_resources()

            except Exception as e:
                logger.error(f"Resource monitor error: {e}")

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
