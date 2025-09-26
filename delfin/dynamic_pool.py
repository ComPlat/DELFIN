"""Dynamic Core Pool Management for optimal resource utilization."""

import threading
import time
import queue
from concurrent.futures import ThreadPoolExecutor, Future
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass
from enum import Enum
from itertools import count

from delfin.common.logging import get_logger

logger = get_logger(__name__)


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
        """Submit a job to the dynamic pool."""
        job.retry_count = 0
        job.next_retry_time = time.time()
        priority_value = job.priority.value
        self._job_queue.put((priority_value, job.next_retry_time, next(self._job_counter), job))
        logger.info(f"Job {job.job_id} queued (cores: {job.cores_min}-{job.cores_optimal}-{job.cores_max})")
        self._signal_state_change()
        return job.job_id

    def _calculate_optimal_allocation(self, job: PoolJob) -> Optional[int]:
        """Calculate optimal core allocation for a job."""
        with self._lock:
            available_cores = self.total_cores - self._allocated_cores
            available_memory = self.total_memory_mb - self._allocated_memory

            # Check if job can run at all
            if (available_cores < job.cores_min or
                available_memory < job.memory_mb or
                len(self._running_jobs) >= self.max_concurrent_jobs):
                return None

            # Calculate allocation based on current load and job requirements
            if len(self._running_jobs) == 0:
                # No other jobs running, give optimal cores
                cores = min(job.cores_optimal, available_cores)
            else:
                # Balance between current jobs and new job
                cores = self._calculate_balanced_allocation(job, available_cores)

            # Ensure we stay within job constraints
            cores = max(job.cores_min, min(cores, job.cores_max, available_cores))

            return cores

    def _calculate_balanced_allocation(self, new_job: PoolJob, available_cores: int) -> int:
        """Calculate balanced allocation considering running jobs."""
        running_count = len(self._running_jobs)

        # Strategy 1: Fair sharing
        fair_share = available_cores // (running_count + 1)

        # Strategy 2: Priority-based allocation
        if new_job.priority == JobPriority.HIGH:
            priority_allocation = min(new_job.cores_optimal, available_cores // 2)
        elif new_job.priority == JobPriority.LOW:
            priority_allocation = max(new_job.cores_min, available_cores // 4)
        else:
            priority_allocation = fair_share

        # Strategy 3: Performance-based allocation
        # Give more cores to jobs that can use them efficiently
        efficiency_factor = min(new_job.cores_optimal / new_job.cores_min, 2.0)
        performance_allocation = int(fair_share * efficiency_factor)

        # Choose the best strategy
        allocation = max(
            new_job.cores_min,
            min(fair_share, priority_allocation, performance_allocation)
        )

        return min(allocation, available_cores)

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
            delay = min(0.2 * (2 ** limited_retry_count), 2.0)
            job.next_retry_time = now + delay
            self._job_queue.put((job.priority.value, job.next_retry_time, next(self._job_counter), job))
            return False

        except queue.Empty:
            return False

    def _start_job(self, job: PoolJob, allocated_cores: int):
        """Start a job with allocated resources."""
        with self._lock:
            # Reserve resources
            self._allocated_cores += allocated_cores
            self._allocated_memory += job.memory_mb

            job.allocated_cores = allocated_cores
            job.actual_start_time = time.time()
            self._running_jobs[job.job_id] = job

            logger.info(f"Starting job {job.job_id} with {allocated_cores} cores "
                       f"({self._allocated_cores}/{self.total_cores} cores used)")

            # Create modified execution function with allocated cores
            def execute_with_cores():
                try:
                    # Add allocated cores to kwargs for ORCA
                    modified_kwargs = job.kwargs.copy()
                    modified_kwargs['cores'] = allocated_cores

                    return job.execute_func(*job.args, **modified_kwargs)
                except Exception as e:
                    logger.error(f"Job {job.job_id} failed: {e}")
                    raise

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

            # Release resources
            self._allocated_cores -= job.allocated_cores
            self._allocated_memory -= job.memory_mb

            duration = time.time() - job.actual_start_time if job.actual_start_time else 0

            try:
                result = future.result()
                self._completed_jobs.append(job_id)
                logger.info(f"Job {job_id} completed in {duration:.1f}s, "
                           f"freed {job.allocated_cores} cores")
            except Exception as e:
                self._failed_jobs.append(job_id)
                logger.error(f"Job {job_id} failed after {duration:.1f}s: {e}")

        # Try to start waiting jobs now that resources are free
        self._try_rebalance_resources()
        self._signal_state_change()

    def _try_rebalance_resources(self):
        """Attempt to rebalance resources and start new jobs."""
        self._schedule_pending_jobs()
        # Consider reallocating cores to running jobs if we have spare capacity
        self._consider_reallocation()

    def _consider_reallocation(self):
        """Consider reallocating cores to running jobs for better utilization."""
        with self._lock:
            if not self._running_jobs:
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

    def shutdown(self):
        """Shutdown the pool gracefully."""
        logger.info("Shutting down dynamic core pool...")
        self._shutdown = True
        self._signal_state_change()
        self._resource_event.set()
        if self._scheduler_thread.is_alive():
            self._scheduler_thread.join(timeout=2)
        self._executor.shutdown(wait=True)
        if self._monitor_thread.is_alive():
            self._monitor_thread.join(timeout=2)

    def _schedule_pending_jobs(self):
        """Start as many queued jobs as resources allow."""
        while True:
            with self._lock:
                if self._job_queue.empty() or len(self._running_jobs) >= self.max_concurrent_jobs:
                    break
            if not self._try_start_next_job():
                break

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


def create_orca_job(job_id: str, inp_file: str, out_file: str,
                   cores_min: int = 2, cores_optimal: int = 8, cores_max: int = 16,
                   priority: JobPriority = JobPriority.NORMAL,
                   estimated_duration: float = 3600.0) -> PoolJob:
    """Create an ORCA job for the dynamic pool."""
    from delfin.orca import run_orca

    def orca_wrapper(inp_file: str, out_file: str, cores: int = 8, timeout: Optional[int] = None):
        """Wrapper to run ORCA with specified cores."""
        # Import here to avoid circular imports
        import tempfile
        import shutil
        from pathlib import Path

        # Modify input file to use allocated cores
        temp_inp = None
        try:
            # Read original input
            with open(inp_file, 'r') as f:
                content = f.read()

            # Replace or add PAL directive
            lines = content.split('\n')
            pal_found = False

            for i, line in enumerate(lines):
                if line.strip().startswith('%pal'):
                    # Replace existing PAL block
                    lines[i] = f'%pal nprocs {cores} end'
                    pal_found = True
                    break

            if not pal_found:
                # Add PAL directive after any existing % blocks
                insert_pos = 0
                for i, line in enumerate(lines):
                    if line.strip().startswith('%') and not line.strip().startswith('%pal'):
                        insert_pos = i + 1
                    elif line.strip() and not line.strip().startswith('%'):
                        break
                lines.insert(insert_pos, f'%pal nprocs {cores} end')

            # Write temporary input file
            temp_inp = inp_file + f'.cores{cores}'
            with open(temp_inp, 'w') as f:
                f.write('\n'.join(lines))

            # Run ORCA with modified input
            run_orca(temp_inp, out_file, timeout=timeout)

            return True

        finally:
            # Cleanup temporary file
            if temp_inp and Path(temp_inp).exists():
                Path(temp_inp).unlink()

    # Estimate memory requirement (rough heuristic)
    memory_mb = cores_optimal * 1000  # 1GB per core baseline

    return PoolJob(
        job_id=job_id,
        cores_min=cores_min,
        cores_optimal=cores_optimal,
        cores_max=cores_max,
        memory_mb=memory_mb,
        priority=priority,
        execute_func=orca_wrapper,
        args=(inp_file, out_file),
        kwargs={'timeout': int(estimated_duration * 1.5)},
        estimated_duration=estimated_duration
    )
