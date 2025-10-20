"""Integration of dynamic pool with OCCUPIER workflow."""

import os
import re
import time
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set

from delfin.common.logging import get_logger
from delfin.dynamic_pool import JobPriority, create_orca_job
from delfin.global_manager import get_global_manager
from delfin.copy_helpers import read_occupier_file
from delfin.imag import run_IMAG
from delfin.orca import run_orca
from delfin.xyz_io import (
    create_s1_optimization_input,
    read_xyz_and_create_input2,
    read_xyz_and_create_input3,
)
from .parallel_classic_manually import (
    WorkflowJob,
    _WorkflowManager,
    _parse_int,
    _update_pal_block,
    _verify_orca_output,
    estimate_parallel_width,
    determine_effective_slots,
    normalize_parallel_token,
)

logger = get_logger(__name__)


class ParallelOccupierManager:
    """Manages parallel execution of OCCUPIER with dynamic resource allocation."""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        global_mgr = get_global_manager()
        if not global_mgr.is_initialized():
            raise RuntimeError(
                "ParallelOccupierManager requires the global job manager to be initialized."
            )

        self.pool = global_mgr.get_pool()
        self.total_cores = max(1, global_mgr.total_cores)
        configured_max_jobs = max(1, _resolve_pal_jobs(config))
        self.max_jobs = max(1, min(configured_max_jobs, self.pool.max_concurrent_jobs))
        self.per_job_share = max(
            1 if self.total_cores == 1 else 2,
            max(1, self.total_cores // max(1, self.max_jobs)),
        )

        # Track if we're running unified (ox+red parallel) to adjust PAL allocation
        self.unified_mode = False
        self.unified_workflow_count = 1

        logger.info(
            "Parallel OCCUPIER manager using GLOBAL SHARED pool (pool_id=%d, %d cores)",
            id(self.pool),
            self.total_cores,
        )

    def execute_parallel_workflows(self, ox_sequence: List[Dict], red_sequence: List[Dict]) -> bool:
        """Execute ox and red workflows in parallel with dynamic resource management."""

        logger.info("Starting parallel OCCUPIER execution: ox_steps + red_steps")

        # Merge both sequences into a unified job list with proper identification
        unified_sequence = []

        if ox_sequence:
            for entry in ox_sequence:
                unified_entry = entry.copy()
                unified_entry['_workflow'] = 'ox'
                unified_entry['_original_index'] = entry['index']
                unified_sequence.append(unified_entry)

        if red_sequence:
            for entry in red_sequence:
                unified_entry = entry.copy()
                unified_entry['_workflow'] = 'red'
                unified_entry['_original_index'] = entry['index']
                unified_sequence.append(unified_entry)

        if not unified_sequence:
            logger.warning("No ox or red sequences to execute")
            return True

        logger.info(
            f"Unified execution: {len(ox_sequence or [])} ox jobs + "
            f"{len(red_sequence or [])} red jobs = {len(unified_sequence)} total jobs"
        )

        # Enable unified mode to split cores between workflows
        self.unified_mode = True
        workflow_count = sum(1 for seq in [ox_sequence, red_sequence] if seq)
        self.unified_workflow_count = workflow_count

        logger.info(
            f"Unified mode enabled: {workflow_count} workflows sharing {self.total_cores} cores"
        )

        # Execute all jobs in a single unified workflow with shared resource pool
        try:
            return self._execute_unified_sequence(unified_sequence)
        finally:
            self.unified_mode = False
            self.unified_workflow_count = 1

    def _execute_unified_sequence(self, unified_sequence: List[Dict]) -> bool:
        """Execute unified ox+red sequence with proper resource sharing."""

        logger.info(f"Starting unified execution with {len(unified_sequence)} jobs")

        # Analyze cross-workflow dependencies
        dependencies = self._analyze_unified_dependencies(unified_sequence)

        # Submit jobs to pool based on dependencies
        submitted_jobs = {}
        completed_jobs = set()
        job_submit_times = {}

        start_time = time.time()
        max_wait_cycles = 1800  # 1 hour timeout
        wait_cycles = 0
        last_progress_time = start_time

        while len(completed_jobs) < len(unified_sequence):
            # Find jobs ready to run
            ready_jobs = []
            for entry in unified_sequence:
                job_key = (entry['_workflow'], entry['_original_index'])
                if (job_key not in submitted_jobs and
                    job_key not in completed_jobs and
                    dependencies[job_key].issubset(completed_jobs)):
                    ready_jobs.append(entry)

            # Submit ready jobs to pool
            for entry in ready_jobs:
                workflow = entry['_workflow']
                idx = entry['_original_index']
                job_key = (workflow, idx)
                job_id = f"{workflow}_steps_job_{idx}"

                inp_file = self._get_workflow_input_filename(workflow, idx)
                out_file = self._get_workflow_output_filename(workflow, idx)

                # Estimate job complexity for resource allocation
                cores_min, cores_opt, cores_max = self._estimate_job_requirements(entry)

                # Create and submit job
                pool_job = create_orca_job(
                    job_id=job_id,
                    inp_file=inp_file,
                    out_file=out_file,
                    cores_min=cores_min,
                    cores_optimal=cores_opt,
                    cores_max=cores_max,
                    priority=JobPriority.NORMAL,
                    estimated_duration=self._estimate_duration(entry)
                )

                self.pool.submit_job(pool_job)
                submitted_jobs[job_key] = job_id
                job_submit_times[job_key] = time.time()

                logger.info(f"Submitted {job_id} to unified pool")

            # Brief pause only on first submission to let jobs start
            if ready_jobs and len(submitted_jobs) == len(ready_jobs):
                time.sleep(0.1)

            # Check for completed jobs
            newly_completed = self._check_completed_unified_jobs(submitted_jobs, completed_jobs)
            completed_jobs.update(newly_completed)

            # Track progress
            if newly_completed or ready_jobs:
                last_progress_time = time.time()
                wait_cycles = 0
            else:
                wait_cycles += 1

            # Deadlock detection
            if not ready_jobs and not newly_completed:
                elapsed_since_progress = time.time() - last_progress_time

                if len(submitted_jobs) > len(completed_jobs):
                    # We have running jobs, just wait
                    if wait_cycles % 30 == 0:
                        running = len(submitted_jobs) - len(completed_jobs)
                        logger.info(
                            f"[unified] Waiting for {running} running jobs to complete "
                            f"({elapsed_since_progress:.0f}s since last progress)"
                        )
                    time.sleep(2)
                elif len(submitted_jobs) < len(unified_sequence):
                    # Deadlock detected
                    pending_keys = set((e['_workflow'], e['_original_index']) for e in unified_sequence) - set(submitted_jobs.keys())
                    logger.error(
                        f"[unified] Deadlock detected: {len(pending_keys)} jobs waiting. "
                        f"Pending: {sorted(pending_keys)}"
                    )
                    for job_key in sorted(pending_keys):
                        deps = dependencies.get(job_key, set())
                        missing_deps = deps - completed_jobs
                        logger.error(f"  Job {job_key} waiting for: {sorted(missing_deps)}")
                    return False
                else:
                    logger.warning("[unified] Unexpected state: breaking loop")
                    break

                if wait_cycles >= max_wait_cycles:
                    logger.error(
                        f"[unified] Timeout after {wait_cycles * 2}s. "
                        f"Completed: {len(completed_jobs)}/{len(unified_sequence)}"
                    )
                    return False

        duration = time.time() - start_time
        logger.info(f"Unified execution completed in {duration:.1f}s")

        return True

    def _analyze_unified_dependencies(self, unified_sequence: List[Dict]) -> Dict[tuple, Set[tuple]]:
        """Analyze dependencies across ox and red workflows."""
        dependencies = {}

        for entry in unified_sequence:
            workflow = entry['_workflow']
            idx = entry['_original_index']
            job_key = (workflow, idx)

            raw_from = entry.get("from", idx - 1)
            dep_indices = self._parse_dependency_field(raw_from)
            dep_indices.discard(idx)

            # Dependencies are within the same workflow
            dep_keys = set((workflow, dep_idx) for dep_idx in dep_indices)
            dependencies[job_key] = dep_keys

        return dependencies

    def _check_completed_unified_jobs(self, submitted_jobs: Dict[tuple, str],
                                     completed_jobs: Set[tuple]) -> Set[tuple]:
        """Check for newly completed jobs in unified execution."""
        newly_completed = set()

        # Use pool feedback first
        finished_ids = set(self.pool.drain_completed_jobs())
        inverse_map = {job_id: key for key, job_id in submitted_jobs.items()}

        for job_id in finished_ids:
            key = inverse_map.get(job_id)
            if key is not None and key not in completed_jobs:
                newly_completed.add(key)
                logger.info(f"Job {job_id} completed (pool notification)")

        # Fallback to file-based detection
        for key, job_id in submitted_jobs.items():
            if key in completed_jobs or key in newly_completed:
                continue
            workflow, idx = key
            out_file = self._get_workflow_output_filename(workflow, idx)
            if self._verify_job_completion(out_file):
                newly_completed.add(key)
                logger.info(f"Job {job_id} completed (output check)")

        return newly_completed

    def _get_workflow_input_filename(self, workflow: str, idx: int) -> str:
        """Get input filename for workflow job."""
        return f"input{'' if idx == 1 else idx}.inp"

    def _get_workflow_output_filename(self, workflow: str, idx: int) -> str:
        """Get output filename for workflow job."""
        return f"output{'' if idx == 1 else idx}.out"

    def _execute_sequence_with_pool(self, workflow_name: str, sequence: List[Dict],
                                   priority: JobPriority) -> bool:
        """Execute a single OCCUPIER sequence using the dynamic pool."""

        logger.info(f"Starting {workflow_name} with {len(sequence)} jobs")

        # Analyze dependencies
        dependencies = self._analyze_dependencies(sequence)

        # Submit jobs to pool based on dependencies
        submitted_jobs = {}
        completed_jobs = set()
        job_submit_times = {}

        start_time = time.time()
        max_wait_cycles = 1800  # 1 hour timeout (1800 * 2s)
        wait_cycles = 0
        last_progress_time = start_time

        while len(completed_jobs) < len(sequence):
            # Find jobs ready to run
            ready_jobs = []
            for entry in sequence:
                idx = entry["index"]
                if (idx not in submitted_jobs and
                    idx not in completed_jobs and
                    dependencies[idx].issubset(completed_jobs)):
                    ready_jobs.append(entry)

            # Submit ready jobs to pool
            for entry in ready_jobs:
                job_id = f"{workflow_name}_job_{entry['index']}"
                inp_file = self._get_input_filename(entry['index'])
                out_file = self._get_output_filename(entry['index'])

                # Estimate job complexity for resource allocation
                cores_min, cores_opt, cores_max = self._estimate_job_requirements(entry)

                # Create and submit job
                pool_job = create_orca_job(
                    job_id=job_id,
                    inp_file=inp_file,
                    out_file=out_file,
                    cores_min=cores_min,
                    cores_optimal=cores_opt,
                    cores_max=cores_max,
                    priority=priority,
                    estimated_duration=self._estimate_duration(entry)
                )

                self.pool.submit_job(pool_job)
                submitted_jobs[entry['index']] = job_id
                job_submit_times[entry['index']] = time.time()

                logger.info(f"Submitted {job_id} to pool")

            # Brief pause only on first submission to let jobs start
            if ready_jobs and len(submitted_jobs) == len(ready_jobs):
                time.sleep(0.1)

            # Check for completed jobs
            newly_completed = self._check_completed_jobs(submitted_jobs, completed_jobs)
            completed_jobs.update(newly_completed)

            # Track progress
            if newly_completed or ready_jobs:
                last_progress_time = time.time()
                wait_cycles = 0
            else:
                wait_cycles += 1

            # Deadlock detection
            if not ready_jobs and not newly_completed:
                elapsed_since_progress = time.time() - last_progress_time

                # Check if we're truly stuck
                if len(submitted_jobs) > len(completed_jobs):
                    # We have running jobs, just wait
                    if wait_cycles % 30 == 0:  # Log every minute
                        running = len(submitted_jobs) - len(completed_jobs)
                        logger.info(
                            f"[{workflow_name}] Waiting for {running} running jobs to complete "
                            f"({elapsed_since_progress:.0f}s since last progress)"
                        )
                    time.sleep(2)
                elif len(submitted_jobs) < len(sequence):
                    # We have unsubmitted jobs but can't submit them (dependency deadlock)
                    pending_indices = set(entry["index"] for entry in sequence) - set(submitted_jobs.keys())
                    logger.error(
                        f"[{workflow_name}] Deadlock detected: {len(pending_indices)} jobs waiting "
                        f"but dependencies not satisfied. Pending: {sorted(pending_indices)}"
                    )
                    # Log dependency info for debugging
                    for idx in sorted(pending_indices):
                        deps = dependencies.get(idx, set())
                        missing_deps = deps - completed_jobs
                        logger.error(
                            f"  Job {idx} waiting for dependencies: {sorted(missing_deps)}"
                        )
                    return False
                else:
                    # All jobs submitted and completed - this shouldn't happen but break anyway
                    logger.warning(f"[{workflow_name}] Unexpected state: breaking loop")
                    break

                # Global timeout
                if wait_cycles >= max_wait_cycles:
                    logger.error(
                        f"[{workflow_name}] Timeout after {wait_cycles * 2}s with no progress. "
                        f"Completed: {len(completed_jobs)}/{len(sequence)}"
                    )
                    return False

        duration = time.time() - start_time
        logger.info(f"{workflow_name} completed in {duration:.1f}s")

        return True

    def _analyze_dependencies(self, sequence: List[Dict]) -> Dict[int, Set[int]]:
        """Analyze dependencies in OCCUPIER sequence."""
        dependencies = {}

        for entry in sequence:
            idx = entry["index"]
            raw_from = entry.get("from", idx - 1)

            dep_indices = self._parse_dependency_field(raw_from)
            dep_indices.discard(idx)
            dependencies[idx] = dep_indices

        return dependencies

    @staticmethod
    def _parse_dependency_field(raw_from: Any) -> Set[int]:
        deps: Set[int] = set()
        if raw_from in (None, "", 0):
            return deps

        def add_value(value: Any) -> None:
            try:
                candidate = int(str(value).strip())
            except (TypeError, ValueError):
                return
            if candidate > 0:
                deps.add(candidate)

        if isinstance(raw_from, (list, tuple, set)):
            for item in raw_from:
                add_value(item)
            return deps

        text = str(raw_from)
        for token in text.replace(";", ",").replace("|", ",").split(","):
            if not token.strip():
                continue
            add_value(token)

        if not deps:
            add_value(raw_from)
        return deps

    def _estimate_job_requirements(self, entry: Dict) -> tuple[int, int, int]:
        """Estimate core requirements for a job."""
        # Adjust total cores if in unified mode (ox+red parallel)
        effective_total_cores = self.total_cores
        if self.unified_mode and self.unified_workflow_count > 1:
            # Split cores between parallel workflows
            effective_total_cores = max(2, self.total_cores // self.unified_workflow_count)
            logger.debug(
                f"Unified mode: reducing effective cores from {self.total_cores} "
                f"to {effective_total_cores} ({self.unified_workflow_count} workflows)"
            )

        # Base requirements
        cores_min = 1 if effective_total_cores == 1 else 2
        burst_capacity = max(cores_min, min(effective_total_cores, self.per_job_share * 2))
        cores_max = burst_capacity

        # Adjust based on multiplicity and job characteristics
        multiplicity = entry.get("m", 1)

        if multiplicity > 3:
            # High-spin jobs might benefit from more cores
            cores_opt = burst_capacity
        else:
            # Standard DFT jobs
            cores_opt = burst_capacity

        # Ensure constraints
        cores_opt = max(cores_min, min(cores_opt, cores_max))

        return cores_min, cores_opt, cores_max

    def _estimate_duration(self, entry: Dict) -> float:
        """Estimate job duration in seconds."""
        base_time = 1800  # 30 minutes base

        # Adjust based on multiplicity
        multiplicity = entry.get("m", 1)
        if multiplicity > 3:
            base_time *= 1.5  # High-spin calculations take longer

        # Adjust based on broken symmetry
        if entry.get("BS"):
            base_time *= 1.3  # BS calculations are more expensive

        return base_time

    def _get_input_filename(self, idx: int) -> str:
        """Get input filename for job index."""
        return f"input{'' if idx == 1 else idx}.inp"

    def _get_output_filename(self, idx: int) -> str:
        """Get output filename for job index."""
        return f"output{'' if idx == 1 else idx}.out"

    def _check_completed_jobs(self, submitted_jobs: Dict[int, str],
                            completed_jobs: Set[int]) -> Set[int]:
        """Check for newly completed jobs."""
        newly_completed = set()

        # Use pool feedback first
        finished_ids = set(self.pool.drain_completed_jobs())
        inverse_map = {job_id: idx for idx, job_id in submitted_jobs.items()}

        for job_id in finished_ids:
            idx = inverse_map.get(job_id)
            if idx is not None and idx not in completed_jobs:
                newly_completed.add(idx)
                logger.info(f"Job {job_id} completed (pool notification)")

        # Fallback to file-based detection for any remaining jobs
        for idx, job_id in submitted_jobs.items():
            if idx in completed_jobs or idx in newly_completed:
                continue
            out_file = self._get_output_filename(idx)
            if self._verify_job_completion(out_file):
                newly_completed.add(idx)
                logger.info(f"Job {job_id} completed (output check)")

        return newly_completed

    def _verify_job_completion(self, out_file: str) -> bool:
        """Verify that an ORCA job completed successfully."""
        try:
            with open(out_file, 'r', errors='ignore') as f:
                content = f.read()
                return "ORCA TERMINATED NORMALLY" in content
        except Exception:
            return False

    def execute_single_sequence(self, sequence: List[Dict], workflow_name: str = "occupier") -> bool:
        """Execute a single OCCUPIER sequence with intelligent parallelization."""

        if len(sequence) <= 1:
            # Single job - run sequentially
            logger.info(f"Running {workflow_name} sequentially (single job)")
            return self._execute_sequential(sequence)

        # Check if parallelization makes sense
        dependencies = self._analyze_dependencies(sequence)
        independent_jobs = sum(1 for deps in dependencies.values() if len(deps) <= 1)  # from=0 or from=1 count as independent
        parallel_potential = sum(1 for level_jobs in self._get_parallel_levels(dependencies).values() if len(level_jobs) >= 2)

        if (independent_jobs >= 2 or parallel_potential >= 1) and self.total_cores >= 4:
            logger.info(f"Running {workflow_name} with dynamic pool parallelization "
                       f"({independent_jobs} independent jobs, {parallel_potential} parallel levels)")
            return self._execute_sequence_with_pool(workflow_name, sequence, JobPriority.NORMAL)
        else:
            logger.info(f"Running {workflow_name} sequentially (insufficient parallelism: "
                       f"{independent_jobs} independent, {parallel_potential} parallel levels, {self.total_cores} cores)")
            return self._execute_sequential(sequence)

    def _get_parallel_levels(self, dependencies: Dict[int, Set[int]]) -> Dict[int, List[int]]:
        """Analyze how many jobs can run in parallel at each dependency level."""
        levels = {}

        def get_level(job_idx: int) -> int:
            if not dependencies[job_idx]:
                return 0
            return max(get_level(dep) for dep in dependencies[job_idx]) + 1

        for job_idx in dependencies:
            level = get_level(job_idx)
            if level not in levels:
                levels[level] = []
            levels[level].append(job_idx)

        return levels

    def _execute_sequential(self, sequence: List[Dict]) -> bool:
        """Fallback sequential execution."""
        from delfin.orca import run_orca

        for entry in sequence:
            idx = entry["index"]
            inp_file = self._get_input_filename(idx)
            out_file = self._get_output_filename(idx)

            logger.info(f"Running sequential job {idx}")
            run_orca(inp_file, out_file)

            if not self._verify_job_completion(out_file):
                logger.error(f"Sequential job {idx} failed")
                return False

        return True

    def get_pool_status(self) -> Dict[str, Any]:
        """Get current status of the dynamic pool."""
        return self.pool.get_status()

    def shutdown(self):
        """Shutdown the parallel manager."""
        logger.debug("Parallel OCCUPIER manager relies on global pool; nothing to shutdown")


@dataclass
class OccupierExecutionContext:
    """Container for OCCUPIER ORCA execution parameters."""

    charge: int
    solvent: str
    metals: List[str]
    main_basisset: str
    metal_basisset: str
    config: Dict[str, Any]
    completed_jobs: Set[str] = field(default_factory=set)
    failed_jobs: Dict[str, str] = field(default_factory=dict)
    skipped_jobs: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class JobDescriptor:
    """Declarative description of an OCCUPIER post-processing job."""

    job_id: str
    description: str
    work: Callable[[int], None]
    produces: Set[str] = field(default_factory=set)
    requires: Set[str] = field(default_factory=set)
    explicit_dependencies: Set[str] = field(default_factory=set)
    preferred_cores: Optional[int] = None


def run_occupier_orca_jobs(context: OccupierExecutionContext, parallel_enabled: bool) -> bool:
    """Execute OCCUPIER post-processing ORCA jobs with optional parallelization."""

    frequency_mode = str(context.config.get('frequency_calculation_OCCUPIER', 'no')).lower()
    if frequency_mode == 'yes':
        logger.info("frequency_calculation_OCCUPIER=yes → skipping ORCA job scheduling")
        return True

    try:
        jobs = _build_occupier_jobs(context)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to prepare OCCUPIER ORCA jobs: %s", exc, exc_info=True)
        return False

    if not jobs:
        logger.info("No OCCUPIER ORCA jobs detected for execution")
        return True

    context.completed_jobs.clear()
    context.failed_jobs.clear()
    context.skipped_jobs.clear()

    pal_jobs_value = _resolve_pal_jobs(context.config)
    parallel_mode = normalize_parallel_token(context.config.get('parallel_workflows', 'auto'))
    width = estimate_parallel_width(jobs)
    requested_parallel = (
        parallel_mode == 'enable'
        or (parallel_mode == 'auto' and width > 1)
    )
    effective_max_jobs = max(1, min(pal_jobs_value, width)) if requested_parallel else 1
    use_parallel = (
        bool(parallel_enabled)
        and requested_parallel
        and pal_jobs_value > 1
        and len(jobs) > 1
        and width > 1
    )

    if use_parallel:
        # Use global pool to ensure coordination with other workflows
        manager = _WorkflowManager(context.config, label="occupier", max_jobs_override=effective_max_jobs)
        try:
            if effective_max_jobs <= 1 and manager.pool.max_concurrent_jobs != 1:
                manager.pool.max_concurrent_jobs = 1
                manager.max_jobs = 1
                manager._sync_parallel_flag()
            for job in jobs:
                manager.add_job(job)
            dynamic_slots = determine_effective_slots(
                manager.total_cores,
                manager._jobs.values(),
                effective_max_jobs,
                len(jobs),
            )
            if dynamic_slots != manager.pool.max_concurrent_jobs:
                logger.info(
                    "[occupier] Adjusting ORCA job slots to %d (width=%d, requested=%d)",
                    dynamic_slots,
                    len(jobs),
                    effective_max_jobs,
                )
                manager.pool.max_concurrent_jobs = dynamic_slots
                manager.max_jobs = dynamic_slots
                manager._sync_parallel_flag()
            manager.run()
            context.completed_jobs.update(manager.completed_jobs)
            context.failed_jobs.update(manager.failed_jobs)
            context.skipped_jobs.update(manager.skipped_jobs)
            if context.failed_jobs or context.skipped_jobs:
                return False
            return True
        except Exception as exc:  # noqa: BLE001
            logger.error("Parallel OCCUPIER ORCA execution failed: %s", exc, exc_info=True)
            try:
                fallback_jobs = _build_occupier_jobs(context)
            except Exception as rebuild_exc:  # noqa: BLE001
                logger.error(
                    "Sequential fallback cannot be prepared after parallel failure: %s",
                    rebuild_exc,
                    exc_info=True,
                )
                return False
            pre_completed = set(getattr(manager, "completed_jobs", set()) or set())
            if pre_completed:
                context.completed_jobs.update(pre_completed)
            failed_map = getattr(manager, "failed_jobs", {}) or {}
            if failed_map:
                context.failed_jobs.update(dict(failed_map))
            skipped_map = getattr(manager, "skipped_jobs", {}) or {}
            if skipped_map:
                context.skipped_jobs.update({key: list(value) for key, value in skipped_map.items()})
            logger.info("Falling back to sequential OCCUPIER ORCA execution")
            return _run_jobs_sequentially(
                fallback_jobs,
                context,
                pal_jobs_value,
                pre_completed=pre_completed,
            )
        finally:
            try:
                manager.shutdown()
            except Exception:  # noqa: BLE001
                logger.debug("Parallel manager shutdown raised", exc_info=True)

    if parallel_enabled and not requested_parallel:
        logger.info(
            "[occupier] Parallel workflows disabled (mode=%s) → running ORCA jobs sequentially",
            parallel_mode,
        )
    elif parallel_enabled and pal_jobs_value <= 1:
        logger.info("[occupier] Parallel execution requested but PAL_JOBS=1 → running sequentially")
    elif len(jobs) <= 1:
        logger.info("[occupier] Single OCCUPIER ORCA job detected → running sequentially")
    elif parallel_enabled and width <= 1:
        logger.info(
            "[occupier] Parallel mode=%s but dependency graph is serial (width=%d) → running sequentially",
            parallel_mode,
            width,
        )

    # Sequential path or fallback after errors
    return _run_jobs_sequentially(jobs, context, pal_jobs_value)


def _run_jobs_sequentially(
    jobs: List[WorkflowJob],
    context: OccupierExecutionContext,
    pal_jobs_value: int,
    *,
    pre_completed: Optional[Set[str]] = None,
) -> bool:
    """Execute OCCUPIER jobs sequentially while respecting PAL limits."""

    total_cores = max(1, _parse_int(context.config.get('PAL'), fallback=1))
    per_job_cores = total_cores
    initial_completed = set(pre_completed or ())
    completed: Set[str] = set(initial_completed)
    pending = {job.job_id: job for job in jobs if job.job_id not in completed}
    failed: Dict[str, str] = {}
    skipped: Dict[str, List[str]] = {}

    context.completed_jobs.clear()
    context.failed_jobs.clear()
    context.skipped_jobs.clear()
    if initial_completed:
        context.completed_jobs.update(initial_completed)

    while pending:
        progressed = False
        for job_id, job in list(pending.items()):
            if not job.dependencies <= completed:
                continue

            allocated = max(job.cores_min, min(job.cores_max, per_job_cores))
            usage_info = f"{job.description}; {allocated}/{total_cores} cores used"
            logger.info(
                "[occupier] Running %s with %d cores (%s)",
                job_id,
                allocated,
                usage_info,
            )
            try:
                job.work(allocated)
            except Exception as exc:  # noqa: BLE001
                failed[job_id] = f"{exc.__class__.__name__}: {exc}"
                pending.pop(job_id, None)
                progressed = True
                continue

            completed.add(job_id)
            pending.pop(job_id)
            progressed = True

        if not progressed:
            unresolved_msgs: List[str] = []
            for job_id, job in list(pending.items()):
                missing = sorted(job.dependencies - completed)
                skipped[job_id] = missing
                if missing:
                    unresolved_msgs.append(f"{job_id} (waiting for {', '.join(missing)})")
                else:
                    unresolved_msgs.append(job_id)
            if unresolved_msgs:
                logger.error("Unresolved OCCUPIER job dependencies: %s", ", ".join(unresolved_msgs))
            pending.clear()
            break

    context.completed_jobs.update(completed)
    context.failed_jobs.update(failed)
    context.skipped_jobs.update(skipped)

    if failed:
        logger.warning(
            "Sequential OCCUPIER execution completed with failures: %s",
            ", ".join(f"{job_id} ({reason})" for job_id, reason in failed.items()),
        )
    if skipped:
        logger.warning(
            "Sequential OCCUPIER execution skipped jobs due to unmet dependencies: %s",
            ", ".join(
                f"{job_id} (missing {', '.join(deps) if deps else 'unknown cause'})"
                for job_id, deps in skipped.items()
            ),
        )

    return not failed and not skipped


def _build_occupier_jobs(context: OccupierExecutionContext) -> List[WorkflowJob]:
    """Create workflow job definitions for OCCUPIER ORCA runs."""

    config = context.config
    jobs: List[WorkflowJob] = []
    descriptors: List[JobDescriptor] = []

    total_cores = max(1, _parse_int(config.get('PAL'), fallback=1))
    pal_jobs_value = _resolve_pal_jobs(config)

    # Check if ox and red will run in parallel - if so, split cores
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    has_ox = len(oxidation_steps) > 0
    has_red = len(reduction_steps) > 0

    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    ox_red_parallel = (has_ox and has_red) and parallel_mode != 'disable'

    if parallel_mode == 'disable':
        pal_jobs_value = 1

    # If ox and red will run in parallel, split cores between them
    effective_total_cores = total_cores
    if ox_red_parallel:
        effective_total_cores = max(2, total_cores // 2)
        logger.info(
            f"[occupier] Ox and red workflows will run in parallel: "
            f"splitting {total_cores} cores → {effective_total_cores} cores per workflow"
        )

    cores_min = 1 if effective_total_cores == 1 else 2
    cores_max = effective_total_cores

    def core_bounds(preferred_opt: Optional[int] = None, job_count_at_level: Optional[int] = None) -> tuple[int, int, int]:
        """Calculate core bounds with awareness of parallel job potential."""
        if pal_jobs_value > 0:
            default_opt = max(cores_min, effective_total_cores // pal_jobs_value)
        else:
            default_opt = max(cores_min, effective_total_cores)

        # If we know multiple jobs will run in parallel, reduce optimal cores for better sharing
        if job_count_at_level and job_count_at_level > 1:
            # Allow more jobs to run simultaneously by reducing per-job cores
            adjusted_opt = max(cores_min, effective_total_cores // min(job_count_at_level, pal_jobs_value))
            default_opt = min(default_opt, adjusted_opt)

        if preferred_opt is None:
            preferred = max(cores_min, min(cores_max, default_opt))
        else:
            preferred = max(cores_min, min(preferred_opt, cores_max))
        return cores_min, preferred, cores_max

    def register_descriptor(descriptor: JobDescriptor) -> None:
        descriptors.append(descriptor)

    def _fallback_multiplicity(folder: str) -> int:
        folder_lower = folder.lower()
        candidates: List[str] = []
        if folder_lower.startswith('initial_'):
            candidates.extend(['multiplicity_0', 'multiplicity'])
        elif folder_lower.startswith('ox_step_'):
            step_match = re.search(r"ox_step_(\d+)_occupier", folder_lower)
            if step_match:
                step = step_match.group(1)
                candidates.append(f'multiplicity_ox{step}')
        elif folder_lower.startswith('red_step_'):
            step_match = re.search(r"red_step_(\d+)_occupier", folder_lower)
            if step_match:
                step = step_match.group(1)
                candidates.append(f'multiplicity_red{step}')

        for key in candidates:
            value = config.get(key)
            if value is None:
                continue
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                logger.warning("[occupier] Cannot parse %s='%s' from CONTROL.txt; ignoring.", key, value)

        logger.warning(
            "[occupier] Missing OCCUPIER multiplicity for %s; defaulting to 1.",
            folder,
        )
        return 1

    def _fallback_additions(folder: str) -> str:
        folder_lower = folder.lower()
        candidates: List[str] = []
        if folder_lower.startswith('initial_'):
            candidates.append('additions_0')
        elif folder_lower.startswith('ox_step_'):
            step_match = re.search(r"ox_step_(\d+)_occupier", folder_lower)
            if step_match:
                step = step_match.group(1)
                candidates.append(f'additions_ox{step}')
        elif folder_lower.startswith('red_step_'):
            step_match = re.search(r"red_step_(\d+)_occupier", folder_lower)
            if step_match:
                step = step_match.group(1)
                candidates.append(f'additions_red{step}')

        for key in candidates:
            value = config.get(key)
            if value is None:
                continue
            if isinstance(value, str):
                stripped = value.strip()
                if not stripped:
                    continue
                if stripped.startswith('%'):
                    return stripped
                digits = re.fullmatch(r"(\d+)\s*,\s*(\d+)", stripped)
                if digits:
                    first, second = digits.groups()
                    return f"%scf BrokenSym {first},{second} end"
                bs_match = re.fullmatch(r"brokensym\s+(\d+)\s*,\s*(\d+)(?:\s+end)?", stripped, re.IGNORECASE)
                if bs_match:
                    first, second = bs_match.groups()
                    return f"%scf BrokenSym {first},{second} end"
                if re.search(r"[A-Za-z]", stripped):
                    logger.debug(
                        "[occupier] Ignoring non-numeric OCCUPIER additions '%s' from CONTROL key %s.",
                        stripped,
                        key,
                    )
                    continue
                return stripped
            if isinstance(value, (list, tuple)):
                tokens = [str(item).strip() for item in value if str(item).strip()]
                if tokens:
                    return f"%scf BrokenSym {','.join(tokens)} end"

        return ""

    def read_occ(folder: str) -> tuple[int, str, Optional[int]]:
        result = read_occupier_file(folder, "OCCUPIER.txt", None, None, None, config)
        if result:
            multiplicity, additions, min_fspe_index = result
            try:
                multiplicity_int = int(multiplicity)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                logger.warning(
                    "[occupier] OCCUPIER multiplicity invalid for %s; falling back to CONTROL settings.",
                    folder,
                )
                multiplicity_int = _fallback_multiplicity(folder)

            if isinstance(additions, str):
                additions_str = additions.strip()
            else:
                additions_str = ""
            return multiplicity_int, additions_str, min_fspe_index

        logger.warning(
            "[occupier] OCCUPIER.txt missing in %s; using CONTROL multiplicity/additions fallback.",
            folder,
        )
        return _fallback_multiplicity(folder), _fallback_additions(folder), None

    solvent = context.solvent
    metals = context.metals
    metal_basis = context.metal_basisset
    main_basis = context.main_basisset
    base_charge = context.charge
    functional = config.get('functional', 'ORCA')

    calc_initial_flag = str(config.get('calc_initial', 'yes')).strip().lower()
    xtb_solvator_enabled = str(config.get('XTB_SOLVATOR', 'no')).strip().lower() == 'yes'
    if calc_initial_flag == 'yes' or xtb_solvator_enabled:
        multiplicity_0, additions_0, _ = read_occ("initial_OCCUPIER")

        if xtb_solvator_enabled:
            solvated_xyz = Path("XTB_SOLVATOR") / "XTB_SOLVATOR.solvator.xyz"
            target_parent_xyz = Path("input_initial_OCCUPIER.xyz")
            if solvated_xyz.exists():
                try:
                    shutil.copyfile(solvated_xyz, target_parent_xyz)
                    logger.info("[occupier] Enforced solvator geometry for %s", target_parent_xyz)
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "[occupier] Could not copy solvator geometry to %s: %s",
                        target_parent_xyz,
                        exc,
                    )

        def run_initial(cores: int,
                        _mult=multiplicity_0,
                        _adds=additions_0) -> None:
            logger.info("[occupier] Preparing initial frequency job")
            read_xyz_and_create_input3(
                "input_initial_OCCUPIER.xyz",
                "initial.inp",
                base_charge,
                _mult,
                solvent,
                metals,
                metal_basis,
                main_basis,
                config,
                _adds,
            )
            _update_pal_block("initial.inp", cores)
            run_orca("initial.inp", "initial.out")
            if not _verify_orca_output("initial.out"):
                raise RuntimeError("ORCA terminated abnormally for initial.out")
            run_IMAG(
                "initial.out",
                "initial",
                base_charge,
                _mult,
                solvent,
                metals,
                config,
                main_basis,
                metal_basis,
                _adds,
            )
            logger.info(
                "%s %s freq & geometry optimization of the initial system complete!",
                functional,
                main_basis,
            )
            initial_xyz = Path("initial.xyz")
            if not initial_xyz.exists():
                source_xyz = Path("input_initial_OCCUPIER.xyz")
                if source_xyz.exists():
                    shutil.copy(source_xyz, initial_xyz)
                else:
                    logger.warning("initial.xyz missing and no backup geometry found")

        register_descriptor(JobDescriptor(
            job_id="occupier_initial",
            description="initial OCCUPIER frequency job",
            work=run_initial,
            produces={"initial.out", "initial.xyz"},
            preferred_cores=None,
        ))

    if str(config.get('absorption_spec', 'no')).strip().lower() == 'yes':
        additions_tddft = config.get('additions_TDDFT', '')

        def run_absorption(cores: int, _adds=additions_tddft) -> None:
            absorption_source = "initial.xyz" if xtb_solvator_enabled else "input_initial_OCCUPIER.xyz"
            read_xyz_and_create_input2(
                absorption_source,
                "absorption_td.inp",
                base_charge,
                1,
                solvent,
                metals,
                config,
                main_basis,
                metal_basis,
                _adds,
            )
            _update_pal_block("absorption_td.inp", cores)
            run_orca("absorption_td.inp", "absorption_spec.out")
            if not _verify_orca_output("absorption_spec.out"):
                raise RuntimeError("ORCA terminated abnormally for absorption_spec.out")
            logger.info("TD-DFT absorption spectra calculation complete!")

        register_descriptor(JobDescriptor(
            job_id="occupier_absorption",
            description="absorption spectrum",
            work=run_absorption,
            produces={"absorption_spec.out"},
            requires={"initial.xyz"} if xtb_solvator_enabled else set(),
        ))

    excitation_flags = str(config.get('excitation', '')).lower()
    emission_enabled = str(config.get('emission_spec', 'no')).strip().lower() == 'yes'
    additions_tddft = config.get('additions_TDDFT', '')
    xyz_initial = "initial.xyz"

    if 't' in excitation_flags and str(config.get('E_00', 'no')).strip().lower() == 'yes':
        def run_t1_state(cores: int, _adds=additions_tddft) -> None:
            if not Path(xyz_initial).exists():
                raise RuntimeError(f"Required geometry '{xyz_initial}' not found")
            read_xyz_and_create_input3(
                xyz_initial,
                "t1_state_opt.inp",
                base_charge,
                3,
                solvent,
                metals,
                metal_basis,
                main_basis,
                config,
                _adds,
            )
            inp_path = Path("t1_state_opt.inp")
            if not inp_path.exists():
                raise RuntimeError("Failed to create t1_state_opt.inp")
            _update_pal_block(str(inp_path), cores)
            run_orca("t1_state_opt.inp", "t1_state_opt.out")
            if not _verify_orca_output("t1_state_opt.out"):
                raise RuntimeError("ORCA terminated abnormally for t1_state_opt.out")
            logger.info(
                "%s %s freq & geometry optimization of T_1 complete!",
                functional,
                main_basis,
            )

        t1_job_id = "occupier_t1_state"
        register_descriptor(JobDescriptor(
            job_id=t1_job_id,
            description="triplet state optimization",
            work=run_t1_state,
            produces={"t1_state_opt.xyz", "t1_state_opt.out"},
            requires={"initial.xyz"},
        ))

        if emission_enabled:
            def run_t1_emission(cores: int, _adds=additions_tddft) -> None:
                read_xyz_and_create_input2(
                    "t1_state_opt.xyz",
                    "emission_t1.inp",
                    base_charge,
                    1,
                    solvent,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    _adds,
                )
                inp_path = Path("emission_t1.inp")
                if not inp_path.exists():
                    raise RuntimeError("Failed to create emission_t1.inp")
                _update_pal_block(str(inp_path), cores)
                run_orca("emission_t1.inp", "emission_t1.out")
                if not _verify_orca_output("emission_t1.out"):
                    raise RuntimeError("ORCA terminated abnormally for emission_t1.out")
                logger.info("TD-DFT T1 emission spectra calculation complete!")

            register_descriptor(JobDescriptor(
                job_id="occupier_t1_emission",
                description="triplet emission spectrum",
                work=run_t1_emission,
                produces={"emission_t1.out"},
                requires={"t1_state_opt.xyz"},
            ))

    if 's' in excitation_flags and str(config.get('E_00', 'no')).strip().lower() == 'yes':
        def run_s1_state(cores: int, _adds=additions_tddft) -> None:
            if not Path(xyz_initial).exists():
                raise RuntimeError(f"Required geometry '{xyz_initial}' not found")
            failed_flag = Path("s1_state_opt.failed")
            if failed_flag.exists():
                try:
                    failed_flag.unlink()
                except Exception:  # noqa: BLE001
                    pass
            create_s1_optimization_input(
                xyz_initial,
                "s1_state_opt.inp",
                base_charge,
                1,
                solvent,
                metals,
                metal_basis,
                main_basis,
                config,
                _adds,
            )
            inp_path = Path("s1_state_opt.inp")
            if not inp_path.exists():
                raise RuntimeError("Failed to create s1_state_opt.inp")
            _update_pal_block(str(inp_path), cores)
            try:
                run_orca("s1_state_opt.inp", "s1_state_opt.out")
                if not _verify_orca_output("s1_state_opt.out"):
                    raise RuntimeError("ORCA terminated abnormally for s1_state_opt.out")
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "[occupier] Skipping singlet state optimization; ORCA failed: %s",
                    exc,
                )
                try:
                    failed_flag.write_text(str(exc))
                except Exception:  # noqa: BLE001
                    pass
                return
            logger.info(
                "%s %s freq & geometry optimization of S_1 complete!",
                functional,
                main_basis,
            )
            if failed_flag.exists():
                try:
                    failed_flag.unlink()
                except Exception:  # noqa: BLE001
                    pass

        s1_job_id = "occupier_s1_state"
        register_descriptor(JobDescriptor(
            job_id=s1_job_id,
            description="singlet state optimization",
            work=run_s1_state,
            produces={"s1_state_opt.xyz", "s1_state_opt.out"},
            requires={"initial.xyz"},
        ))

        if emission_enabled:
            def run_s1_emission(cores: int, _adds=additions_tddft) -> None:
                failed_flag = Path("s1_state_opt.failed")
                if failed_flag.exists():
                    logger.info(
                        "[occupier] Skipping singlet emission; singlet optimization failed (see %s).",
                        failed_flag,
                    )
                    return
                if not Path("s1_state_opt.xyz").exists():
                    logger.info(
                        "[occupier] Skipping singlet emission; missing geometry 's1_state_opt.xyz'.",
                    )
                    return
                read_xyz_and_create_input2(
                    "s1_state_opt.xyz",
                    "emission_s1.inp",
                    base_charge,
                    1,
                    solvent,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    _adds,
                )
                inp_path = Path("emission_s1.inp")
                if not inp_path.exists():
                    raise RuntimeError("Failed to create emission_s1.inp")
                _update_pal_block(str(inp_path), cores)
                run_orca("emission_s1.inp", "emission_s1.out")
                if not _verify_orca_output("emission_s1.out"):
                    raise RuntimeError("ORCA terminated abnormally for emission_s1.out")
                logger.info("TD-DFT S1 emission spectra calculation complete!")

            register_descriptor(JobDescriptor(
                job_id="occupier_s1_emission",
                description="singlet emission spectrum",
                work=run_s1_emission,
                produces={"emission_s1.out"},
                requires={"s1_state_opt.xyz"},
            ))

    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    for step in oxidation_steps:
        folder = f"ox_step_{step}_OCCUPIER"
        multiplicity_step, additions_step, _ = read_occ(folder)
        if xtb_solvator_enabled:
            if step == 1:
                primary_geom = Path("initial.xyz")
                requires = {"initial.out"}
            else:
                primary_geom = Path(f"ox_step_{step - 1}.xyz")
                requires = {f"ox_step_{step - 1}.out"}
        else:
            primary_geom = Path(f"input_ox_step_{step}_OCCUPIER.xyz")
            requires = set()

        if primary_geom.exists() or xtb_solvator_enabled:
            xyz_source_path = primary_geom
        else:
            fallback_geom = Path(f"input_ox_step_{step}_OCCUPIER.xyz")
            if fallback_geom.exists():
                if primary_geom != fallback_geom:
                    logger.warning(
                        "[occupier] Primary oxidation geometry %s missing; using OCCUPIER fallback %s",
                        primary_geom,
                        fallback_geom,
                    )
                xyz_source_path = fallback_geom
            else:
                logger.warning(
                    "[occupier] No geometry found for oxidation step %d; proceeding with %s",
                    step,
                    primary_geom,
                )
                xyz_source_path = primary_geom
        xyz_source = str(xyz_source_path)
        inp_path = f"ox_step_{step}.inp"
        out_path = f"ox_step_{step}.out"
        step_charge = base_charge + step

        def make_oxidation_work(idx: int, mult: int, adds: str,
                                xyz_path: str, inp: str, out: str,
                                charge_value: int) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                read_xyz_and_create_input3(
                    xyz_path,
                    inp,
                    charge_value,
                    mult,
                    solvent,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    adds,
                )
                inp_path = Path(inp)
                if not inp_path.exists():
                    raise RuntimeError(f"Failed to create {inp}")
                _update_pal_block(str(inp_path), cores)
                run_orca(inp, out)
                if not _verify_orca_output(out):
                    raise RuntimeError(f"ORCA terminated abnormally for {out}")
                logger.info(
                    "%s %s freq & geometry optimization cation (step %d) complete!",
                    functional,
                    main_basis,
                    idx,
                )

            return _work

        register_descriptor(JobDescriptor(
            job_id=f"occupier_ox_{step}",
            description=f"oxidation step {step}",
            work=make_oxidation_work(step, multiplicity_step, additions_step, xyz_source, inp_path, out_path, step_charge),
            produces={out_path, f"ox_step_{step}.xyz"},
            requires=requires,
        ))

    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    for step in reduction_steps:
        folder = f"red_step_{step}_OCCUPIER"
        multiplicity_step, additions_step, _ = read_occ(folder)
        if xtb_solvator_enabled:
            if step == 1:
                primary_geom = Path("initial.xyz")
                requires = {"initial.out"}
            else:
                primary_geom = Path(f"red_step_{step - 1}.xyz")
                requires = {f"red_step_{step - 1}.out"}
        else:
            primary_geom = Path(f"input_red_step_{step}_OCCUPIER.xyz")
            requires = set()

        if primary_geom.exists() or xtb_solvator_enabled:
            xyz_source_path = primary_geom
        else:
            fallback_geom = Path(f"input_red_step_{step}_OCCUPIER.xyz")
            if fallback_geom.exists():
                if primary_geom != fallback_geom:
                    logger.warning(
                        "[occupier] Primary reduction geometry %s missing; using OCCUPIER fallback %s",
                        primary_geom,
                        fallback_geom,
                    )
                xyz_source_path = fallback_geom
            else:
                logger.warning(
                    "[occupier] No geometry found for reduction step %d; proceeding with %s",
                    step,
                    primary_geom,
                )
                xyz_source_path = primary_geom
        xyz_source = str(xyz_source_path)
        inp_path = f"red_step_{step}.inp"
        out_path = f"red_step_{step}.out"
        step_charge = base_charge - step

        def make_reduction_work(idx: int, mult: int, adds: str,
                                xyz_path: str, inp: str, out: str,
                                charge_value: int) -> Callable[[int], None]:
            def _work(cores: int) -> None:
                read_xyz_and_create_input3(
                    xyz_path,
                    inp,
                    charge_value,
                    mult,
                    solvent,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    adds,
                )
                inp_path = Path(inp)
                if not inp_path.exists():
                    raise RuntimeError(f"Failed to create {inp}")
                _update_pal_block(str(inp_path), cores)
                run_orca(inp, out)
                if not _verify_orca_output(out):
                    raise RuntimeError(f"ORCA terminated abnormally for {out}")
                logger.info(
                    "%s %s freq & geometry optimization anion (step %d) complete!",
                    functional,
                    main_basis,
                    idx,
                )

            return _work

        register_descriptor(JobDescriptor(
            job_id=f"occupier_red_{step}",
            description=f"reduction step {step}",
            work=make_reduction_work(step, multiplicity_step, additions_step, xyz_source, inp_path, out_path, step_charge),
            produces={out_path, f"red_step_{step}.xyz"},
            requires=requires,
        ))

    # Resolve implicit dependencies based on produced artifacts
    produced_by: Dict[str, str] = {}
    for descriptor in descriptors:
        for artifact in descriptor.produces:
            produced_by.setdefault(artifact, descriptor.job_id)

    # Build dependency graph
    job_deps: Dict[str, Set[str]] = {}
    for descriptor in descriptors:
        dependencies: Set[str] = set(descriptor.explicit_dependencies)
        for requirement in descriptor.requires:
            producer = produced_by.get(requirement)
            if producer and producer != descriptor.job_id:
                dependencies.add(producer)
        job_deps[descriptor.job_id] = dependencies

    # Calculate dependency levels for better parallelization
    def get_dependency_level(job_id: str, memo: Dict[str, int]) -> int:
        """Get the dependency level of a job (0 = no deps, 1 = depends on level 0, etc.)."""
        if job_id in memo:
            return memo[job_id]
        deps = job_deps.get(job_id, set())
        if not deps:
            memo[job_id] = 0
            return 0
        level = max(get_dependency_level(dep, memo) for dep in deps) + 1
        memo[job_id] = level
        return level

    level_memo: Dict[str, int] = {}
    job_levels: Dict[str, int] = {}
    for descriptor in descriptors:
        job_levels[descriptor.job_id] = get_dependency_level(descriptor.job_id, level_memo)

    # Count jobs at each level for better core allocation
    levels_count: Dict[int, int] = {}
    for level in job_levels.values():
        levels_count[level] = levels_count.get(level, 0) + 1

    # Build WorkflowJob objects with optimized core allocation
    for descriptor in descriptors:
        dependencies = job_deps[descriptor.job_id]
        job_level = job_levels[descriptor.job_id]
        parallel_jobs_at_level = levels_count.get(job_level, 1)

        # Use parallel job count to optimize core allocation
        cores_min_v, cores_opt_v, cores_max_v = core_bounds(
            descriptor.preferred_cores,
            job_count_at_level=parallel_jobs_at_level if parallel_jobs_at_level > 1 else None
        )

        jobs.append(
            WorkflowJob(
                job_id=descriptor.job_id,
                work=descriptor.work,
                description=descriptor.description,
                dependencies=dependencies,
                cores_min=cores_min_v,
                cores_optimal=cores_opt_v,
                cores_max=cores_max_v,
            )
        )

    _log_job_plan_with_levels(descriptors, job_levels, levels_count)
    return jobs


def _resolve_pal_jobs(config: Dict[str, Any]) -> int:
    value = config.get('pal_jobs')
    parsed = _parse_int(value, fallback=0)
    if parsed <= 0:
        total = max(1, _parse_int(config.get('PAL'), fallback=1))
        return max(1, min(4, max(1, total // 2)))
    return parsed


def _log_job_plan_with_levels(
    descriptors: List[JobDescriptor],
    job_levels: Dict[str, int],
    levels_count: Dict[int, int]
) -> None:
    """Log job plan with dependency levels for parallelization analysis."""
    logger.info("Planned OCCUPIER ORCA jobs (%d total):", len(descriptors))

    # Group jobs by level
    jobs_by_level: Dict[int, List[JobDescriptor]] = {}
    for descriptor in descriptors:
        level = job_levels.get(descriptor.job_id, 0)
        if level not in jobs_by_level:
            jobs_by_level[level] = []
        jobs_by_level[level].append(descriptor)

    # Log summary of parallelization potential
    max_parallel = max(levels_count.values()) if levels_count else 0
    logger.info(
        "Parallelization potential: %d levels, max %d jobs in parallel",
        len(levels_count),
        max_parallel
    )

    # Log jobs grouped by level
    for level in sorted(jobs_by_level.keys()):
        job_list = jobs_by_level[level]
        logger.info("  Level %d (%d jobs can run in parallel):", level, len(job_list))
        for descriptor in job_list:
            deps = sorted(descriptor.explicit_dependencies | descriptor.requires)
            produces = sorted(descriptor.produces)
            logger.info(
                "    - %s: %s | deps=%s | outputs=%s",
                descriptor.job_id,
                descriptor.description,
                deps or ['none'],
                produces or ['none'],
            )


def _parse_step_list(raw_steps: Any) -> List[int]:
    if not raw_steps:
        return []
    tokens: List[str]
    if isinstance(raw_steps, str):
        cleaned = raw_steps.replace(';', ',')
        tokens = [token.strip() for token in cleaned.split(',')]
    else:
        tokens = []
        for item in raw_steps:
            tokens.extend(str(item).split(','))
    result: Set[int] = set()
    for token in tokens:
        if not token:
            continue
        try:
            value = int(token)
        except ValueError:
            continue
        if value >= 1:
            result.add(value)
    return sorted(result)


def should_use_parallel_occupier(config: Dict[str, Any]) -> bool:
    """Determine if parallel OCCUPIER execution would be beneficial."""
    total_cores = config.get('PAL', 1)

    # Enable parallel execution if we have sufficient resources
    # Lowered threshold - even 4 cores can benefit from parallelization
    return total_cores >= 4
