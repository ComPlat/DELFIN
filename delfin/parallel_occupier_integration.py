"""Integration of dynamic pool with OCCUPIER workflow."""

import time
from typing import Dict, List, Set, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

from delfin.common.logging import get_logger
from delfin.dynamic_pool import DynamicCorePool, PoolJob, JobPriority, create_orca_job

logger = get_logger(__name__)


class ParallelOccupierManager:
    """Manages parallel execution of OCCUPIER with dynamic resource allocation."""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.total_cores = config.get('PAL', 1)
        self.total_memory = config.get('maxcore', 1000) * self.total_cores

        # Create dynamic pool
        self.pool = DynamicCorePool(
            total_cores=self.total_cores,
            total_memory_mb=self.total_memory,
            max_jobs=min(4, self.total_cores // 2)  # Max 4 concurrent jobs
        )

        logger.info(f"Parallel OCCUPIER manager initialized with {self.total_cores} cores")

    def execute_parallel_workflows(self, ox_sequence: List[Dict], red_sequence: List[Dict]) -> bool:
        """Execute ox and red workflows in parallel with dynamic resource management."""

        logger.info("Starting parallel OCCUPIER execution: ox_steps + red_steps")

        # Create workflow executors
        workflows = []

        if ox_sequence:
            workflows.append(('ox_steps', ox_sequence, JobPriority.NORMAL))

        if red_sequence:
            workflows.append(('red_steps', red_sequence, JobPriority.NORMAL))

        if not workflows:
            logger.warning("No ox or red sequences to execute")
            return True

        # Execute workflows in parallel
        with ThreadPoolExecutor(max_workers=len(workflows)) as executor:
            futures = []

            for workflow_name, sequence, priority in workflows:
                future = executor.submit(
                    self._execute_sequence_with_pool,
                    workflow_name, sequence, priority
                )
                futures.append((workflow_name, future))

            # Wait for all workflows to complete
            all_success = True
            for workflow_name, future in futures:
                try:
                    success = future.result()
                    if success:
                        logger.info(f"Workflow {workflow_name} completed successfully")
                    else:
                        logger.error(f"Workflow {workflow_name} failed")
                        all_success = False
                except Exception as e:
                    logger.error(f"Workflow {workflow_name} raised exception: {e}")
                    all_success = False

        return all_success

    def _execute_sequence_with_pool(self, workflow_name: str, sequence: List[Dict],
                                   priority: JobPriority) -> bool:
        """Execute a single OCCUPIER sequence using the dynamic pool."""

        logger.info(f"Starting {workflow_name} with {len(sequence)} jobs")

        # Analyze dependencies
        dependencies = self._analyze_dependencies(sequence)

        # Submit jobs to pool based on dependencies
        submitted_jobs = {}
        completed_jobs = set()

        start_time = time.time()

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

                logger.info(f"Submitted {job_id} to pool")

            # Wait for some jobs to complete
            if ready_jobs:
                time.sleep(1)  # Brief pause to let jobs start

            # Check for completed jobs
            newly_completed = self._check_completed_jobs(submitted_jobs, completed_jobs)
            completed_jobs.update(newly_completed)

            # Avoid busy waiting
            if not ready_jobs and not newly_completed:
                time.sleep(2)

        duration = time.time() - start_time
        logger.info(f"{workflow_name} completed in {duration:.1f}s")

        return True

    def _analyze_dependencies(self, sequence: List[Dict]) -> Dict[int, Set[int]]:
        """Analyze dependencies in OCCUPIER sequence."""
        dependencies = {}

        for entry in sequence:
            idx = entry["index"]
            raw_from = entry.get("from", idx - 1)

            try:
                from_idx = int(raw_from)
                if from_idx == idx or from_idx <= 0:
                    dependencies[idx] = set()
                else:
                    dependencies[idx] = {from_idx}
            except (ValueError, TypeError):
                dependencies[idx] = set()

        return dependencies

    def _estimate_job_requirements(self, entry: Dict) -> tuple[int, int, int]:
        """Estimate core requirements for a job."""
        # Base requirements
        cores_min = 2
        cores_max = min(16, self.total_cores)

        # Adjust based on multiplicity and job characteristics
        multiplicity = entry.get("m", 1)

        if multiplicity > 3:
            # High-spin jobs might benefit from more cores
            cores_opt = min(12, self.total_cores // 2)
        else:
            # Standard DFT jobs
            cores_opt = min(8, self.total_cores // 3)

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

        for idx, job_id in submitted_jobs.items():
            if idx not in completed_jobs:
                out_file = self._get_output_filename(idx)
                if self._verify_job_completion(out_file):
                    newly_completed.add(idx)
                    logger.info(f"Job {job_id} completed")

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
        logger.info("Shutting down parallel OCCUPIER manager")
        self.pool.shutdown()


def should_use_parallel_occupier(config: Dict[str, Any]) -> bool:
    """Determine if parallel OCCUPIER execution would be beneficial."""
    total_cores = config.get('PAL', 1)

    # Enable parallel execution if we have sufficient resources
    # Lowered threshold - even 4 cores can benefit from parallelization
    return total_cores >= 4