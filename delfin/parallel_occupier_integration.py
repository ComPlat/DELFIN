"""Integration of dynamic pool with OCCUPIER workflow."""

import time
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

from delfin.common.logging import get_logger
from delfin.dynamic_pool import DynamicCorePool, PoolJob, JobPriority, create_orca_job
from delfin.copy_helpers import read_occupier_file
from delfin.imag import run_IMAG
from delfin.orca import run_orca
from delfin.xyz_io import (
    read_xyz_and_create_input2,
    read_xyz_and_create_input3,
    read_xyz_and_create_input4,
)
from .parallel_classic import (
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
        self.total_cores = config.get('PAL', 1)
        self.max_jobs = max(1, _resolve_pal_jobs(config))
        self.total_memory = config.get('maxcore', 1000) * self.total_cores
        self.per_job_share = max(1 if self.total_cores == 1 else 2,
                                 max(1, self.total_cores // max(1, self.max_jobs)))

        # Create dynamic pool
        self.pool = DynamicCorePool(
            total_cores=self.total_cores,
            total_memory_mb=self.total_memory,
            max_jobs=min(self.max_jobs, max(1, self.total_cores))
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
        # Base requirements
        cores_min = 1 if self.total_cores == 1 else 2
        burst_capacity = max(cores_min, min(self.total_cores, self.per_job_share * 2))
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


@dataclass
class OccupierExecutionContext:
    """Container for OCCUPIER ORCA execution parameters."""

    charge: int
    solvent: str
    metals: List[str]
    main_basisset: str
    metal_basisset: str
    config: Dict[str, Any]


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
        # First, run the initial job separately (may trigger IMAG)
        initial_job = None
        other_jobs = []
        for job in jobs:
            if job.job_id == "occupier_initial":
                initial_job = job
            else:
                other_jobs.append(job)

        # Run initial job first
        if initial_job:
            logger.info("[occupier] Running initial job first before parallel execution")
            initial_job.work(max(1, _parse_int(context.config.get('PAL'), fallback=1)))

        # Rebuild remaining jobs after potential IMAG modifications
        try:
            updated_jobs = _build_occupier_jobs(context)
            # Filter out the initial job since it's already done
            remaining_jobs = [job for job in updated_jobs if job.job_id != "occupier_initial"]

            if not remaining_jobs:
                logger.info("No remaining jobs after initial completion")
                return True

        except Exception as rebuild_exc:
            logger.warning("Failed to rebuild jobs after initial: %s", rebuild_exc)
            remaining_jobs = other_jobs  # fallback to original jobs

        # Run remaining jobs in parallel
        manager = _WorkflowManager(context.config, label="occupier", max_jobs_override=effective_max_jobs)
        try:
            for job in remaining_jobs:
                manager.add_job(job)
            dynamic_slots = determine_effective_slots(
                manager.total_cores,
                manager._jobs.values(),
                effective_max_jobs,
                len(remaining_jobs),
            )
            if dynamic_slots != manager.pool.max_concurrent_jobs:
                logger.info(
                    "[occupier] Adjusting ORCA job slots to %d (width=%d, requested=%d)",
                    dynamic_slots,
                    len(remaining_jobs),
                    effective_max_jobs,
                )
                manager.pool.max_concurrent_jobs = dynamic_slots
                manager.max_jobs = dynamic_slots
            manager.run()
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
            logger.info("Falling back to sequential OCCUPIER ORCA execution")
            return _run_jobs_sequentially(fallback_jobs, context, pal_jobs_value)
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


def _run_jobs_sequentially(jobs: List[WorkflowJob], context: OccupierExecutionContext,
                           pal_jobs_value: int) -> bool:
    """Execute OCCUPIER jobs sequentially while respecting PAL limits."""

    total_cores = max(1, _parse_int(context.config.get('PAL'), fallback=1))
    per_job_cores = total_cores
    pending = {job.job_id: job for job in jobs}
    completed: Set[str] = set()

    while pending:
        progressed = False
        for job_id, job in list(pending.items()):
            if not job.dependencies <= completed:
                continue

            allocated = max(job.cores_min, min(job.cores_max, per_job_cores))
            logger.info(
                "[occupier] Running %s with %d cores (%s)",
                job_id,
                allocated,
                job.description,
            )
            try:
                job.work(allocated)
            except Exception as exc:  # noqa: BLE001
                logger.error("Sequential OCCUPIER job %s failed: %s", job_id, exc, exc_info=True)
                return False

            completed.add(job_id)
            pending.pop(job_id)
            progressed = True

        if not progressed:
            logger.error("Unresolved OCCUPIER job dependencies: %s", ", ".join(sorted(pending)))
            return False

    return True


def _build_occupier_jobs(context: OccupierExecutionContext) -> List[WorkflowJob]:
    """Create workflow job definitions for OCCUPIER ORCA runs."""

    config = context.config
    jobs: List[WorkflowJob] = []
    descriptors: List[JobDescriptor] = []

    total_cores = max(1, _parse_int(config.get('PAL'), fallback=1))
    pal_jobs_value = _resolve_pal_jobs(config)
    cores_min = 1 if total_cores == 1 else 2
    cores_max = total_cores

    def core_bounds(preferred_opt: Optional[int] = None) -> tuple[int, int, int]:
        if pal_jobs_value > 0:
            default_opt = max(cores_min, total_cores // pal_jobs_value)
        else:
            default_opt = max(cores_min, total_cores)
        if preferred_opt is None:
            preferred = max(cores_min, min(cores_max, default_opt))
        else:
            preferred = max(cores_min, min(preferred_opt, cores_max))
        return cores_min, preferred, cores_max

    def register_descriptor(descriptor: JobDescriptor) -> None:
        descriptors.append(descriptor)

    def read_occ(folder: str) -> tuple[int, str, Optional[int]]:
        result = read_occupier_file(folder, "OCCUPIER.txt", None, None, None, config)
        if not result:
            raise RuntimeError(f"read_occupier_file failed for '{folder}'")
        multiplicity, additions, min_fspe_index = result
        return int(multiplicity), additions, min_fspe_index

    solvent = context.solvent
    metals = context.metals
    metal_basis = context.metal_basisset
    main_basis = context.main_basisset
    base_charge = context.charge
    functional = config.get('functional', 'ORCA')

    if str(config.get('calc_initial', 'yes')).strip().lower() == 'yes':
        multiplicity_0, additions_0, _ = read_occ("initial_OCCUPIER")

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
            read_xyz_and_create_input2(
                "input_initial_OCCUPIER.xyz",
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
            read_xyz_and_create_input4(
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
            run_orca("s1_state_opt.inp", "s1_state_opt.out")
            if not _verify_orca_output("s1_state_opt.out"):
                raise RuntimeError("ORCA terminated abnormally for s1_state_opt.out")
            logger.info(
                "%s %s freq & geometry optimization of S_1 complete!",
                functional,
                main_basis,
            )

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
        xyz_source = f"input_ox_step_{step}_OCCUPIER.xyz"
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
            produces={out_path},
        ))

    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    for step in reduction_steps:
        folder = f"red_step_{step}_OCCUPIER"
        multiplicity_step, additions_step, _ = read_occ(folder)
        xyz_source = f"input_red_step_{step}_OCCUPIER.xyz"
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
            produces={out_path},
        ))

    # Resolve implicit dependencies based on produced artifacts
    produced_by: Dict[str, str] = {}
    for descriptor in descriptors:
        for artifact in descriptor.produces:
            produced_by.setdefault(artifact, descriptor.job_id)

    for descriptor in descriptors:
        dependencies: Set[str] = set(descriptor.explicit_dependencies)
        for requirement in descriptor.requires:
            producer = produced_by.get(requirement)
            if producer and producer != descriptor.job_id:
                dependencies.add(producer)

        cores_min_v, cores_opt_v, cores_max_v = core_bounds(descriptor.preferred_cores)
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

    _log_job_plan(descriptors)
    return jobs


def _resolve_pal_jobs(config: Dict[str, Any]) -> int:
    value = config.get('pal_jobs')
    parsed = _parse_int(value, fallback=0)
    if parsed <= 0:
        total = max(1, _parse_int(config.get('PAL'), fallback=1))
        return max(1, min(4, max(1, total // 2)))
    return parsed


def _log_job_plan(descriptors: List[JobDescriptor]) -> None:
    logger.info("Planned OCCUPIER ORCA jobs (%d total):", len(descriptors))
    for descriptor in descriptors:
        deps = sorted(descriptor.explicit_dependencies | descriptor.requires)
        produces = sorted(descriptor.produces)
        logger.info(
            "  - %s: %s | deps=%s | outputs=%s",
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
