"""Integration of dynamic pool with OCCUPIER workflow."""

import os
import re
import time
import shutil
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Set

from delfin.common.logging import get_logger
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
    _add_moinp_block,
    _verify_orca_output,
    estimate_parallel_width,
    determine_effective_slots,
    normalize_parallel_token,
)

if TYPE_CHECKING:
    from .global_scheduler import GlobalOrcaScheduler

logger = get_logger(__name__)


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


def run_occupier_orca_jobs(
    context: OccupierExecutionContext,
    parallel_enabled: bool,
    *,
    scheduler: Optional["GlobalOrcaScheduler"] = None,
    jobs: Optional[List[WorkflowJob]] = None,
) -> bool:
    """Execute OCCUPIER post-processing ORCA jobs with optional parallelization."""

    frequency_mode = str(context.config.get('frequency_calculation_OCCUPIER', 'no')).lower()
    if frequency_mode == 'yes':
        logger.info("frequency_calculation_OCCUPIER=yes → skipping ORCA job scheduling")
        return True

    if jobs is None:
        try:
            jobs = build_occupier_jobs(context)
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to prepare OCCUPIER ORCA jobs: %s", exc, exc_info=True)
            return False

    if not jobs:
        logger.info("No OCCUPIER ORCA jobs detected for execution")
        return True

    context.completed_jobs.clear()
    context.failed_jobs.clear()
    context.skipped_jobs.clear()

    if scheduler is not None:
        scheduler.add_jobs(jobs)
        result = scheduler.run()
        context.completed_jobs.update(result.completed)
        context.failed_jobs.update(result.failed)
        context.skipped_jobs.update(result.skipped)
        return result.success

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
                fallback_jobs = jobs if jobs is not None else build_occupier_jobs(context)
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


def build_occupier_jobs(
    context: OccupierExecutionContext,
    *,
    planning_only: bool = False,
    include_auxiliary: bool = True,
) -> List[WorkflowJob]:
    """Create workflow job definitions for OCCUPIER ORCA runs."""

    config = context.config
    jobs: List[WorkflowJob] = []
    descriptors: List[JobDescriptor] = []
    occ_results: Dict[str, Dict[str, Any]] = config.setdefault('_occ_results_runtime', {})

    total_cores = max(1, _parse_int(config.get('PAL'), fallback=1))
    pal_jobs_value = _resolve_pal_jobs(config)

    # Determine whether oxidation and reduction flows may run side-by-side
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    has_ox = len(oxidation_steps) > 0
    has_red = len(reduction_steps) > 0

    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    ox_red_parallel = (has_ox and has_red) and parallel_mode != 'disable'

    if parallel_mode == 'disable':
        pal_jobs_value = 1

    max_allocatable = total_cores
    cores_min = 1 if max_allocatable == 1 else 2

    # Suggested share if both workflows run concurrently
    workflow_parallel_share = max_allocatable
    if ox_red_parallel:
        workflow_parallel_share = max(cores_min, max_allocatable // 2)
        logger.info(
            f"[occupier] Oxidation and reduction may run in parallel – "
            f"target share ≈ {workflow_parallel_share}/{max_allocatable} cores per workflow"
        )

    def _preferred_share(job_count: Optional[int]) -> int:
        """Heuristic for an optimal core share per job."""
        share = workflow_parallel_share

        if job_count and job_count > 1:
            divisor = max(1, job_count)
            share = min(share, max_allocatable // divisor if max_allocatable >= divisor else cores_min)
            if pal_jobs_value > 0:
                pal_div = max(1, pal_jobs_value)
                share = min(share, max_allocatable // pal_div if max_allocatable >= pal_div else cores_min)
        else:
            # Single job → allow full PAL
            share = max_allocatable

        return max(cores_min, min(max_allocatable, share))

    def core_bounds(preferred_opt: Optional[int] = None,
                    job_count_at_level: Optional[int] = None) -> tuple[int, int, int]:
        """Calculate core bounds with awareness of parallel job potential."""
        default_opt = _preferred_share(job_count_at_level)
        if preferred_opt is not None:
            preferred = max(cores_min, min(preferred_opt, max_allocatable))
        else:
            preferred = default_opt
        return cores_min, preferred, max_allocatable

    def register_descriptor(descriptor: JobDescriptor) -> None:
        descriptors.append(descriptor)

    def _control_multiplicity(step_type: str, step: Optional[int] = None) -> int:
        if step_type == "initial":
            keys = ["multiplicity_0", "multiplicity"]
        elif step_type == "ox":
            keys = [f"multiplicity_ox{step}"] if step else []
        else:
            keys = [f"multiplicity_red{step}"] if step else []
        for key in keys:
            value = config.get(key)
            if value is None:
                continue
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                logger.debug("[occupier] Cannot parse %s='%s' from CONTROL.txt.", key, value)
        return 1

    def _control_additions(step_type: str, step: Optional[int] = None) -> str:
        if step_type == "initial":
            keys = ["additions_0"]
        elif step_type == "ox":
            keys = [f"additions_ox{step}"] if step else []
        else:
            keys = [f"additions_red{step}"] if step else []
        for key in keys:
            value = config.get(key)
            if value is None:
                continue
            if isinstance(value, str):
                stripped = value.strip()
                if not stripped:
                    continue
                if stripped.startswith("%"):
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
                    logger.debug("[occupier] Ignoring non-numeric OCCUPIER additions '%s' from CONTROL key %s.", stripped, key)
                    continue
                return stripped
            if isinstance(value, (list, tuple)):
                tokens = [str(item).strip() for item in value if str(item).strip()]
                if tokens:
                    return f"%scf BrokenSym {','.join(tokens)} end"
        return ""

    def read_occ_from_control(step_type: str, step: Optional[int] = None) -> tuple[int, str, Optional[int]]:
        return _control_multiplicity(step_type, step), _control_additions(step_type, step), None

    def read_occ(folder: str, step_type: str, step: Optional[int]) -> tuple[int, str, Optional[int]]:
        folder_path = Path(folder)
        report_path = folder_path / "OCCUPIER.txt"
        cached = occ_results.get(folder)
        if cached:
            cached_mult = cached.get("multiplicity")
            cached_adds = cached.get("additions", "")
            cached_index = cached.get("preferred_index")
            if cached_mult is not None:
                return cached_mult, cached_adds, cached_index
        if not folder_path.is_dir() or not report_path.is_file():
            raise RuntimeError(
                f"Required OCCUPIER results for '{folder}' not available (missing OCCUPIER.txt)."
            )

        result = read_occupier_file(folder, "OCCUPIER.txt", None, None, None, config, verbose=False)
        if not result:
            raise RuntimeError(
                f"Unable to read OCCUPIER results for '{folder}'."
            )

        multiplicity, additions, min_fspe_index, _gbw_path = result
        try:
            multiplicity_int = int(multiplicity)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            raise RuntimeError(
                f"Preferred multiplicity missing or invalid in OCCUPIER results for '{folder}'."
            ) from None

        additions_str = additions.strip() if isinstance(additions, str) else ""
        occ_results[folder] = {
            "multiplicity": multiplicity_int,
            "additions": additions_str,
            "preferred_index": min_fspe_index,
        }
        return multiplicity_int, additions_str, min_fspe_index

    solvent = context.solvent
    metals = context.metals
    metal_basis = context.metal_basisset
    main_basis = context.main_basisset
    base_charge = context.charge
    functional = config.get('functional', 'ORCA')

    # Cache OCCUPIER outcomes (multiplicity/additions/index) for reuse by post-jobs
    occ_results: Dict[str, Dict[str, Any]] = {}
    config['_occ_results_runtime'] = occ_results

    calc_initial_flag = str(config.get('calc_initial', 'yes')).strip().lower()
    xtb_solvator_enabled = str(config.get('XTB_SOLVATOR', 'no')).strip().lower() == 'yes'
    if calc_initial_flag == 'yes' or xtb_solvator_enabled:
        multiplicity_0, additions_0, _ = read_occ("initial_OCCUPIER", "initial", None)

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
            mult_val = _mult
            adds_val = _adds
            try:
                dyn_mult, dyn_adds, _ = read_occ("initial_OCCUPIER", "initial", None)
                if dyn_mult:
                    mult_val = dyn_mult
                if isinstance(dyn_adds, str):
                    adds_val = dyn_adds
            except Exception:  # noqa: BLE001
                logger.debug("[occupier] Using fallback multiplicity/additions for initial job", exc_info=True)

            read_xyz_and_create_input3(
                "input_initial_OCCUPIER.xyz",
                "initial.inp",
                base_charge,
                mult_val,
                solvent,
                metals,
                metal_basis,
                main_basis,
                config,
                adds_val,
            )
            _update_pal_block("initial.inp", cores)

            # Add %moinp block to reuse OCCUPIER wavefunction
            gbw_initial = Path("input_initial_OCCUPIER.gbw")
            if gbw_initial.exists():
                _add_moinp_block("initial.inp", str(gbw_initial))
                logger.info("[occupier_initial] Using GBW from OCCUPIER: %s", gbw_initial)

            if not run_orca("initial.inp", "initial.out"):
                raise RuntimeError("ORCA terminated abnormally for initial.out")
            run_IMAG(
                "initial.out",
                "initial",
                base_charge,
                mult_val,
                solvent,
                metals,
                config,
                main_basis,
                metal_basis,
                adds_val,
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

    if include_auxiliary and str(config.get('absorption_spec', 'no')).strip().lower() == 'yes':
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
            if not run_orca("absorption_td.inp", "absorption_spec.out"):
                raise RuntimeError("ORCA terminated abnormally for absorption_spec.out")
            logger.info("TD-DFT absorption spectra calculation complete!")

        absorption_requires: Set[str] = {"initial.out"}
        if xtb_solvator_enabled:
            absorption_requires.add("initial.xyz")

        absorption_explicit: Set[str] = {"occupier_initial"}
        register_descriptor(JobDescriptor(
            job_id="occupier_absorption",
            description="absorption spectrum",
            work=run_absorption,
            produces={"absorption_spec.out"},
            requires=absorption_requires,
            explicit_dependencies=absorption_explicit,
        ))

    excitation_flags = str(config.get('excitation', '')).lower()
    emission_enabled = str(config.get('emission_spec', 'no')).strip().lower() == 'yes'
    additions_tddft = config.get('additions_TDDFT', '')
    xyz_initial = "initial.xyz"

    if include_auxiliary and 't' in excitation_flags and str(config.get('E_00', 'no')).strip().lower() == 'yes':
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
            if not run_orca("t1_state_opt.inp", "t1_state_opt.out"):
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
            requires={"initial.xyz", "initial.out"},
            explicit_dependencies={"occupier_initial"},
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
                if not run_orca("emission_t1.inp", "emission_t1.out"):
                    raise RuntimeError("ORCA terminated abnormally for emission_t1.out")
                logger.info("TD-DFT T1 emission spectra calculation complete!")

            register_descriptor(JobDescriptor(
                job_id="occupier_t1_emission",
                description="triplet emission spectrum",
                work=run_t1_emission,
                produces={"emission_t1.out"},
                requires={"t1_state_opt.xyz", "t1_state_opt.out"},
                explicit_dependencies={t1_job_id},
            ))

    if include_auxiliary and 's' in excitation_flags and str(config.get('E_00', 'no')).strip().lower() == 'yes':
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
                if not run_orca("s1_state_opt.inp", "s1_state_opt.out"):
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
            requires={"initial.xyz", "initial.out"},
            explicit_dependencies={"occupier_initial"},
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
                if not run_orca("emission_s1.inp", "emission_s1.out"):
                    raise RuntimeError("ORCA terminated abnormally for emission_s1.out")
                logger.info("TD-DFT S1 emission spectra calculation complete!")

            register_descriptor(JobDescriptor(
                job_id="occupier_s1_emission",
                description="singlet emission spectrum",
                work=run_s1_emission,
                produces={"emission_s1.out"},
                requires={"s1_state_opt.xyz", "s1_state_opt.out"},
                explicit_dependencies={s1_job_id},
            ))

    initial_job_enabled = calc_initial_flag == 'yes'
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    for step in oxidation_steps:
        folder = f"ox_step_{step}_OCCUPIER"
        multiplicity_step, additions_step, _ = read_occ(folder, "ox", step)
        if step == 1:
            requires: Set[str] = set()
            explicit_deps: Set[str] = set()
            if initial_job_enabled:
                requires.add("initial.out")
                explicit_deps.add("occupier_initial")
            if xtb_solvator_enabled:
                requires.add("initial.xyz")
        else:
            requires = {f"ox_step_{step - 1}.out"}
            explicit_deps = {f"occupier_ox_{step - 1}"}

        if xtb_solvator_enabled:
            primary_geom = Path("initial.xyz") if step == 1 else Path(f"ox_step_{step - 1}.xyz")
        else:
            primary_geom = Path(f"input_ox_step_{step}_OCCUPIER.xyz")

        if primary_geom.exists() or xtb_solvator_enabled:
            xyz_source_path = primary_geom
        else:
            fallback_geom = Path(f"input_ox_step_{step}_OCCUPIER.xyz")
            if fallback_geom.exists():
                if primary_geom != fallback_geom and not planning_only:
                    logger.warning(
                        "[occupier] Primary oxidation geometry %s missing; using OCCUPIER fallback %s",
                        primary_geom,
                        fallback_geom,
                    )
                xyz_source_path = fallback_geom
            else:
                if not planning_only:
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
                dyn_mult = mult
                dyn_adds = adds
                try:
                    refreshed_mult, refreshed_adds, _ = read_occ(f"ox_step_{idx}_OCCUPIER", "ox", idx)
                    if refreshed_mult:
                        dyn_mult = refreshed_mult
                    if isinstance(refreshed_adds, str):
                        dyn_adds = refreshed_adds
                except Exception:  # noqa: BLE001
                    logger.debug("[occupier] Using fallback multiplicity/additions for ox_step_%d", idx, exc_info=True)

                read_xyz_and_create_input3(
                    xyz_path,
                    inp,
                    charge_value,
                    dyn_mult,
                    solvent,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    dyn_adds,
                )
                inp_path = Path(inp)
                if not inp_path.exists():
                    raise RuntimeError(f"Failed to create {inp}")
                _update_pal_block(str(inp_path), cores)

                # Add %moinp block to reuse OCCUPIER wavefunction
                gbw_ox = Path(f"input_ox_step_{idx}_OCCUPIER.gbw")
                if gbw_ox.exists():
                    _add_moinp_block(str(inp_path), str(gbw_ox))
                    logger.info("[occupier_ox%d] Using GBW from OCCUPIER: %s", idx, gbw_ox)

                if not run_orca(inp, out):
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
            explicit_dependencies=explicit_deps,
        ))

    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    for step in reduction_steps:
        folder = f"red_step_{step}_OCCUPIER"
        multiplicity_step, additions_step, _ = read_occ(folder, "red", step)
        if step == 1:
            requires: Set[str] = set()
            explicit_deps: Set[str] = set()
            if initial_job_enabled:
                requires.add("initial.out")
                explicit_deps.add("occupier_initial")
            if xtb_solvator_enabled:
                requires.add("initial.xyz")
        else:
            requires = {f"red_step_{step - 1}.out"}
            explicit_deps = {f"occupier_red_{step - 1}"}

        if xtb_solvator_enabled:
            primary_geom = Path("initial.xyz") if step == 1 else Path(f"red_step_{step - 1}.xyz")
        else:
            primary_geom = Path(f"input_red_step_{step}_OCCUPIER.xyz")

        if primary_geom.exists() or xtb_solvator_enabled:
            xyz_source_path = primary_geom
        else:
            fallback_geom = Path(f"input_red_step_{step}_OCCUPIER.xyz")
            if fallback_geom.exists():
                if primary_geom != fallback_geom and not planning_only:
                    logger.warning(
                        "[occupier] Primary reduction geometry %s missing; using OCCUPIER fallback %s",
                        primary_geom,
                        fallback_geom,
                    )
                xyz_source_path = fallback_geom
            else:
                if not planning_only:
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
                dyn_mult = mult
                dyn_adds = adds
                try:
                    refreshed_mult, refreshed_adds, _ = read_occ(f"red_step_{idx}_OCCUPIER", "red", idx)
                    if refreshed_mult:
                        dyn_mult = refreshed_mult
                    if isinstance(refreshed_adds, str):
                        dyn_adds = refreshed_adds
                except Exception:  # noqa: BLE001
                    logger.debug("[occupier] Using fallback multiplicity/additions for red_step_%d", idx, exc_info=True)

                read_xyz_and_create_input3(
                    xyz_path,
                    inp,
                    charge_value,
                    dyn_mult,
                    solvent,
                    metals,
                    metal_basis,
                    main_basis,
                    config,
                    dyn_adds,
                )
                inp_path = Path(inp)
                if not inp_path.exists():
                    raise RuntimeError(f"Failed to create {inp}")
                _update_pal_block(str(inp_path), cores)

                # Add %moinp block to reuse OCCUPIER wavefunction
                gbw_red = Path(f"input_red_step_{idx}_OCCUPIER.gbw")
                if gbw_red.exists():
                    _add_moinp_block(str(inp_path), str(gbw_red))
                    logger.info("[occupier_red%d] Using GBW from OCCUPIER: %s", idx, gbw_red)

                if not run_orca(inp, out):
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
            explicit_dependencies=explicit_deps,
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


def log_orca_job_plan(label: str, jobs: List[WorkflowJob]) -> None:
    header = f"[{label}] ORCA job plan ({len(jobs)} jobs):"
    logger.info(header)
    print(header)

    # Build adjacency for levels and dependents
    job_map = {job.job_id: job for job in jobs}
    dependents: Dict[str, Set[str]] = {job.job_id: set() for job in jobs}
    for job in jobs:
        for dep in job.dependencies:
            if dep in dependents:
                dependents[dep].add(job.job_id)

    # topological levels
    levels: Dict[str, int] = {}

    def _level(job_id: str) -> int:
        if job_id not in job_map:
            return 0
        if job_id in levels:
            return levels[job_id]
        job = job_map[job_id]
        if not job.dependencies:
            levels[job_id] = 0
            return 0
        lvl = 1 + max((_level(dep) for dep in job.dependencies if dep in job_map), default=0)
        levels[job_id] = lvl
        return lvl

    for jid in job_map:
        _level(jid)

    for job in jobs:
        deps = ", ".join(sorted(job.dependencies)) if job.dependencies else "none"
        desc = job.description or "no description"
        level = levels.get(job.job_id, 0)
        outs = ", ".join(sorted(dependents.get(job.job_id, ()))) or "none"
        line = f"  - {job.job_id} → {desc} | level: {level} | deps: {deps} | unlocks: {outs}"
        logger.info(line)
        print(line)


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


def build_occupier_process_jobs(config: Dict[str, Any]) -> List[WorkflowJob]:
    """Build scheduler jobs for OCCUPIER process execution (initial, ox_*, red_*).

    These jobs prepare folders and run run_OCCUPIER() as scheduler-managed tasks,
    enabling dynamic core allocation across all OCCUPIER steps.

    NOTE: Initial job is ALWAYS included to ensure dependencies work correctly.
    The recalc logic inside OCCUPIER determines if it actually runs or skips.

    Args:
        config: DELFIN configuration dict

    Returns:
        List of WorkflowJob objects for scheduler execution (always includes initial)
    """
    from delfin.copy_helpers import prepare_occ_folder_only_setup
    from delfin.thread_safe_helpers import prepare_occ_folder_2_only_setup
    from delfin.occupier import run_OCCUPIER
    import os

    jobs: List[WorkflowJob] = []
    original_cwd = Path.cwd()
    occ_results: Dict[str, Dict[str, Any]] = config.setdefault('_occ_results_runtime', {})

    # Pre-calculate total number of parallel jobs to optimize core allocation
    total_cores = max(1, _parse_int(config.get('PAL'), fallback=1))
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    reduction_steps = _parse_step_list(config.get('reduction_steps'))

    # Calculate how many jobs might run in parallel
    # Level 0: initial (1 job)
    # Level 1: ox_step_1 and red_step_1 (up to 2 jobs in parallel)
    # Level 2: ox_step_2 and red_step_2 (up to 2 jobs in parallel)
    # etc.
    max_parallel_at_any_level = 1
    if oxidation_steps and reduction_steps:
        # Both ox and red can run in parallel at each level
        max_parallel_at_any_level = 2

    # Smart core allocation strategy:
    # The key insight is that OCCUPIER processes don't need ALL cores because:
    # 1. They often run sequentially due to dependencies (red_1 → red_2 → red_3)
    # 2. Post-processing ORCA jobs (initial.inp, red_step_1.inp, etc.) could run
    #    in parallel while later OCCUPIER processes are still running
    # 3. With 64 cores, we can allocate e.g. 48 cores to OCCUPIER and reserve
    #    16 cores for post-processing, achieving better overall throughput

    # Strategy: Reserve ~25-30% of cores for potential parallel post-processing
    # unless we have very few cores (< 16) or parallel ox/red workflows
    if max_parallel_at_any_level > 1 and total_cores >= 8:
        # Both ox and red: split cores between them
        cores_optimal_per_job = max(4, total_cores // max_parallel_at_any_level)
    elif total_cores >= 32:
        # Sequential OCCUPIER but enough cores to enable parallel post-processing
        # Use ~70-75% of cores for OCCUPIER, reserve rest for post-processing
        cores_optimal_per_job = max(16, int(total_cores * 0.75))
    else:
        # Too few cores to benefit from reservation - use all cores
        cores_optimal_per_job = total_cores

    logger.info(
        "[occupier_all] Core allocation strategy: %d cores total, "
        "%d cores optimal per OCCUPIER process (max_parallel=%d)",
        total_cores,
        cores_optimal_per_job,
        max_parallel_at_any_level,
    )

    def _record_occ_result(folder_name: str, folder_path: Path) -> None:
        try:
            result = read_occupier_file(
                str(folder_path),
                "OCCUPIER.txt",
                None,
                None,
                None,
                config,
                verbose=False,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("[%s] Failed to process OCCUPIER output: %s", folder_name, exc)
            return

        if not result:
            logger.warning("[%s] OCCUPIER.txt missing or invalid; keeping fallback settings", folder_name)
            return

        raw_mult, raw_adds, preferred_index, gbw_path = result
        try:
            mult_int = int(raw_mult) if raw_mult is not None else None
        except (TypeError, ValueError):
            mult_int = None

        additions_str = raw_adds.strip() if isinstance(raw_adds, str) else ""
        occ_results[folder_name] = {
            "multiplicity": mult_int,
            "additions": additions_str,
            "preferred_index": preferred_index,
            "gbw_path": str(gbw_path) if gbw_path else None,
        }

        log_suffix = f", gbw={gbw_path}" if gbw_path else ""
        logger.info(
            "[%s] Propagated preferred OCCUPIER geometry (index=%s%s)",
            folder_name,
            preferred_index,
            log_suffix,
        )

    def make_occupier_job(job_id: str, folder_name: str, charge_delta: int,
                         source_folder: Optional[str] = None,
                         dependencies: Optional[Set[str]] = None) -> WorkflowJob:
        """Create a WorkflowJob that prepares folder and runs OCCUPIER."""

        def work(cores: int) -> None:
            # Prepare folder (without running OCCUPIER)
            if source_folder is None:
                # Initial OCCUPIER
                folder_path = prepare_occ_folder_only_setup(folder_name, charge_delta, parent_dir=original_cwd)
            else:
                # ox/red step
                folder_path = prepare_occ_folder_2_only_setup(
                    folder_name, source_folder, charge_delta, config, original_cwd
                )

            if folder_path is None:
                raise RuntimeError(f"Failed to prepare folder {folder_name}")

            # Set environment variables for global core management
            global_cfg = {
                'PAL': cores,
                'maxcore': int(config.get('maxcore', 1000) or 1000),
            }
            pal_jobs_raw = config.get('pal_jobs')
            if pal_jobs_raw not in (None, ''):
                try:
                    global_cfg['pal_jobs'] = int(pal_jobs_raw)
                except Exception:  # noqa: BLE001
                    pass

            import json
            import sys
            import subprocess

            # Run OCCUPIER as subprocess to avoid CWD race conditions
            # Each OCCUPIER runs in its own process with its own CWD
            child_env = os.environ.copy()
            child_env['DELFIN_CHILD_GLOBAL_MANAGER'] = json.dumps(global_cfg)
            child_env['DELFIN_SUBPROCESS'] = '1'  # Flag to indicate subprocess mode

            cmd = [
                sys.executable,
                "-c",
                (
                    "from delfin.common.logging import configure_logging; "
                    "configure_logging(); "
                    "from delfin.global_manager import bootstrap_global_manager_from_env; "
                    "bootstrap_global_manager_from_env(); "
                    "import delfin.occupier as _occ; _occ.run_OCCUPIER()"
                ),
            ]

            log_prefix = f"[{folder_name}]"
            separator = "-" * (len(log_prefix) + 18)
            print(separator)
            print(f"{log_prefix} OCCUPIER start")
            print(separator)

            try:
                result = subprocess.run(
                    cmd,
                    cwd=folder_path,
                    check=False,
                    capture_output=True,
                    text=True,
                    env=child_env,
                )
            except Exception as e:
                raise RuntimeError(f"Failed to launch OCCUPIER in {folder_path}: {e}")

            # Print output
            if result.stdout:
                for line in result.stdout.splitlines():
                    print(f"{log_prefix} {line}")
            if result.stderr:
                for line in result.stderr.splitlines():
                    print(f"{log_prefix} [stderr] {line}")

            # Save logs
            try:
                (folder_path / "occupier_stdout.log").write_text(result.stdout or "")
                (folder_path / "occupier_stderr.log").write_text(result.stderr or "")
            except Exception:
                pass

            if result.returncode != 0:
                print(f"{log_prefix} OCCUPIER failed (exit={result.returncode})")
                print(separator)
                raise RuntimeError(f"OCCUPIER process in {folder_path} exited with code {result.returncode}")

            print(f"{log_prefix} OCCUPIER completed")
            print(separator)
            print()

            _record_occ_result(folder_name, folder_path)

        # Core bounds - use the pre-calculated cores_optimal_per_job
        cores_min = 1 if total_cores == 1 else 2
        cores_max = total_cores

        return WorkflowJob(
            job_id=job_id,
            work=work,
            description=f"OCCUPIER process for {folder_name}",
            dependencies=dependencies or set(),
            cores_min=cores_min,
            cores_optimal=cores_optimal_per_job,
            cores_max=cores_max,
        )

    # Always include initial OCCUPIER job (recalc logic in OCCUPIER decides if it runs)
    # This ensures dependencies work correctly even when initial already exists
    jobs.append(make_occupier_job(
        job_id="occ_proc_initial",
        folder_name="initial_OCCUPIER",
        charge_delta=0,
        source_folder=None,
        dependencies=set(),
    ))

    # Oxidation steps
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    for step in oxidation_steps:
        folder_name = f"ox_step_{step}_OCCUPIER"
        source_folder = "initial_OCCUPIER" if step == 1 else f"ox_step_{step-1}_OCCUPIER"
        deps = {"occ_proc_initial"} if step == 1 else {f"occ_proc_ox_{step-1}"}

        jobs.append(make_occupier_job(
            job_id=f"occ_proc_ox_{step}",
            folder_name=folder_name,
            charge_delta=step,
            source_folder=source_folder,
            dependencies=deps,
        ))

    # Reduction steps
    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    for step in reduction_steps:
        folder_name = f"red_step_{step}_OCCUPIER"
        source_folder = "initial_OCCUPIER" if step == 1 else f"red_step_{step-1}_OCCUPIER"
        deps = {"occ_proc_initial"} if step == 1 else {f"occ_proc_red_{step-1}"}

        jobs.append(make_occupier_job(
            job_id=f"occ_proc_red_{step}",
            folder_name=folder_name,
            charge_delta=-step,
            source_folder=source_folder,
            dependencies=deps,
        ))

    logger.info(
        "Built %d OCCUPIER process jobs (initial + ox=%d + red=%d)",
        len(jobs),
        len(oxidation_steps),
        len(reduction_steps),
    )

    return jobs


def build_combined_occupier_and_postprocessing_jobs(config: Dict[str, Any]) -> List[WorkflowJob]:
    """Build BOTH OCCUPIER process jobs AND post-processing ORCA jobs in one scheduler.

    This enables true parallelization: while red_step_3_OCCUPIER runs, initial.inp
    post-processing can run in parallel, maximizing core utilization.

    The dependency structure:
    - occ_proc_initial → occupier_initial (post-processing)
    - occ_proc_red_1 (depends on occ_proc_initial) → occupier_red_1
    - occ_proc_red_2 (depends on occ_proc_red_1) → occupier_red_2
    - etc.

    This way, while occ_proc_red_2 runs, occupier_initial can run in parallel.

    Args:
        config: DELFIN configuration dict

    Returns:
        Combined list of OCCUPIER process + post-processing jobs
    """
    # Reset staged scheduling helpers each run
    config['_occ_post_planned'] = set()
    config.pop('_post_attach_callback', None)

    # First, build OCCUPIER process jobs
    occupier_process_jobs = build_occupier_process_jobs(config)

    # Check if frequency calculation is done within OCCUPIER
    # If yes, skip post-processing ORCA jobs (they're already done inside OCCUPIER)
    frequency_mode = str(config.get('frequency_calculation_OCCUPIER', 'no')).lower()
    if frequency_mode == 'yes':
        logger.info(
            "[combined] frequency_calculation_OCCUPIER=yes → post-processing is done "
            "within OCCUPIER processes; returning OCCUPIER jobs only"
        )
        return occupier_process_jobs

    # Build post-processing ORCA jobs with dependencies on OCCUPIER processes
    # We need to create an OccupierExecutionContext (will be filled during execution)
    solvent = config.get('solvent', '')
    metals = config.get('metals', [])
    main_basisset = config.get('main_basisset', 'def2-SVP')
    metal_basisset = config.get('metal_basisset', 'def2-TZVP')
    charge = int(config.get('charge', 0))

    context = OccupierExecutionContext(
        charge=charge,
        solvent=solvent,
        metals=metals,
        main_basisset=main_basisset,
        metal_basisset=metal_basisset,
        config=config,
    )

    # Build post-processing jobs (but don't execute them yet)
    postprocessing_jobs = None
    try:
        postprocessing_jobs = build_occupier_jobs(context, planning_only=True, include_auxiliary=True)
    except RuntimeError as exc:
        exc_text = str(exc)
        missing_initial = "initial_OCCUPIER" in exc_text and "missing OCCUPIER.txt" in exc_text
        if not missing_initial:
            logger.warning(
                "[combined] Could not build post-processing jobs: %s; "
                "falling back to OCCUPIER-only execution",
                exc,
                exc_info=True,
            )
            return occupier_process_jobs

    # If planning failed because initial OCCUPIER hasn't produced results yet,
    # we'll attach a callback to build the post-processing jobs after occ_proc_initial completes.
    if postprocessing_jobs is None:
        logger.info("[combined] Deferring post-processing job generation until occ_proc_initial completes.")

        planned_stages: Set[str] = config.setdefault('_occ_post_planned', set())

        def _schedule_subset(manager: _WorkflowManager, stage: str, *, step: Optional[int] = None) -> None:
            key = f"{stage}:{step}" if step is not None else stage
            if key in planned_stages:
                return

            logger.info("[combined] Scheduling post-processing stage %s", key)

            staged_config = dict(config)

            if stage == "initial":
                staged_config['oxidation_steps'] = ''
                staged_config['reduction_steps'] = ''
            elif stage == "ox":
                staged_config['oxidation_steps'] = str(step)
                staged_config['reduction_steps'] = ''
                staged_config['calc_initial'] = 'no'
            elif stage == "red":
                staged_config['oxidation_steps'] = ''
                staged_config['reduction_steps'] = str(step)
                staged_config['calc_initial'] = 'no'
            else:
                return

            staged_context = OccupierExecutionContext(
                charge=context.charge,
                solvent=context.solvent,
                metals=context.metals,
                main_basisset=context.main_basisset,
                metal_basisset=context.metal_basisset,
                config=staged_config,
            )

            try:
                new_jobs = build_occupier_jobs(
                    staged_context,
                    planning_only=True,
                    include_auxiliary=True,
                )
            except Exception as build_exc:  # noqa: BLE001
                logger.warning(
                    "[combined] Could not build post-processing jobs for stage %s: %s",
                    key,
                    build_exc,
                    exc_info=True,
                )
                return

            def _is_stage_job(job_id: str) -> bool:
                if not job_id.startswith("occupier_"):
                    return False
                if stage == "initial":
                    return not (
                        job_id.startswith("occupier_ox_")
                        or job_id.startswith("occupier_red_")
                    )
                if stage == "ox" and step is not None:
                    return job_id == f"occupier_ox_{step}"
                if stage == "red" and step is not None:
                    return job_id == f"occupier_red_{step}"
                return False

            filtered_jobs = [job for job in new_jobs if _is_stage_job(job.job_id)]

            candidate_jobs = [job.job_id for job in filtered_jobs]
            logger.info("[combined] Candidate jobs for stage %s: %s", key, candidate_jobs or "<none>")
            added_any = False
            for job in filtered_jobs:
                try:
                    logger.info(
                        "[combined] Attempting to register job %s (deps=%s)",
                        job.job_id,
                        sorted(job.dependencies),
                    )
                    manager.add_job(job)
                    logger.info(
                        "[combined] Registered post-processing job %s for stage %s",
                        job.job_id,
                        key,
                    )
                    added_any = True
                except ValueError:
                    logger.info("[combined] Skip duplicate job %s for stage %s", job.job_id, key)
                    continue
                except Exception as register_exc:  # noqa: BLE001
                    logger.warning(
                        "[combined] Failed to register job %s for stage %s: %s",
                        job.job_id,
                        key,
                        register_exc,
                        exc_info=True,
                    )
                    continue
            if added_any:
                logger.info("[combined] Enqueued post-processing jobs for stage %s", key)
                manager.reschedule_pending()
                planned_stages.add(key)
            else:
                logger.warning("[combined] No post-processing jobs enqueued for stage %s", key)

        def _attach_postprocessing(manager: _WorkflowManager) -> None:
            if not manager:
                return

            def on_occ_initial_complete(job_id: str) -> None:
                if job_id != "occ_proc_initial":
                    return
                manager.unregister_completion_listener(on_occ_initial_complete)
                _schedule_subset(manager, "initial")

                def on_followup_complete(dep_job_id: str) -> None:
                    if dep_job_id == "occ_proc_initial":
                        return
                    if dep_job_id.startswith("occ_proc_ox_"):
                        try:
                            step_val = int(dep_job_id.rsplit('_', 1)[-1])
                        except ValueError:
                            return
                        _schedule_subset(manager, "ox", step=step_val)
                    elif dep_job_id.startswith("occ_proc_red_"):
                        try:
                            step_val = int(dep_job_id.rsplit('_', 1)[-1])
                        except ValueError:
                            return
                        _schedule_subset(manager, "red", step=step_val)

                manager.register_completion_listener(on_followup_complete)

            manager.register_completion_listener(on_occ_initial_complete)

        config['_post_attach_callback'] = _attach_postprocessing
        return occupier_process_jobs

    # Build a mapping of what each OCCUPIER process "produces" (in terms of files)
    # This allows post-processing jobs to find their dependencies
    occupier_produces = {
        "occ_proc_initial": {"input_initial_OCCUPIER.xyz", "input_initial_OCCUPIER.gbw", "initial.xyz"},
        "occupier_initial": {"initial.out", "initial.xyz"},  # Post-processing outputs
        "occupier_absorption": {"absorption_spec.out"},
        "occupier_t1_state": {"t1_state_opt.xyz", "t1_state_opt.out"},
        "occupier_t1_emission": {"emission_t1.out"},
        "occupier_s1_state": {"s1_state_opt.xyz", "s1_state_opt.out"},
        "occupier_s1_emission": {"emission_s1.out"},
    }

    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    reduction_steps = _parse_step_list(config.get('reduction_steps'))

    for step in oxidation_steps:
        occupier_produces[f"occ_proc_ox_{step}"] = {
            f"input_ox_step_{step}_OCCUPIER.xyz",
            f"input_ox_step_{step}_OCCUPIER.gbw",
            f"ox_step_{step}.xyz",
        }
        occupier_produces[f"occupier_ox_{step}"] = {f"ox_step_{step}.out", f"ox_step_{step}.xyz"}

    for step in reduction_steps:
        occupier_produces[f"occ_proc_red_{step}"] = {
            f"input_red_step_{step}_OCCUPIER.xyz",
            f"input_red_step_{step}_OCCUPIER.gbw",
            f"red_step_{step}.xyz",
        }
        occupier_produces[f"occupier_red_{step}"] = {f"red_step_{step}.out", f"red_step_{step}.xyz"}

    # Build reverse mapping: file → job that produces it
    produced_by: Dict[str, str] = {}
    for job_id, files in occupier_produces.items():
        for file in files:
            produced_by.setdefault(file, job_id)

    # Also add post-processing job products to the mapping
    for job in postprocessing_jobs:
        # Post-processing jobs already have their produces set
        # We need to extract them from the job description
        # Since WorkflowJob doesn't have a produces field, we'll infer from job_id
        pass

    # Update dependencies for post-processing jobs based on file requirements
    for job in postprocessing_jobs:
        # The job already has dependencies based on file requirements
        # We need to map those file requirements to OCCUPIER process dependencies

        # Create a new set of dependencies
        new_deps = set(job.dependencies)

        # Check if this is a post-processing job that needs its OCCUPIER counterpart
        if job.job_id.startswith("occupier_"):
            # Map to corresponding OCCUPIER process
            if job.job_id == "occupier_initial":
                new_deps.add("occ_proc_initial")
            elif job.job_id == "occupier_absorption":
                new_deps.add("occ_proc_initial")
            elif job.job_id == "occupier_t1_state":
                new_deps.add("occ_proc_initial")
            elif job.job_id == "occupier_t1_emission":
                new_deps.add("occ_proc_initial")
            elif job.job_id == "occupier_s1_state":
                new_deps.add("occ_proc_initial")
            elif job.job_id == "occupier_s1_emission":
                new_deps.add("occ_proc_initial")
            elif job.job_id.startswith("occupier_ox_"):
                step = job.job_id.replace("occupier_ox_", "")
                new_deps.add(f"occ_proc_ox_{step}")
            elif job.job_id.startswith("occupier_red_"):
                step = job.job_id.replace("occupier_red_", "")
                new_deps.add(f"occ_proc_red_{step}")

        # Update the job's dependencies
        job.dependencies = new_deps
        logger.debug(
            "[combined] Updated dependencies for %s: %s",
            job.job_id,
            sorted(new_deps),
        )

    # Combine both lists
    combined_jobs = occupier_process_jobs + postprocessing_jobs

    logger.info(
        "[combined] Built %d total jobs (%d OCCUPIER processes + %d post-processing)",
        len(combined_jobs),
        len(occupier_process_jobs),
        len(postprocessing_jobs),
    )

    return combined_jobs
