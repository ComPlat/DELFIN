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
        if folder_path.is_dir() and report_path.is_file():
            result = read_occupier_file(folder, "OCCUPIER.txt", None, None, None, config)
        else:
            result = None
        if result:
            multiplicity, additions, min_fspe_index, _gbw_path = result
            try:
                multiplicity_int = int(multiplicity)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                logger.debug("[occupier] OCCUPIER multiplicity invalid for %s; using CONTROL fallback.", folder)
                return read_occ_from_control(step_type, step)
            additions_str = additions.strip() if isinstance(additions, str) else ""
            return multiplicity_int, additions_str, min_fspe_index
        return read_occ_from_control(step_type, step)

    solvent = context.solvent
    metals = context.metals
    metal_basis = context.metal_basisset
    main_basis = context.main_basisset
    base_charge = context.charge
    functional = config.get('functional', 'ORCA')

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

            # Add %moinp block to reuse OCCUPIER wavefunction
            gbw_initial = Path("input_initial_OCCUPIER.gbw")
            if gbw_initial.exists():
                _add_moinp_block("initial.inp", str(gbw_initial))
                logger.info("[occupier_initial] Using GBW from OCCUPIER: %s", gbw_initial)

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
        multiplicity_step, additions_step, _ = read_occ(folder, "ox", step)
        if step == 1:
            requires: Set[str] = set()
            if xtb_solvator_enabled:
                requires.add("initial.xyz")
        else:
            requires = {f"ox_step_{step - 1}.out"}

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

                # Add %moinp block to reuse OCCUPIER wavefunction
                gbw_ox = Path(f"input_ox_step_{idx}_OCCUPIER.gbw")
                if gbw_ox.exists():
                    _add_moinp_block(str(inp_path), str(gbw_ox))
                    logger.info("[occupier_ox%d] Using GBW from OCCUPIER: %s", idx, gbw_ox)

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
        multiplicity_step, additions_step, _ = read_occ(folder, "red", step)
        if step == 1:
            requires: Set[str] = set()
            if xtb_solvator_enabled:
                requires.add("initial.xyz")
        else:
            requires = {f"red_step_{step - 1}.out"}

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

                # Add %moinp block to reuse OCCUPIER wavefunction
                gbw_red = Path(f"input_red_step_{idx}_OCCUPIER.gbw")
                if gbw_red.exists():
                    _add_moinp_block(str(inp_path), str(gbw_red))
                    logger.info("[occupier_red%d] Using GBW from OCCUPIER: %s", idx, gbw_red)

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

        # Core bounds
        total_cores = max(1, _parse_int(config.get('PAL'), fallback=1))
        cores_min = 1 if total_cores == 1 else 2
        cores_optimal = total_cores  # Start with full allocation
        cores_max = total_cores

        return WorkflowJob(
            job_id=job_id,
            work=work,
            description=f"OCCUPIER process for {folder_name}",
            dependencies=dependencies or set(),
            cores_min=cores_min,
            cores_optimal=cores_optimal,
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
