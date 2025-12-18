"""Integration of dynamic pool with OCCUPIER workflow."""

import os
import re
import shutil
import threading
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Set

from delfin.common.logging import get_logger
from delfin.copy_helpers import read_occupier_file, copy_preferred_files_with_names
from delfin.occupier_sequences import resolve_sequences_for_delta
from delfin.imag import run_IMAG
from delfin.orca import run_orca
from delfin.occupier_flat_extraction import _cwd_lock
from delfin.xyz_io import read_xyz_and_create_input3
from .parallel_classic_manually import (
    WorkflowJob,
    _WorkflowManager,
    _parse_int,
    _update_pal_block,
    _add_moinp_block,
    estimate_parallel_width,
    determine_effective_slots,
    normalize_parallel_token,
)

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


@dataclass
class _StageSpec:
    folder_name: str
    charge_delta: int
    source_folder: Optional[str]
    stage_prefix: str
    ensure_setup: Callable[[], Path]
    folder_path: Path
    stage_charge: int
    depends_on: Set[str] = field(default_factory=set)


def _resolve_workspace_root(base: Optional[Path] = None) -> Path:
    """Return project root even if current path points inside *_OCCUPIER folders."""
    current = Path(base) if base else Path.cwd()
    while current.name.endswith("_OCCUPIER"):
        current = current.parent
    return current


class _AutoStageController:
    """Feeds OCCUPIER FoB jobs to the scheduler stage-by-stage in auto mode."""

    def __init__(
        self,
        *,
        stage_specs: List[_StageSpec],
        build_stage_fn: Callable[[_StageSpec], tuple[List[WorkflowJob], Optional[str]]],
        post_job_builder: Callable[[], List[WorkflowJob]],
        logger_label: str = "occupier_flat",
    ) -> None:
        self._specs = {spec.stage_prefix: spec for spec in stage_specs}
        self._deps = {spec.stage_prefix: set(spec.depends_on) for spec in stage_specs}
        self._started: Set[str] = set()
        self._completed: Set[str] = set()
        self._pending_callbacks: Dict[str, List[str]] = {}
        self._build_stage = build_stage_fn
        self._build_post_jobs = post_job_builder
        self._manager: Optional[_WorkflowManager] = None
        self._lock = threading.RLock()
        self._best_by_prefix: Dict[str, Optional[str]] = {}
        self._best_by_job_id: Dict[str, str] = {}
        self._post_jobs_enqueued = False
        self._active = True
        self._label = logger_label

    def attach(self, manager: _WorkflowManager) -> None:
        with self._lock:
            self._manager = manager
            manager.register_completion_listener(self._handle_completion)

    def bootstrap(self) -> List[WorkflowJob]:
        with self._lock:
            return self._start_ready_specs(via_manager=False)

    def _handle_completion(self, job_id: str) -> None:
        with self._lock:
            if not self._active:
                return
            prefix = self._best_by_job_id.pop(job_id, None)
            if not prefix:
                return
            self._completed.add(prefix)
            self._start_ready_specs(via_manager=True)
            if len(self._completed) == len(self._specs):
                self._enqueue_post_jobs()

    def _start_ready_specs(self, *, via_manager: bool) -> List[WorkflowJob]:
        ready_jobs: List[WorkflowJob] = []
        for prefix, spec in self._specs.items():
            if prefix in self._started:
                continue
            deps = self._deps.get(prefix, set())
            if any(dep not in self._completed for dep in deps):
                continue
            jobs, best_id = self._build_stage(spec)
            self._started.add(prefix)
            self._best_by_prefix[prefix] = best_id
            if best_id:
                self._best_by_job_id[best_id] = prefix
            else:
                self._completed.add(prefix)

            if jobs:
                if via_manager and self._manager is not None:
                    logger.info(
                        "[%s] Auto OCCUPIER: queued stage %s",
                        self._label,
                        spec.folder_name,
                    )
                    for job in jobs:
                        self._manager.add_job(job)
                else:
                    ready_jobs.extend(jobs)
        if via_manager and not ready_jobs and self._manager and len(self._started) == len(self._specs):
            self._enqueue_post_jobs()
        return ready_jobs

    def _enqueue_post_jobs(self) -> None:
        if self._post_jobs_enqueued:
            return
        self._post_jobs_enqueued = True
        post_jobs = self._build_post_jobs()
        if post_jobs and self._manager:
            logger.info(
                "[%s] Auto OCCUPIER: adding %d post-processing jobs",
                self._label,
                len(post_jobs),
            )
            for job in post_jobs:
                self._manager.add_job(job)
        if self._manager:
            self._manager.unregister_completion_listener(self._handle_completion)
        self._active = False


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
    workspace_root = _resolve_workspace_root()
    base_charge = _parse_int(context.charge, fallback=0)

    neutral_electrons = config.get("_neutral_electrons")
    try:
        neutral_electrons = int(neutral_electrons)
    except (TypeError, ValueError):
        neutral_electrons = None

    multiplicity_guess = config.get("_multiplicity_guess")
    try:
        multiplicity_guess = int(multiplicity_guess)
    except (TypeError, ValueError):
        multiplicity_guess = None

    if neutral_electrons is None or multiplicity_guess is None:
        try:
            from delfin.utils import calculate_total_electrons_txt

            control_path = workspace_root / "CONTROL.txt"
            electron_info = calculate_total_electrons_txt(str(control_path))
        except Exception:  # noqa: BLE001
            electron_info = None

        if electron_info:
            if neutral_electrons is None and electron_info[0] is not None:
                neutral_electrons = int(electron_info[0])
                config["_neutral_electrons"] = neutral_electrons
            if multiplicity_guess is None and electron_info[1] is not None:
                multiplicity_guess = int(electron_info[1])
                config["_multiplicity_guess"] = multiplicity_guess
            logger.debug(
                "[occupier] Derived electron info from CONTROL: neutral=%s, guess=%s",
                electron_info[0],
                electron_info[1],
            )
        else:
            logger.debug(
                "[occupier] Could not derive electron info from CONTROL at %s",
                control_path,
            )

    @contextmanager
    def _use_workspace() -> Iterator[None]:
        prev_cwd = os.getcwd()
        try:
            os.chdir(workspace_root)
            yield
        finally:
            try:
                os.chdir(prev_cwd)
            except Exception:  # noqa: BLE001
                logger.debug(
                    "[occupier] Failed to restore working directory to %s",
                    prev_cwd,
                    exc_info=True,
                )

    def _wrap_with_workspace(func: Callable[[int], None]) -> Callable[[int], None]:
        def _runner(cores: int) -> None:
            prev_cwd = os.getcwd()
            try:
                os.chdir(workspace_root)
                func(cores)
            finally:
                try:
                    os.chdir(prev_cwd)
                except Exception:  # noqa: BLE001
                    logger.debug(
                        "[occupier] Failed to restore working directory to %s",
                        prev_cwd,
                        exc_info=True,
                    )
        return _runner
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

    def _stage_delta(step_type: str, step: Optional[int]) -> int:
        if step_type == "ox" and step:
            return step
        if step_type == "red" and step:
            return -step
        return 0

    def _parity_aligned_default(delta: int) -> int:
        if neutral_electrons is not None:
            actual = neutral_electrons - (base_charge + delta)
            return 1 if (actual % 2 == 0) else 2
        if multiplicity_guess is not None:
            base_even = (multiplicity_guess % 2 == 1)
            stage_even = base_even if (delta % 2 == 0) else not base_even
            return 1 if stage_even else 2
        # Assume neutral species is even-electron as last resort
        assumed_even = (delta % 2 == 0)
        return 1 if assumed_even else 2

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
        delta = _stage_delta(step_type, step)
        derived = _parity_aligned_default(delta)
        logger.info(
            "[occupier] Fallback multiplicity m=%d for %s step=%s (delta=%+d, neutral=%s, guess=%s)",
            derived,
            step_type,
            step or 0,
            delta,
            neutral_electrons,
            multiplicity_guess,
        )
        return derived

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

        preferred_override = None
        override_map = config.get("_occ_preferred_override")
        if isinstance(override_map, dict):
            preferred_override = override_map.get(folder)
            if preferred_override:
                logger.info(
                    "[occupier] Applying Preferred Index override for %s: %s",
                    folder,
                    preferred_override,
                )

        result = read_occupier_file(
            folder,
            "OCCUPIER.txt",
            None,
            None,
            None,
            config,
            verbose=False,
            preferred_index_override=preferred_override,
        )
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
    occ_results: Dict[str, Dict[str, Any]] = config.setdefault('_occ_results_runtime', {})

    calc_initial_flag = str(config.get('calc_initial', 'yes')).strip().lower()
    xtb_solvator_enabled = str(config.get('XTB_SOLVATOR', 'no')).strip().lower() == 'yes'
    if calc_initial_flag == 'yes' or xtb_solvator_enabled:
        try:
            multiplicity_0, additions_0, _ = read_occ("initial_OCCUPIER", "initial", None)
        except RuntimeError:
            multiplicity_0, additions_0, _ = read_occ_from_control("initial", None)

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
            with _use_workspace():
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

                geom_source: Path = workspace_root / "input_initial_OCCUPIER.xyz"
                if xtb_solvator_enabled:
                    solv_geom = workspace_root / "XTB_SOLVATOR" / "XTB_SOLVATOR.solvator.xyz"
                    if solv_geom.exists():
                        geom_source = solv_geom
                        logger.info("[occupier_initial] Using solvator geometry from %s", geom_source)
                    else:
                        logger.warning(
                            "[occupier_initial] Expected solvator geometry %s missing; falling back to %s",
                            solv_geom,
                            geom_source,
                        )

                inp_initial = workspace_root / "initial.inp"
                out_initial = workspace_root / "initial.out"

                # In recalc mode, only regenerate .inp if this job is being forced
                if _should_regenerate_inp(inp_initial, out_initial, config):
                    read_xyz_and_create_input3(
                        str(geom_source),
                        str(inp_initial),
                        base_charge,
                        mult_val,
                        solvent,
                        metals,
                        metal_basis,
                        main_basis,
                        config,
                        adds_val,
                    )
                    _update_pal_block(str(inp_initial), cores)

                    # Add %moinp block to reuse OCCUPIER wavefunction
                    gbw_initial = workspace_root / "input_initial_OCCUPIER.gbw"
                    if not xtb_solvator_enabled and gbw_initial.exists():
                        _add_moinp_block(str(inp_initial), str(gbw_initial))
                        logger.info("[occupier_initial] Using GBW from OCCUPIER: %s", gbw_initial)
                    elif xtb_solvator_enabled:
                        logger.debug("[occupier_initial] Skipping OCCUPIER GBW reuse because XTB_SOLVATOR is enabled")

                if not run_orca(str(inp_initial), str(out_initial)):
                    raise RuntimeError("ORCA terminated abnormally for initial.out")
                run_IMAG(
                    str(workspace_root / "initial.out"),
                    "initial",
                    base_charge,
                    mult_val,
                    solvent,
                    metals,
                    config,
                    main_basis,
                    metal_basis,
                    adds_val,
                    source_input=str(workspace_root / "initial.inp"),
                )
                logger.info(
                    "%s %s freq & geometry optimization of the initial system complete!",
                    functional,
                    main_basis,
                )
                initial_xyz = workspace_root / "initial.xyz"
                if not initial_xyz.exists():
                    if geom_source.exists():
                        shutil.copy(geom_source, initial_xyz)
                    else:
                        logger.warning("initial.xyz missing and no backup geometry found")

        register_descriptor(JobDescriptor(
            job_id="occupier_initial",
            description="initial OCCUPIER frequency job",
            work=run_initial,
            produces={"initial.out", "initial.xyz"},
            preferred_cores=None,
        ))

    excitation_flags = str(config.get('excitation', '')).lower()
    # Note: Excited state calculations (E_00, S1, T1) are now handled by the ESD module
    # Set ESD_modul=yes and states=[S0,S1,T1] in CONTROL.txt to calculate excited states

    initial_job_enabled = calc_initial_flag == 'yes'
    oxidation_steps = _parse_step_list(config.get('oxidation_steps'))
    for step in oxidation_steps:
        folder = f"ox_step_{step}_OCCUPIER"
        try:
            multiplicity_step, additions_step, _ = read_occ(folder, "ox", step)
        except RuntimeError:
            multiplicity_step, additions_step, _ = read_occ_from_control("ox", step)
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

        inp_path = f"ox_step_{step}.inp"
        out_path = f"ox_step_{step}.out"
        step_charge = base_charge + step

        def make_oxidation_work(idx: int, mult: int, adds: str,
                                inp: str, out: str,
                                charge_value: int,
                                use_solvator: bool) -> Callable[[int], None]:
            # Use workspace_root instead of os.getcwd() to avoid race conditions
            root_dir = Path(workspace_root).resolve()

            def _work(cores: int) -> None:
                # Ensure we're in the root directory for this job
                import os
                original_cwd = None
                with _cwd_lock:
                    original_cwd = os.getcwd()
                    try:
                        os.chdir(root_dir)
                    except Exception as e:  # noqa: BLE001
                        logger.warning("[occupier_ox%d] Could not change to root_dir %s: %s", idx, root_dir, e)

                try:
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

                    # Determine geometry paths using absolute paths from root_dir
                    if use_solvator:
                        primary_geom = root_dir / ("initial.xyz" if idx == 1 else f"ox_step_{idx - 1}.xyz")
                        fallback_geom = root_dir / f"input_ox_step_{idx}_OCCUPIER.xyz"
                    else:
                        occ_geom = root_dir / f"input_ox_step_{idx}_OCCUPIER.xyz"
                        run_geom = root_dir / f"ox_step_{idx}.xyz"
                        primary_geom = run_geom if run_geom.exists() else occ_geom
                        fallback_geom = occ_geom

                    geom_path = primary_geom if primary_geom.exists() else fallback_geom
                    if not geom_path.exists() and _recalc_enabled():
                        restored = _restore_occ_geometry("ox", idx)
                        if restored:
                            geom_path = restored
                    if not geom_path.exists():
                        import os
                        cwd = os.getcwd()
                        abs_primary = primary_geom.resolve() if primary_geom else None
                        abs_fallback = fallback_geom.resolve() if fallback_geom else None
                        logger.error(
                            "[occupier] Geometry not found! CWD=%s, primary=%s (abs=%s, exists=%s), fallback=%s (abs=%s, exists=%s)",
                            cwd, primary_geom, abs_primary, primary_geom.exists() if primary_geom else False,
                            fallback_geom, abs_fallback, fallback_geom.exists() if fallback_geom else False
                        )
                        raise FileNotFoundError(
                            f"Geometry for oxidation step {idx} not found "
                            f"(checked {primary_geom} and {fallback_geom})"
                        )

                    if primary_geom != fallback_geom and primary_geom != geom_path:
                        logger.warning(
                            "[occupier] Primary oxidation geometry %s missing; using fallback %s",
                            primary_geom,
                            geom_path,
                        )
    
                    # Use absolute path for inp to avoid race conditions with CWD changes
                    inp_abs = root_dir / inp
                    out_abs = root_dir / out

                    # In recalc mode, only regenerate .inp if this job is being forced
                    if _should_regenerate_inp(inp_abs, out_abs, config):
                        read_xyz_and_create_input3(
                            str(geom_path),
                            str(inp_abs),
                            charge_value,
                            dyn_mult,
                            solvent,
                            metals,
                            metal_basis,
                            main_basis,
                            config,
                            dyn_adds,
                        )
                        if not inp_abs.exists():
                            raise RuntimeError(f"Failed to create {inp_abs}")
                        _update_pal_block(str(inp_abs), cores)

                        # Add %moinp block to reuse OCCUPIER wavefunction
                        gbw_ox = root_dir / f"input_ox_step_{idx}_OCCUPIER.gbw"
                        if not xtb_solvator_enabled and gbw_ox.exists():
                            _add_moinp_block(str(inp_abs), str(gbw_ox))
                            logger.info("[occupier_ox%d] Using GBW from OCCUPIER: %s", idx, gbw_ox)
                        elif xtb_solvator_enabled and gbw_ox.exists():
                            logger.debug("[occupier_ox%d] Skipping OCCUPIER GBW reuse because XTB_SOLVATOR is enabled", idx)

                    if not run_orca(str(inp_abs), str(out_abs)):
                        raise RuntimeError(f"ORCA terminated abnormally for {out_abs}")
                    run_IMAG(
                        str(out_abs),
                        str(root_dir / f"ox_step_{idx}"),
                        step_charge,
                        multiplicity_step,
                        solvent,
                        metals,
                        config,
                        main_basis,
                        metal_basis,
                        additions_step,
                        step_name=f"ox_step_{idx}",
                        source_input=str(inp_abs),
                    )
                    logger.info(
                        "%s %s freq & geometry optimization cation (step %d) complete!",
                        functional,
                        main_basis,
                        idx,
                    )
                finally:
                    # Restore original working directory
                    if original_cwd is not None:
                        with _cwd_lock:
                            try:
                                os.chdir(original_cwd)
                            except Exception:  # noqa: BLE001
                                pass

            return _work

        register_descriptor(JobDescriptor(
            job_id=f"occupier_ox_{step}",
            description=f"oxidation step {step}",
            work=make_oxidation_work(step, multiplicity_step, additions_step, inp_path, out_path, step_charge, xtb_solvator_enabled),
            produces={out_path, f"ox_step_{step}.xyz"},
            requires=requires,
            explicit_dependencies=explicit_deps,
        ))

    reduction_steps = _parse_step_list(config.get('reduction_steps'))
    for step in reduction_steps:
        folder = f"red_step_{step}_OCCUPIER"
        try:
            multiplicity_step, additions_step, _ = read_occ(folder, "red", step)
        except RuntimeError:
            multiplicity_step, additions_step, _ = read_occ_from_control("red", step)
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

        inp_path = f"red_step_{step}.inp"
        out_path = f"red_step_{step}.out"
        step_charge = base_charge - step

        def make_reduction_work(idx: int, mult: int, adds: str,
                                 inp: str, out: str,
                                 charge_value: int,
                                 use_solvator: bool) -> Callable[[int], None]:
            # Use workspace_root instead of os.getcwd() to avoid race conditions
            root_dir = Path(workspace_root).resolve()

            def _work(cores: int) -> None:
                # Ensure we're in the root directory for this job
                import os
                original_cwd = None
                with _cwd_lock:
                    original_cwd = os.getcwd()
                    try:
                        os.chdir(root_dir)
                    except Exception as e:  # noqa: BLE001
                        logger.warning("[occupier_red%d] Could not change to root_dir %s: %s", idx, root_dir, e)

                try:
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

                # Determine geometry paths using absolute paths from root_dir
                    if use_solvator:
                        primary_geom = root_dir / ("initial.xyz" if idx == 1 else f"red_step_{idx - 1}.xyz")
                        fallback_geom = root_dir / f"input_red_step_{idx}_OCCUPIER.xyz"
                    else:
                        occ_geom = root_dir / f"input_red_step_{idx}_OCCUPIER.xyz"
                        run_geom = root_dir / f"red_step_{idx}.xyz"
                        primary_geom = run_geom if run_geom.exists() else occ_geom
                        fallback_geom = occ_geom
    
                    geom_path = primary_geom if primary_geom.exists() else fallback_geom
                    if not geom_path.exists() and _recalc_enabled():
                        restored = _restore_occ_geometry("red", idx)
                        if restored:
                            geom_path = restored
                    if not geom_path.exists():
                        import os
                        cwd = os.getcwd()
                        abs_primary = primary_geom.resolve() if primary_geom else None
                        abs_fallback = fallback_geom.resolve() if fallback_geom else None
                        logger.error(
                            "[occupier] Geometry not found! CWD=%s, primary=%s (abs=%s, exists=%s), fallback=%s (abs=%s, exists=%s)",
                            cwd, primary_geom, abs_primary, primary_geom.exists() if primary_geom else False,
                            fallback_geom, abs_fallback, fallback_geom.exists() if fallback_geom else False
                        )
                        raise FileNotFoundError(
                            f"Geometry for reduction step {idx} not found "
                            f"(checked {primary_geom} and {fallback_geom})"
                        )
    
                    if primary_geom != fallback_geom and primary_geom != geom_path:
                        logger.warning(
                            "[occupier] Primary reduction geometry %s missing; using fallback %s",
                            primary_geom,
                            geom_path,
                        )
    
                    # Use absolute path for inp to avoid race conditions with CWD changes
                    inp_abs = root_dir / inp
                    out_abs = root_dir / out

                    # In recalc mode, only regenerate .inp if this job is being forced
                    if _should_regenerate_inp(inp_abs, out_abs, config):
                        read_xyz_and_create_input3(
                            str(geom_path),
                            str(inp_abs),
                            charge_value,
                            dyn_mult,
                            solvent,
                            metals,
                            metal_basis,
                            main_basis,
                            config,
                            dyn_adds,
                        )
                        if not inp_abs.exists():
                            raise RuntimeError(f"Failed to create {inp_abs}")
                        _update_pal_block(str(inp_abs), cores)

                        # Add %moinp block to reuse OCCUPIER wavefunction
                        gbw_red = root_dir / f"input_red_step_{idx}_OCCUPIER.gbw"
                        if not xtb_solvator_enabled and gbw_red.exists():
                            _add_moinp_block(str(inp_abs), str(gbw_red))
                            logger.info("[occupier_red%d] Using GBW from OCCUPIER: %s", idx, gbw_red)
                        elif xtb_solvator_enabled and gbw_red.exists():
                            logger.debug("[occupier_red%d] Skipping OCCUPIER GBW reuse because XTB_SOLVATOR is enabled", idx)

                    if not run_orca(str(inp_abs), str(out_abs)):
                        raise RuntimeError(f"ORCA terminated abnormally for {out_abs}")
                    run_IMAG(
                        str(out_abs),
                        str(root_dir / f"red_step_{idx}"),
                        step_charge,
                        multiplicity_step,
                        solvent,
                        metals,
                        config,
                        main_basis,
                        metal_basis,
                        additions_step,
                        step_name=f"red_step_{idx}",
                        source_input=str(inp_abs),
                )
                    logger.info(
                        "%s %s freq & geometry optimization anion (step %d) complete!",
                        functional,
                        main_basis,
                        idx,
                    )
                finally:
                    # Restore original working directory
                    if original_cwd is not None:
                        with _cwd_lock:
                            try:
                                os.chdir(original_cwd)
                            except Exception:  # noqa: BLE001
                                pass

            return _work

        register_descriptor(JobDescriptor(
            job_id=f"occupier_red_{step}",
            description=f"reduction step {step}",
            work=make_reduction_work(step, multiplicity_step, additions_step, inp_path, out_path, step_charge, xtb_solvator_enabled),
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

        # Progressive core scaling for long-running jobs (>5 minutes estimated)
        # Jobs that historically take longer get boosted core allocation
        is_long_job = False
        long_job_patterns = [
            "s1_state",      # Singlet state optimizations are typically slow
            "t1_state",      # Triplet state optimizations
            "red_1_fob",     # Reduction step 1 FoBs (especially m=4)
            "red_2_fob",     # Reduction step 2 FoBs
        ]

        for pattern in long_job_patterns:
            if pattern in descriptor.job_id:
                is_long_job = True
                break

        if is_long_job and parallel_jobs_at_level <= 1:
            # If this long job runs alone, give it more cores
            cores_opt_boosted = min(cores_max_v, int(cores_opt_v * 1.5))
            if cores_opt_boosted > cores_opt_v:
                logger.info(
                    "[occupier] Progressive scaling: %s boosted from %d to %d cores (long-running job)",
                    descriptor.job_id, cores_opt_v, cores_opt_boosted
                )
                cores_opt_v = cores_opt_boosted

        jobs.append(
            WorkflowJob(
                job_id=descriptor.job_id,
                work=_wrap_with_workspace(descriptor.work),
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
    if isinstance(raw_steps, int):
        return [raw_steps] if raw_steps > 0 else []
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


def build_flat_occupier_fob_jobs(config: Dict[str, Any]) -> List[WorkflowJob]:
    """Build OCCUPIER FoB jobs for initial and red/ox stages without nested managers."""

    from delfin.copy_helpers import prepare_occ_folder_only_setup
    from delfin.thread_safe_helpers import prepare_occ_folder_2_only_setup
    from delfin.occupier_flat_extraction import _create_occupier_fob_jobs, _update_runtime_cache
    from delfin.occupier import run_OCCUPIER
    from delfin.utils import calculate_total_electrons_txt, search_transition_metals, set_main_basisset

    jobs: List[WorkflowJob] = []
    stage_completion: Dict[str, str] = {}
    original_cwd = Path.cwd()
    total_cores = max(1, _parse_int(config.get("PAL"), fallback=1))
    oxidation_steps = _parse_step_list(config.get("oxidation_steps"))
    reduction_steps = _parse_step_list(config.get("reduction_steps"))
    occ_results: Dict[str, Dict[str, Any]] = config.setdefault("_occ_results_runtime", {})

    control_path = original_cwd / "CONTROL.txt"
    electron_info = calculate_total_electrons_txt(str(control_path))
    neutral_electrons = electron_info[0] if electron_info else None
    multiplicity_guess = electron_info[1] if electron_info else None
    base_charge = _parse_int(config.get("charge"), fallback=0)

    input_entry = str(config.get("input_file", "input.txt"))
    input_path = (original_cwd / input_entry).resolve()
    try:
        metals = search_transition_metals(str(input_path))
    except Exception:
        metals = []
    main_basisset, metal_basisset = set_main_basisset(metals, config)
    solvent = str(config.get("solvent", ""))

    def compute_is_even(delta: int) -> Optional[bool]:
        if neutral_electrons is not None:
            actual = neutral_electrons - (base_charge + delta)
            return (actual % 2) == 0
        if multiplicity_guess is not None:
            base_even = (multiplicity_guess % 2 == 1)
            return base_even if (delta % 2 == 0) else not base_even
        return None

    def select_sequence(delta: int) -> tuple[str, List[Dict[str, Any]]]:
        parity = compute_is_even(delta)
        seq_bundle = resolve_sequences_for_delta(config, delta)
        candidate_keys: List[str] = []
        if parity is True:
            candidate_keys.extend(["initial_sequence", "even_seq"])
        elif parity is False:
            candidate_keys.extend(["initial_sequence", "odd_seq"])
        else:
            candidate_keys.append("initial_sequence")
        candidate_keys.extend(["even_seq", "odd_seq"])
        seen: Set[str] = set()
        for key in candidate_keys:
            if not key or key in seen:
                continue
            seen.add(key)
            if key in ("even_seq", "odd_seq"):
                seq = seq_bundle.get(key, [])
                if seq:
                    return key, seq
                continue
            seq = config.get(key) or []
            if seq:
                try:
                    return key, copy.deepcopy(seq)
                except Exception:
                    return key, list(seq)
        return "even_seq", seq_bundle.get("even_seq", [])

    fallback_cwd_lock = threading.RLock()

    def make_setup(folder_name: str, charge_delta: int, source_folder: Optional[str]) -> Callable[[], Path]:
        setup_lock = threading.Lock()
        state: Dict[str, Any] = {"done": False, "path": None, "error": None}

        def _ensure() -> Path:
            if state["done"]:
                if state["error"]:
                    raise state["error"]
                return state["path"]

            with setup_lock:
                if state["done"]:
                    if state["error"]:
                        raise state["error"]
                    return state["path"]
                try:
                    if source_folder is None:
                        path = prepare_occ_folder_only_setup(folder_name, charge_delta, parent_dir=original_cwd)
                    else:
                        path = prepare_occ_folder_2_only_setup(
                            folder_name,
                            source_folder,
                            charge_delta,
                            config,
                            original_cwd,
                        )
                except Exception as exc:  # noqa: BLE001
                    error = exc if isinstance(exc, RuntimeError) else RuntimeError(str(exc))
                    state["error"] = error
                    state["done"] = True
                    raise error

                if path is None:
                    error = RuntimeError(f"Failed to prepare folder {folder_name}")
                    state["error"] = error
                    state["done"] = True
                    raise error

                state["path"] = path
                state["done"] = True
                return path

        return _ensure

    def _make_spec(folder: str, delta: int, source: Optional[str], prefix: str) -> _StageSpec:
        ensure_setup = make_setup(folder, delta, source)
        depends_on: Set[str] = set()
        if source:
            source_prefix = source.replace("_OCCUPIER", "").replace("_step_", "_")
            if source_prefix:
                depends_on.add(source_prefix)
        return _StageSpec(
            folder_name=folder,
            charge_delta=delta,
            source_folder=source,
            stage_prefix=prefix,
            ensure_setup=ensure_setup,
            folder_path=original_cwd / folder,
            stage_charge=base_charge + delta,
            depends_on=depends_on,
        )

    stage_specs: List[_StageSpec] = [_make_spec("initial_OCCUPIER", 0, None, "initial")]
    stage_specs.extend(
        _make_spec(
            f"red_step_{step}_OCCUPIER",
            -step,
            "initial_OCCUPIER" if step == 1 else f"red_step_{step-1}_OCCUPIER",
            f"red_{step}",
        )
        for step in reduction_steps
    )
    stage_specs.extend(
        _make_spec(
            f"ox_step_{step}_OCCUPIER",
            step,
            "initial_OCCUPIER" if step == 1 else f"ox_step_{step-1}_OCCUPIER",
            f"ox_{step}",
        )
        for step in oxidation_steps
    )

    def make_fallback_job(
        folder_name: str,
        stage_prefix: str,
        ensure: Callable[[], Path],
        dependencies: Set[str],
        *,
        working_dir: Optional[Path] = None,
    ) -> WorkflowJob:
        def _work(cores: int) -> None:
            folder_dir = ensure()
            with fallback_cwd_lock:
                prev_cwd = os.getcwd()
                try:
                    os.chdir(folder_dir)
                    logger.info("[%s] Fallback OCCUPIER execution with %d cores", folder_name, cores)
                    run_OCCUPIER()
                finally:
                    os.chdir(prev_cwd)
            _update_runtime_cache(folder_name, folder_dir, config, occ_results)

        workdir = (working_dir or Path(folder_name)).resolve()

        return WorkflowJob(
            job_id=f"{stage_prefix}_fallback",
            work=_work,
            description=f"{folder_name} (fallback OCCUPIER run)",
            dependencies=dependencies,
            cores_min=max(1, min(total_cores, 2)),
            cores_optimal=max(2, min(total_cores, total_cores // 2 or 2)),
            cores_max=total_cores,
            working_dir=workdir,
        )

    def build_stage(spec: _StageSpec) -> tuple[List[WorkflowJob], Optional[str]]:
        dependencies: Set[str] = set()
        if spec.source_folder:
            source_prefix = spec.source_folder.replace("_OCCUPIER", "").replace("_step_", "_")
            source_token = stage_completion.get(source_prefix)
            if source_token:
                dependencies.add(source_token)
            else:
                logger.debug(
                    "[occupier_flat] Missing completion token for %s while preparing %s; "
                    "stage dependencies may be incomplete",
                    source_prefix,
                    spec.stage_prefix,
                )

        sequence_key, sequence = select_sequence(spec.charge_delta)

        if not sequence:
            logger.warning(
                "[occupier_flat] No sequence entries for %s – running OCCUPIER sequentially",
                spec.folder_name,
            )
            fallback_job = make_fallback_job(
                spec.folder_name,
                spec.stage_prefix,
                spec.ensure_setup,
                dependencies=set(dependencies),
                working_dir=spec.folder_path,
            )
            stage_completion[spec.stage_prefix] = fallback_job.job_id
            if spec.stage_prefix == "initial":
                stage_completion["initial"] = fallback_job.job_id
            return [fallback_job], fallback_job.job_id

        inner_config = dict(config)
        inner_config['PAL'] = total_cores

        inner_jobs, best_job_id = _create_occupier_fob_jobs(
            folder_name=spec.folder_name,
            folder_path=spec.folder_path,
            stage_prefix=spec.stage_prefix,
            sequence=sequence,
            sequence_label=sequence_key,
            total_cores=total_cores,
            global_config=inner_config,
            ensure_setup=spec.ensure_setup,
            source_folder=spec.source_folder,
            metals=metals,
            metal_basisset=metal_basisset,
            main_basisset=main_basisset,
            solvent=solvent,
            occ_results=occ_results,
            stage_charge=spec.stage_charge,
        )

        if dependencies:
            for job in inner_jobs:
                job.dependencies = set(job.dependencies) | set(dependencies)

        stage_completion[spec.stage_prefix] = best_job_id
        if spec.stage_prefix == "initial":
            stage_completion["initial"] = best_job_id
        return inner_jobs, best_job_id

    context = OccupierExecutionContext(
        charge=base_charge,
        solvent=solvent,
        metals=metals,
        main_basisset=main_basisset,
        metal_basisset=metal_basisset,
        config=config,
    )

    def build_post_processing_jobs() -> List[WorkflowJob]:
        post_jobs = build_occupier_jobs(context, planning_only=True, include_auxiliary=True)

        def _map_occ_dependency(dep: str) -> Optional[str]:
            if dep == "occ_proc_initial":
                return stage_completion.get("initial")
            match = re.match(r"occ_proc_(ox|red)_(\d+)", dep)
            if match:
                stage_type, step = match.groups()
                return stage_completion.get(f"{stage_type}_{step}")
            return dep

        filtered: List[WorkflowJob] = []
        for job in post_jobs:
            remapped: Set[str] = set()
            skip = False
            for dep in job.dependencies:
                if dep.startswith("occ_proc_"):
                    mapped = _map_occ_dependency(dep)
                    if not mapped or mapped.startswith("occ_proc_"):
                        logger.debug(
                            "[occupier_flat] Skipping post-processing job %s due to missing dependency %s",
                            job.job_id,
                            dep,
                        )
                        skip = True
                        break
                    remapped.add(mapped)
                else:
                    remapped.add(dep)
            if skip:
                continue
            if job.job_id == "occupier_initial":
                initial_token = stage_completion.get("initial")
                if initial_token:
                    remapped.add(initial_token)
            elif job.job_id.startswith("occupier_ox_"):
                step_id = job.job_id.split("_")[2]
                token = stage_completion.get(f"ox_{step_id}")
                if token:
                    remapped.add(token)
            elif job.job_id.startswith("occupier_red_"):
                step_id = job.job_id.split("_")[2]
                token = stage_completion.get(f"red_{step_id}")
                if token:
                    remapped.add(token)
            job.dependencies = remapped
            filtered.append(job)
        return filtered

    sequential_auto = str(config.get("OCCUPIER_method", "auto")).strip().lower() == "auto"
    if stage_specs and sequential_auto:
        logger.info("[occupier_flat] Auto OCCUPIER → enabling dependency-driven scheduling")
        controller = _AutoStageController(
            stage_specs=stage_specs,
            build_stage_fn=build_stage,
            post_job_builder=build_post_processing_jobs,
        )
        initial_jobs = controller.bootstrap()
        config['_post_attach_callback'] = controller.attach
        jobs.extend(initial_jobs)
    else:
        for spec in stage_specs:
            stage_jobs, _ = build_stage(spec)
            jobs.extend(stage_jobs)

        post_jobs = build_post_processing_jobs()
        if post_jobs:
            logger.info(
                "[occupier_flat] Appending %d post-processing jobs to OCCUPIER plan",
                len(post_jobs),
            )
            jobs.extend(post_jobs)

    logger.info(
        "[occupier_flat] Prepared %d OCCUPIER jobs (FoBs + post-processing) across %d stages (ox=%d, red=%d)",
        len(jobs),
        len(stage_specs),
        len(oxidation_steps),
        len(reduction_steps),
    )
    return jobs


def build_occupier_process_jobs(config: Dict[str, Any]) -> List[WorkflowJob]:
    """Build scheduler jobs for OCCUPIER process execution (initial, ox_*, red_*).

    DEPRECATED: This creates nested managers which causes deadlocks.
    Use build_flat_occupier_fob_jobs() instead!

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

    # Strategy: Use cores_optimal as a target, but allow scheduler to allocate
    # more cores (up to cores_max) when fewer jobs are running.
    # This ensures ORCA processes can use all available cores when possible,
    # while still allowing parallelism when multiple workflows are ready.
    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))

    if parallel_mode == 'disable':
        cores_optimal_per_job = total_cores
    elif max_parallel_at_any_level > 1 and total_cores >= 8:
        # Both ox and red can run in parallel: target 50% split as optimal
        # but cores_max=total_cores allows scheduler to give all cores when only 1 runs
        cores_optimal_per_job = max(6, total_cores // max_parallel_at_any_level)
    elif total_cores >= 32:
        # Sequential OCCUPIER but enough cores to enable parallel post-processing
        # Use ~70-75% of cores for OCCUPIER as optimal, reserve rest for post-processing
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

            # Run OCCUPIER in-process with temporary CWD change
            # This allows all OCCUPIER folders to share the same global pool!
            # Thread-safe CWD changes are handled by a lock.
            log_prefix = f"[{folder_name}]"
            separator = "-" * (len(log_prefix) + 18)
            logger.info("%s", separator)
            logger.info("%s OCCUPIER start", log_prefix)
            logger.info("%s", separator)

            # Thread-safe CWD change
            import threading
            _cwd_lock = getattr(build_occupier_process_jobs, '_cwd_lock', None)
            if _cwd_lock is None:
                _cwd_lock = threading.RLock()
                build_occupier_process_jobs._cwd_lock = _cwd_lock

            old_cwd = os.getcwd()
            try:
                with _cwd_lock:
                    os.chdir(folder_path)

                    # run_OCCUPIER() reads CONTROL.txt from current directory
                    # It expects no parameters!
                    # The global manager is already initialized and will handle core allocation
                    from delfin.occupier import run_OCCUPIER
                    run_OCCUPIER()

                with _cwd_lock:
                    os.chdir(old_cwd)

                logger.info("%s OCCUPIER completed", log_prefix)
                logger.info("%s", separator)
                logger.info("")

            except Exception as e:
                with _cwd_lock:
                    try:
                        os.chdir(old_cwd)
                    except:
                        pass
                logger.error("%s OCCUPIER failed: %s", log_prefix, e)
                logger.info("%s", separator)
                raise RuntimeError(f"OCCUPIER failed in {folder_path}: {e}")

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
            working_dir=folder_path.resolve(),
        )

    imag_enabled = str(config.get('IMAG', 'no')).strip().lower() == 'yes'

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
        missing_folder = None
        if "missing OCCUPIER.txt" in exc_text:
            match = re.search(r"'([^']+_OCCUPIER)'", exc_text)
            if match:
                missing_folder = match.group(1)
        if missing_folder:
            deferred_set = config.setdefault('_occ_post_missing', set())
            deferred_set.add(missing_folder)
            logger.info(
                "[combined] Post-processing for %s deferred until OCCUPIER results are regenerated.",
                missing_folder,
            )
        else:
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
        missing_targets = sorted(config.get('_occ_post_missing', set()))
        if missing_targets:
            logger.info(
                "[combined] Deferring post-processing job generation; pending OCCUPIER results for: %s",
                ", ".join(missing_targets),
            )
        else:
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

    # Guard OCCUPIER processes until the IMAG refinement they rely on has finished,
    # if requested via CONTROL settings.
    imag_enabled = str(config.get('IMAG', 'no')).strip().lower() == 'yes'
    if imag_enabled:
        try:
            imag_option = int(config.get('IMAG_option', 2) or 2)
        except (TypeError, ValueError):
            imag_option = 2
        if imag_option != 2:
            logger.debug(
                "[combined] IMAG_option=%s → skipping additional OCCUPIER dependencies",
                imag_option,
            )
        else:
            imag_scope = str(config.get('IMAG_scope', 'initial')).strip().lower()
            job_lookup = {job.job_id: job for job in combined_jobs}

            def _attach_imag_dependency(target_id: str, dependency_id: str) -> None:
                target_job = job_lookup.get(target_id)
                dependency_job = job_lookup.get(dependency_id)
                if not target_job or not dependency_job:
                    return
                if dependency_id in target_job.dependencies:
                    return
                target_job.dependencies.add(dependency_id)
                logger.debug(
                    "[combined] Added IMAG dependency: %s waits for %s",
                    target_id,
                    dependency_id,
                )

            if imag_scope in {"initial", "all"}:
                _attach_imag_dependency("occ_proc_red_1", "occupier_initial")
                _attach_imag_dependency("occ_proc_ox_1", "occupier_initial")

            if imag_scope == "all":
                for step in oxidation_steps:
                    if step <= 1:
                        continue
                    prev_job = f"occupier_ox_{step - 1}"
                    current_proc = f"occ_proc_ox_{step}"
                    _attach_imag_dependency(current_proc, prev_job)
                for step in reduction_steps:
                    if step <= 1:
                        continue
                    prev_job = f"occupier_red_{step - 1}"
                    current_proc = f"occ_proc_red_{step}"
                    _attach_imag_dependency(current_proc, prev_job)

    logger.info(
        "[combined] Built %d total jobs (%d OCCUPIER processes + %d post-processing)",
        len(combined_jobs),
        len(occupier_process_jobs),
        len(postprocessing_jobs),
    )

    return combined_jobs
def _recalc_enabled() -> bool:
    return str(os.environ.get("DELFIN_RECALC", "0")).lower() in ("1", "true", "yes", "on")


def _should_regenerate_inp(inp_path: Path, out_path: Path, config: Dict[str, Any]) -> bool:
    """Check if .inp file should be regenerated in recalc mode.

    Returns:
        True if .inp should be regenerated, False if existing .inp should be preserved.
    """
    if not _recalc_enabled():
        return True  # Always regenerate in normal mode

    # In recalc mode, only regenerate if this specific job is being forced
    if not inp_path.exists():
        return True  # Must create it

    # Check if this output is in the force list
    force_outputs = config.get("_recalc_force_outputs", set())
    if not force_outputs:
        # No force list - preserve existing .inp files
        logger.debug("[recalc] No force list; preserving existing %s", inp_path.name)
        return False

    # Check if this specific output is forced
    out_resolved = out_path.resolve() if out_path.exists() else out_path
    is_forced = any(
        out_resolved == (p.resolve() if hasattr(p, 'resolve') else Path(p).resolve())
        for p in force_outputs
    )

    if is_forced:
        logger.info("[recalc] Regenerating %s (output %s is forced)", inp_path.name, out_path.name)
        return True

    # Not forced - preserve existing .inp
    logger.debug("[recalc] Preserving existing %s (not forced)", inp_path.name)
    return False


def _restore_occ_geometry(stage: str, idx: int) -> Optional[Path]:
    """Rebuild OCCUPIER geometry/GBW for the requested stage from OCCUPIER outputs.

    Returns path to the restored fallback geometry (input_*_OCCUPIER.xyz).
    """
    stage = stage.lower()
    folder = Path(f"{stage}_step_{idx}_OCCUPIER")
    if not folder.is_dir():
        return None

    dest_dir = folder.parent
    dest_input_name = f"input_{stage}_step_{idx}_OCCUPIER.xyz"
    try:
        _, dest_input_path, preferred_index = copy_preferred_files_with_names(
            folder_name=str(folder),
            dest_output_filename=f"occ_{stage}_step_{idx}_preferred.out",
            dest_input_filename=dest_input_name,
            dest_dir=str(dest_dir),
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "[recalc] Failed to reconstruct OCCUPIER geometry for %s_step_%d: %s",
            stage,
            idx,
            exc,
        )
        return None

    # Try to copy the corresponding GBW if available (best-effort)
    try:
        pref_idx = int(preferred_index) if preferred_index is not None else None
    except (TypeError, ValueError):
        pref_idx = None

    if pref_idx:
        gbw_source_candidates = [
            folder / ("input.gbw" if pref_idx == 1 else f"input{pref_idx}.gbw"),
            folder / ("input.gbw_hs" if pref_idx == 1 else f"input{pref_idx}.gbw_hs"),
            folder / ("input.gbw_bs" if pref_idx == 1 else f"input{pref_idx}.gbw_bs"),
        ]
        dest_gbw = dest_dir / f"input_{stage}_step_{idx}_OCCUPIER.gbw"
        for candidate in gbw_source_candidates:
            if candidate.exists():
                try:
                    shutil.copy(candidate, dest_gbw)
                except Exception as gbw_exc:  # noqa: BLE001
                    logger.warning(
                        "[recalc] Failed to copy GBW %s → %s: %s",
                        candidate,
                        dest_gbw,
                        gbw_exc,
                    )
                break

    restored_path = dest_dir / dest_input_name
    if restored_path.exists():
        # Also refresh the run-level XYZ (ox_step_n.xyz / red_step_n.xyz) so primary_geom exists
        run_xyz = dest_dir / f"{stage}_step_{idx}.xyz"
        try:
            shutil.copy(restored_path, run_xyz)
        except Exception as run_exc:  # noqa: BLE001
            logger.warning(
                "[recalc] Failed to refresh %s from %s: %s",
                run_xyz,
                restored_path,
                run_exc,
            )
        logger.info("[recalc] Restored OCCUPIER geometry for %s_step_%d from %s", stage, idx, folder)
        return restored_path
    return None
