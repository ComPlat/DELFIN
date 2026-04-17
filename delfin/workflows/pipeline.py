"""Workflow execution helpers for CLI orchestration."""

from __future__ import annotations

import difflib
import json
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.copy_helpers import copy_if_exists, read_occupier_file
from delfin.workflows.scheduling.manager import get_global_manager
from delfin.workflows.engine.scheduler import GlobalOrcaScheduler
from delfin.workflows.engine.occupier import OccupierExecutionContext, run_occupier_orca_jobs
from delfin.workflows.engine.classic import execute_classic_workflows, execute_manually_workflows, normalize_parallel_token, WorkflowRunResult
from delfin.xtb_crest import XTB, XTB_GOAT, XTB_SOLVATOR, run_crest_workflow
from delfin.cli_calculations import calculate_redox_potentials, select_final_potentials
from delfin.energies import find_gibbs_energy
from delfin.esd_module import execute_esd_jobs as execute_esd_module, parse_esd_config
from delfin.esd_results import collect_esd_results, ESDSummary

logger = get_logger(__name__)


def _resolve_pipeline_job_limits(*, requested_cores: int, requested_maxcore: int) -> tuple[int, int]:
    manager = get_global_manager()
    try:
        if manager.is_initialized():
            return manager.resolve_job_resources(
                requested_cores=requested_cores,
                requested_maxcore=requested_maxcore,
            )
    except Exception:
        logger.debug("Could not resolve global job limits; using requested pipeline limits", exc_info=True)
    return max(1, int(requested_cores)), max(256, int(requested_maxcore))




@dataclass
class FileBundle:
    """Collection of frequently accessed file names for workflows."""

    xyz_initial: str = "initial.xyz"
    xyz_red1: str = "red_step_1.xyz"
    xyz_red2: str = "red_step_2.xyz"
    xyz_red3: str = "red_step_3.xyz"
    xyz_ox1: str = "ox_step_1.xyz"
    xyz_ox2: str = "ox_step_2.xyz"
    xyz_ox3: str = "ox_step_3.xyz"
    output_initial: str = "initial.inp"
    output_absorption: str = "absorption_td.inp"
    output_t1: str = "t1_state_opt.inp"
    output_s1: str = "s1_state_opt.inp"
    output_emission: str = "emission_td.inp"
    output_ox1: str = "ox_step_1.inp"
    output_ox2: str = "ox_step_2.inp"
    output_ox3: str = "ox_step_3.inp"
    output_red1: str = "red_step_1.inp"
    output_red2: str = "red_step_2.inp"
    output_red3: str = "red_step_3.inp"


@dataclass
class PipelineContext:
    """Aggregated state that downstream workflow helpers rely on."""

    config: Dict[str, Any]
    control_file_path: Path
    input_file: str
    charge: int
    PAL: int
    multiplicity: int
    solvent: str
    metals: List[str]
    main_basisset: str
    metal_basisset: str
    number_explicit_solv_molecules: int
    total_electrons_txt: int
    start_time: float
    name: str
    file_bundle: FileBundle = field(default_factory=FileBundle)
    extra: Dict[str, Any] = field(default_factory=dict)

    def clone_with(self, **updates: Any) -> "PipelineContext":
        data = {
            'config': self.config,
            'control_file_path': self.control_file_path,
            'input_file': self.input_file,
            'charge': self.charge,
            'PAL': self.PAL,
            'multiplicity': self.multiplicity,
            'solvent': self.solvent,
            'metals': self.metals,
            'main_basisset': self.main_basisset,
            'metal_basisset': self.metal_basisset,
            'number_explicit_solv_molecules': self.number_explicit_solv_molecules,
            'total_electrons_txt': self.total_electrons_txt,
            'start_time': self.start_time,
            'name': self.name,
            'file_bundle': self.file_bundle,
            'extra': dict(self.extra),
        }
        data.update(updates)
        return PipelineContext(**data)


# ---------------------------------------------------------------------------
# OCCUPIER helpers
# ---------------------------------------------------------------------------


# DEPRECATED FUNCTIONS REMOVED:
# - _execute_oxidation_workflow (replaced by build_occupier_process_jobs)
# - _execute_reduction_workflow (replaced by build_occupier_process_jobs)
# - _should_parallelize (no longer needed with scheduler-based execution)
# - _execute_parallel_workflows (replaced by inline scheduler calls)
# - _execute_sequential_workflows (replaced by inline scheduler calls)
# - _run_occ_workflows (replaced by inline scheduler calls with proper dependency handling)


def run_occuper_phase(ctx: PipelineContext) -> bool:
    """Execute OCCUPIER-specific preparation and post-processing."""

    config = ctx.config
    multiplicity = ctx.multiplicity
    charge = ctx.charge

    _override_skip = _skip_preprocessing_for_override(config, ctx.control_file_path.parent)
    if config['XTB_OPT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_OPT: --occupier-override active, upstream geometry reused.")
        else:
            XTB(multiplicity, charge, config)

    if config['XTB_GOAT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_GOAT: --occupier-override active, upstream geometry reused.")
        elif _skip_xtb_goat_after_guppy(config):
            logger.info("Skipping XTB_GOAT: GUPPY already provided GOAT-refined winner geometry.")
        else:
            XTB_GOAT(multiplicity, charge, config)

    if config['CREST'] == "yes":
        run_crest_workflow(ctx.PAL, ctx.solvent, charge, multiplicity, ctx.config.get('input_file'))

    metals_list = list(ctx.metals) if isinstance(ctx.metals, (list, tuple, set)) else ([ctx.metals] if ctx.metals else [])

    if config['XTB_SOLVATOR'] == "no":
        # Run ALL OCCUPIER jobs (initial + ox/red) + post-processing in ONE scheduler run
        if "yes" in config.get("calc_initial", ""):
            print("\nOCCUPIER for the initial system:\n")

            from delfin.workflows.engine.occupier import build_flat_occupier_fob_jobs
            from delfin.workflows.engine.classic import _WorkflowManager

            # Build ALL OCCUPIER FoBs as flat top-level jobs
            # This avoids nested managers and deadlocks!
            all_jobs = build_flat_occupier_fob_jobs(config)
            sc_embedded = False
            if str(config.get("stability_constant", "no")).strip().lower() == "yes":
                from delfin.stability_constant import (
                    build_stability_constant_plan,
                    build_stability_reaction_plan,
                )

                try:
                    sc_mode = str(config.get("stability_constant_mode", "auto")).strip().lower() or "auto"
                    if sc_mode == "reaction":
                        sc_plan = build_stability_reaction_plan(
                            ctx,
                            initial_completion_dependency=config.get("_occ_initial_energy_job"),
                        )
                    else:
                        sc_plan = build_stability_constant_plan(
                            ctx,
                            initial_completion_dependency=config.get("_occ_initial_energy_job"),
                        )
                    all_jobs.extend(sc_plan.jobs)
                    sc_embedded = True
                    logger.info(
                        "[SC] Embedded thermodynamics jobs into the OCCUPIER shared scheduler "
                        "(mode=%s, initial dependency: %s).",
                        sc_mode,
                        config.get("_occ_initial_energy_job") or "none",
                    )
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "[SC] Could not embed thermodynamics workflow into OCCUPIER run: %s",
                        exc,
                        exc_info=True,
                    )

            if all_jobs:
                manager = _WorkflowManager(config, label="occupier_all")
                try:
                    # Check if sequential execution is requested
                    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
                    if parallel_mode == 'disable':
                        logger.info("[occupier_all] parallel_workflows=no → enforcing sequential execution")
                        manager.enforce_sequential_allocation()
                        if manager.pool.max_concurrent_jobs != 1:
                            manager.pool.max_concurrent_jobs = 1
                            manager.max_jobs = 1
                            manager._sync_parallel_flag()

                    for job in all_jobs:
                        manager.add_job(job)

                    manager.run()

                    if sc_embedded and "sc_postprocess" in manager.completed_jobs:
                        ctx.extra["stability_constant_embedded_complete"] = True
                        logger.info("[SC] Embedded thermodynamics workflow completed within OCCUPIER run.")
                    elif sc_embedded:
                        logger.warning(
                            "[SC] Embedded stability constant workflow did not complete inside OCCUPIER run; "
                            "the standalone SC phase remains available as fallback."
                        )

                    if manager.failed_jobs:
                        logger.warning("OCCUPIER workflows completed with failures (continuing): %s", list(manager.failed_jobs.keys()))
                finally:
                    manager.shutdown()

            logger.info("All OCCUPIER workflows completed successfully")

            # Post-processing jobs are scheduled within the combined OCCUPIER run
            config['_used_combined_occupier'] = True
        else:
            # No initial calculation, but still run ox/red if configured
            from delfin.workflows.engine.occupier import build_flat_occupier_fob_jobs
            from delfin.workflows.engine.classic import _WorkflowManager

            # Build all jobs (OCCUPIER FoBs + post-processing) with flat architecture
            jobs = build_flat_occupier_fob_jobs(config)

            if jobs:
                manager = _WorkflowManager(config, label="occupier_workflows")
                try:
                    # Check if sequential execution is requested
                    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
                    if parallel_mode == 'disable':
                        logger.info("[occupier_workflows] parallel_workflows=no → enforcing sequential execution")
                        manager.enforce_sequential_allocation()
                        if manager.pool.max_concurrent_jobs != 1:
                            manager.pool.max_concurrent_jobs = 1
                            manager.max_jobs = 1
                            manager._sync_parallel_flag()

                    for job in jobs:
                        manager.add_job(job)

                    manager.run()

                    if manager.failed_jobs:
                        logger.warning("OCCUPIER workflows completed with failures (continuing): %s", list(manager.failed_jobs.keys()))
                finally:
                    manager.shutdown()

            # Mark that we used combined execution
            config['_used_combined_occupier'] = True

        if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == "yes":
            multiplicity_0, broken_sym_0, _, gbw_initial = read_occupier_file(
                "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
            )
            ctx.extra['multiplicity_0'] = multiplicity_0
            ctx.extra['broken_sym_0'] = broken_sym_0
            ctx.extra['ground_broken_sym'] = broken_sym_0
            ctx.extra['gbw_initial'] = gbw_initial

            copy_if_exists("./initial_OCCUPIER", "initial.out", "initial.xyz")
            for step in (1, 2, 3):
                copy_if_exists(f"./ox_step_{step}_OCCUPIER", f"ox_step_{step}.out", f"ox_step_{step}.xyz")
                copy_if_exists(f"./red_step_{step}_OCCUPIER", f"red_step_{step}.out", f"red_step_{step}.xyz")

    else:  # XTB_SOLVATOR == "yes"
        calc_initial_flag = str(config.get("calc_initial", "")).strip().lower()
        initial_requested = "yes" in calc_initial_flag
        initial_folder = Path("initial_OCCUPIER")
        initial_report = initial_folder / "OCCUPIER.txt"

        initial_rerun = False
        need_occ_workflows = False

        if initial_requested or not initial_report.exists():
            print("\nOCCUPIER for the initial system:\n")
            initial_rerun = True
            need_occ_workflows = True
        else:
            logger.info(
                "Reusing existing OCCUPIER results in %s (calc_initial=%s)",
                initial_folder,
                config.get("calc_initial"),
            )

        if config.get("oxidation_steps", "").strip() or config.get("reduction_steps", "").strip():
            def _extract_steps(raw: str) -> List[int]:
                if not raw:
                    return []
                return [int(token) for token in re.findall(r"\d+", str(raw)) if token.strip()]

            # Check if we need to run ox/red workflows
            if not need_occ_workflows:
                for step in _extract_steps(config.get("oxidation_steps", "")):
                    if not (Path(f"ox_step_{step}_OCCUPIER") / "OCCUPIER.txt").exists():
                        need_occ_workflows = True
                        break

            if not need_occ_workflows:
                for step in _extract_steps(config.get("reduction_steps", "")):
                    if not (Path(f"red_step_{step}_OCCUPIER") / "OCCUPIER.txt").exists():
                        need_occ_workflows = True
                        break

            if need_occ_workflows:
                logger.info("Running OCCUPIER workflows (initial + ox/red) prior to solvation")

                from delfin.workflows.engine.occupier import build_flat_occupier_fob_jobs
                from delfin.workflows.engine.classic import _WorkflowManager

                # Build ALL OCCUPIER FoBs as flat top-level jobs
                # This avoids nested managers and deadlocks!
                all_jobs = build_flat_occupier_fob_jobs(config)

                if all_jobs:
                    manager = _WorkflowManager(config, label="occupier_all")
                    try:
                        # Check if sequential execution is requested
                        parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
                        if parallel_mode == 'disable':
                            logger.info("[occupier_all] parallel_workflows=no → enforcing sequential execution")
                            manager.enforce_sequential_allocation()
                            if manager.pool.max_concurrent_jobs != 1:
                                manager.pool.max_concurrent_jobs = 1
                                manager.max_jobs = 1
                                manager._sync_parallel_flag()

                        for job in all_jobs:
                            manager.add_job(job)

                        manager.run()

                        if manager.failed_jobs:
                            logger.warning("OCCUPIER workflows completed with failures (continuing): %s", list(manager.failed_jobs.keys()))
                    finally:
                        manager.shutdown()

                # Mark that we used combined execution
                config['_used_combined_occupier'] = True
            else:
                logger.info("Reusing existing OCCUPIER oxidation/reduction workflows")

        multiplicity_0, broken_sym_0, _, gbw_initial = read_occupier_file(
            "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
        )
        ctx.extra['multiplicity_0'] = multiplicity_0
        ctx.extra['broken_sym_0'] = broken_sym_0
        ctx.extra['ground_broken_sym'] = broken_sym_0
        ctx.extra['gbw_initial'] = gbw_initial

        preferred_parent_xyz = Path("input_initial_OCCUPIER.xyz")
        if not preferred_parent_xyz.exists():
            logger.warning(
                "Preferred OCCUPIER geometry %s missing; falling back to start.txt for solvator run.",
                preferred_parent_xyz,
            )
            solvator_source = Path("start.txt")
        else:
            solvator_source = preferred_parent_xyz

        XTB_SOLVATOR(
            str(solvator_source.resolve()),
            multiplicity_0,
            charge,
            ctx.solvent,
            ctx.number_explicit_solv_molecules,
            config,
        )

        solvated_xyz = Path("XTB_SOLVATOR") / "XTB_SOLVATOR.solvator.xyz"
        target_parent_xyz = Path("input_initial_OCCUPIER.xyz")
        if solvated_xyz.exists():
            try:
                shutil.copyfile(solvated_xyz, target_parent_xyz)
                logger.info("Propagated solvated geometry to %s", target_parent_xyz)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "Failed to update %s with solvated coordinates: %s",
                    target_parent_xyz,
                    exc,
                )
        else:
            logger.warning(
                "XTB_SOLVATOR completed but %s is missing; OCCUPIER workflows will reuse unsolvated geometry.",
                solvated_xyz,
            )

    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    parallel_enabled = parallel_mode != 'disable'
    metals_list = list(ctx.metals) if isinstance(ctx.metals, (list, tuple, set)) else ([ctx.metals] if ctx.metals else [])

    # Check if we already ran post-processing in combined mode
    used_combined = config.get('_used_combined_occupier', False)

    if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() != "yes" and not used_combined:
        # Only run separate post-processing if we didn't use combined execution
        logger.info("[pipeline] Running separate post-processing ORCA jobs")

        occ_context = OccupierExecutionContext(
            charge=charge,
            solvent=ctx.solvent,
            metals=metals_list,
            main_basisset=ctx.main_basisset,
            metal_basisset=ctx.metal_basisset,
            config=config,
        )

        scheduler = GlobalOrcaScheduler(config, label="occupier")
        try:
            occ_success = run_occupier_orca_jobs(
                occ_context,
                parallel_enabled,
                scheduler=scheduler,
            )
        finally:
            scheduler.shutdown()
        if not occ_success:
            if occ_context.failed_jobs or occ_context.skipped_jobs:
                failed_desc = ", ".join(
                    f"{job_id} ({reason})" for job_id, reason in occ_context.failed_jobs.items()
                ) or "none"
                skipped_desc = ", ".join(
                    f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
                    for job_id, deps in occ_context.skipped_jobs.items()
                ) or "none"
                logger.warning(
                    "OCCUPIER post-processing encountered failures (continuing run). Failed jobs: %s | Skipped jobs: %s",
                    failed_desc,
                    skipped_desc,
                )
            else:
                logger.error("OCCUPIER post-processing failed; aborting run")
                return False
    elif used_combined:
        logger.info("[pipeline] Skipping separate post-processing (already done in combined mode)")

    return True


# ---------------------------------------------------------------------------
# Classic / Manually workflows
# ---------------------------------------------------------------------------


def run_classic_phase(ctx: PipelineContext) -> Dict[str, Any]:
    config = ctx.config
    multiplicity = ctx.multiplicity
    charge = ctx.charge

    _override_skip = _skip_preprocessing_for_override(config, ctx.control_file_path.parent)
    if config['XTB_OPT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_OPT: --occupier-override active, upstream geometry reused.")
        else:
            XTB(multiplicity, charge, config)

    if config['XTB_GOAT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_GOAT: --occupier-override active, upstream geometry reused.")
        elif _skip_xtb_goat_after_guppy(config):
            logger.info("Skipping XTB_GOAT: GUPPY already provided GOAT-refined winner geometry.")
        else:
            XTB_GOAT(multiplicity, charge, config)

    if config['CREST'] == "yes":
        run_crest_workflow(ctx.PAL, ctx.solvent, charge, multiplicity, ctx.config.get('input_file'))

    if config['XTB_SOLVATOR'] == "yes":
        XTB_SOLVATOR(
            ctx.config.get('input_file') or 'start.txt',
            multiplicity,
            charge,
            ctx.solvent,
            ctx.number_explicit_solv_molecules,
            config,
        )

    ground_multiplicity = multiplicity
    classic_kwargs = {
        'total_electrons_txt': ctx.total_electrons_txt,
        'xyz_file': ctx.file_bundle.xyz_initial,
        'xyz_file2': ctx.file_bundle.xyz_red1,
        'xyz_file3': ctx.file_bundle.xyz_red2,
        'xyz_file4': ctx.file_bundle.xyz_ox1,
        'xyz_file8': ctx.file_bundle.xyz_ox2,
        'output_file5': ctx.file_bundle.output_ox1,
        'output_file9': ctx.file_bundle.output_ox2,
        'output_file10': ctx.file_bundle.output_ox3,
        'output_file6': ctx.file_bundle.output_red1,
        'output_file7': ctx.file_bundle.output_red2,
        'output_file8': ctx.file_bundle.output_red3,
        'solvent': ctx.solvent,
        'metals': ctx.metals,
        'metal_basisset': ctx.metal_basisset,
        'main_basisset': ctx.main_basisset,
        'broken_sym': "",
        'input_file_path': ctx.input_file,
        'output_initial': ctx.file_bundle.output_initial,
        'ground_multiplicity': ground_multiplicity,
        'include_excited_jobs': True,
    }

    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    allow_parallel = parallel_mode != 'disable'
    mode_label = "parallel" if allow_parallel else "sequential"
    logger.info("[classic] Dispatching workflows to scheduler (%s mode)", mode_label)

    # Check if ESD is enabled - we need to know this before creating the scheduler
    from delfin.esd_module import parse_emission_rates, parse_esd_config
    esd_enabled, states, iscs, ics = parse_esd_config(config)
    emission_rates = parse_emission_rates(config)

    # If ESD is enabled, skip initial job (S0 from ESD replaces it)
    if esd_enabled and (states or iscs or ics or emission_rates):
        original_calc_initial = config.get('calc_initial', 'yes')
        if str(original_calc_initial).strip().lower() == 'yes':
            logger.info("ESD module enabled: skipping initial job (using S0 from ESD instead)")
            config['calc_initial'] = 'no'

    scheduler = GlobalOrcaScheduler(config, label="classic")
    try:
        # Add ESD jobs to scheduler FIRST (before execute_classic_workflows which calls run())
        # This allows ESD to run in parallel with ox/red steps after initial completes
        if esd_enabled and (states or iscs or ics or emission_rates):
            from delfin.esd_module import add_esd_jobs_to_scheduler
            add_esd_jobs_to_scheduler(
                scheduler,
                config,
                charge=charge,
                solvent=ctx.solvent,
                metals=ctx.metals if isinstance(ctx.metals, list) else [ctx.metals] if ctx.metals else [],
                main_basisset=ctx.main_basisset,
                metal_basisset=ctx.metal_basisset,
                dependency_job_id="classic_initial",
            )

            # Store ESD configuration in context for later use
            ctx.extra['esd_added_to_classic'] = esd_enabled
            ctx.extra['esd_states'] = states
            ctx.extra['esd_iscs'] = iscs
            ctx.extra['esd_ics'] = ics
            ctx.extra['esd_emission_rates'] = sorted(emission_rates)
        else:
            esd_enabled = False

        # Add classic jobs to scheduler and run all jobs together
        # execute_classic_workflows will call scheduler.run() internally
        result = execute_classic_workflows(
            config,
            allow_parallel=allow_parallel,
            scheduler=scheduler,
            **classic_kwargs,
        )

    finally:
        scheduler.shutdown()

    if not result.success:
        failed_desc = ", ".join(
            f"{job_id} ({reason})" for job_id, reason in result.failed.items()
        ) or "none"
        skipped_desc = ", ".join(
            f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
            for job_id, deps in result.skipped.items()
        ) or "none"
        logger.warning(
            "Classic workflows completed with issues; continuing. Failed jobs: %s | Skipped jobs: %s",
            failed_desc,
            skipped_desc,
        )

    ctx.extra['classic_result'] = result
    ctx.extra['ground_multiplicity'] = ground_multiplicity

    # If ESD jobs were added, store ESD results separately for downstream processing
    if esd_enabled:
        # Extract ESD-specific results from combined result
        esd_completed = {job_id for job_id in result.completed if job_id.startswith('esd_')}
        esd_failed = {job_id: reason for job_id, reason in result.failed.items() if job_id.startswith('esd_')}
        esd_skipped = {job_id: deps for job_id, deps in result.skipped.items() if job_id.startswith('esd_')}

        ctx.extra['esd_result'] = WorkflowRunResult(
            completed=esd_completed,
            failed=esd_failed,
            skipped=esd_skipped,
        )

        if esd_failed or esd_skipped:
            logger.warning(
                "ESD jobs completed with issues in classic phase. Failed: %s | Skipped: %s",
                ", ".join(f"{jid} ({reason})" for jid, reason in esd_failed.items()) or "none",
                ", ".join(f"{jid} (deps: {', '.join(deps)})" for jid, deps in esd_skipped.items()) or "none",
            )
        elif esd_completed:
            logger.info("ESD jobs completed successfully in parallel with classic phase")
    return {'result': result, 'ground_multiplicity': ground_multiplicity}


def run_manual_phase(ctx: PipelineContext) -> Dict[str, Any]:
    config = ctx.config
    multiplicity = config.get('multiplicity_0') or ctx.multiplicity

    _override_skip = _skip_preprocessing_for_override(config, ctx.control_file_path.parent)
    if config['XTB_OPT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_OPT: --occupier-override active, upstream geometry reused.")
        else:
            XTB(multiplicity, ctx.charge, config)

    if config['XTB_GOAT'] == "yes":
        if _override_skip:
            logger.info("[recalc] Skipping XTB_GOAT: --occupier-override active, upstream geometry reused.")
        elif _skip_xtb_goat_after_guppy(config):
            logger.info("Skipping XTB_GOAT: GUPPY already provided GOAT-refined winner geometry.")
        else:
            XTB_GOAT(multiplicity, ctx.charge, config)

    if config['CREST'] == "yes":
        run_crest_workflow(ctx.PAL, ctx.solvent, ctx.charge, multiplicity, ctx.config.get('input_file'))

    if config['XTB_SOLVATOR'] == "yes":
        XTB_SOLVATOR(
            ctx.config.get('input_file') or 'start.txt',
            multiplicity,
            ctx.charge,
            ctx.solvent,
            ctx.number_explicit_solv_molecules,
            config,
        )

    manual_kwargs = {
        'total_electrons_txt': ctx.total_electrons_txt,
        'xyz_file': ctx.file_bundle.xyz_initial,
        'xyz_file2': ctx.file_bundle.xyz_red1,
        'xyz_file3': ctx.file_bundle.xyz_red2,
        'xyz_file4': ctx.file_bundle.xyz_ox1,
        'xyz_file8': ctx.file_bundle.xyz_ox2,
        'output_file5': ctx.file_bundle.output_ox1,
        'output_file9': ctx.file_bundle.output_ox2,
        'output_file10': ctx.file_bundle.output_ox3,
        'output_file6': ctx.file_bundle.output_red1,
        'output_file7': ctx.file_bundle.output_red2,
        'output_file8': ctx.file_bundle.output_red3,
        'solvent': ctx.solvent,
        'metals': ctx.metals,
        'metal_basisset': ctx.metal_basisset,
        'main_basisset': ctx.main_basisset,
        'broken_sym': ctx.extra.get('ground_broken_sym', ""),
        'input_file_path': ctx.input_file,
        'output_initial': ctx.file_bundle.output_initial,
        'ground_multiplicity': config.get('multiplicity_0', 1),
        'ground_broken_sym': ctx.extra.get('ground_broken_sym', ""),
        'include_excited_jobs': True,
    }

    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    allow_parallel = parallel_mode != 'disable'
    mode_label = "parallel" if allow_parallel else "sequential"
    logger.info("[manually] Dispatching workflows to scheduler (%s mode)", mode_label)
    scheduler = GlobalOrcaScheduler(config, label="manually")
    try:
        result = execute_manually_workflows(
            config,
            allow_parallel=allow_parallel,
            scheduler=scheduler,
            **manual_kwargs,
        )
    finally:
        scheduler.shutdown()

    if not result.success:
        failed_desc = ", ".join(
            f"{job_id} ({reason})" for job_id, reason in result.failed.items()
        ) or "none"
        skipped_desc = ", ".join(
            f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
            for job_id, deps in result.skipped.items()
        ) or "none"
        logger.info(
            "Manual workflows completed with issues; continuing. Failed jobs: %s | Skipped jobs: %s",
            failed_desc,
            skipped_desc,
        )

    ctx.extra['manual_result'] = result
    ctx.multiplicity = int(config.get('multiplicity_0', ctx.multiplicity))
    return {'result': result}


# ---------------------------------------------------------------------------
# ESD (Excited State Dynamics) Module
# ---------------------------------------------------------------------------


def run_esd_phase(ctx: PipelineContext) -> bool:
    """Execute ESD calculations and store the outcome in the context."""
    config = ctx.config

    result = execute_esd_module(
        config=config,
        charge=ctx.charge,
        solvent=ctx.solvent,
        metals=ctx.metals if isinstance(ctx.metals, list) else [ctx.metals] if ctx.metals else [],
        main_basisset=ctx.main_basisset,
        metal_basisset=ctx.metal_basisset,
    )

    ctx.extra['esd_result'] = result

    # Generate UV-Vis spectrum report if ESD calculations were successful
    try:
        from delfin.reporting.uv_vis_report import generate_all_esd_uv_vis_reports
        esd_dir = Path("ESD").resolve()
        if esd_dir.exists():
            logger.info("Generating UV-Vis spectrum report from ESD results")
            generate_all_esd_uv_vis_reports(esd_dir)
    except Exception as e:
        logger.warning(f"Failed to generate UV-Vis spectrum report: {e}")

    if not result.success:
        failed_desc = ", ".join(
            f"{job_id} ({reason})" for job_id, reason in result.failed.items()
        ) or "none"
        skipped_desc = ", ".join(
            f"{job_id} (missing {', '.join(deps) if deps else 'unknown'})"
            for job_id, deps in result.skipped.items()
        ) or "none"
        logger.warning(
            "ESD module completed with issues; continuing. Failed jobs: %s | Skipped jobs: %s",
            failed_desc,
            skipped_desc,
        )
        return True

    logger.info("ESD module completed successfully")
    return True


# ---------------------------------------------------------------------------
# xTB Hyperpolarizability Phase
# ---------------------------------------------------------------------------


def run_hyperpol_xtb_phase(ctx: PipelineContext) -> bool:
    """Execute xTB-based hyperpolarizability calculation."""
    config = ctx.config
    if str(config.get('hyperpol_xTB', 'no')).strip().lower() != 'yes':
        return True

    from delfin.hyperpol import run_single_hyperpol_workflow
    from delfin.tadf_xtb import WorkflowEntry

    xyz_filename = str(config.get('hyperpol_xTB_xyz', 'start.txt')).strip()
    xyz_path = ctx.control_file_path.parent / xyz_filename
    if not xyz_path.exists():
        logger.error("hyperpol_xTB: source file '%s' not found — skipping", xyz_path)
        return False

    preopt = str(config.get('hyperpol_xTB_preopt', 'none')).strip().lower()
    engine = str(config.get('hyperpol_xTB_engine', 'std2')).strip().lower()
    use_bfw = str(config.get('hyperpol_xTB_bfw', 'no')).strip().lower() == 'yes'
    import math
    raw_wl = str(config.get('hyperpol_xTB_wavelengths', '')).strip()
    if raw_wl and raw_wl.lower() not in ('', 'none', 'static'):
        wavelengths = [float(w.strip()) for w in raw_wl.split(',') if w.strip()]
    else:
        wavelengths = [math.inf]  # static only
    energy_window = float(config.get('hyperpol_xTB_energy_window', 15.0))
    maxcore = int(str(config.get('maxcore', 1000)).strip())
    resolved_cores, resolved_maxcore = _resolve_pipeline_job_limits(
        requested_cores=ctx.PAL,
        requested_maxcore=maxcore,
    )
    label = ctx.name or "mol"
    workdir = ctx.control_file_path.parent / "hyperpol_xtb"

    if xyz_filename.endswith('.xyz'):
        entry = WorkflowEntry(label=label, xyz_path=str(xyz_path))
    else:
        entry = WorkflowEntry(label=label, xyz_text=xyz_path.read_text(encoding="utf-8", errors="ignore"))

    try:
        result = run_single_hyperpol_workflow(
            entry,
            charge=ctx.charge,
            multiplicity=ctx.multiplicity,
            engine=engine,
            preopt=preopt,
            wavelengths_nm=wavelengths,
            energy_window_ev=energy_window,
            cores=resolved_cores,
            maxcore=resolved_maxcore,
            workdir=workdir,
            use_bfw=use_bfw,
        )
        ctx.extra['hyperpol_xtb_result'] = result

        # Write JSON summary
        import json, math as _math
        from dataclasses import asdict
        json_path = workdir / "hyperpol_xtb_summary.json"
        json_payload = {
            "hyperpol_xtb": {
                "settings": {
                    "engine": result.engine,
                    "preopt": result.preopt,
                    "charge": result.charge,
                    "multiplicity": result.multiplicity,
                    "energy_window_ev": result.energy_window_ev,
                    "bfw": use_bfw,
                    "requested_wavelengths_nm": result.requested_wavelengths_nm,
                },
                "dipole_moment": {
                    "x_au": result.dipole_x_au,
                    "y_au": result.dipole_y_au,
                    "z_au": result.dipole_z_au,
                    "total_debye": result.dipole_total_debye,
                },
                "response_points": result.response_points,
                "files": {
                    "initial_xyz": result.initial_xyz,
                    "selected_xyz": result.selected_xyz,
                    "xtb4stda_output": result.xtb4stda_output,
                    "response_output": result.response_output,
                    "beta_hrs_file": result.beta_hrs_file,
                    "beta_tensor_file": result.beta_tensor_file,
                },
            }
        }
        json_path.write_text(json.dumps(json_payload, indent=2, default=str), encoding="utf-8")
        logger.info("hyperpol_xTB JSON written to %s", json_path)

        # Write human-readable TXT report (DELFIN.txt style)
        def _fmt_beta(label, point, au_key):
            v_au = point.get(au_key)
            if v_au is None:
                return f"  {label:<36s} = n/a"
            v_esu = point.get(au_key.replace("_au", "_esu"))
            v_esu_30 = point.get(au_key.replace("_au", "_esu_30"))
            return (
                f"  {label:<36s} = {float(v_au):>12.3f} au"
                f"    ({float(v_esu):.3e} esu; {float(v_esu_30):.3f} x10^-30 esu)"
            )

        sep = "-" * 72
        lines = [
            sep,
            "DELFIN — xTB Hyperpolarizability (sTD-DFT-xTB)",
            sep,
            "",
            "Calculation settings:",
            f"  Label          = {result.label}",
            f"  Engine         = {result.engine}",
            f"  Preopt         = {result.preopt}",
            f"  Charge         = {result.charge}",
            f"  Multiplicity   = {result.multiplicity}",
            f"  Energy window  = {result.energy_window_ev} eV",
            f"  BFW            = {'yes' if use_bfw else 'no'}",
            "",
            sep,
            "Dipole moment:",
            f"  mu_x           = {result.dipole_x_au:.6f} au" if result.dipole_x_au is not None else "  mu_x           = n/a",
            f"  mu_y           = {result.dipole_y_au:.6f} au" if result.dipole_y_au is not None else "  mu_y           = n/a",
            f"  mu_z           = {result.dipole_z_au:.6f} au" if result.dipole_z_au is not None else "  mu_z           = n/a",
            f"  |mu|           = {result.dipole_total_debye:.4f} Debye" if result.dipole_total_debye is not None else "  |mu|           = n/a",
            "",
        ]
        for point in result.response_points:
            wl = point.get("wavelength_nm")
            is_static = point.get("is_static", False)
            tag = "static" if is_static else f"{float(wl):.1f} nm"
            lines.append(sep)
            lines.append(f"First hyperpolarizability ({tag}):")
            lines.append("")
            lines.append(_fmt_beta(f"beta_HRS({tag})", point, "beta_hrs_au"))
            if point.get("beta_zzz_au") is not None:
                lines.append(_fmt_beta(f"beta_zzz({tag})", point, "beta_zzz_au"))
            if point.get("beta_zzz_aligned_au") is not None:
                lines.append(_fmt_beta(f"|beta_zzz| aligned to dipole({tag})", point, "beta_zzz_aligned_au"))
            for axis in ("x", "y", "z"):
                key = f"beta_vec_{axis}_au"
                if point.get(key) is not None:
                    lines.append(_fmt_beta(f"beta_vec_{axis}({tag})", point, key))
            if point.get("dr") is not None:
                lines.append(f"  {'Depolarization ratio (DR)':<36s} = {float(point['dr']):>12.3f}")
            if point.get("beta_sq_zzz_au2") is not None:
                lines.append(f"  {'<beta_zzz^2>':<36s} = {float(point['beta_sq_zzz_au2']):>12.3f} au^2")
            if point.get("beta_sq_xzz_au2") is not None:
                lines.append(f"  {'<beta_xzz^2>':<36s} = {float(point['beta_sq_xzz_au2']):>12.3f} au^2")
            if point.get("kleinman_score_raw_rotation") is not None:
                lines.append(f"  {'Kleinman symmetry score':<36s} = {float(point['kleinman_score_raw_rotation']):>12.6f}")
            lines.append("")
        lines.append(sep)
        lines.append(f"Working directory: {result.workdir}")
        lines.append(sep)
        txt_path = workdir / "hyperpol_xtb_summary.txt"
        txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        logger.info("hyperpol_xTB TXT report written to %s", txt_path)

        logger.info("hyperpol_xTB phase completed in %s", workdir)
        return True
    except Exception as exc:
        logger.error("hyperpol_xTB phase failed: %s", exc, exc_info=True)
        return False


# ---------------------------------------------------------------------------
# xTB TADF Screening Phase
# ---------------------------------------------------------------------------


def run_tadf_xtb_phase(ctx: PipelineContext) -> bool:
    """Execute xTB-based TADF screening calculation."""
    config = ctx.config
    if str(config.get('tadf_xTB', 'no')).strip().lower() != 'yes':
        return True

    from delfin.tadf_xtb import WorkflowEntry, run_single_tadf_xtb

    xyz_filename = str(config.get('tadf_xTB_xyz', 'start.txt')).strip()
    xyz_path = ctx.control_file_path.parent / xyz_filename
    if not xyz_path.exists():
        logger.error("tadf_xTB: source file '%s' not found — skipping", xyz_path)
        return False

    preopt = str(config.get('tadf_xTB_preopt', 'none')).strip().lower()
    excited_method = str(config.get('tadf_xTB_excited_method', 'stda')).strip().lower()
    use_bfw = str(config.get('tadf_xTB_bfw', 'no')).strip().lower() == 'yes'
    energy_window = float(config.get('tadf_xTB_energy_window', 10.0))
    run_t1_opt = str(config.get('tadf_xTB_run_t1_opt', 'yes')).strip().lower() == 'yes'
    xtb_method = str(config.get('xTB_method', 'XTB2')).strip()
    maxcore = int(str(config.get('maxcore', 1000)).strip())
    resolved_cores, resolved_maxcore = _resolve_pipeline_job_limits(
        requested_cores=ctx.PAL,
        requested_maxcore=maxcore,
    )
    label = ctx.name or "mol"
    workdir = ctx.control_file_path.parent / "tadf_xtb"

    use_crest = (preopt == "crest")
    use_goat = (preopt == "goat")
    optimize_s0 = (preopt == "xtb")

    if xyz_filename.endswith('.xyz'):
        entry = WorkflowEntry(label=label, xyz_path=str(xyz_path))
    else:
        entry = WorkflowEntry(label=label, xyz_text=xyz_path.read_text(encoding="utf-8", errors="ignore"))

    try:
        result = run_single_tadf_xtb(
            entry,
            charge=ctx.charge,
            multiplicity=ctx.multiplicity,
            xtb_method=xtb_method,
            excited_method=excited_method,
            energy_window=energy_window,
            cores=resolved_cores,
            maxcore=resolved_maxcore,
            workdir=workdir,
            use_crest=use_crest,
            use_goat=use_goat,
            run_t1_opt=run_t1_opt,
            t1_multiplicity=3,
            optimize_s0=optimize_s0,
            use_bfw=use_bfw,
        )
        ctx.extra['tadf_xtb_result'] = result

        # Write structured JSON summary
        import json
        json_path = workdir / "tadf_xtb_summary.json"

        def _opt_f(v, d=6):
            return round(float(v), d) if v is not None else None

        json_payload = {
            "tadf_xtb": {
                "settings": {
                    "label": result.label,
                    "excited_method": result.excited_method,
                    "xtb_method": result.xtb_method,
                    "charge": result.charge,
                    "multiplicity": result.multiplicity,
                    "energy_window_ev": energy_window,
                    "preopt": preopt,
                    "bfw": use_bfw,
                    "run_t1_opt": run_t1_opt,
                },
                "excited_states": {
                    "s1_ev": _opt_f(result.s1_ev, 4),
                    "s1_nm": _opt_f(result.s1_nm, 1),
                    "s1_f_osc": _opt_f(result.s1_f_osc, 5),
                    "t1_ev": _opt_f(result.t1_ev, 4),
                    "t1_nm": _opt_f(result.t1_nm, 1),
                    "delta_est_vertical_ev": _opt_f(result.delta_est_ev, 4),
                    "first_allowed_singlet": {
                        "state": result.first_allowed_singlet_state,
                        "ev": _opt_f(result.first_allowed_singlet_ev, 4),
                        "nm": _opt_f(result.first_allowed_singlet_nm, 1),
                        "f_osc": _opt_f(result.first_allowed_singlet_f_osc, 5),
                    } if result.first_allowed_singlet_ev is not None else None,
                    "brightest_singlet": {
                        "state": result.brightest_singlet_state,
                        "ev": _opt_f(result.brightest_singlet_ev, 4),
                        "nm": _opt_f(result.brightest_singlet_nm, 1),
                        "f_osc": _opt_f(result.brightest_singlet_f_osc, 5),
                    } if result.brightest_singlet_ev is not None else None,
                },
                "energetics": {
                    "s0_xtb_energy_eh": _opt_f(result.s0_xtb_energy_eh, 8),
                    "t1_vertical_xtb_ev": _opt_f(result.t1_vertical_xtb_ev, 4),
                    "t1_adiabatic_ev": _opt_f(result.t1_adiabatic_ev, 4),
                    "t1_relaxed_ev": _opt_f(result.t1_relaxed_ev, 4),
                    "t1_relaxation_ev": _opt_f(result.t1_relaxation_ev, 4),
                    "s1_relaxed_est_ev": _opt_f(result.s1_relaxed_est_ev, 4),
                },
                "photophysics": {
                    "lambda_abs_nm": _opt_f(result.lambda_abs_nm, 1),
                    "lambda_em_s1_nm": _opt_f(result.lambda_em_s1_nm, 1),
                    "lambda_pl_est_nm": _opt_f(result.lambda_pl_est_nm, 1),
                    "stokes_shift_est_ev": _opt_f(result.stokes_shift_est_ev, 4),
                    "delta_est_relaxed_est_ev": _opt_f(result.delta_est_relaxed_est_ev, 4),
                },
                "files": {
                    "workdir": str(result.workdir),
                    "s0_xyz": result.s0_xyz,
                    "t1_opt_xyz": result.t1_opt_xyz,
                    "crest_xyz": result.crest_xyz,
                    "goat_xyz": result.goat_xyz,
                },
            }
        }
        json_path.write_text(json.dumps(json_payload, indent=2, default=str), encoding="utf-8")
        logger.info("tadf_xTB JSON written to %s", json_path)

        # Write human-readable TXT report (DELFIN.txt style)
        def _fmt(v, d=3):
            return f"{v:.{d}f}" if v is not None else "n/a"

        sep = "-" * 72
        txt_path = workdir / "tadf_xtb_summary.txt"
        lines = [
            sep,
            "DELFIN — xTB TADF Screening (sTD-DFT-xTB)",
            sep,
            "",
            "Calculation settings:",
            f"  Label             = {result.label}",
            f"  Excited method    = {result.excited_method}",
            f"  xTB method        = {result.xtb_method}",
            f"  Charge            = {result.charge}",
            f"  Multiplicity      = {result.multiplicity}",
            f"  Preopt            = {preopt}",
            f"  BFW               = {'yes' if use_bfw else 'no'}",
            f"  Energy window     = {energy_window} eV",
            f"  T1 optimization   = {'yes' if run_t1_opt else 'no'}",
            "",
            sep,
            "Vertical excited states (at S0 geometry):",
            f"  S1                = {_fmt(result.s1_ev, 4):>10s} eV    ({_fmt(result.s1_nm, 1):>8s} nm,  f = {_fmt(result.s1_f_osc, 5)})",
            f"  T1                = {_fmt(result.t1_ev, 4):>10s} eV    ({_fmt(result.t1_nm, 1):>8s} nm)",
            f"  Delta_E(S-T) vert = {_fmt(result.delta_est_ev, 4):>10s} eV",
            "",
        ]

        if result.first_allowed_singlet_ev is not None:
            lines.append(
                f"  First allowed     = S{result.first_allowed_singlet_state:<3d}"
                f" {_fmt(result.first_allowed_singlet_ev, 4):>10s} eV"
                f"    ({_fmt(result.first_allowed_singlet_nm, 1):>8s} nm,"
                f"  f = {_fmt(result.first_allowed_singlet_f_osc, 5)})"
            )
        if result.brightest_singlet_ev is not None:
            lines.append(
                f"  Brightest singlet = S{result.brightest_singlet_state:<3d}"
                f" {_fmt(result.brightest_singlet_ev, 4):>10s} eV"
                f"    ({_fmt(result.brightest_singlet_nm, 1):>8s} nm,"
                f"  f = {_fmt(result.brightest_singlet_f_osc, 5)})"
            )
        if result.first_allowed_singlet_ev is not None or result.brightest_singlet_ev is not None:
            lines.append("")

        lines.append(sep)
        lines.append("Total energies:")
        if result.s0_xtb_energy_eh is not None:
            lines.append(f"  E(S0, xTB)        = {_fmt(result.s0_xtb_energy_eh, 8):>16s} Eh")
        if result.t1_xtb_energy_eh is not None:
            lines.append(f"  E(T1, xTB)        = {_fmt(result.t1_xtb_energy_eh, 8):>16s} Eh")
        lines.append("")

        lines.append(sep)
        lines.append("Adiabatic & relaxed energetics:")
        if result.t1_vertical_xtb_ev is not None:
            lines.append(f"  E(T1, vertical)   = {_fmt(result.t1_vertical_xtb_ev, 4):>10s} eV")
        if result.t1_adiabatic_ev is not None:
            lines.append(f"  E(T1, adiabatic)  = {_fmt(result.t1_adiabatic_ev, 4):>10s} eV")
        if result.t1_relaxed_ev is not None:
            lines.append(f"  E(T1, relaxed)    = {_fmt(result.t1_relaxed_ev, 4):>10s} eV")
        if result.t1_relaxation_ev is not None:
            lines.append(f"  T1 relaxation     = {_fmt(result.t1_relaxation_ev, 4):>10s} eV")
        if result.s1_relaxed_est_ev is not None:
            lines.append(f"  E(S1, rel. est.)  = {_fmt(result.s1_relaxed_est_ev, 4):>10s} eV")
        if result.delta_est_relaxed_est_ev is not None:
            lines.append(f"  Delta_E(S-T) rel. = {_fmt(result.delta_est_relaxed_est_ev, 4):>10s} eV")
        lines.append("")

        lines.append(sep)
        lines.append("Photophysical properties:")
        lines.append(f"  lambda_abs(S0->S1)  = {_fmt(result.lambda_abs_nm, 1):>8s} nm")
        if result.lambda_em_s1_nm is not None:
            lines.append(f"  lambda_em(S1)       = {_fmt(result.lambda_em_s1_nm, 1):>8s} nm")
        lines.append(f"  lambda_PL(est.)     = {_fmt(result.lambda_pl_est_nm, 1):>8s} nm")
        if result.stokes_shift_est_ev is not None:
            lines.append(f"  Stokes shift (est.) = {_fmt(result.stokes_shift_est_ev, 4):>10s} eV")
        lines.append("")

        lines.append(sep)
        lines.append(f"Working directory: {result.workdir}")
        lines.append(sep)

        txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        logger.info("tadf_xTB TXT report written to %s", txt_path)

        logger.info("tadf_xTB phase completed in %s", workdir)
        return True
    except Exception as exc:
        logger.error("tadf_xTB phase failed: %s", exc, exc_info=True)
        return False


# ---------------------------------------------------------------------------
# Post-processing and reporting
# ---------------------------------------------------------------------------


def collect_gibbs_energies(ctx: PipelineContext) -> Dict[str, Optional[float]]:
    """Collect Gibbs free energies for all redox states in the pipeline run.

    Thin wrapper around the canonical helper in ``delfin.energies``. Both
    this function and ``delfin.reporting.delfin_collector.collect_gibbs_energies``
    delegate to ``collect_gibbs_energies_from_dir`` so the state-file
    mapping and fallback rules live in exactly one place.
    """
    from delfin.energies import collect_gibbs_energies_from_dir

    esd_enabled, _, _, _ = parse_esd_config(ctx.config)
    return collect_gibbs_energies_from_dir(
        ctx.control_file_path.parent,
        esd_enabled=esd_enabled,
        occupier_fallback=True,
        log_progress=True,
    )


@dataclass
class SummaryResults:
    E_ox: Optional[float]
    E_ox_2: Optional[float]
    E_ox_3: Optional[float]
    E_red: Optional[float]
    E_red_2: Optional[float]
    E_red_3: Optional[float]
    E_00_t1: Optional[float]
    E_00_s1: Optional[float]
    multiplicity: int
    duration: float
    esd_summary: Optional[ESDSummary]


def compute_summary(ctx: PipelineContext, E_ref: float) -> SummaryResults:
    # Arbeitsverzeichnis: Ordner, in dem die CONTROL-Datei liegt
    working_dir = ctx.control_file_path.parent

    energies = collect_gibbs_energies(ctx)

    missing_potential_inputs = {
        '0': ('initial.out', ['E_ox', 'E_red', 'E_ox_2', 'E_red_2', 'E_ox_3', 'E_red_3']),
        '+1': ('ox_step_1.out', ['E_ox', 'E_ox_2']),
        '+2': ('ox_step_2.out', ['E_ox_2', 'E_ox_3']),
        '+3': ('ox_step_3.out', ['E_ox_3']),
        '-1': ('red_step_1.out', ['E_red', 'E_red_2']),
        '-2': ('red_step_2.out', ['E_red_2', 'E_red_3']),
        '-3': ('red_step_3.out', ['E_red_3']),
    }

    for key, (filename, potentials) in missing_potential_inputs.items():
        value = energies.get(key)
        file_path = working_dir / filename
        if value is None and file_path.exists():
            logger.info(
                "Skipping potentials %s (Gibbs data unavailable in %s)",
                ", ".join(potentials),
                filename,
            )

    m1_avg, m2_step, m3_mix, use_flags = calculate_redox_potentials(ctx.config, energies, E_ref)
    E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3 = select_final_potentials(m1_avg, m2_step, m3_mix, use_flags)

    # E_00 calculation from ESD results (if ESD module is enabled)
    E_00_t1 = None
    E_00_s1 = None
    E_00_t2 = None
    E_00_s2 = None

    esd_summary: Optional[ESDSummary] = None
    esd_enabled, esd_states, esd_iscs, esd_ics = parse_esd_config(ctx.config)
    if esd_enabled:
        esd_dir = working_dir / "ESD"
        esd_summary = collect_esd_results(esd_dir, esd_states, esd_iscs, esd_ics, config=ctx.config)

        # Calculate E_00 energies from ESD state results
        # E_00 = [E(excited) - E(S0)] + [ZPE(excited) - ZPE(S0)]
        if esd_summary and esd_summary.states:
            s0_data = esd_summary.states.get("S0")
            if s0_data and s0_data.fspe is not None and s0_data.zpe is not None:
                # E_00 for T1
                t1_data = esd_summary.states.get("T1")
                if t1_data and t1_data.fspe is not None and t1_data.zpe is not None:
                    E_00_t1 = ((t1_data.fspe - s0_data.fspe) + (t1_data.zpe - s0_data.zpe)) * 27.211386245988
                    logger.info("E_00(T1) from ESD (eV): %s", E_00_t1)

                # E_00 for S1
                s1_data = esd_summary.states.get("S1")
                if s1_data and s1_data.fspe is not None and s1_data.zpe is not None:
                    E_00_s1 = ((s1_data.fspe - s0_data.fspe) + (s1_data.zpe - s0_data.zpe)) * 27.211386245988
                    logger.info("E_00(S1) from ESD (eV): %s", E_00_s1)

                # E_00 for T2
                t2_data = esd_summary.states.get("T2")
                if t2_data and t2_data.fspe is not None and t2_data.zpe is not None:
                    E_00_t2 = ((t2_data.fspe - s0_data.fspe) + (t2_data.zpe - s0_data.zpe)) * 27.211386245988
                    logger.info("E_00(T2) from ESD (eV): %s", E_00_t2)

                # E_00 for S2
                s2_data = esd_summary.states.get("S2")
                if s2_data and s2_data.fspe is not None and s2_data.zpe is not None:
                    E_00_s2 = ((s2_data.fspe - s0_data.fspe) + (s2_data.zpe - s0_data.zpe)) * 27.211386245988
                    logger.info("E_00(S2) from ESD (eV): %s", E_00_s2)

    duration = time.time() - ctx.start_time
    return SummaryResults(
        E_ox=E_ox,
        E_ox_2=E_ox_2,
        E_ox_3=E_ox_3,
        E_red=E_red,
        E_red_2=E_red_2,
        E_red_3=E_red_3,
        E_00_t1=E_00_t1,
        E_00_s1=E_00_s1,
        multiplicity=ctx.multiplicity,
        duration=duration,
        esd_summary=esd_summary,
    )



def interpret_method_alias(raw_method: str) -> Tuple[str, Optional[str]]:
    method_aliases = {
        'classic': 'classic',
        'manual': 'manually',
        'manually': 'manually',
        'occupier': 'OCCUPIER',
        'occ': 'OCCUPIER',
        'occuper': 'OCCUPIER',
    }

    canonical_method = method_aliases.get(raw_method.lower())
    if canonical_method is None:
        suggestions = difflib.get_close_matches(raw_method, method_aliases.keys(), n=1)
        return raw_method, suggestions[0] if suggestions else None
    return canonical_method, None


def _is_truthy_token(value: Any) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1", "on"}


def _skip_xtb_goat_after_guppy(config: Dict[str, Any]) -> bool:
    return _is_truthy_token(config.get("_guppy_goat_completed", "no"))


def _skip_preprocessing_for_override(config: Dict[str, Any], workdir: Path) -> bool:
    """Return True when an --occupier-override is active and upstream artefacts
    (start.txt + initial ORCA output with OK marker) already exist. In that
    case GUPPY/XTB/XTB_GOAT must not be rerun — the override only affects
    downstream stages.
    """
    if not config.get("_occ_preferred_override"):
        return False
    start_txt = workdir / "start.txt"
    initial_out = workdir / "initial.out"
    if not (start_txt.exists() and initial_out.exists()):
        return False
    try:
        from delfin import smart_recalc
        if not smart_recalc.has_ok_marker(initial_out):
            return False
    except Exception:
        return False
    return True


def _run_guppy_for_smiles(smiles: str, start_path: Path, config: Dict[str, Any]) -> None:
    """Run GUPPY sampling for a SMILES string and write best geometry to start_path.

    Used when ``smiles_converter=GUPPY`` is selected in CONTROL.txt. After
    completion, start_path contains the lowest-energy
    XTB2-optimised geometry in DELFIN coordinate format (no XYZ header),
    ready for XTB_OPT / XTB_GOAT / GOAT / subsequent ORCA steps.
    """
    import os
    from delfin.guppy_sampling import run_sampling

    workdir = start_path.parent / "GUPPY"
    guppy_input = start_path.parent / "guppy_input.txt"
    guppy_output = workdir / "GUPPY_try.xyz"
    best_coord = workdir / "best_coordniation.xyz"
    guppy_settings_path = start_path.parent / "guppy_settings.json"

    guppy_input.write_text(smiles + "\n", encoding="utf-8")

    pal_raw = config.get('PAL') or os.environ.get('SLURM_CPUS_PER_TASK') or '40'
    maxcore_raw = config.get('maxcore') or os.environ.get('DELFIN_MAXCORE') or '6000'
    try:
        pal = int(str(pal_raw).strip())
    except (ValueError, TypeError):
        pal = 40
    try:
        maxcore = int(str(maxcore_raw).strip())
    except (ValueError, TypeError):
        maxcore = 6000

    try:
        runs = int(str(config.get('GUPPY_RUNS', 20)).strip())
    except (ValueError, TypeError):
        runs = 20
    if runs <= 0:
        runs = 20

    try:
        parallel_jobs = int(str(config.get('GUPPY_PARALLEL_JOBS', 4)).strip())
    except (ValueError, TypeError):
        parallel_jobs = 4
    if parallel_jobs <= 0:
        parallel_jobs = 4

    try:
        seed = int(str(config.get('GUPPY_SEED', 31)).strip())
    except (ValueError, TypeError):
        seed = 31

    method = str(config.get('xTB_method') or 'XTB2').strip() or 'XTB2'
    try:
        goat_topk = int(str(config.get('GUPPY_GOAT', 0)).strip())
    except (ValueError, TypeError):
        goat_topk = 0
    goat_topk = max(0, min(3, goat_topk))
    goat_parallel_jobs = parallel_jobs
    config['_guppy_goat_completed'] = 'no'

    guppy_settings = {
        'source': 'CONTROL.txt',
        'mode': 'guppy',
        'smiles': smiles,
        'runs': runs,
        'pal': pal,
        'maxcore': maxcore,
        'parallel_jobs': parallel_jobs,
        'goat_topk': goat_topk,
        'goat_parallel_jobs': goat_parallel_jobs,
        'seed': seed,
        'method': method,
        'input_file': guppy_input.name,
        'output_file': guppy_output.name,
        'winner_file': best_coord.name,
        'start_file': start_path.name,
        'cli_command': (
            f"python -m delfin.guppy_sampling {guppy_input.name} "
            f"--runs {runs} --pal {pal} --maxcore {maxcore} "
            f"--parallel-jobs {parallel_jobs} --goat-topk {goat_topk} "
            f"--goat-parallel-jobs {goat_parallel_jobs} --seed {seed} "
            f"--method {method} --output {guppy_output.name} --workdir {workdir.name}"
        ),
    }
    guppy_settings_path.write_text(json.dumps(guppy_settings, indent=2), encoding="utf-8")

    logger.info("smiles_converter=GUPPY: starting GUPPY sampling for SMILES → %s", start_path.name)
    ret = run_sampling(
        input_file=guppy_input,
        runs=runs,
        charge=None,
        pal=pal,
        maxcore=maxcore,
        parallel_jobs=parallel_jobs,
        method=method,
        output_file=guppy_output,
        workdir=workdir,
        seed=seed,
        allow_partial=True,
        goat_topk=goat_topk,
        goat_parallel_jobs=goat_parallel_jobs,
    )
    if ret != 0:
        raise RuntimeError("GUPPY sampling returned non-zero status")

    summary_path = workdir / "guppy_goat_summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if str(summary.get("winner_source", "")).strip().lower() == "goat":
                config['_guppy_goat_completed'] = 'yes'
                logger.info("GUPPY GOAT refinement selected final winner; downstream XTB_GOAT will be skipped.")
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to parse GUPPY GOAT summary (%s): %s", summary_path, exc)

    if not best_coord.exists():
        raise RuntimeError(
            "GUPPY sampling did not produce best_coordniation.xyz - "
            "check GUPPY/ subdirectory for details."
        )

    # Strip 2-line XYZ header (natoms + comment) → DELFIN coordinate format
    lines = best_coord.read_text(encoding="utf-8").splitlines()
    coord_lines = [ln for ln in lines[2:] if ln.strip()]
    start_path.write_text("\n".join(coord_lines) + "\n", encoding="utf-8")
    logger.info("GUPPY: best geometry written to %s (%d atoms)", start_path.name, len(coord_lines))


def _resolve_smiles_converter(config: Dict[str, Any]) -> str:
    """Return the effective SMILES converter with legacy fallback."""
    raw_value = str(config.get('smiles_converter', '') or '').strip().upper()
    if raw_value in {'QUICK', 'NORMAL', 'GUPPY', 'ARCHITECTOR'}:
        return raw_value
    if str(config.get('GUPPY', 'no')).strip().lower() == 'yes':
        return 'GUPPY'
    return 'NORMAL'


def _write_smiles_xyz_to_start(start_path: Path, xyz_content: str) -> None:
    """Write DELFIN-format coordinates to start.txt."""
    start_path.write_text(xyz_content, encoding='utf-8')
    coord_count = len([line for line in xyz_content.splitlines() if line.strip()])
    logger.info("Converted SMILES to XYZ coordinates in %s (%d atoms)", start_path.name, coord_count)


def normalize_input_file(config: Dict[str, Any], control_path: Path) -> str:
    input_entry = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    entry_path = Path(input_entry)
    if entry_path.is_absolute():
        input_path = resolve_path(entry_path)
    else:
        input_path = resolve_path(control_path.parent / entry_path)

    # Check if input contains SMILES
    from delfin.smiles_converter import (
        is_smiles_string,
        smiles_to_xyz,
        smiles_to_xyz_architector,
        smiles_to_xyz_quick,
    )

    is_smiles = False
    try:
        content = input_path.read_text(encoding='utf-8', errors='ignore')
        is_smiles = is_smiles_string(content)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not check for SMILES in '%s': %s", input_path, exc)

    if input_path.suffix.lower() == '.xyz':
        target = input_path.with_suffix('.txt')
        from delfin.define import convert_xyz_to_input_txt

        convert_xyz_to_input_txt(str(input_path), str(target))
        result_path = target
    else:
        result_path = input_path

    start_path = result_path.parent / 'start.txt'

    # Handle SMILES conversion: write XYZ directly to start.txt, keep input.txt unchanged
    if is_smiles:
        smiles_line = None
        for line in content.split('\n'):
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('*'):
                smiles_line = line
                break

        if smiles_line:
            logger.info("Detected SMILES in %s: %s", input_path.name, smiles_line)
            converter = _resolve_smiles_converter(config)
            logger.info("Using smiles_converter=%s for %s", converter, input_path.name)

            if converter == 'GUPPY':
                if _skip_preprocessing_for_override(config, start_path.parent):
                    logger.info(
                        "[recalc] Skipping GUPPY sampling: --occupier-override active and "
                        "start.txt + initial.out already exist (override is downstream-only)."
                    )
                    config['_guppy_goat_completed'] = 'yes'
                else:
                    try:
                        _run_guppy_for_smiles(smiles_line, start_path, config)
                    except Exception as exc:  # noqa: BLE001
                        logger.error("GUPPY sampling failed: %s", exc)
                        raise ValueError(f"GUPPY sampling failed: {exc}") from exc
            elif converter == 'QUICK':
                xyz_content, error = smiles_to_xyz_quick(smiles_line)

                if error:
                    logger.error("QUICK SMILES conversion failed: %s", error)
                    raise ValueError(f"QUICK SMILES conversion failed: {error}")

                try:
                    _write_smiles_xyz_to_start(start_path, xyz_content)
                except Exception as exc:  # noqa: BLE001
                    logger.error("Could not write QUICK-converted coordinates to '%s': %s", start_path, exc)
                    raise ValueError(f"Could not write QUICK-converted coordinates: {exc}") from exc
            elif converter == 'ARCHITECTOR':
                xyz_content, error = smiles_to_xyz_architector(smiles_line)

                if error:
                    logger.error("ARCHITECTOR SMILES conversion failed: %s", error)
                    raise ValueError(f"ARCHITECTOR SMILES conversion failed: {error}")

                try:
                    _write_smiles_xyz_to_start(start_path, xyz_content)
                except Exception as exc:  # noqa: BLE001
                    logger.error("Could not write ARCHITECTOR coordinates to '%s': %s", start_path, exc)
                    raise ValueError(f"Could not write ARCHITECTOR coordinates: {exc}") from exc
            else:
                # NORMAL: use the robust DELFIN converter first, then fall back
                # to the quick single-conformer path if the full conversion fails.
                xyz_content, error = smiles_to_xyz(smiles_line)

                if error:
                    logger.warning("NORMAL SMILES conversion failed, trying quick fallback: %s", error)
                    xyz_content, error = smiles_to_xyz_quick(smiles_line)

                if error:
                    logger.error("NORMAL SMILES conversion failed after quick fallback: %s", error)
                    raise ValueError(f"NORMAL SMILES conversion failed after quick fallback: {error}")

                try:
                    _write_smiles_xyz_to_start(start_path, xyz_content)
                except Exception as exc:  # noqa: BLE001
                    logger.error("Could not write converted coordinates to '%s': %s", start_path, exc)
                    raise ValueError(f"Could not write converted coordinates: {exc}") from exc
        else:
            logger.warning("SMILES format detected but no valid SMILES string found")
            # Fallback to normal copy
            try:
                if result_path.resolve() != start_path.resolve():
                    shutil.copyfile(result_path, start_path)
                elif not start_path.exists():
                    shutil.copyfile(result_path, start_path)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Could not create geometry backup '%s': %s", start_path, exc)
    else:
        # Normal case: copy input.txt to start.txt
        try:
            if result_path.resolve() != start_path.resolve():
                shutil.copyfile(result_path, start_path)
            elif not start_path.exists():
                shutil.copyfile(result_path, start_path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not create geometry backup '%s': %s", start_path, exc)

    work_path = start_path if start_path.exists() else result_path

    config.setdefault('input_file_original', str(result_path))
    config['input_file_backup'] = str(start_path)
    config['input_file'] = str(work_path)
    return str(work_path)
