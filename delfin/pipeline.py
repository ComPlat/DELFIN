"""Workflow execution helpers for CLI orchestration."""

from __future__ import annotations

import difflib
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.copy_helpers import copy_if_exists, read_occupier_file
from delfin.global_scheduler import GlobalOrcaScheduler
from delfin.parallel_occupier import OccupierExecutionContext, run_occupier_orca_jobs
from delfin.parallel_classic_manually import execute_classic_workflows, execute_manually_workflows, normalize_parallel_token, WorkflowRunResult
from delfin.xtb_crest import XTB, XTB_GOAT, XTB_SOLVATOR, run_crest_workflow
from delfin.cli_calculations import calculate_redox_potentials, select_final_potentials
from delfin.energies import find_gibbs_energy
from delfin.esd_module import run_esd_phase as execute_esd_module, parse_esd_config
from delfin.esd_results import collect_esd_results, ESDSummary

logger = get_logger(__name__)




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

    if config['XTB_OPT'] == "yes":
        XTB(multiplicity, charge, config)

    if config['XTB_GOAT'] == "yes":
        XTB_GOAT(multiplicity, charge, config)

    if config['CREST'] == "yes":
        run_crest_workflow(ctx.PAL, ctx.solvent, charge, multiplicity, ctx.config.get('input_file'))

    metals_list = list(ctx.metals) if isinstance(ctx.metals, (list, tuple, set)) else ([ctx.metals] if ctx.metals else [])

    if config['XTB_SOLVATOR'] == "no":
        # Run ALL OCCUPIER jobs (initial + ox/red) + post-processing in ONE scheduler run
        if "yes" in config.get("calc_initial", ""):
            print("\nOCCUPIER for the initial system:\n")

            from delfin.parallel_occupier import build_flat_occupier_fob_jobs
            from delfin.parallel_classic_manually import _WorkflowManager

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

            logger.info("All OCCUPIER workflows completed successfully")

            # Post-processing jobs are scheduled within the combined OCCUPIER run
            config['_used_combined_occupier'] = True
        else:
            # No initial calculation, but still run ox/red if configured
            from delfin.parallel_occupier import build_flat_occupier_fob_jobs
            from delfin.parallel_classic_manually import _WorkflowManager

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
            multiplicity_0, additions_0, _, gbw_initial = read_occupier_file(
                "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
            )
            ctx.extra['multiplicity_0'] = multiplicity_0
            ctx.extra['additions_0'] = additions_0
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

                from delfin.parallel_occupier import build_flat_occupier_fob_jobs
                from delfin.parallel_classic_manually import _WorkflowManager

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

        multiplicity_0, additions_0, _, gbw_initial = read_occupier_file(
            "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
        )
        ctx.extra['multiplicity_0'] = multiplicity_0
        ctx.extra['additions_0'] = additions_0
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

    if config['XTB_OPT'] == "yes":
        XTB(multiplicity, charge, config)

    if config['XTB_GOAT'] == "yes":
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
        'additions': "",
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

    if config['XTB_OPT'] == "yes":
        XTB(multiplicity, ctx.charge, config)

    if config['XTB_GOAT'] == "yes":
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
        'additions': ctx.extra.get('ground_additions', ""),
        'input_file_path': ctx.input_file,
        'output_initial': ctx.file_bundle.output_initial,
        'ground_multiplicity': config.get('multiplicity_0', 1),
        'ground_additions': ctx.extra.get('ground_additions', ""),
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
# Post-processing and reporting
# ---------------------------------------------------------------------------


def collect_gibbs_energies(ctx: PipelineContext) -> Dict[str, Optional[float]]:
    bundle = ctx.file_bundle
    file_map = {
        '0': 'initial.out',
        '+1': 'ox_step_1.out',
        '+2': 'ox_step_2.out',
        '+3': 'ox_step_3.out',
        '-1': 'red_step_1.out',
        '-2': 'red_step_2.out',
        '-3': 'red_step_3.out',
    }

    energies: Dict[str, Optional[float]] = {}

    # Determine working directory from control file path
    working_dir = ctx.control_file_path.parent

    for key, filename in file_map.items():
        path = working_dir / filename

        # Special handling for initial.out: use ESD/S0.out if initial.out doesn't exist and ESD is enabled
        if key == '0' and not path.exists():
            esd_enabled, _, _, _ = parse_esd_config(ctx.config)
            if esd_enabled:
                esd_s0_path = working_dir / "ESD" / "S0.out"
                if esd_s0_path.exists():
                    path = esd_s0_path
                    logger.info("Using ESD/S0.out for ground state Gibbs energy (initial.out not found)")

        if not path.exists():
            energies[key] = None
            continue

        value = find_gibbs_energy(str(path))
        energies[key] = value
        if value is not None:
            logger.info("Free Gibbs Free Energy %s (H): %s", key, value)
        else:
            logger.info(
                "Skipping Gibbs energy for state %s (data unavailable in %s)",
                key,
                path.name,
            )
    return energies


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


def normalize_input_file(config: Dict[str, Any], control_path: Path) -> str:
    input_entry = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    entry_path = Path(input_entry)
    if entry_path.is_absolute():
        input_path = resolve_path(entry_path)
    else:
        input_path = resolve_path(control_path.parent / entry_path)

    # Check if input contains SMILES
    from .smiles_converter import is_smiles_string, smiles_to_xyz

    is_smiles = False
    try:
        content = input_path.read_text(encoding='utf-8', errors='ignore')
        is_smiles = is_smiles_string(content)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not check for SMILES in '%s': %s", input_path, exc)

    if input_path.suffix.lower() == '.xyz':
        target = input_path.with_suffix('.txt')
        from .define import convert_xyz_to_input_txt

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
            xyz_content, error = smiles_to_xyz(smiles_line)

            if error:
                logger.error("SMILES conversion failed: %s", error)
                raise ValueError(f"SMILES conversion failed: {error}")

            # Write converted XYZ to start.txt
            try:
                start_path.write_text(xyz_content, encoding='utf-8')
                logger.info("Converted SMILES to XYZ coordinates in %s (input.txt unchanged)", start_path.name)
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
