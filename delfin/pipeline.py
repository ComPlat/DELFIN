"""Workflow execution helpers for CLI orchestration."""

from __future__ import annotations

import difflib
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path
from delfin.copy_helpers import copy_if_exists, prepare_occ_folder, read_occupier_file
from delfin.parallel_occupier import OccupierExecutionContext, run_occupier_orca_jobs
from delfin.parallel_classic_manually import execute_classic_workflows, execute_manually_workflows, normalize_parallel_token
from delfin.xtb_crest import XTB, XTB_GOAT, XTB_SOLVATOR, run_crest_workflow
from delfin.cli_calculations import calculate_redox_potentials, select_final_potentials
from delfin.thread_safe_helpers import prepare_occ_folder_2_threadsafe
from delfin.energies import find_gibbs_energy, find_ZPE, find_electronic_energy

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


def _execute_oxidation_workflow(config: Dict[str, Any], pal_override: Optional[int] = None) -> bool:
    """Execute oxidation steps workflow."""

    logger.info("Starting oxidation workflow")
    original_cwd = Path.cwd()

    try:
        if "1" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the first oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "ox_step_1_OCCUPIER",
                source_occ_folder="initial_OCCUPIER",
                charge_delta=+1,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare ox_step_1_OCCUPIER")
                return False

        if "2" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the second oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "ox_step_2_OCCUPIER",
                source_occ_folder="ox_step_1_OCCUPIER",
                charge_delta=+2,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare ox_step_2_OCCUPIER")
                return False

        if "3" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the third oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "ox_step_3_OCCUPIER",
                source_occ_folder="ox_step_2_OCCUPIER",
                charge_delta=+3,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare ox_step_3_OCCUPIER")
                return False

        logger.info("Oxidation workflow completed")
        return True

    except Exception as exc:  # noqa: BLE001
        logger.error("Oxidation workflow failed: %s", exc)
        return False


def _execute_reduction_workflow(config: Dict[str, Any], pal_override: Optional[int] = None) -> bool:
    """Execute reduction steps workflow."""

    logger.info("Starting reduction workflow")
    original_cwd = Path.cwd()

    try:
        if "1" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the first reduction step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "red_step_1_OCCUPIER",
                source_occ_folder="initial_OCCUPIER",
                charge_delta=-1,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare red_step_1_OCCUPIER")
                return False

        if "2" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the second reduction step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "red_step_2_OCCUPIER",
                source_occ_folder="red_step_1_OCCUPIER",
                charge_delta=-2,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare red_step_2_OCCUPIER")
                return False

        if "3" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the third reduction step:\n")
            success = prepare_occ_folder_2_threadsafe(
                "red_step_3_OCCUPIER",
                source_occ_folder="red_step_2_OCCUPIER",
                charge_delta=-3,
                config=config,
                original_cwd=original_cwd,
                pal_override=pal_override,
            )
            if not success:
                logger.error("Failed to prepare red_step_3_OCCUPIER")
                return False

        logger.info("Reduction workflow completed")
        return True

    except Exception as exc:  # noqa: BLE001
        logger.error("Reduction workflow failed: %s", exc)
        return False


def _should_parallelize(config: Dict[str, Any]) -> Tuple[bool, bool, bool]:
    parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
    has_ox = bool(config.get("oxidation_steps", "").strip())
    has_red = bool(config.get("reduction_steps", "").strip())
    can_parallel = has_ox and has_red
    enable_parallel = (
        parallel_mode == 'enable'
        or (parallel_mode == 'auto' and can_parallel)
    )
    return enable_parallel, has_ox, has_red


def _execute_parallel_workflows(config: Dict[str, Any]) -> bool:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    total_pal = max(1, int(config.get('PAL', 12)))
    _, has_ox, has_red = _should_parallelize(config)
    reduced_pal = max(2, total_pal // 2) if has_ox and has_red else None

    if reduced_pal:
        def _format_banner(orig_pal: int, split_pal: int) -> str:
            box_width = 65
            inner_width = box_width - 2

            def _box_line(text: str = "") -> str:
                return f"║{text.ljust(inner_width)}║"

            def _label(label: str, value: str) -> str:
                return f"{label:<29}: {value}".ljust(inner_width)

            header = "PARALLEL WORKFLOW EXECUTION".center(inner_width)
            body_lines = [
                header,
                _label("Original PAL", str(orig_pal)),
                _label("Reduced PAL per workflow", str(split_pal)),
                _label("Reason", "ox and red running in parallel"),
            ]

            top = "╔" + "═" * (box_width - 2) + "╗"
            bottom = "╚" + "═" * (box_width - 2) + "╝"
            return "\n".join([top, *[_box_line(line) for line in body_lines], bottom])

        print(_format_banner(total_pal, reduced_pal))
    else:
        logger.info("Single workflow execution - using full PAL=%s", total_pal)

    tasks: List[Tuple[str, Any]] = []
    if config.get("oxidation_steps", ""):
        tasks.append(("oxidation", (_execute_oxidation_workflow, reduced_pal)))
    if config.get("reduction_steps", ""):
        tasks.append(("reduction", (_execute_reduction_workflow, reduced_pal)))

    all_success = True
    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        futures = {
            executor.submit(func, config, pal_override): name
            for name, (func, pal_override) in tasks
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                success = future.result()
            except Exception as exc:  # noqa: BLE001
                logger.error("Parallel OCCUPIER %s workflow raised exception: %s", name, exc)
                all_success = False
            else:
                if success:
                    logger.info("OCCUPIER %s workflow completed successfully", name)
                else:
                    logger.error("OCCUPIER %s workflow failed", name)
                    all_success = False

    if all_success:
        logger.info("All parallel workflows completed successfully")
    else:
        logger.error("Some parallel workflows failed")

    return all_success


def _execute_sequential_workflows(config: Dict[str, Any]) -> bool:
    all_success = True
    total_pal = max(1, int(config.get('PAL', 12)))
    logger.info("Sequential workflow execution - each workflow uses full PAL=%s", total_pal)

    if config.get("oxidation_steps", ""):
        all_success &= _execute_oxidation_workflow(config, pal_override=None)
    if config.get("reduction_steps", ""):
        all_success &= _execute_reduction_workflow(config, pal_override=None)

    return all_success


def _run_occ_workflows(config: Dict[str, Any]) -> None:
    needs_ox_red = config.get("oxidation_steps", "").strip() or config.get("reduction_steps", "").strip()
    if not needs_ox_red:
        logger.info("No oxidation or reduction steps configured")
        return

    should_parallel, has_ox, has_red = _should_parallelize(config)
    if should_parallel and has_ox and has_red:
        _execute_parallel_workflows(config)
    else:
        _execute_sequential_workflows(config)


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

    if config['XTB_SOLVATOR'] == "no":
        if "yes" in config.get("calc_initial", ""):
            print("\nOCCUPIER for the initial system:\n")
            prepare_occ_folder("initial_OCCUPIER", charge_delta=0)
            logger.info("Initial OCCUPIER completed successfully")

        _run_occ_workflows(config)

        if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == "yes":
            multiplicity_0, additions_0, _ = read_occupier_file(
                "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
            )
            ctx.extra['multiplicity_0'] = multiplicity_0
            ctx.extra['additions_0'] = additions_0

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
        if initial_requested or not initial_report.exists():
            print("\nOCCUPIER for the initial system:\n")
            prepare_occ_folder("initial_OCCUPIER", charge_delta=0)
            initial_rerun = True
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

            need_occ_workflows = initial_rerun

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
                logger.info("Preparing OCCUPIER oxidation/reduction workflows prior to solvation")
                _run_occ_workflows(config)
            else:
                logger.info("Reusing existing OCCUPIER oxidation/reduction workflows")

        multiplicity_0, additions_0, _ = read_occupier_file(
            "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
        )
        ctx.extra['multiplicity_0'] = multiplicity_0
        ctx.extra['additions_0'] = additions_0

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

    if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() != "yes":
        occ_context = OccupierExecutionContext(
            charge=charge,
            solvent=ctx.solvent,
            metals=metals_list,
            main_basisset=ctx.main_basisset,
            metal_basisset=ctx.metal_basisset,
            config=config,
        )
        occ_success = run_occupier_orca_jobs(occ_context, parallel_enabled)
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
    result = execute_classic_workflows(config, allow_parallel=allow_parallel, **classic_kwargs)

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
    result = execute_manually_workflows(config, allow_parallel=allow_parallel, **manual_kwargs)

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

    for key, filename in file_map.items():
        path = Path(filename)
        if not path.exists():
            energies[key] = None
            continue

        value = find_gibbs_energy(filename)
        energies[key] = value
        if value is not None:
            logger.info("Free Gibbs Free Energy %s (H): %s", key, value)
        else:
            logger.info(
                "Skipping Gibbs energy for state %s (data unavailable in %s)",
                key,
                filename,
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


def compute_summary(ctx: PipelineContext, E_ref: float) -> SummaryResults:
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
        if value is None and Path(filename).exists():
            logger.info(
                "Skipping potentials %s (Gibbs data unavailable in %s)",
                ", ".join(potentials),
                filename,
            )

    m1_avg, m2_step, m3_mix, use_flags = calculate_redox_potentials(ctx.config, energies, E_ref)
    E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3 = select_final_potentials(m1_avg, m2_step, m3_mix, use_flags)

    ZPE_S0 = find_ZPE('initial.out')
    ZPE_T1 = find_ZPE('t1_state_opt.out') if Path('t1_state_opt.out').exists() else None
    ZPE_S1 = find_ZPE('s1_state_opt.out') if Path('s1_state_opt.out').exists() else None
    E_0 = find_electronic_energy('initial.out')
    E_T1 = find_electronic_energy('t1_state_opt.out') if Path('t1_state_opt.out').exists() else None
    E_S1 = find_electronic_energy('s1_state_opt.out') if Path('s1_state_opt.out').exists() else None

    E_00_t1 = None
    E_00_s1 = None
    def _missing_components(components: dict[str, Optional[float]]) -> tuple[list[str], list[str]]:
        missing = []
        missing_files = []
        for label, value in components.items():
            if value is not None:
                continue
            if '(' in label and label.endswith(')'):
                filename = label[label.find('(') + 1:-1]
                if not Path(filename).exists():
                    missing_files.append(filename)
                    continue
            missing.append(label)
        return missing, missing_files

    if ctx.config['E_00'] == "yes":
        excitation_flags = ctx.config.get("excitation", "")

        if "t" in excitation_flags:
            requirements = {
                "ZPE(initial.out)": ZPE_S0,
                "ZPE(t1_state_opt.out)": ZPE_T1,
                "Energy(initial.out)": E_0,
                "Energy(t1_state_opt.out)": E_T1,
            }
            missing, missing_files = _missing_components(requirements)
            if not missing and not missing_files:
                E_00_t1 = ((E_T1 - E_0) + (ZPE_T1 - ZPE_S0)) * 27.211386245988
                logger.info("E_00_t (eV): %s", E_00_t1)
            else:
                if missing:
                    logger.info(
                        "Skipping E_00_t calculation (data unavailable: %s)",
                        ", ".join(missing),
                    )

        if "s" in excitation_flags:
            requirements = {
                "ZPE(initial.out)": ZPE_S0,
                "ZPE(s1_state_opt.out)": ZPE_S1,
                "Energy(initial.out)": E_0,
                "Energy(s1_state_opt.out)": E_S1,
            }
            missing, missing_files = _missing_components(requirements)
            if not missing and not missing_files:
                E_00_s1 = ((E_S1 - E_0) + (ZPE_S1 - ZPE_S0)) * 27.211386245988
                logger.info("E_00_s (eV): %s", E_00_s1)
            else:
                if missing:
                    logger.info(
                        "Skipping E_00_s calculation (data unavailable: %s)",
                        ", ".join(missing),
                    )

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
    if input_path.suffix.lower() == '.xyz':
        target = input_path.with_suffix('.txt')
        from .define import convert_xyz_to_input_txt

        convert_xyz_to_input_txt(str(input_path), str(target))
        result_path = target
    else:
        result_path = input_path

    start_path = result_path.parent / 'start.txt'
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
