"""ESD (Excited State Dynamics) module for DELFIN.

This module handles calculations of electronic states (S0, S1, T1, T2)
and their transitions (ISCs and ICs) in a separate ESD directory.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set
import threading

from delfin.common.logging import get_logger
from delfin.esd_input_generator import (
    create_ic_input,
    create_isc_input,
    create_state_input,
)
from delfin.orca import run_orca
from delfin.parallel_classic_manually import (
    WorkflowJob,
    WorkflowRunResult,
    _WorkflowManager,
)

logger = get_logger(__name__)

# Thread lock for input file generation to avoid race conditions
_input_generation_lock = threading.Lock()


def _resolve_esd_pal(config: Dict[str, Any], total_cores: int) -> int:
    """Resolve PAL (core count) for ESD transition jobs."""
    raw = config.get('ESD_PAL')
    try:
        pal_val = int(str(raw).strip()) if raw is not None else None
    except (TypeError, ValueError):
        pal_val = None
    if pal_val is None or pal_val <= 0:
        base = config.get('PAL', total_cores)
        try:
            pal_val = int(str(base).strip())
        except (TypeError, ValueError):
            pal_val = total_cores
    pal_val = max(1, min(pal_val, total_cores))
    return pal_val


def parse_esd_config(config: Dict[str, Any]) -> tuple[bool, List[str], List[str], List[str]]:
    """Parse ESD module configuration from control file.

    Args:
        config: Configuration dictionary

    Returns:
        Tuple of (esd_enabled, states, iscs, ics)
    """
    esd_enabled = str(config.get('ESD_modul', 'no')).strip().lower() == 'yes'

    # Handle both list and string formats
    states_raw = config.get('states', [])
    if isinstance(states_raw, list):
        states = [s.strip().upper() for s in states_raw if s.strip()]
    else:
        states = [s.strip().upper() for s in str(states_raw).split(',') if s.strip()]

    iscs_raw = config.get('ISCs', [])
    if isinstance(iscs_raw, list):
        iscs = [isc.strip() for isc in iscs_raw if isc.strip()]
    else:
        iscs = [isc.strip() for isc in str(iscs_raw).split(',') if isc.strip()]

    ics_raw = config.get('ICs', [])
    if isinstance(ics_raw, list):
        ics = [ic.strip() for ic in ics_raw if ic.strip()]
    else:
        ics = [ic.strip() for ic in str(ics_raw).split(',') if ic.strip()]

    logger.info(f"ESD config: enabled={esd_enabled}, states={states}, ISCs={iscs}, ICs={ics}")

    return esd_enabled, states, iscs, ics


def setup_esd_directory(esd_dir: Path, states: List[str]) -> None:
    """Set up ESD working directory and copy initial files if needed.

    Args:
        esd_dir: Path to ESD directory
        states: List of requested states
    """
    esd_dir.mkdir(exist_ok=True)
    logger.info(f"ESD directory created: {esd_dir}")

    # If S0 is in states and initial files exist, copy them to ESD directory
    if "S0" in states:
        files_to_copy = {
            "initial.out": "S0.out",
            "initial.gbw": "S0.gbw",
        }

        for src_name, dst_name in files_to_copy.items():
            src = Path(src_name)
            dst = esd_dir / dst_name

            if src.exists() and not dst.exists():
                try:
                    shutil.copy2(src, dst)
                    logger.info(f"Copied {src} → {dst}")
                except Exception as exc:
                    logger.warning(f"Failed to copy {src} to {dst}: {exc}")


def _populate_state_jobs(
    manager: _WorkflowManager,
    states: List[str],
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Add state calculation jobs to workflow manager.

    Args:
        manager: Workflow manager
        states: List of states to calculate
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary
    """
    # Define dependencies between states
    state_deps: Dict[str, Set[str]] = {
        "S0": set(),
        "S1": {"esd_S0"},
        "T1": {"esd_S0"},
        "T2": {"esd_T1"},
    }

    for state in states:
        state_upper = state.upper()

        # Check if this state should be calculated
        # S0 is special: only skip if both .out and .hess already exist
        if state_upper == "S0":
            s0_out = esd_dir / "S0.out"
            s0_hess = esd_dir / "S0.hess"
            if s0_out.exists() and s0_hess.exists():
                logger.info(f"S0 calculation skipped (S0.out/.hess exist in {esd_dir})")
                # Mark as completed so dependencies can proceed
                manager._completed.add("esd_S0")
                continue

        deps = state_deps.get(state_upper, set())

        def make_state_work(st: str) -> Callable[[int], None]:
            """Create work function for state calculation."""
            def work(cores: int) -> None:
                st_upper = st.upper()

                # Generate input file (thread-safe)
                with _input_generation_lock:
                    input_file = create_state_input(
                        st_upper,
                        esd_dir,
                        charge,
                        solvent,
                        metals,
                        main_basisset,
                        metal_basisset,
                        config,
                    )

                # Convert to absolute path before any chdir operations
                abs_input = Path(input_file).resolve()

                # Update PAL in input file (use absolute path)
                _update_pal_block(str(abs_input), cores)

                # Run ORCA in ESD directory
                output_file = esd_dir / f"{st_upper}.out"
                hess_file = esd_dir / f"{st_upper}.hess"

                logger.info(f"Running ORCA for state {st_upper} in {esd_dir}")

                # Run ORCA with absolute paths (no chdir needed)
                abs_output = output_file.resolve()

                # If Hessian missing but old output exists (e.g., recalc mode), force rerun
                if not hess_file.exists() and abs_output.exists():
                    try:
                        abs_output.unlink()
                    except Exception:  # noqa: BLE001
                        logger.debug("[esd] Could not remove stale output %s", abs_output, exc_info=True)

                # Change to ESD directory for calculation (esd_dir is already absolute)
                import os
                original_dir = os.getcwd()
                try:
                    os.chdir(esd_dir)
                    if not run_orca(abs_input.name, abs_output.name):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for {st_upper} state"
                        )
                finally:
                    os.chdir(original_dir)

                logger.info(f"State {st_upper} calculation completed")

            return work

        cores_min, cores_opt, cores_max = manager.derive_core_bounds()
        manager.add_job(
            WorkflowJob(
                job_id=f"esd_{state_upper}",
                work=make_state_work(state_upper),
                description=f"ESD {state_upper} optimization",
                dependencies=deps,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )


def _populate_isc_jobs(
    manager: _WorkflowManager,
    iscs: List[str],
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Add ISC calculation jobs to workflow manager.

    Args:
        manager: Workflow manager
        iscs: List of ISC transitions (e.g., ["S1>T1", "T1>S1"])
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary
    """
    for isc in iscs:
        initial_state, final_state = isc.split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()

        # ISC depends on both initial and final states
        deps = {f"esd_{initial_state}", f"esd_{final_state}"}

        job_id = f"esd_isc_{initial_state}_{final_state}"

        def make_isc_work(isc_pair: str) -> Callable[[int], None]:
            """Create work function for ISC calculation."""
            def work(cores: int) -> None:
                # Generate input file
                input_file = create_isc_input(
                    isc_pair,
                    esd_dir,
                    charge,
                    solvent,
                    metals,
                    main_basisset,
                    metal_basisset,
                    config,
                )

                # Convert to absolute path before any chdir operations
                abs_input = Path(input_file).resolve()

                # Update PAL in input file (use absolute path)
                _update_pal_block(str(abs_input), cores)

                # Determine output file name
                init_st, fin_st = isc_pair.split(">")
                init_st = init_st.strip().upper()
                fin_st = fin_st.strip().upper()
                output_file = esd_dir / f"{init_st}_{fin_st}_ISC.out"

                logger.info(f"Running ORCA for ISC {isc_pair} in {esd_dir}")

                # Run ORCA in ESD directory (esd_dir is already absolute)
                import os
                original_dir = os.getcwd()
                try:
                    os.chdir(esd_dir)
                    scratch_token = Path("scratch") / f"ISC_{init_st}_{fin_st}"
                    if not run_orca(
                        abs_input.name,
                        Path(output_file).name,
                        scratch_subdir=scratch_token,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for ISC {isc_pair}"
                        )
                finally:
                    os.chdir(original_dir)

                logger.info(f"ISC {isc_pair} calculation completed")

            return work

        pal_value = _resolve_esd_pal(config, manager.total_cores)
        cores_min = cores_opt = cores_max = pal_value

        manager.add_job(
            WorkflowJob(
                job_id=job_id,
                work=make_isc_work(isc),
                description=f"ISC {initial_state}→{final_state}",
                dependencies=deps,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )


def _populate_ic_jobs(
    manager: _WorkflowManager,
    ics: List[str],
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Add IC calculation jobs to workflow manager.

    Args:
        manager: Workflow manager
        ics: List of IC transitions (e.g., ["S1>S0", "T1>T2"])
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary
    """
    for ic in ics:
        initial_state, final_state = ic.split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()

        # IC depends on both initial and final states
        deps = {f"esd_{initial_state}", f"esd_{final_state}"}

        job_id = f"esd_ic_{initial_state}_{final_state}"

        def make_ic_work(ic_pair: str) -> Callable[[int], None]:
            """Create work function for IC calculation."""
            def work(cores: int) -> None:
                # Generate input file
                input_file = create_ic_input(
                    ic_pair,
                    esd_dir,
                    charge,
                    solvent,
                    metals,
                    main_basisset,
                    metal_basisset,
                    config,
                )

                # Convert to absolute path before any chdir operations
                abs_input = Path(input_file).resolve()

                # Update PAL in input file (use absolute path)
                _update_pal_block(str(abs_input), cores)

                # Determine output file name
                init_st, fin_st = ic_pair.split(">")
                init_st = init_st.strip().upper()
                fin_st = fin_st.strip().upper()
                output_file = esd_dir / f"{init_st}_{fin_st}_IC.out"

                logger.info(f"Running ORCA for IC {ic_pair} in {esd_dir}")

                # Run ORCA in ESD directory (esd_dir is already absolute)
                import os
                original_dir = os.getcwd()
                try:
                    os.chdir(esd_dir)
                    scratch_token = Path("scratch") / f"IC_{init_st}_{fin_st}"
                    if not run_orca(
                        abs_input.name,
                        Path(output_file).name,
                        scratch_subdir=scratch_token,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for IC {ic_pair}"
                        )
                finally:
                    os.chdir(original_dir)

                logger.info(f"IC {ic_pair} calculation completed")

            return work

        pal_value = _resolve_esd_pal(config, manager.total_cores)
        cores_min = cores_opt = cores_max = pal_value

        manager.add_job(
            WorkflowJob(
                job_id=job_id,
                work=make_ic_work(ic),
                description=f"IC {initial_state}→{final_state}",
                dependencies=deps,
                cores_min=cores_min,
                cores_optimal=cores_opt,
                cores_max=cores_max,
            )
        )


def _update_pal_block(input_path: str, cores: int) -> None:
    """Update %pal block in ORCA input file with given core count.

    Args:
        input_path: Path to input file
        cores: Number of cores to use
    """
    # Wait briefly for filesystem sync (parallel write issues)
    import time
    input_file_obj = Path(input_path)
    max_wait = 2  # seconds
    wait_step = 0.05
    elapsed = 0
    while not input_file_obj.exists() and elapsed < max_wait:
        time.sleep(wait_step)
        elapsed += wait_step

    try:
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as stream:
            lines = stream.readlines()
    except FileNotFoundError as exc:
        raise RuntimeError(f"Input file '{input_path}' missing") from exc

    pal_line = f"%pal nprocs {cores} end\n"
    replaced = False

    for idx, line in enumerate(lines):
        if line.strip().startswith('%pal'):
            lines[idx] = pal_line
            replaced = True
            break

    if not replaced:
        # Insert after other % blocks
        insert_idx = 0
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('%') and not stripped.startswith('%pal'):
                insert_idx = idx + 1
            elif stripped and not stripped.startswith('%'):
                break
        lines.insert(insert_idx, pal_line)

    with open(input_path, 'w', encoding='utf-8') as stream:
        stream.writelines(lines)


def run_esd_phase(
    config: Dict[str, Any],
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
) -> WorkflowRunResult:
    """Execute ESD module calculations.

    This is the main entry point for the ESD module. It:
    1. Parses ESD configuration
    2. Sets up ESD directory
    3. Schedules state, ISC, and IC calculations
    4. Executes all jobs in parallel

    Args:
        config: Configuration dictionary
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set

    Returns:
        WorkflowRunResult with completed/failed/skipped jobs
    """
    esd_enabled, states, iscs, ics = parse_esd_config(config)

    if not esd_enabled:
        logger.info("ESD module disabled (ESD_modul=no)")
        return WorkflowRunResult()

    if not states and not iscs and not ics:
        logger.info("ESD module enabled but no states/ISCs/ICs configured")
        return WorkflowRunResult()

    logger.info("Starting ESD module calculations")

    # Setup ESD directory (use absolute path to avoid chdir issues in parallel jobs)
    esd_dir = Path("ESD").resolve()
    setup_esd_directory(esd_dir, states)

    # Create workflow manager
    manager = _WorkflowManager(config, label="esd")

    try:
        # Populate jobs
        if states:
            _populate_state_jobs(
                manager,
                states,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )

        if iscs:
            _populate_isc_jobs(
                manager,
                iscs,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )

        if ics:
            _populate_ic_jobs(
                manager,
                ics,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )

        if not manager.has_jobs():
            logger.info("No ESD jobs to execute")
            return WorkflowRunResult()

        # Run all jobs
        logger.info(f"Executing {len(manager._jobs)} ESD jobs")
        manager.run()

        # Build result
        result = WorkflowRunResult(
            completed=set(manager.completed_jobs),
            failed=dict(manager.failed_jobs),
            skipped={
                job_id: list(deps) for job_id, deps in manager.skipped_jobs.items()
            },
        )

        if result.failed:
            logger.warning(
                f"ESD module completed with {len(result.failed)} failed jobs"
            )
        elif result.skipped:
            logger.warning(
                f"ESD module completed with {len(result.skipped)} skipped jobs"
            )
        else:
            logger.info("ESD module completed successfully")

        return result

    except Exception as exc:
        logger.error(f"ESD module failed: {exc}")
        result = WorkflowRunResult(
            completed=set(getattr(manager, 'completed_jobs', set())),
            failed=dict(getattr(manager, 'failed_jobs', {}) or {}),
            skipped={
                job_id: list(deps)
                for job_id, deps in (getattr(manager, 'skipped_jobs', {}) or {}).items()
            },
        )
        result.failed.setdefault('esd_error', f"{exc.__class__.__name__}: {exc}")
        return result

    finally:
        manager.shutdown()
