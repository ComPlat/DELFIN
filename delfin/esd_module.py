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
    create_fluor_input,
    create_ic_input,
    create_isc_input,
    create_phosp_input,
    create_state_input,
    _format_ms_suffix,
)
from delfin.orca import run_orca_with_intelligent_recovery
from delfin.parallel_classic_manually import (
    WorkflowJob,
    WorkflowRunResult,
    _WorkflowManager,
)

logger = get_logger(__name__)

# Thread lock for input file generation to avoid race conditions
_input_generation_lock = threading.Lock()


def _run_orca_esd(
    inp_file: Path | str,
    out_file: Path | str,
    *,
    working_dir: Path,
    config: Dict[str, Any],
    scratch_subdir: Optional[Path] = None,
) -> bool:
    """Run ORCA for ESD with intelligent recovery if enabled."""
    return run_orca_with_intelligent_recovery(
        str(inp_file),
        str(out_file),
        scratch_subdir=scratch_subdir,
        working_dir=working_dir,
        config=config,
    )


def _convert_start_to_xyz(start_txt: Path, output_xyz: Path) -> None:
    """Convert start.txt (DELFIN format) to XYZ format with header.

    DELFIN format (start.txt): element x y z (no header)
    XYZ format: atom_count, comment line, then element x y z

    Args:
        start_txt: Path to start.txt file
        output_xyz: Path to output XYZ file
    """
    # Read start.txt
    with open(start_txt, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Count atoms
    atom_count = len(lines)

    # Write XYZ format
    with open(output_xyz, 'w', encoding='utf-8') as f:
        f.write(f"{atom_count}\n")
        f.write("Converted from start.txt\n")
        for line in lines:
            f.write(f"{line}\n")


def _run_xtb_goat_if_needed(config: Dict[str, Any], charge: int) -> None:
    """Run xTB optimization and/or GOAT if enabled in config.

    This is called before ESD calculations to ensure we use optimized geometries.
    Updates start.txt in-place with optimized coordinates.

    Args:
        config: DELFIN configuration dictionary
        charge: Molecular charge
    """
    from delfin.xtb_crest import XTB, XTB_GOAT

    # Determine multiplicity (for ESD, typically singlet S0)
    multiplicity = 1  # Default for ground state S0

    if config.get('XTB_OPT') == 'yes':
        logger.info("Running xTB optimization before ESD calculations")
        try:
            XTB(multiplicity, charge, config)
        except Exception as exc:
            logger.warning(f"xTB optimization failed: {exc}")

    if config.get('XTB_GOAT') == 'yes':
        logger.info("Running GOAT optimization before ESD calculations")
        try:
            XTB_GOAT(multiplicity, charge, config)
        except Exception as exc:
            logger.warning(f"GOAT optimization failed: {exc}")


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

    # Ensure ground-state S0 is present when ESD module is enabled
    # S0 is always calculated as minimum when ESD_modul=yes
    if esd_enabled:
        states_set = {s.upper() for s in states}
        if "S0" not in states_set:
            states = ["S0"] + states

    # logger.info(f"ESD config: enabled={esd_enabled}, states={states}, ISCs={iscs}, ICs={ics}")

    return esd_enabled, states, iscs, ics


def parse_emission_rates(config: Dict[str, Any]) -> Set[str]:
    """Parse CONTROL emission_rates into a set like {'f','p'}.

    Accepted separators: comma, semicolon, whitespace. Case-insensitive.
    """
    raw = config.get("emission_rates", "")
    if raw is None:
        return set()
    if isinstance(raw, list):
        tokens = [str(x) for x in raw]
    else:
        tokens = str(raw).replace(";", ",").replace(" ", ",").split(",")
    rates = {t.strip().lower() for t in tokens if str(t).strip()}
    return {r for r in rates if r in {"f", "p"}}


def setup_esd_directory(esd_dir: Path, states: List[str]) -> None:
    """Set up ESD working directory (creates directory only, no file copies).

    Args:
        esd_dir: Path to ESD directory
        states: List of requested states

    Note:
        File copying (initial.xyz → S0.xyz, etc.) is handled by individual jobs
        to ensure they run AFTER classic_initial completes.
    """
    esd_dir.mkdir(exist_ok=True)
    logger.info(f"ESD directory created: {esd_dir}")


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
    esd_modus = str(config.get("ESD_modus", "TDDFT")).strip().lower()
    # If pipe-separated options (e.g., "TDDFT|deltaSCF"), take first as default
    if "|" in esd_modus:
        esd_modus = esd_modus.split("|")[0].strip()

    # Define dependencies between states
    # In TDDFT mode: all excited states depend only on S0 (ground state wavefunction)
    # In deltaSCF mode: sequential dependencies for orbital optimization
    state_deps: Dict[str, Set[str]] = {
        "S0": set(),
        "S1": {"esd_S0"},
        "S2": {"esd_S0"} if esd_modus == "tddft" else {"esd_S1"},
        "S3": {"esd_S0"} if esd_modus == "tddft" else {"esd_S2"},
        "S4": {"esd_S0"} if esd_modus == "tddft" else {"esd_S3"},
        "S5": {"esd_S0"} if esd_modus == "tddft" else {"esd_S4"},
        "S6": {"esd_S0"} if esd_modus == "tddft" else {"esd_S5"},
        "T1": {"esd_S0"},
        "T2": {"esd_S0"} if esd_modus == "tddft" else {"esd_T1"},
        "T3": {"esd_S0"} if esd_modus == "tddft" else {"esd_T2"},
        "T4": {"esd_S0"} if esd_modus == "tddft" else {"esd_T3"},
        "T5": {"esd_S0"} if esd_modus == "tddft" else {"esd_T4"},
        "T6": {"esd_S0"} if esd_modus == "tddft" else {"esd_T5"},
    }

    for state in states:
        state_upper = state.upper()

        # Check if this state should be calculated
        # In normal mode: skip if both .out and .hess exist
        # In recalc mode: let the recalc wrapper decide (it checks TERMINATED NORMALLY)
        import os
        recalc_mode = os.environ.get("DELFIN_RECALC", "0") == "1"

        if not recalc_mode:
            state_out = esd_dir / f"{state_upper}.out"

            # Check if output is complete (TERMINATED NORMALLY)
            is_complete = False
            if state_out.exists():
                try:
                    with state_out.open("r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                        is_complete = "ORCA TERMINATED NORMALLY" in content
                except Exception:
                    pass

            if is_complete:
                logger.info(f"{state_upper} calculation skipped (ORCA TERMINATED NORMALLY in {state_out.name})")
                # Mark as completed so dependencies can proceed
                manager._completed.add(f"esd_{state_upper}")
                continue

        deps = set(state_deps.get(state_upper, set()))
        # In TDDFT mode, all excited states read S0.xyz; enforce dependency on S0 unless already present
        if esd_modus == "tddft" and state_upper != "S0" and "esd_S0" not in deps:
            deps.add("esd_S0")

        def make_state_work(st: str) -> Callable[[int], None]:
            """Create work function for state calculation."""
            def work(cores: int) -> None:
                st_upper = st.upper()

                # For S0: copy optimized initial.xyz to S0.xyz OR convert start.txt
                if st_upper == "S0":
                    initial_xyz = Path("initial.xyz")
                    start_txt = Path("start.txt")
                    s0_xyz = esd_dir / "S0.xyz"

                    if not s0_xyz.exists():
                        if initial_xyz.exists():
                            # Prefer optimized initial.xyz if available
                            import shutil
                            shutil.copy2(initial_xyz, s0_xyz)
                            logger.info(f"Copied optimized {initial_xyz} → {s0_xyz}")
                        elif start_txt.exists():
                            # Fallback: convert start.txt to S0.xyz
                            _convert_start_to_xyz(start_txt, s0_xyz)
                            logger.info(f"Converted {start_txt} → {s0_xyz} (no initial.xyz found)")
                        else:
                            raise RuntimeError(
                                f"Cannot create S0.xyz: neither {initial_xyz} nor {start_txt} exist."
                            )

                    # deltaSCF optimization: reuse initial.* files if inputs are identical
                    if esd_modus == "deltascf":
                        initial_out = Path("initial.out")
                        initial_gbw = Path("initial.gbw")
                        initial_hess = Path("initial.hess")

                        if initial_out.exists() and initial_gbw.exists() and initial_hess.exists():
                            # Copy initial files to S0 files to skip redundant Opt+Freq calculation
                            try:
                                s0_out = esd_dir / "S0.out"
                                s0_gbw = esd_dir / "S0.gbw"
                                s0_hess = esd_dir / "S0.hess"
                                shutil.copy2(initial_out, s0_out)
                                shutil.copy2(initial_gbw, s0_gbw)
                                shutil.copy2(initial_hess, s0_hess)
                                logger.info(f"[deltaSCF optimization] Reused initial.{{out,gbw,hess}} → S0.* (skipping redundant Opt+Freq)")

                                # Still need to run TDDFT check job for state identification
                                # Create a separate TDDFT-only input
                                tddft_input = esd_dir / "S0_TDDFT.inp"
                                _create_s0_tddft_check_input(tddft_input, esd_dir, charge, solvent, metals, main_basisset, metal_basisset, config)
                                _update_pal_block(str(tddft_input), cores)

                                # Run only the TDDFT check
                                tddft_output = esd_dir / "S0_TDDFT.out"
                                logger.info(f"Running TDDFT check for S0 state identification in {esd_dir}")
                                if not _run_orca_esd(
                                    tddft_input,
                                    tddft_output,
                                    working_dir=esd_dir,
                                    config=config,
                                ):
                                    raise RuntimeError("ORCA terminated abnormally for S0 TDDFT check")

                                logger.info(f"State S0 calculation completed (Opt+Freq reused, TDDFT check executed)")
                                return
                            except Exception as exc:
                                logger.warning(f"deltaSCF optimization failed: {exc}, running full S0 calculation instead")

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

                # Check if hybrid1 mode and if we need to run two-step calculation
                mode = str(config.get("ESD_modus", "TDDFT")).strip().lower()
                if "|" in mode:
                    mode = mode.split("|")[0].strip()

                is_hybrid1 = mode == "hybrid1"
                # S0 and T1 are lowest states in their spin manifolds - stable without two-step
                needs_two_step = st_upper not in ("S0", "T1")

                if is_hybrid1 and needs_two_step:
                    # Hybrid1 mode for higher excited states (S1+, T2+): run two sequential ORCA calculations
                    # Step 1: TDDFT optimization (first_TDDFT)
                    first_input = esd_dir / f"{st_upper}_first_TDDFT.inp"
                    first_output = esd_dir / f"{st_upper}_first_TDDFT.out"

                    if not first_input.exists():
                        raise RuntimeError(f"Missing hybrid1 first step input: {first_input}")

                    abs_first_input = first_input.resolve()
                    abs_first_output = first_output.resolve()

                    # Update PAL for first step
                    _update_pal_block(str(abs_first_input), cores)

                    logger.info(f"Running ORCA hybrid1 step 1 (TDDFT) for {st_upper} in {esd_dir}")

                    if not _run_orca_esd(
                        abs_first_input,
                        abs_first_output,
                        working_dir=esd_dir,
                        config=config,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for {st_upper} hybrid1 step 1 (TDDFT)"
                        )

                    logger.info(f"Hybrid1 step 1 (TDDFT) completed for {st_upper}")

                    # Step 2: deltaSCF optimization (second)
                    # input_file already points to S1_second.inp from create_state_input
                    abs_input = Path(input_file).resolve()
                    _update_pal_block(str(abs_input), cores)

                    # Use consistent output naming: S1_second.inp -> S1_second.out
                    output_file = abs_input.with_suffix('.out')
                    abs_output = output_file.resolve()

                    logger.info(f"Running ORCA hybrid1 step 2 (deltaSCF) for {st_upper} in {esd_dir}")

                    if not _run_orca_esd(
                        abs_input,
                        abs_output,
                        working_dir=esd_dir,
                        config=config,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for {st_upper} hybrid1 step 2 (deltaSCF)"
                        )

                    logger.info(f"Hybrid1 step 2 (deltaSCF) completed for {st_upper}")
                else:
                    # Standard single-step calculation (TDDFT or deltaSCF)
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

                    # Run ORCA in ESD directory using working_dir parameter (no os.chdir needed)
                    if not _run_orca_esd(
                        abs_input,
                        abs_output,
                        working_dir=esd_dir,
                        config=config,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for {st_upper} state"
                        )

                logger.info(f"State {st_upper} calculation completed")

                # If S0 completed and initial.xyz doesn't exist, copy S0.xyz to initial.xyz
                # This allows ox/red jobs to use the ESD S0 geometry
                if st_upper == "S0":
                    initial_xyz = Path("initial.xyz")
                    s0_xyz = esd_dir / "S0.xyz"
                    if not initial_xyz.exists() and s0_xyz.exists():
                        import shutil
                        shutil.copy2(s0_xyz, initial_xyz)
                        logger.info(f"Copied {s0_xyz} → {initial_xyz} for ox/red jobs")

            return work

        def make_tddft_check_work(st: str) -> Callable[[int], None]:
            """Create work function for TDDFT check on optimized singlet states."""
            def work(cores: int) -> None:
                st_upper = st.upper()
                tddft_input = esd_dir / f"{st_upper}_TDDFT.inp"
                tddft_output = esd_dir / f"{st_upper}_TDDFT.out"

                if tddft_output.exists():
                    try:
                        with tddft_output.open("r", encoding="utf-8", errors="replace") as f:
                            content = f.read()
                        if "ORCA TERMINATED NORMALLY" in content:
                            logger.info(f"{st_upper} TDDFT check skipped (ORCA TERMINATED NORMALLY in {tddft_output.name})")
                            return
                    except Exception:
                        pass

                if not tddft_input.exists():
                    raise RuntimeError(f"Missing TDDFT check input: {tddft_input}")

                _update_pal_block(str(tddft_input.resolve()), cores)

                if not _run_orca_esd(
                    tddft_input.resolve(),
                    tddft_output.resolve(),
                    working_dir=esd_dir,
                    config=config,
                ):
                    raise RuntimeError(
                        f"ORCA terminated abnormally for {st_upper} TDDFT check"
                    )

                logger.info(f"TDDFT check for {st_upper} completed")

            return work

        # Allow ESD jobs to run in parallel with oxidation/reduction
        # cores_min based on total_cores / max_jobs (e.g., 72/6=12, 12/4=3)
        min_cores = max(2, manager.total_cores // manager.max_jobs)
        half_cores = max(min_cores, manager.total_cores // 2)

        manager.add_job(
            WorkflowJob(
                job_id=f"esd_{state_upper}",
                work=make_state_work(state_upper),
                description=f"ESD {state_upper} optimization",
                dependencies=deps,
                cores_min=min_cores,
                cores_optimal=half_cores,
                cores_max=manager.total_cores,
            )
        )

        if esd_modus == "deltascf" and state_upper.startswith("S") and state_upper != "S0":
            manager.add_job(
                WorkflowJob(
                    job_id=f"esd_{state_upper}_tddft",
                    work=make_tddft_check_work(state_upper),
                    description=f"TDDFT check {state_upper}",
                    dependencies={f"esd_{state_upper}"},
                    cores_min=min_cores,
                    cores_optimal=half_cores,
                    cores_max=manager.total_cores,
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
    # Parse TROOTSSL values from config
    trootssl_raw = config.get('TROOTSSL', '0')

    # Handle list/tuple directly
    if isinstance(trootssl_raw, (list, tuple)):
        trootssl_values = [int(x) for x in trootssl_raw]
    else:
        # Handle string format
        trootssl_str = str(trootssl_raw).strip()
        # Remove brackets if present (e.g., "['-1', '0', '1']" or "[-1, 0, 1]")
        trootssl_str = trootssl_str.strip('[]')
        if ',' in trootssl_str:
            # Split and clean each value (remove quotes)
            trootssl_values = [int(x.strip().strip("'\"")) for x in trootssl_str.split(',')]
        else:
            trootssl_values = [int(trootssl_str.strip("'\""))]

    for isc in iscs:
        initial_state, final_state = isc.split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()

        # ISC depends on both initial and final states
        deps = {f"esd_{initial_state}", f"esd_{final_state}"}

        # Create a separate job for each TROOTSSL value
        for trootssl in trootssl_values:
            ms_suffix = _format_ms_suffix(trootssl)
            job_id = f"esd_isc_{initial_state}_{final_state}_{ms_suffix}"

            def make_isc_work(isc_pair: str, trootssl_val: int) -> Callable[[int], None]:
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
                        trootssl=trootssl_val,
                    )

                    # Convert to absolute path before any chdir operations
                    abs_input = Path(input_file).resolve()

                    # Update PAL in input file (use absolute path)
                    _update_pal_block(str(abs_input), cores)

                    # Determine output file name (with TROOTSSL suffix)
                    init_st, fin_st = isc_pair.split(">")
                    init_st = init_st.strip().upper()
                    fin_st = fin_st.strip().upper()
                    ms_sfx = _format_ms_suffix(trootssl_val)
                    output_file = esd_dir / f"{init_st}_{fin_st}_ISC_{ms_sfx}.out"

                    logger.info(f"Running ORCA for ISC {isc_pair} (Ms={trootssl_val}) in {esd_dir}")

                    # Run ORCA in ESD directory using working_dir parameter (no os.chdir needed)
                    scratch_token = Path("scratch") / f"ISC_{init_st}_{fin_st}_{ms_sfx}"
                    if not _run_orca_esd(
                        abs_input,
                        output_file.resolve(),
                        scratch_subdir=scratch_token,
                        working_dir=esd_dir,
                        config=config,
                    ):
                        raise RuntimeError(
                            f"ORCA terminated abnormally for ISC {isc_pair} (Ms={trootssl_val})"
                        )

                    logger.info(f"ISC {isc_pair} (Ms={trootssl_val}) calculation completed")

                return work

            # Allow ISC jobs to run in parallel with other jobs
            # cores_min based on total_cores / max_jobs (e.g., 72/6=12, 12/4=3)
            min_cores = max(2, manager.total_cores // manager.max_jobs)
            half_cores = max(min_cores, manager.total_cores // 2)

            # Format description with sign
            ms_str = f"{trootssl:+d}" if trootssl != 0 else "0"
            manager.add_job(
                WorkflowJob(
                    job_id=job_id,
                    work=make_isc_work(isc, trootssl),
                    description=f"ISC {initial_state}→{final_state} (Ms={ms_str})",
                    dependencies=deps,
                    cores_min=min_cores,
                    cores_optimal=half_cores,
                    cores_max=manager.total_cores,
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
    unsupported_ics: list[str] = []

    for ic in ics:
        initial_state, final_state = ic.split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()

        # ORCA IC support: Transitions to S1, T1, or S0 (fluorescence only)
        # Examples: S2→S1, S3→S1, S1→S0 (fluorescence), T2→T1, T3→T1
        # NOT allowed: S2→S0, S3→S0 (only S1→S0 for fluorescence)

        # Allow Sn→S1 (IC to S1) or S1→S0 (fluorescence)
        allowed_sn_to_s1 = initial_state.startswith("S") and final_state == "S1"
        allowed_fluorescence = initial_state == "S1" and final_state == "S0"

        # Allow Tn→T1 (IC to T1, excluding T1→T1)
        allowed_tn_to_t1 = initial_state.startswith("T") and initial_state != "T1" and final_state == "T1"

        if not (allowed_sn_to_s1 or allowed_fluorescence or allowed_tn_to_t1):
            unsupported_ics.append(ic)
            logger.warning(
                "Skipping IC %s: ORCA IC support is limited to Sn→S1, S1→S0 (fluorescence), or Tn→T1; calculation not scheduled.",
                ic,
            )
            continue

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

                # Run ORCA in ESD directory using working_dir parameter (no os.chdir needed)
                scratch_token = Path("scratch") / f"IC_{init_st}_{fin_st}"
                if not _run_orca_esd(
                    abs_input,
                    output_file.resolve(),
                    scratch_subdir=scratch_token,
                    working_dir=esd_dir,
                    config=config,
                ):
                    raise RuntimeError(
                        f"ORCA terminated abnormally for IC {ic_pair}"
                    )

                logger.info(f"IC {ic_pair} calculation completed")

            return work

        # Allow IC jobs to run in parallel with other jobs
        # cores_min based on total_cores / max_jobs (e.g., 72/6=12, 12/4=3)
        min_cores = max(2, manager.total_cores // manager.max_jobs)
        half_cores = max(min_cores, manager.total_cores // 2)

        manager.add_job(
            WorkflowJob(
                job_id=job_id,
                work=make_ic_work(ic),
                description=f"IC {initial_state}→{final_state}",
                dependencies=deps,
                cores_min=min_cores,
                cores_optimal=half_cores,
                cores_max=manager.total_cores,
            )
        )

    if unsupported_ics:
        logger.info(
            "Unsupported IC transitions skipped: %s",
            ", ".join(unsupported_ics),
        )


def _populate_fluor_jobs(
    manager: _WorkflowManager,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Add fluorescence (S1→S0) ESD(FLUOR) job when requested via emission_rates=f."""
    deps = {"esd_S0", "esd_S1"}
    job_id = "esd_fluor_S1_S0"

    def work(cores: int) -> None:
        input_file = create_fluor_input(
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
            initial_state="S1",
            final_state="S0",
        )

        abs_input = Path(input_file).resolve()
        _update_pal_block(str(abs_input), cores)

        output_file = esd_dir / "S1_S0_FLUOR.out"
        logger.info("Running ORCA for fluorescence (S1→S0) in %s", esd_dir)

        scratch_token = Path("scratch") / "FLUOR_S1_S0"
        if not _run_orca_esd(
            abs_input,
            output_file.resolve(),
            scratch_subdir=scratch_token,
            working_dir=esd_dir,
            config=config,
        ):
            raise RuntimeError("ORCA terminated abnormally for fluorescence (S1→S0)")

        logger.info("Fluorescence (S1→S0) calculation completed")

    # cores_min based on total_cores / max_jobs (e.g., 72/6=12, 12/4=3)
    min_cores = max(2, manager.total_cores // manager.max_jobs)
    half_cores = max(min_cores, manager.total_cores // 2)
    manager.add_job(
        WorkflowJob(
            job_id=job_id,
            work=work,
            description="Fluorescence S1→S0",
            dependencies=deps,
            cores_min=min_cores,
            cores_optimal=half_cores,
            cores_max=manager.total_cores,
        )
    )


def _populate_phosp_jobs(
    manager: _WorkflowManager,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Add phosphorescence (T1→S0) ESD(PHOSP) job when requested via emission_rates=p.

    ORCA phosphorescence requires SOC and typically three IROOT jobs (1..3) after SOC splitting.
    """
    deps = {"esd_S0", "esd_T1"}
    job_id = "esd_phosp_T1_S0"

    def work(cores: int) -> None:
        input_file = create_phosp_input(
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
            initial_state="T1",
            final_state="S0",
        )

        abs_input = Path(input_file).resolve()
        _update_pal_block(str(abs_input), cores)

        output_file = esd_dir / "T1_S0_PHOSP.out"
        logger.info("Running ORCA for phosphorescence (T1→S0) in %s", esd_dir)

        scratch_token = Path("scratch") / "PHOSP_T1_S0"
        if not _run_orca_esd(
            abs_input,
            output_file.resolve(),
            scratch_subdir=scratch_token,
            working_dir=esd_dir,
            config=config,
        ):
            raise RuntimeError("ORCA terminated abnormally for phosphorescence (T1→S0)")

        logger.info("Phosphorescence (T1→S0) calculation completed")

    # cores_min based on total_cores / max_jobs (e.g., 72/6=12, 12/4=3)
    min_cores = max(2, manager.total_cores // manager.max_jobs)
    half_cores = max(min_cores, manager.total_cores // 2)
    manager.add_job(
        WorkflowJob(
            job_id=job_id,
            work=work,
            description="Phosphorescence T1→S0",
            dependencies=deps,
            cores_min=min_cores,
            cores_optimal=half_cores,
            cores_max=manager.total_cores,
        )
    )


def _create_s0_tddft_check_input(
    output_path: Path,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> None:
    """Create standalone TDDFT check input for S0 state identification.

    This creates a single-point TDDFT calculation on the optimized S0 geometry
    for identifying excited state energies.
    """
    from delfin.xyz_io import select_rel_and_aux

    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    implicit_solvation = config.get('implicit_solvation_model', '')
    pal = config.get("PAL", 12)
    maxcore = config.get("maxcore", 6000)
    nroots = config.get('ESD_nroots', 15)
    tda_flag = str(config.get('TDA', 'FALSE')).upper()
    esd_maxdim = config.get('ESD_maxdim', None)
    maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))

    _, aux_jk, _ = select_rel_and_aux(metals, config)

    # Build solvation keyword
    from delfin.esd_input_generator import _build_solvation_keyword
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)

    # TDDFT keyword line (RKS for vertical excitations from S0)
    tddft_keywords = [
        functional,
        "RKS",
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
    ]
    if solvation_kw:
        tddft_keywords.append(solvation_kw)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("! " + " ".join(tddft_keywords) + "\n")
        f.write('%base "S0_TDDFT"\n')
        f.write(f"%pal nprocs {pal} end\n")
        f.write(f"%maxcore {maxcore}\n")
        f.write("\n%tddft\n")
        f.write(f"  nroots {nroots}\n")
        f.write(f"  maxdim {maxdim}\n")
        f.write(f"  tda {tda_flag}\n")
        tddft_maxiter = config.get('ESD_TDDFT_maxiter')
        if tddft_maxiter is not None:
            f.write(f"  maxiter {tddft_maxiter}\n")
        f.write("  triplets true\n")
        f.write("end\n")
        f.write("\n")
        f.write(f"* xyzfile {charge} 1 S0.xyz\n")

    logger.info(f"Created standalone TDDFT check input: {output_path}")


def _update_pal_block(input_path: str, cores: int) -> None:
    """Update ALL %pal blocks in ORCA input file with given core count.

    For input files with multiple $new_job sections, this updates every %pal block
    to ensure all jobs use the same number of cores.

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
    replaced_count = 0

    # Replace ALL %pal blocks (not just the first one)
    for idx, line in enumerate(lines):
        if line.strip().startswith('%pal'):
            lines[idx] = pal_line
            replaced_count += 1

    if replaced_count == 0:
        # Insert after other % blocks (first occurrence)
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


def add_esd_jobs_to_scheduler(
    scheduler: Any,
    config: Dict[str, Any],
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    *,
    dependency_job_id: str = "classic_initial",
) -> tuple[bool, List[str], List[str], List[str]]:
    """Add ESD jobs to an existing scheduler (e.g., classic phase scheduler).

    This allows ESD jobs to run in parallel with oxidation/reduction jobs
    after the initial calculation completes, instead of waiting for the entire
    classic phase to finish.

    Args:
        scheduler: GlobalOrcaScheduler to add jobs to
        config: Configuration dictionary
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        dependency_job_id: Job ID that ESD jobs should depend on (default: "classic_initial")

    Returns:
        Tuple of (esd_enabled, states, iscs, ics) for tracking what was added
    """
    esd_enabled, states, iscs, ics = parse_esd_config(config)
    emission_rates = parse_emission_rates(config)

    if not esd_enabled:
        logger.debug("ESD module disabled (ESD_modul=no)")
        return esd_enabled, states, iscs, ics

    if not states and not iscs and not ics and not emission_rates:
        logger.debug("ESD module enabled but no states/ISCs/ICs configured")
        return esd_enabled, states, iscs, ics

    logger.info("Adding ESD jobs to scheduler (will run after %s)", dependency_job_id)

    # Ensure required states exist when emission rates are requested
    states_effective = list(states)
    if "f" in emission_rates:
        state_set = {s.upper() for s in states_effective}
        for required in ("S0", "S1"):
            if required not in state_set:
                states_effective.append(required)
                state_set.add(required)
    if "p" in emission_rates:
        state_set = {s.upper() for s in states_effective}
        for required in ("S0", "T1"):
            if required not in state_set:
                states_effective.append(required)
                state_set.add(required)

    # Setup ESD directory (use absolute path to avoid chdir issues in parallel jobs)
    esd_dir = Path("ESD").resolve()
    setup_esd_directory(esd_dir, states_effective)

    # Get the scheduler's internal manager
    manager = scheduler.manager

    # Temporarily store original dependencies and modify them
    # We need to inject the dependency_job_id into state dependencies
    original_state_deps = {
        "S0": set(),
        "S1": {"esd_S0"},
        "S2": {"esd_S1"},
        "S3": {"esd_S2"},
        "S4": {"esd_S3"},
        "S5": {"esd_S4"},
        "S6": {"esd_S5"},
        "T1": {"esd_S0"},
        "T2": {"esd_T1"},
        "T3": {"esd_T2"},
        "T4": {"esd_T3"},
        "T5": {"esd_T4"},
        "T6": {"esd_T5"},
    }

    # Populate jobs (they will be added to the shared scheduler manager)
    if states_effective:
        _populate_state_jobs(
            manager,
            states_effective,
            esd_dir,
            charge,
            solvent,
            metals,
            main_basisset,
            metal_basisset,
            config,
        )

        # Add dependency on classic_initial only if it will actually run
        # When ESD module is enabled, calc_initial is automatically set to 'no' to avoid redundancy
        # In that case, S0 should start directly from start.txt without depending on classic_initial
        if "esd_S0" in manager._jobs:
            calc_initial = str(config.get('calc_initial', 'no')).strip().lower() == 'yes'

            if calc_initial:
                # classic_initial will run -> add dependency
                manager._jobs["esd_S0"].dependencies.add(dependency_job_id)
                logger.debug("Added dependency: esd_S0 depends on %s", dependency_job_id)
            else:
                # No classic_initial -> esd_S0 will use start.txt directly
                logger.info("Skipping dependency on %s (calc_initial=no)", dependency_job_id)
                logger.info("esd_S0 will use start.txt or existing initial.xyz directly")

    # Check if frequency calculations are enabled (required for ISC/IC)
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')

    if iscs and esd_frequency_enabled:
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
    elif iscs and not esd_frequency_enabled:
        logger.info("ISC calculations skipped (ESD_frequency=no - ZPE not available)")

    if ics and esd_frequency_enabled:
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
    elif ics and not esd_frequency_enabled:
        logger.info("IC calculations skipped (ESD_frequency=no - ZPE not available)")

    if "f" in emission_rates and esd_frequency_enabled:
        _populate_fluor_jobs(
            manager,
            esd_dir,
            charge,
            solvent,
            metals,
            main_basisset,
            metal_basisset,
            config,
        )
    elif "f" in emission_rates and not esd_frequency_enabled:
        logger.info("Fluorescence calculations skipped (ESD_frequency=no - Hessians not available)")

    if "p" in emission_rates and esd_frequency_enabled:
        _populate_phosp_jobs(
            manager,
            esd_dir,
            charge,
            solvent,
            metals,
            main_basisset,
            metal_basisset,
            config,
        )
    elif "p" in emission_rates and not esd_frequency_enabled:
        logger.info("Phosphorescence calculations skipped (ESD_frequency=no - Hessians not available)")

    logger.info("Added %d ESD jobs to scheduler", len([j for j in manager._jobs if j.startswith("esd_")]))

    return esd_enabled, states, iscs, ics


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
    emission_rates = parse_emission_rates(config)

    if not esd_enabled:
        logger.info("ESD module disabled (ESD_modul=no)")
        return WorkflowRunResult()

    if not states and not iscs and not ics and not emission_rates:
        logger.info("ESD module enabled but no states/ISCs/ICs configured")
        return WorkflowRunResult()

    logger.info("Starting ESD module calculations")

    # Run xTB/GOAT optimization if enabled (before ESD calculations)
    _run_xtb_goat_if_needed(config, charge)

    # Setup ESD directory (use absolute path to avoid chdir issues in parallel jobs)
    esd_dir = Path("ESD").resolve()
    states_effective = list(states)
    if "f" in emission_rates:
        state_set = {s.upper() for s in states_effective}
        for required in ("S0", "S1"):
            if required not in state_set:
                states_effective.append(required)
                state_set.add(required)
    if "p" in emission_rates:
        state_set = {s.upper() for s in states_effective}
        for required in ("S0", "T1"):
            if required not in state_set:
                states_effective.append(required)
                state_set.add(required)
    setup_esd_directory(esd_dir, states_effective)

    # Create workflow manager
    manager = _WorkflowManager(config, label="esd")

    # Check if frequency calculations are enabled (required for ISC/IC)
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')

    try:
        # Populate jobs
        if states_effective:
            _populate_state_jobs(
                manager,
                states_effective,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )

        if iscs and esd_frequency_enabled:
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
        elif iscs and not esd_frequency_enabled:
            logger.info("ISC calculations skipped (ESD_frequency=no - ZPE not available)")

        if ics and esd_frequency_enabled:
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
        elif ics and not esd_frequency_enabled:
            logger.info("IC calculations skipped (ESD_frequency=no - ZPE not available)")

        if "f" in emission_rates and esd_frequency_enabled:
            _populate_fluor_jobs(
                manager,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )
        elif "f" in emission_rates and not esd_frequency_enabled:
            logger.info("Fluorescence calculations skipped (ESD_frequency=no - Hessians not available)")

        if "p" in emission_rates and esd_frequency_enabled:
            _populate_phosp_jobs(
                manager,
                esd_dir,
                charge,
                solvent,
                metals,
                main_basisset,
                metal_basisset,
                config,
            )
        elif "p" in emission_rates and not esd_frequency_enabled:
            logger.info("Phosphorescence calculations skipped (ESD_frequency=no - Hessians not available)")

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
