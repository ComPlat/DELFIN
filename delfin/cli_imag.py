"""
cli_imag.py

Handles `delfin --imag` mode:
- Runs IMAG elimination on existing .out/.hess files
- Automatically generates report with updated potentials
"""

from pathlib import Path
from typing import Dict, Any

from delfin.common.logging import get_logger
from delfin.global_manager import get_global_manager
from delfin.global_scheduler import GlobalOrcaScheduler
from delfin.imag import run_IMAG
from delfin.parallel_classic_manually import WorkflowJob
from delfin.utils import search_transition_metals, set_main_basisset

logger = get_logger(__name__)


def _extract_charge_mult_from_inp(inp_file: Path) -> tuple[int, int]:
    """Extract charge and multiplicity from ORCA .inp file.

    Returns:
        (charge, multiplicity)
    """
    try:
        with open(inp_file, 'r') as f:
            for line in f:
                if line.startswith('* xyz'):
                    parts = line.split()
                    if len(parts) >= 4:
                        charge = int(parts[2])
                        mult = int(parts[3])
                        return (charge, mult)
    except Exception as e:
        logger.warning("Could not extract charge/mult from %s: %s", inp_file, e)
    return (0, 1)  # fallback


def run_imag_mode(config: Dict[str, Any], control_file_path: Path) -> int:
    """
    Run IMAG elimination on existing output files, then generate report.

    Args:
        config: Parsed CONTROL.txt configuration
        control_file_path: Path to CONTROL.txt

    Returns:
        Exit code (0 = success, non-zero = error)
    """
    logger.info("=" * 70)
    logger.info("IMAG ELIMINATION MODE")
    logger.info("=" * 70)

    # Check IMAG is enabled
    if str(config.get("IMAG", "no")).lower() != "yes":
        logger.error("IMAG=yes must be set in CONTROL.txt to use --imag mode")
        return 2

    imag_scope = str(config.get("IMAG_scope", "initial")).lower()
    allow_imaginary_freq = config.get("allow_imaginary_freq", 0)

    logger.info("IMAG_scope: %s", imag_scope)
    logger.info("allow_imaginary_freq: %s", allow_imaginary_freq)
    logger.info("=" * 70)

    # Get basic settings
    charge = int(config.get("charge", 0))
    solvent = config.get("solvent", "")

    # Get input file to detect metals
    input_file = config.get("input_file", "input.txt")
    input_path = control_file_path.parent / input_file

    if not input_path.exists():
        logger.error("Input file not found: %s", input_path)
        return 2

    # Detect metals
    metals = search_transition_metals(str(input_path))
    if metals:
        logger.info("Found transition metals: %s", ", ".join(metals))
    else:
        logger.info("No transition metals found")

    # Get basis sets
    main_basisset, metal_basisset = set_main_basisset(metals, config)

    # Determine which steps to process
    calc_initial = str(config.get("calc_initial", "yes")).lower() == "yes"
    ox_steps_raw = str(config.get("oxidation_steps", "")).strip()
    red_steps_raw = str(config.get("reduction_steps", "")).strip()

    ox_steps = [int(x) for x in ox_steps_raw.split(",") if x.strip().isdigit()] if ox_steps_raw else []
    red_steps = [int(x) for x in red_steps_raw.split(",") if x.strip().isdigit()] if red_steps_raw else []

    # Build list of steps to process
    steps_to_process = []

    if calc_initial and imag_scope in ["initial", "all"]:
        steps_to_process.append({
            "name": "initial",
            "inp_file": "initial.inp",
            "out_file": "initial.out",
            "hess_file": "initial.hess",
            "additions": str(config.get("additions_0", "")).strip(),
        })

    if imag_scope == "all":
        # Add oxidation steps
        for idx in ox_steps:
            add_key = f"additions_ox{idx}"
            steps_to_process.append({
                "name": f"ox_step_{idx}",
                "inp_file": f"ox_step_{idx}.inp",
                "out_file": f"ox_step_{idx}.out",
                "hess_file": f"ox_step_{idx}.hess",
                "additions": str(config.get(add_key, "")).strip(),
            })

        # Add reduction steps
        for idx in red_steps:
            add_key = f"additions_red{idx}"
            steps_to_process.append({
                "name": f"red_step_{idx}",
                "inp_file": f"red_step_{idx}.inp",
                "out_file": f"red_step_{idx}.out",
                "hess_file": f"red_step_{idx}.hess",
                "additions": str(config.get(add_key, "")).strip(),
            })

    if not steps_to_process:
        logger.error("No steps to process! Check CONTROL.txt settings.")
        return 2

    logger.info("Will process %d step(s):", len(steps_to_process))
    for step in steps_to_process:
        logger.info("  - %s", step["name"])

    # Process each step
    success_count = 0
    workspace = control_file_path.parent
    metals_list = list(metals) if metals else []

    # Ensure the global job manager is initialized so we can attach to the shared pool
    try:
        get_global_manager().initialize(config)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not initialize global job manager for IMAG mode: %s", exc)

    try:
        scheduler = GlobalOrcaScheduler(config, label="imag_elimination")
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to create global scheduler for IMAG mode: %s", exc)
        return 1

    jobs_registered = 0
    manager = scheduler.manager

    try:
        for step in steps_to_process:
            logger.info("")
            logger.info("=" * 70)
            logger.info("Processing: %s", step["name"])
            logger.info("=" * 70)

            inp_file = workspace / step["inp_file"]
            out_file = workspace / step["out_file"]
            hess_file = workspace / step["hess_file"]

            # Check files exist
            if not inp_file.exists():
                logger.warning("Skipping %s: %s not found", step["name"], inp_file)
                continue

            if not out_file.exists():
                logger.warning("Skipping %s: %s not found", step["name"], out_file)
                continue

            if not hess_file.exists():
                logger.warning("Skipping %s: %s not found", step["name"], hess_file)
                continue

            # Extract charge and multiplicity from .inp file
            step_charge, step_mult = _extract_charge_mult_from_inp(inp_file)
            logger.info("Charge: %d, Multiplicity: %d (from %s)", step_charge, step_mult, inp_file.name)

            additions = step["additions"]

            hess_base = str(hess_file).replace(".hess", "")
            job_id = f"imag_{step['name']}"
            additions_payload = additions

            cores_min, cores_opt, cores_max = scheduler.manager.derive_core_bounds(
                hint=f"IMAG {step['name']}"
            )

            def make_work(step_name: str,
                          out_path: Path,
                          hess: str,
                          step_charge_val: int,
                          step_mult_val: int,
                          additions_val,
                          metals_val: list[str]):
                def _work(cores: int) -> None:
                    logger.info("[IMAG] %s → start (cores=%d)", step_name, cores)
                    run_IMAG(
                        str(out_path),
                        hess,
                        step_charge_val,
                        step_mult_val,
                        solvent,
                        metals_val,
                        config,
                        main_basisset,
                        metal_basisset,
                        additions_val,
                        step_name=step_name,
                        source_input=str(inp_file),
                        pal_override=cores,
                        maxcore_override=None,
                    )
                    logger.info("✅ %s completed successfully", step_name)

                return _work

            work_fn = make_work(
                step["name"],
                out_file,
                hess_base,
                step_charge,
                step_mult,
                additions_payload,
                metals_list.copy(),
            )

            scheduler.add_job(
                WorkflowJob(
                    job_id=job_id,
                    work=work_fn,
                    description=f"IMAG elimination for {step['name']}",
                    dependencies=set(),
                    cores_min=cores_min,
                    cores_optimal=cores_opt,
                    cores_max=cores_max,
                )
            )
            jobs_registered += 1

        if jobs_registered:
            try:
                result = scheduler.run()
            except Exception as exc:  # noqa: BLE001
                logger.error("IMAG scheduler execution failed: %s", exc, exc_info=True)
                result = None
            else:
                success_count = len(result.completed)
                if result.failed:
                    for job_id, reason in result.failed.items():
                        logger.error("❌ IMAG job %s failed: %s", job_id, reason)
                if result.skipped:
                    for job_id, deps in result.skipped.items():
                        dep_desc = ", ".join(deps) if deps else "unknown dependencies"
                        logger.warning("⚠️ IMAG job %s skipped (missing %s)", job_id, dep_desc)
        else:
            logger.warning("No IMAG jobs were scheduled; nothing to execute.")

    finally:
        scheduler.shutdown()

    logger.info("")
    logger.info("=" * 70)
    logger.info("IMAG elimination completed: %d/%d steps processed", success_count, len(steps_to_process))
    logger.info("=" * 70)

    # Now run report mode to recalculate potentials
    if success_count > 0:
        logger.info("")
        logger.info("=" * 70)
        logger.info("Generating report with updated potentials...")
        logger.info("=" * 70)

        from .cli_report import run_report_mode

        return run_report_mode(config)
    else:
        logger.error("No IMAG eliminations were successful. Skipping report generation.")
        return 1
