# cli_report.py
# Functionality for --report flag: Recompute potentials from existing outputs

import os
from pathlib import Path
from typing import Dict, Any, Optional

from delfin.common.logging import get_logger
from delfin.energies import find_gibbs_energy
from delfin.config import get_E_ref
from delfin.cli_calculations import calculate_redox_potentials, select_final_potentials
from delfin.reporting.delfin_reports import generate_summary_report_DELFIN

logger = get_logger(__name__)


def extract_energies_from_outputs(config: Dict[str, Any]) -> Dict[str, Optional[float]]:
    """Extract Gibbs energies from existing ORCA output files.

    Args:
        config: Configuration dictionary

    Returns:
        Dictionary mapping charge states to Gibbs energies
    """
    energies = {}

    # Determine which files to check based on config
    calc_initial = str(config.get('calc_initial', 'yes')).lower() == 'yes'
    oxidation_steps = config.get('oxidation_steps', '')
    reduction_steps = config.get('reduction_steps', '')

    # Parse steps
    ox_steps = [int(s.strip()) for s in str(oxidation_steps).split(',') if s.strip().isdigit()]
    red_steps = [int(s.strip()) for s in str(reduction_steps).split(',') if s.strip().isdigit()]

    # Initial state
    if calc_initial:
        initial_out = 'initial.out'
        if os.path.exists(initial_out):
            g = find_gibbs_energy(initial_out)
            if g is not None:
                energies['0'] = g
                logger.info(f"Extracted G(initial) = {g:.6f} Eh from {initial_out}")
            else:
                logger.warning(f"Could not extract Gibbs energy from {initial_out}")
        else:
            logger.warning(f"File {initial_out} not found")

    # Oxidation steps
    for step in ox_steps:
        ox_out = f'ox_step_{step}.out'
        if os.path.exists(ox_out):
            g = find_gibbs_energy(ox_out)
            if g is not None:
                energies[f'+{step}'] = g
                logger.info(f"Extracted G(+{step}) = {g:.6f} Eh from {ox_out}")
            else:
                logger.warning(f"Could not extract Gibbs energy from {ox_out}")
        else:
            logger.warning(f"File {ox_out} not found")

    # Reduction steps
    for step in red_steps:
        red_out = f'red_step_{step}.out'
        if os.path.exists(red_out):
            g = find_gibbs_energy(red_out)
            if g is not None:
                energies[f'-{step}'] = g
                logger.info(f"Extracted G(-{step}) = {g:.6f} Eh from {red_out}")
            else:
                logger.warning(f"Could not extract Gibbs energy from {red_out}")
        else:
            logger.warning(f"File {red_out} not found")

    return energies


def run_report_mode(config: Dict[str, Any]) -> int:
    """Run report mode: recompute potentials from existing outputs.

    Args:
        config: Configuration dictionary

    Returns:
        Exit code (0 for success, non-zero for error)
    """
    logger.info("="*70)
    logger.info("DELFIN --report MODE")
    logger.info("Recomputing redox potentials from existing output files")
    logger.info("="*70)

    # Extract energies from existing output files
    free_gibbs_energies = extract_energies_from_outputs(config)

    if not free_gibbs_energies:
        logger.error("No Gibbs energies found in output files!")
        logger.error("Make sure you have run DELFIN calculations first.")
        return 1

    logger.info(f"Found {len(free_gibbs_energies)} Gibbs energies")

    # Get E_ref (will use ONIOM-adjusted value if applicable)
    E_ref = get_E_ref(config)
    logger.info(f"Using E_ref = {E_ref:.3f} V")

    # Calculate redox potentials
    m1_avg, m2_step, m3_mix, use_flags = calculate_redox_potentials(
        config, free_gibbs_energies, E_ref
    )

    # Select final potentials
    E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3 = select_final_potentials(
        m1_avg, m2_step, m3_mix, use_flags
    )

    logger.info("="*70)
    logger.info("RECOMPUTED REDOX POTENTIALS")
    logger.info("="*70)
    if E_ox is not None:
        logger.info(f"E_ox   = {E_ox:+7.3f} V vs. Fc+/Fc")
    if E_ox_2 is not None:
        logger.info(f"E_ox_2 = {E_ox_2:+7.3f} V vs. Fc+/Fc")
    if E_ox_3 is not None:
        logger.info(f"E_ox_3 = {E_ox_3:+7.3f} V vs. Fc+/Fc")
    if E_red is not None:
        logger.info(f"E_red  = {E_red:+7.3f} V vs. Fc+/Fc")
    if E_red_2 is not None:
        logger.info(f"E_red_2= {E_red_2:+7.3f} V vs. Fc+/Fc")
    if E_red_3 is not None:
        logger.info(f"E_red_3= {E_red_3:+7.3f} V vs. Fc+/Fc")
    logger.info("="*70)

    # Generate DELFIN.txt report
    # Extract additional info from config
    charge = config.get('charge', 0)
    multiplicity = config.get('multiplicity_0', 1)
    solvent = config.get('solvent', '')
    NAME = config.get('NAME', '')
    main_basisset = config.get('main_basisset', '')
    metal_basisset = config.get('metal_basisset', '')

    # Read metals from initial.out if available
    metals = []
    if os.path.exists('initial.out'):
        try:
            with open('initial.out', 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read(5000)
                # Simple heuristic: look for metal atoms in coordinates
                metal_symbols = ['Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mn', 'Cr', 'V', 'Ti', 'Ru', 'Rh', 'Pd', 'Ag', 'Os', 'Ir', 'Pt', 'Au']
                for metal in metal_symbols:
                    if f' {metal} ' in content or f' {metal}\t' in content:
                        if metal not in metals:
                            metals.append(metal)
        except Exception:
            pass

    # Generate report
    logger.info("Generating DELFIN.txt report...")

    try:
        generate_summary_report_DELFIN(
            charge=charge,
            multiplicity=multiplicity,
            solvent=solvent,
            E_ox=E_ox,
            E_ox_2=E_ox_2,
            E_ox_3=E_ox_3,
            E_red=E_red,
            E_red_2=E_red_2,
            E_red_3=E_red_3,
            E_00_t1=None,
            E_00_s1=None,
            metals=metals,
            metal_basisset=metal_basisset,
            NAME=NAME,
            main_basisset=main_basisset,
            config=config,
            duration=0.0,  # No duration for report mode
            E_ref=E_ref,
            esd_summary=None,
            output_dir=Path('.')
        )
        logger.info("DELFIN.txt updated successfully!")
    except Exception as e:
        logger.error(f"Failed to generate DELFIN.txt: {e}")
        return 1

    return 0
