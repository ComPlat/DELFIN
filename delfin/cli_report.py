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
from delfin.copy_helpers import extract_preferred_spin
from delfin.utils import _TM_LIST

logger = get_logger(__name__)


def _extract_energy_for_state(
    state_label: str,
    out_filename: str,
    energies: Dict[str, Optional[float]]
) -> None:
    """Extract Gibbs energy for a specific state and update energies dict.

    Args:
        state_label: Label for the state (e.g., '0', '+1', '-1')
        out_filename: Output file to extract energy from
        energies: Dictionary to update with extracted energy
    """
    if os.path.exists(out_filename):
        g = find_gibbs_energy(out_filename)
        if g is not None:
            energies[state_label] = g
            logger.info(f"Extracted G({state_label}) = {g:.6f} Eh from {out_filename}")
        else:
            logger.warning(f"Could not extract Gibbs energy from {out_filename}")
    else:
        logger.warning(f"File {out_filename} not found")


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
        _extract_energy_for_state('0', 'initial.out', energies)

    # Oxidation steps
    for step in ox_steps:
        _extract_energy_for_state(f'+{step}', f'ox_step_{step}.out', energies)

    # Reduction steps
    for step in red_steps:
        _extract_energy_for_state(f'-{step}', f'red_step_{step}.out', energies)

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
                for metal in _TM_LIST:
                    if f' {metal} ' in content or f' {metal}\t' in content:
                        if metal not in metals:
                            metals.append(metal)
        except Exception:
            pass

    # Refresh multiplicity label from OCCUPIER, if available
    _apply_occuper_spin_metadata(config, Path('.'))

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


def _apply_occuper_spin_metadata(config: Dict[str, Any], workspace: Path) -> None:
    occ_folder = workspace / "initial_OCCUPIER"
    mult, bs = extract_preferred_spin(occ_folder)
    if mult is None:
        return
    config['multiplicity_0'] = mult
    if bs:
        config['_multiplicity_display'] = f"{mult} (BS {bs})"
    else:
        config['_multiplicity_display'] = str(mult)
