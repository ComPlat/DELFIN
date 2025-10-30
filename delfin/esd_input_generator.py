"""ORCA input file generator for ESD module (excited state dynamics).

This module generates ORCA input files for:
- Electronic states (S0, S1, T1, T2)
- Intersystem crossings (ISCs)
- Internal conversions (ICs)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# Conversion factor: Hartree to cm^-1
HARTREE_TO_CM1 = 219474.63


def calculate_dele_cm1(state1_file: str, state2_file: str) -> Optional[float]:
    """Calculate adiabatic energy difference (DELE) between two states.

    DELE = E(initial_state) - E(final_state)
    Both energies evaluated at their respective optimized geometries.

    Args:
        state1_file: Path to initial state .out file
        state2_file: Path to final state .out file

    Returns:
        DELE in cm^-1, or None if energies cannot be extracted
    """
    from delfin.energies import find_electronic_energy
    from pathlib import Path

    # Check if files exist
    if not Path(state1_file).exists() or not Path(state2_file).exists():
        logger.warning(f"Cannot calculate DELE: missing {state1_file} or {state2_file}")
        return None

    # Extract electronic energies
    e1 = find_electronic_energy(state1_file)
    e2 = find_electronic_energy(state2_file)

    if e1 is None or e2 is None:
        logger.warning(f"Cannot calculate DELE: failed to extract energies from outputs")
        return None

    # Calculate DELE in cm^-1
    dele_hartree = e1 - e2
    dele_cm1 = dele_hartree * HARTREE_TO_CM1

    logger.info(f"Calculated DELE: {dele_cm1:.2f} cm⁻¹ ({e1:.6f} - {e2:.6f} Eh)")

    return dele_cm1


def create_state_input(
    state: str,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> str:
    """Generate ORCA input file for electronic state calculation.

    Args:
        state: State identifier (S0, S1, T1, T2)
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary

    Returns:
        Path to generated input file
    """
    state_upper = state.upper()
    input_file = esd_dir / f"{state_upper}.inp"

    # Determine multiplicity based on state type
    # Singlet states (S0, S1, etc.): M = 1
    # Triplet states (T1, T2, etc.): M = 3
    if state_upper.startswith('T'):
        multiplicity = 3  # Triplet states
    else:
        multiplicity = 1  # Singlet states

    # Determine source geometry
    if state_upper == "S0":
        xyz_file = "initial.xyz"
        moinp_gbw = None
        use_deltascf = False
    elif state_upper == "S1":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "T1":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = False
    elif state_upper == "T2":
        xyz_file = "S0.xyz"
        moinp_gbw = "T1.gbw"
        use_deltascf = True
    else:
        raise ValueError(f"Unknown state: {state}")

    # Functional and basis set
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')

    # Solvation
    implicit_solvation = config.get('implicit_solvation_model', 'CPCM')

    # Geometry optimization token from CONTROL (fall back to OPT)
    geom_token = str(config.get('geom_opt', 'OPT')).strip()

    # Build simple keyword line
    # S0 is closed-shell (RKS), all other states need UKS
    scf_type = "RKS" if state_upper == "S0" else "UKS"

    keywords = [
        functional,
        scf_type,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
        f"{implicit_solvation}({solvent})",
    ]

    if geom_token:
        keywords.append(geom_token)
    keywords.append("FREQ")

    if use_deltascf:
        keywords.append("deltaSCF")

    if moinp_gbw:
        keywords.append("NODIIS")
        keywords.append("MOREAD")

    simple_line = "! " + " ".join(keywords)

    # Blocks
    blocks = []

    # Base block
    blocks.append(f'%base "{state_upper}"')

    # MO input
    if moinp_gbw:
        blocks.append(f'%moinp "{moinp_gbw}"')

    # PAL
    pal = config.get('PAL', 12)
    blocks.append(f"%pal nprocs {pal} end")

    # Maxcore
    maxcore = config.get('maxcore', 6000)
    blocks.append(f"%maxcore {maxcore}")

    # SCF settings for deltaSCF
    if use_deltascf:
        domom = str(config.get('deltaSCF_DOMOM', 'true')).lower()  # Changed default to true
        pmom = str(config.get('deltaSCF_PMOM', 'true')).lower()
        keepinitialref = str(config.get('deltaSCF_keepinitialref', 'true')).lower()
        soscfhessup = config.get('deltaSCF_SOSCFHESSUP', 'LSR1')  # Changed to LSR1 (better for excited states)

        scf_block = [
            "%scf",
            f"  DOMOM {domom}",
            f"  pmom {pmom}",
            f"  keepinitialref {keepinitialref}",
        ]

        # State-specific orbital configurations
        if state_upper == "S1":
            scf_block.extend([
                "  alphaconf 0,1",
                "  betaconf 0",
            ])
        elif state_upper == "T2":
            scf_block.extend([
                "  alphaconf 0,1",
                "  betaconf 0",
            ])

        scf_block.append(f"  SOSCFHESSUP {soscfhessup}")
        scf_block.append("end")
        blocks.append("\n".join(scf_block))

    # Geometry - read from input.txt or xyz file
    if xyz_file == "initial.xyz":
        # For S0: read from input.txt in main directory
        xyz_path = Path("input.txt")
        skip_lines = 0  # input.txt has no header
    else:
        # For S1, T1, T2: read from ESD directory (XYZ format with header)
        xyz_path = esd_dir / xyz_file
        skip_lines = 2  # Skip atom count and comment line

    # Read coordinates
    try:
        with open(xyz_path, 'r', encoding='utf-8') as f:
            all_lines = f.readlines()
            coord_lines = all_lines[skip_lines:]  # Skip header if needed
    except FileNotFoundError:
        logger.error(f"Coordinate file not found: {xyz_path}")
        raise

    # Write input file
    with open(input_file, 'w', encoding='utf-8') as f:
        f.write(simple_line + "\n")
        for block in blocks:
            f.write(block + "\n")
        f.write("\n")
        f.write(f"* xyz {charge} {multiplicity}\n")
        for line in coord_lines:
            f.write(line)
        f.write("*\n")

        # Add TDDFT check job for S0 to identify excited states
        if state_upper == "S0":
            f.write("\n")
            f.write("#==========================================\n")
            f.write("# TDDFT Check: Identify excited states\n")
            f.write("#==========================================\n")
            f.write("\n")
            f.write("$new_job\n")

            # TDDFT keyword line (RKS for vertical excitations from S0)
            tddft_keywords = [
                functional,
                "RKS",
                main_basisset,
                disp_corr,
                ri_jkx,
                aux_jk,
                f"{implicit_solvation}({solvent})",
            ]
            f.write("! " + " ".join(tddft_keywords) + "\n")

            # Base block for TDDFT check
            f.write(f'%base "S0_TDDFT_check"\n')

            # PAL and maxcore
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # TDDFT block for both singlets and triplets
            nroots = config.get('NROOTS', 15)
            tda_flag = str(config.get('TDA', 'FALSE')).upper()
            maxdim = max(5, int(nroots / 2))  # Davidson subspace dimension
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            f.write("  triplets true\n")
            f.write("end\n")

            # Geometry reference
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 S0.xyz\n")

            logger.info(f"Added TDDFT check job to S0 input for state identification")

    logger.info(f"Created ESD state input: {input_file}")
    return str(input_file)


def create_isc_input(
    isc_pair: str,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> str:
    """Generate ORCA input file for intersystem crossing (ISC) calculation.

    Args:
        isc_pair: ISC transition (e.g., "S1>T1")
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary

    Returns:
        Path to generated input file
    """
    initial_state, final_state = isc_pair.split(">")
    initial_state = initial_state.strip().upper()
    final_state = final_state.strip().upper()

    job_name = f"{initial_state}_{final_state}_ISC"
    input_file = esd_dir / f"{job_name}.inp"

    # Determine source geometry (use optimized geometry of initial state)
    xyz_file = f"{initial_state}.xyz"

    # Calculate adiabatic energy difference (DELE) for ISC
    # DELE = E(initial) - E(final) in cm^-1
    dele = calculate_dele_cm1(
        str(esd_dir / f"{initial_state}.out"),
        str(esd_dir / f"{final_state}.out"),
    )

    # Build input
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')
    implicit_solvation = config.get('implicit_solvation_model', 'CPCM')

    # Simple keyword line (no RKS/UKS flag - let ORCA decide based on multiplicity)
    keywords = [
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
        f"{implicit_solvation}({solvent})",
        "ESD(ISC)",
    ]

    simple_line = "! " + " ".join(keywords)

    # Blocks
    blocks = []

    # Base
    blocks.append(f'%base "{job_name}"')

    # TDDFT block (aligned with reference layout)
    nroots = config.get('ESD_ISC_NROOTS', config.get('NROOTS', 10))  # Increased default to 10
    trootssl = str(config.get('TROOTSSL', '0')).strip()
    dosoc_flag = "TRUE"
    tddft_block = [
        f"%TDDFT  NROOTS  {int(nroots):>2}",
        "        SROOT   1",
        "        TROOT   1",
        f"        TROOTSSL {trootssl}",
        f"        DOSOC   {dosoc_flag}",
        "END",
    ]
    blocks.append("\n".join(tddft_block))

    # ESD block
    temperature = config.get('temperature', 298.15)
    doht_flag = str(config.get('DOHT', 'TRUE')).upper()
    esd_block = [
        "%ESD",
        f'  ISCISHESS       "{initial_state}.hess"',
        f'  ISCFSHESS       "{final_state}.hess"',
        "  USEJ            TRUE",
        f"  DOHT            {doht_flag}",
        f"  TEMP            {temperature}",
    ]
    if dele is not None:
        esd_block.append(f"  DELE            {int(dele)}")
    esd_block.append("END")
    blocks.append("\n".join(esd_block))

    # PAL and maxcore
    pal = config.get('PAL', 12)
    maxcore = config.get('maxcore', 6000)
    blocks.append(f"%pal nprocs {pal} end")
    blocks.append(f"%maxcore {maxcore}")

    # Geometry - read coordinates (XYZ format with header)
    xyz_path = esd_dir / xyz_file
    try:
        with open(xyz_path, 'r', encoding='utf-8') as f:
            all_lines = f.readlines()
            coord_lines = all_lines[2:]  # Skip atom count and comment line
    except FileNotFoundError:
        logger.error(f"Coordinate file not found: {xyz_path}")
        raise

    # Write input file
    with open(input_file, 'w', encoding='utf-8') as f:
        f.write(simple_line + "\n")
        for block in blocks:
            f.write(block + "\n")
        f.write("\n")
        f.write(f"* xyz {charge} 1\n")
        for line in coord_lines:
            f.write(line)
        f.write("*\n")

    logger.info(f"Created ISC input: {input_file}")
    return str(input_file)


def create_ic_input(
    ic_pair: str,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> str:
    """Generate ORCA input file for internal conversion (IC) calculation.

    Args:
        ic_pair: IC transition (e.g., "S1>S0")
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
        config: Configuration dictionary

    Returns:
        Path to generated input file
    """
    initial_state, final_state = ic_pair.split(">")
    initial_state = initial_state.strip().upper()
    final_state = final_state.strip().upper()

    job_name = f"{initial_state}_{final_state}_IC"
    input_file = esd_dir / f"{job_name}.inp"

    # Determine source geometry (use optimized geometry of initial state)
    xyz_file = f"{initial_state}.xyz"

    # Build input (same as ISC but labeled as IC)
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')
    implicit_solvation = config.get('implicit_solvation_model', 'CPCM')

    # Simple keyword line (no RKS/UKS flag - let ORCA decide based on multiplicity)
    keywords = [
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
        f"{implicit_solvation}({solvent})",
        "ESD(IC)",
    ]

    simple_line = "! " + " ".join(keywords)

    # Blocks
    blocks = []

    # Base
    blocks.append(f'%base "{job_name}"')

    # TDDFT block tailored for IC calculations
    nroots = config.get('ESD_IC_NROOTS', config.get('NROOTS', 10))  # Increased default to 10
    iroot = config.get('IROOT', 1)
    tda_flag = str(config.get('TDA', 'FALSE')).upper()
    nacme_flag = str(config.get('NACME', 'TRUE')).upper()
    etf_flag = str(config.get('ETF', 'TRUE')).upper()
    tddft_block = [
        "%TDDFT",
        f"  TDA      {tda_flag}",
        f"  NROOTS   {nroots}",
        f"  IROOT    {iroot}",
        f"  NACME    {nacme_flag}",
        f"  ETF      {etf_flag}",
        "END",
    ]
    blocks.append("\n".join(tddft_block))

    # ESD block
    temperature = config.get('temperature', 298.15)
    esd_block = [
        "%ESD",
        f'  GSHESSIAN       "{initial_state}.hess"',
        f'  ESHESSIAN       "{final_state}.hess"',
        "  USEJ            TRUE",
        f"  TEMP            {temperature}",
    ]
    esd_block.append("END")
    blocks.append("\n".join(esd_block))

    # PAL and maxcore
    pal = config.get('PAL', 12)
    maxcore = config.get('maxcore', 6000)
    blocks.append(f"%pal nprocs {pal} end")
    blocks.append(f"%maxcore {maxcore}")

    # Geometry - read coordinates (XYZ format with header)
    xyz_path = esd_dir / xyz_file
    try:
        with open(xyz_path, 'r', encoding='utf-8') as f:
            all_lines = f.readlines()
            coord_lines = all_lines[2:]  # Skip atom count and comment line
    except FileNotFoundError:
        logger.error(f"Coordinate file not found: {xyz_path}")
        raise

    # Write input file
    with open(input_file, 'w', encoding='utf-8') as f:
        f.write(simple_line + "\n")
        for block in blocks:
            f.write(block + "\n")
        f.write("\n")
        f.write(f"* xyz {charge} 1\n")
        for line in coord_lines:
            f.write(line)
        f.write("*\n")

    logger.info(f"Created IC input: {input_file}")
    return str(input_file)
