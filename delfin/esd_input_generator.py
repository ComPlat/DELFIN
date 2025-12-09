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
from delfin.common.orca_blocks import resolve_maxiter, collect_output_blocks

logger = get_logger(__name__)

# Conversion factor: Hartree to cm^-1
HARTREE_TO_CM1 = 219474.63


def _format_ms_suffix(trootssl: int) -> str:
    """Format TROOTSSL value as ms suffix (e.g., -1 -> 'msm1', 0 -> 'ms0', 1 -> 'msp1')."""
    if trootssl < 0:
        return f"msm{abs(trootssl)}"
    elif trootssl > 0:
        return f"msp{trootssl}"
    else:
        return "ms0"


def _parse_state_root(label: str) -> tuple[str, int]:
    """Return (state_type, root_index) from labels like 'S1', 'T2'. Defaults to 1 on parse errors."""
    if not label:
        return "", 1
    label = label.strip().upper()
    state_type = label[0]
    try:
        root = int(label[1:]) if len(label) > 1 else 1
    except Exception:
        root = 1
    return state_type, root


def _resolve_tddft_maxiter(config: Dict[str, Any]) -> Optional[int]:
    """Prefer an ESD-specific TDDFT maxiter override, then fall back to global."""
    esd_override = resolve_maxiter(config, key="ESD_TDDFT_maxiter")
    if esd_override is not None:
        return esd_override
    return resolve_maxiter(config, key="TDDFT_maxiter")


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
    """Generate ORCA input for a state, respecting ESD_modus (deltaSCF|TDDFT)."""
    mode = str(config.get("ESD_modus", "TDDFT")).strip().lower()
    # If pipe-separated options (e.g., "TDDFT|deltaSCF"), take first as default
    if "|" in mode:
        mode = mode.split("|")[0].strip()
    if mode == "tddft":
        return _create_state_input_tddft(
            state=state,
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
        )
    return _create_state_input_delta_scf(
        state=state,
        esd_dir=esd_dir,
        charge=charge,
        solvent=solvent,
        metals=metals,
        main_basisset=main_basisset,
        metal_basisset=metal_basisset,
        config=config,
    )


def _create_state_input_delta_scf(
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
    elif state_upper == "S2":
        xyz_file = "S1.xyz"
        moinp_gbw = "S1.gbw"
        use_deltascf = True
    elif state_upper == "T1":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = False
    elif state_upper == "T2":
        xyz_file = "T1.xyz"
        moinp_gbw = "T1.gbw"
        use_deltascf = True
    elif state_upper == "T3":
        xyz_file = "S0.xyz"
        moinp_gbw = "T2.gbw"
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
    geom_token_raw = config.get('geom_opt', 'OPT')
    geom_token = str(geom_token_raw).strip() or "OPT"

    # Check if frequency calculations are enabled for ESD
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')

    # Frequency calculation type from CONTROL (FREQ or numFREQ)
    freq_type = str(config.get('freq_type', 'FREQ')).strip().upper()
    if freq_type not in ('FREQ', 'NUMFREQ'):
        freq_type = 'FREQ'

    # Initial guess from CONTROL (e.g., PModel)
    initial_guess = (str(config.get("initial_guess", "")).split() or [""])[0]

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

    # Only add FREQ/numFREQ if frequency calculations are enabled
    if esd_frequency_enabled:
        keywords.append(freq_type)

    if use_deltascf:
        keywords.append("deltaSCF")

    if moinp_gbw:
        keywords.append("NODIIS")
        keywords.append("MOREAD")

    if initial_guess:
        keywords.append(initial_guess)

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

    # Optional TDDFT iteration limit for follow-up TDDFT checks
    tddft_maxiter = _resolve_tddft_maxiter(config)

    # Optional output blocks (e.g., print_MOs)
    blocks.extend(collect_output_blocks(config, allow=True))

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

    # Geometry - read from start.txt or xyz file
    if xyz_file == "initial.xyz":
        # Prefer optimized initial.xyz; fallback to start.txt
        if Path("initial.xyz").exists():
            xyz_path = Path("initial.xyz")
            skip_lines = 2  # initial.xyz has header
        else:
            xyz_path = Path("start.txt")
            skip_lines = 0  # start.txt has no header
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
            f.write(f'%base "S0_TDDFT"\n')

            # PAL and maxcore
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # TDDFT block for both singlets and triplets
            nroots = config.get('ESD_nroots', 15)
            tda_flag = str(config.get('TDA', 'FALSE')).upper()
            # Use ESD_maxdim if set, otherwise default to nroots/2 (min 5)
            esd_maxdim = config.get('ESD_maxdim', None)
            maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            if tddft_maxiter is not None:
                f.write(f"  maxiter {tddft_maxiter}\n")
            f.write("  triplets true\n")
            f.write("end\n")

            # Geometry reference
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 S0.xyz\n")

            logger.info(f"Added TDDFT check job to S0 input for state identification")

    logger.info(f"Created ESD state input: {input_file}")
    return str(input_file)


def _create_state_input_tddft(
    state: str,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> str:
    """Generate ORCA input for TDDFT-driven ESD mode."""
    state_upper = state.upper()
    input_file = esd_dir / f"{state_upper}.inp"

    functional = config.get("functional", "PBE0")
    disp_corr = config.get("disp_corr", "D4")
    ri_jkx = config.get("ri_jkx", "RIJCOSX")
    aux_jk = config.get("aux_jk", "def2/J")
    implicit_solvation = config.get("implicit_solvation_model", "CPCM")
    geom_token_raw = config.get("geom_opt", "OPT")
    geom_token = str(geom_token_raw).strip() or "OPT"
    pal = config.get("PAL", 12)
    maxcore = config.get("maxcore", 6000)
    nroots = config.get("ESD_nroots", 15)
    tda_flag = str(config.get("ESD_TDA", config.get("TDA", "FALSE"))).upper()
    # Use ESD_maxdim if set, otherwise default to nroots/2 (min 5)
    esd_maxdim = config.get("ESD_maxdim", None)
    maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
    tddft_maxiter = _resolve_tddft_maxiter(config)
    followiroot = str(config.get("ESD_followiroot", "true")).lower() in ("true", "yes", "1", "on")
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')
    output_blocks = collect_output_blocks(config, allow=True)

    def _join_keywords(parts: List[str]) -> str:
        """Join keyword fragments while skipping empty entries."""
        return " ".join(str(p) for p in parts if str(p).strip())

    def _build_keywords(scf_type: str, with_freq: bool = True) -> List[str]:
        """Build keyword list with optional frequency calculation."""
        kw = [
            functional,
            scf_type,
            main_basisset,
            disp_corr,
            ri_jkx,
            aux_jk,
            f"{implicit_solvation}({solvent})",
            geom_token,
        ]
        if with_freq and esd_frequency_enabled:
            kw.append("numFREQ")
        return kw

    # Coordinate source
    if state_upper == "S0":
        if Path("initial.xyz").exists():
            xyz_path = Path("initial.xyz")
            skip_lines = 2
        else:
            xyz_path = Path("start.txt")
            skip_lines = 0
    elif state_upper == "S1":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "S2":
        xyz_path = esd_dir / "S1.xyz"
        skip_lines = 2
    elif state_upper == "S3":
        xyz_path = esd_dir / "S2.xyz"
        skip_lines = 2
    elif state_upper == "S4":
        xyz_path = esd_dir / "S3.xyz"
        skip_lines = 2
    elif state_upper == "S5":
        xyz_path = esd_dir / "S4.xyz"
        skip_lines = 2
    elif state_upper == "S6":
        xyz_path = esd_dir / "S5.xyz"
        skip_lines = 2
    elif state_upper == "T1":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T2":
        xyz_path = esd_dir / "T1.xyz"
        skip_lines = 2
    elif state_upper == "T3":
        xyz_path = esd_dir / "T2.xyz"
        skip_lines = 2
    elif state_upper == "T4":
        xyz_path = esd_dir / "T3.xyz"
        skip_lines = 2
    elif state_upper == "T5":
        xyz_path = esd_dir / "T4.xyz"
        skip_lines = 2
    elif state_upper == "T6":
        xyz_path = esd_dir / "T5.xyz"
        skip_lines = 2
    else:
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2

    try:
        with open(xyz_path, "r", encoding="utf-8") as f:
            all_lines = f.readlines()
            coord_lines = all_lines[skip_lines:]
    except FileNotFoundError:
        logger.error(f"Coordinate file not found: {xyz_path}")
        raise

    def _write_tddft_block(
        fh,
        iroot: Optional[int] = None,
        irootmult: Optional[str] = None,
        *,
        triplets: bool = False,
    ) -> None:
        fh.write("%tddft\n")
        fh.write(f"  nroots {nroots}\n")
        fh.write(f"  maxdim {maxdim}\n")
        fh.write(f"  tda {tda_flag}\n")
        if tddft_maxiter is not None:
            fh.write(f"  maxiter {tddft_maxiter}\n")

        if triplets:
            fh.write("  triplets true\n")
        if iroot is not None:
            fh.write(f"  iroot {iroot}\n")
        if irootmult:
            fh.write(f"  irootmult {irootmult}\n")
        if iroot is not None and followiroot:
            fh.write("  followiroot true\n")
        fh.write("end\n")

    def _write_output_blocks(fh) -> None:
        for block in output_blocks:
            fh.write(block if block.endswith("\n") else block + "\n")

    with open(input_file, "w", encoding="utf-8") as f:
        if state_upper == "S0":
            keywords = [
                functional,
                "RKS",
                main_basisset,
                disp_corr,
                ri_jkx,
                aux_jk,
                f"{implicit_solvation}({solvent})",
            ]
            if geom_token:
                keywords.append(geom_token)
            if esd_frequency_enabled:
                keywords.append("numFREQ")
            f.write("! " + " ".join(keywords) + "\n")
            f.write('%base "S0"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n\n")

            f.write("$new_job\n")
            f.write("! " + " ".join(keywords[:-2]) + "\n")
            f.write('%base "S0_TDDFT"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, triplets=True)
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 S0.xyz\n")
        elif state_upper == "S1":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S1"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            f.write("\n")
            _write_tddft_block(f, iroot=1, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n\n")
        elif state_upper == "T1":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T1"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=1, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "T2":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T2"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=2, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "S2":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S2"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=2, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "T3":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T3"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=3, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "S3":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S3"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=3, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "S4":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S4"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=4, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "S5":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S5"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=5, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "S6":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "S6"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=6, irootmult="singlet")
            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "T4":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T4"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=4, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "T5":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T5"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=5, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        elif state_upper == "T6":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T6"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=6, irootmult="triplet")
            f.write(f"\n* xyz {charge} 1\n")  # Multiplicity 1 with irootmult=triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")
        else:
            raise ValueError(f"Unknown state: {state}")

    logger.info(f"Created ESD TDDFT state input: {input_file}")
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
    trootssl: Optional[int] = None,
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
        trootssl: Triplet sublevel (-1, 0, or 1). If None, uses config['TROOTSSL']

    Returns:
        Path to generated input file
    """
    initial_state, final_state = isc_pair.split(">")
    initial_state = initial_state.strip().upper()
    final_state = final_state.strip().upper()
    init_type, init_root = _parse_state_root(initial_state)
    final_type, final_root = _parse_state_root(final_state)

    # Determine TROOTSSL value for this calculation
    if trootssl is None:
        trootssl_str = str(config.get('TROOTSSL', '0')).strip()
    else:
        trootssl_str = str(trootssl)

    # Generate job name with TROOTSSL suffix
    trootssl_int = int(trootssl_str)
    ms_suffix = _format_ms_suffix(trootssl_int)
    job_name = f"{initial_state}_{final_state}_ISC_{ms_suffix}"
    input_file = esd_dir / f"{job_name}.inp"

    # Determine source geometry (use optimized geometry of FINAL state, per ORCA manual)
    xyz_file = f"{final_state}.xyz"
    # Use restricted (closed-shell) reference for SOC; keep multiplicity 1 to avoid UKS
    final_multiplicity = 1

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

    # Simple keyword line (restricted reference)
    keywords = [
        "RKS",
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

    # Map roots to correct spin manifolds based on states, per ORCA ESD docs
    s_root = init_root if init_type == "S" else (final_root if final_type == "S" else 1)
    t_root = init_root if init_type == "T" else (final_root if final_type == "T" else 1)

    dosoc_flag = "TRUE"
    tddft_maxiter = _resolve_tddft_maxiter(config)
    tddft_block = [
        f"%TDDFT  NROOTS  {int(nroots):>2}",
        f"        SROOT   {int(s_root)}",
        f"        TROOT   {int(t_root)}",
        f"        TROOTSSL {trootssl_str}",
        f"        DOSOC   {dosoc_flag}",
    ]
    if tddft_maxiter is not None:
        tddft_block.append(f"        maxiter {tddft_maxiter}")
    tddft_block.append(
        "END",
    )
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
        f.write(f"* xyz {charge} {final_multiplicity}\n")
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
    init_type, init_root = _parse_state_root(initial_state)
    final_type, _ = _parse_state_root(final_state)

    job_name = f"{initial_state}_{final_state}_IC"
    input_file = esd_dir / f"{job_name}.inp"

    # Determine source geometry (use lower-state geometry for IC, per ORCA manual)
    # For S1>S0: use S0.xyz; for Tn>T1: use T1.xyz
    xyz_file = f"{final_state}.xyz"
    # Multiplicity follows final state (triplet -> 3, singlet -> 1)
    final_type, _ = _parse_state_root(final_state)
    final_multiplicity = 3 if final_type == "T" else 1

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

    # Calculate IROOT: For Tn->T1 IC, T1 is the SCF ground state (multiplicity 3)
    # and Tn is the (n-1)-th excited state above T1
    # For Sn->S0 IC, S0 is the SCF ground state and Sn is the n-th excited state
    if final_type == "T" and final_state == "T1":
        # Triplet IC: T2->T1 uses IROOT=1, T3->T1 uses IROOT=2, etc.
        iroot = config.get('IROOT', init_root - 1)
    else:
        # Singlet IC: S1->S0 uses IROOT=1, S2->S0 uses IROOT=2, etc.
        iroot = config.get('IROOT', init_root)

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
    tddft_maxiter = _resolve_tddft_maxiter(config)
    if tddft_maxiter is not None:
        tddft_block.insert(-1, f"  maxiter  {tddft_maxiter}")
    blocks.append("\n".join(tddft_block))

    # ESD block
    # For IC: GSHESSIAN = ground state (final), ESHESSIAN = excited state (initial)
    # Example: S1>S0 IC → GSHESSIAN=S0.hess, ESHESSIAN=S1.hess
    temperature = config.get('temperature', 298.15)
    esd_block = [
        "%ESD",
        f'  GSHESSIAN       "{final_state}.hess"',
        f'  ESHESSIAN       "{initial_state}.hess"',
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
        f.write(f"* xyz {charge} {final_multiplicity}\n")
        for line in coord_lines:
            f.write(line)
        f.write("*\n")

    logger.info(f"Created IC input: {input_file}")
    return str(input_file)
