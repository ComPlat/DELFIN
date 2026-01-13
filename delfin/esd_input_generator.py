"""ORCA input file generator for ESD module (excited state dynamics).

This module generates ORCA input files for:
- Electronic states (S0, S1, T1, T2)
- Intersystem crossings (ISCs)
- Internal conversions (ICs)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from delfin.common.logging import get_logger
from delfin.common.orca_blocks import resolve_maxiter, collect_output_blocks

logger = get_logger(__name__)

# Conversion factor: Hartree to cm^-1
HARTREE_TO_CM1 = 219474.63


def _resolve_state_filename(state: str, extension: str, esd_mode: str) -> str:
    """Resolve the correct filename for a state based on ESD mode.

    In hybrid1 mode:
    - S0: S0.{extension}
    - T1: T1.{extension} (special case - single step)
    - S1, S2, T2, etc:
      - .out files: {state}_second.out
      - .xyz, .hess, .gbw files: {state}_second_deltaSCF.{extension}

    In tddft mode:
    - All states: {state}.{extension}

    Args:
        state: State label (e.g., 'S1', 'T2')
        extension: File extension without dot (e.g., 'out', 'xyz', 'gbw', 'hess')
        esd_mode: ESD mode ('hybrid1' or 'tddft')

    Returns:
        Filename like 'S1_second.out' or 'S1_second_deltaSCF.xyz'
    """
    state_upper = state.upper()

    # In tddft mode, always use simple naming
    if esd_mode != 'hybrid1':
        return f"{state_upper}.{extension}"

    # In hybrid1 mode, S0 and T1 use simple naming
    if state_upper in ('S0', 'T1'):
        return f"{state_upper}.{extension}"

    # All other states in hybrid1 use _second suffix
    # .out files: S1_second.out
    # .xyz, .hess, .gbw: S1_second_deltaSCF.xyz
    if extension == 'out':
        return f"{state_upper}_second.{extension}"
    else:
        return f"{state_upper}_second_deltaSCF.{extension}"


def _parse_tddft_and_derive_deltascf(s0_out_path: Path, state: str) -> Optional[tuple[str, str]]:
    """Parse TD-DFT results from S0.out and derive ALPHACONF/BETACONF for deltaSCF.

    Args:
        s0_out_path: Path to S0.out file containing TD-DFT results
        state: Target state (e.g., 'S1', 'S2', 'T1', 'T2')

    Returns:
        Tuple of (alphaconf_string, betaconf_string) or None if parsing fails

    Example:
        For S2 with dominant excitation 117a -> 119a (HOMO-1 -> LUMO):
        Returns ("0,1,1", "0")
    """
    if not s0_out_path.exists():
        logger.warning(f"S0.out not found at {s0_out_path}, cannot derive deltaSCF config")
        return None

    import re

    # Parse state type and number (e.g., 'S2' -> 'S', 2)
    match = re.match(r'([ST])(\d+)', state.upper())
    if not match:
        return None
    state_type, state_num = match.group(1), int(match.group(2))

    # Read S0.out and find HOMO orbital number
    try:
        content = s0_out_path.read_text(encoding='utf-8', errors='ignore')
    except Exception as exc:
        logger.warning(f"Failed to read {s0_out_path}: {exc}")
        return None

    # Find HOMO orbital number (last occupied orbital before JOB 2)
    # Look for orbital energies before "JOB NUMBER  2"
    job2_match = re.search(r'\$+\s*JOB NUMBER\s+2\s+\$+', content)
    if job2_match:
        content_before_job2 = content[:job2_match.start()]
    else:
        content_before_job2 = content

    # Find last ORBITAL ENERGIES section before JOB 2
    # Support both RHF/RKS (no spin labels) and UHF/UKS (SPIN UP/DOWN labels)
    orbital_sections = list(re.finditer(r'ORBITAL ENERGIES\s*\n-+', content_before_job2))
    if not orbital_sections:
        logger.warning("Could not find ORBITAL ENERGIES section in S0.out")
        return None

    last_orbital_section = orbital_sections[-1]
    orbital_text = content_before_job2[last_orbital_section.end():][:500000]  # Take next 500000 chars (supports systems with 1000+ orbitals)

    # Find HOMO (last orbital with OCC = 2.0000 or 1.0000)
    homo_num = None
    for line in orbital_text.split('\n'):
        orbital_match = re.match(r'\s*(\d+)\s+(2\.0000|1\.0000)\s+', line)
        if orbital_match:
            homo_num = int(orbital_match.group(1))

    if homo_num is None:
        logger.warning("Could not determine HOMO orbital number from S0.out")
        return None

    lumo_num = homo_num + 1
    logger.info(f"Detected HOMO = orbital {homo_num}, LUMO = orbital {lumo_num}")

    # Find TD-DFT excitations for the requested state
    # Note: TD-DFT results may be in a different JOB than ORBITAL ENERGIES
    if state_type == 'S':
        section_header = r'TD-DFT EXCITED STATES \(SINGLETS\)'
    else:  # T
        section_header = r'TD-DFT EXCITED STATES \(TRIPLETS\)'

    # Search in full content (not just before JOB 2)
    tddft_match = re.search(section_header, content)
    if not tddft_match:
        logger.warning(f"Could not find {section_header} section in S0.out")
        return None

    tddft_text = content[tddft_match.end():][:50000]  # Take next 50000 chars for TD-DFT data

    # Find STATE N where N = state_num
    # Match STATE line, then capture all excitations until next STATE or empty line
    state_pattern = rf'STATE\s+{state_num}:.*?\n(.*?)(?:\n\s*\n|STATE|$)'
    state_match = re.search(state_pattern, tddft_text, re.DOTALL)

    if not state_match:
        logger.warning(f"Could not find STATE {state_num} in TD-DFT results")
        return None

    excitations_text = state_match.group(1)

    # Parse excitations and find dominant one (highest weight)
    excitations = []
    for line in excitations_text.strip().split('\n'):
        exc_match = re.match(r'\s*(\d+)([ab])\s+->\s+(\d+)([ab])\s+:\s+([\d.]+)', line.strip())
        if exc_match:
            from_orb = int(exc_match.group(1))
            from_spin = exc_match.group(2)
            to_orb = int(exc_match.group(3))
            to_spin = exc_match.group(4)
            weight = float(exc_match.group(5))
            excitations.append((from_orb, from_spin, to_orb, to_spin, weight))

    if not excitations:
        logger.warning(f"No excitations found for STATE {state_num}")
        return None

    # Get dominant excitation (highest weight)
    dominant = max(excitations, key=lambda x: x[4])
    from_orb, from_spin, to_orb, to_spin, weight = dominant

    logger.info(f"{state}: Dominant excitation {from_orb}{from_spin} -> {to_orb}{to_spin} (weight={weight:.4f})")

    # Calculate ALPHACONF/BETACONF based on excitation
    # n = homo_num - from_orb  # How many orbitals below HOMO
    # m = to_orb - lumo_num    # How many orbitals above LUMO
    #
    # For Singlets (S1, S2, S3, ...): Standard single excitation
    #   HOMO-n → LUMO+m: alphaconf [0]*(m+1) + [1]*(n+1), betaconf "0"
    #   Examples:
    #     HOMO→LUMO (n=0, m=0): alphaconf "0,1", betaconf "0"
    #     HOMO-1→LUMO (n=1, m=0): alphaconf "0,1,1", betaconf "0"
    #
    # For Triplets (T1, T2, T3, ...): Two electrons with parallel spin
    #   HOMO-n → LUMO+m:
    #     ALPHACONF: "1,1" if m=0, else [0]*m + [1]
    #     BETACONF: [0] + [1]*n
    #   Examples:
    #     T1 HOMO→LUMO (n=0, m=0): alphaconf "1,1", betaconf "0"
    #     T2 HOMO→LUMO+1 (n=0, m=1): alphaconf "0,1", betaconf "0"
    #     HOMO-3→LUMO+2 (n=3, m=2): alphaconf "0,0,1", betaconf "0,1,1,1"

    n = homo_num - from_orb  # How many orbitals below HOMO
    m = to_orb - lumo_num    # How many orbitals above LUMO

    # Special handling for ALL Triplet states (T1, T2, T3, ...)
    if state_type == 'T':
        # Triplet configuration for HOMO-n → LUMO+m excitation:
        # ALPHACONF: "1,1" if m=0, else [0]*m + [1]
        # BETACONF: [0] + [1]*n
        # Examples:
        #   T1 (HOMO→LUMO, n=0, m=0): alphaconf "1,1", betaconf "0"
        #   T2 (HOMO→LUMO+1, n=0, m=1): alphaconf "0,1", betaconf "0"
        #   HOMO-3→LUMO+2 (n=3, m=2): alphaconf "0,0,1", betaconf "0,1,1,1"
        if m == 0:
            alphaconf = "1,1"
        else:
            alphaconf_list = [0] * m + [1]
            alphaconf = ','.join(map(str, alphaconf_list))

        betaconf_list = [0] + [1] * n
        betaconf = ','.join(map(str, betaconf_list))
        logger.info(f"{state}: Triplet config (HOMO-{n}→LUMO+{m}) - ALPHACONF {alphaconf}, BETACONF {betaconf}")
    elif from_spin == 'a':
        # Alpha excitation (standard for Singlets)
        alphaconf_list = [0] * (m + 1) + [1] * (n + 1)
        alphaconf = ','.join(map(str, alphaconf_list))
        betaconf = "0"
        logger.info(f"{state}: Singlet config - ALPHACONF {alphaconf}, BETACONF {betaconf}")
    else:
        # Beta excitation (fallback for other cases)
        alphaconf = "0,1"  # Still promote alpha
        betaconf_list = [0] * (m + 1) + [1] * (n + 1)
        betaconf = ','.join(map(str, betaconf_list))
        logger.info(f"{state}: Beta excitation - ALPHACONF {alphaconf}, BETACONF {betaconf}")

    return (alphaconf, betaconf)


def _build_solvation_keyword(implicit_solvation_model: str, solvent: str) -> str:
    """Build solvation keyword string, returning empty string if model is not set.

    Args:
        implicit_solvation_model: Solvation model (e.g., 'CPCM', 'SMD')
        solvent: Solvent name (e.g., 'water', 'DMF')

    Returns:
        Formatted solvation keyword (e.g., 'CPCM(water)') or empty string
    """
    model = str(implicit_solvation_model).strip()
    if not model:
        return ""
    if solvent and str(solvent).strip():
        return f"{model}({solvent})"
    return model


def _format_ms_suffix(trootssl: int) -> str:
    """Format TROOTSSL value as ms suffix (e.g., -1 -> 'msm1', 0 -> 'ms0', 1 -> 'msp1')."""
    if trootssl < 0:
        return f"msm{abs(trootssl)}"
    elif trootssl > 0:
        return f"msp{trootssl}"
    else:
        return "ms0"


def generate_elprop_block(config: Dict[str, Any]) -> str:
    """Generate %elprop block based on individual property settings.

    Args:
        config: Configuration dictionary containing elprop_* parameters

    Returns:
        Formatted %elprop block string, or empty string if no properties are specified

    Example config values:
        elprop_Dipole = true
        elprop_Polar = true
        elprop_PolarVelocity = true
        elprop_PolarDipQuad = true
        elprop_Hyperpol = true
    """
    def _parse_bool(value: Any, default: bool = False) -> bool:
        """Parse boolean value from config."""
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        s = str(value).strip().lower()
        if s in ('true', 'yes', '1', 'on'):
            return True
        if s in ('false', 'no', '0', 'off'):
            return False
        return default

    # Main properties in order
    main_properties = [
        ('Dipole', 'elprop_Dipole'),
        ('Quadrupole', 'elprop_Quadrupole'),
        ('Polar', 'elprop_Polar'),
        ('Hyperpol', 'elprop_Hyperpol'),
    ]

    # Sub-properties for Polar
    polar_sub_properties = [
        ('PolarVelocity', 'elprop_PolarVelocity'),
        ('PolarDipQuad', 'elprop_PolarDipQuad'),
        ('PolarQuadQuad', 'elprop_PolarQuadQuad'),
    ]

    # Collect enabled properties
    enabled_properties = []
    polar_enabled = False

    for prop_name, config_key in main_properties:
        value = config.get(config_key)
        if _parse_bool(value, default=False):
            enabled_properties.append(prop_name)
            if prop_name == 'Polar':
                polar_enabled = True

    # If no properties enabled, return empty string
    if not enabled_properties:
        return ""

    # Build the %elprop block
    lines = ["%elprop"]

    for prop_name in enabled_properties:
        lines.append(f"  {prop_name} true")

        # Add Polar sub-properties if Polar is enabled
        if prop_name == 'Polar' and polar_enabled:
            for sub_name, sub_key in polar_sub_properties:
                sub_value = config.get(sub_key)
                if _parse_bool(sub_value, default=False):
                    lines.append(f"  {sub_name} true")

    lines.append("end")

    return "\n".join(lines)


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
    """Prefer TDDFT_maxiter, fall back to ESD_TDDFT_maxiter (legacy), then global TDDFT_maxiter."""
    # Try new naming first
    tddft_maxiter = resolve_maxiter(config, key="TDDFT_maxiter")
    if tddft_maxiter is not None:
        return tddft_maxiter
    # Fall back to legacy ESD_TDDFT_maxiter for backwards compatibility
    esd_override = resolve_maxiter(config, key="ESD_TDDFT_maxiter")
    if esd_override is not None:
        return esd_override
    return None


def _get_tddft_param(config: Dict[str, Any], param_name: str, default: Any) -> Any:
    """Get TDDFT parameter with fallback to legacy ESD_ naming.

    Args:
        config: Configuration dictionary
        param_name: Parameter name without prefix (e.g., 'nroots', 'maxdim', 'TDA')
        default: Default value if not found

    Returns:
        Parameter value from TDDFT_{param_name} or ESD_{param_name} (legacy)
    """
    # Try new TDDFT_ prefix first
    new_key = f"TDDFT_{param_name}"
    if new_key in config:
        return config[new_key]

    # Fall back to legacy ESD_ prefix for backwards compatibility
    legacy_key = f"ESD_{param_name}"
    return config.get(legacy_key, default)


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


def _parse_step_set(raw: Any) -> Set[int]:
    """Parse oxidation/reduction step spec into a set of ints."""
    if raw is None:
        return set()
    if isinstance(raw, (int, float)):
        return {int(raw)}
    tokens = str(raw).replace(';', ',').replace('[', '').replace(']', '').split(',')
    steps: Set[int] = set()
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        try:
            steps.add(int(tok))
        except Exception:  # noqa: BLE001
            continue
    return steps


def _parse_iroot_spec(raw: Any, *, default: List[int]) -> List[int]:
    """Parse an IROOT specification from CONTROL.

    Semantics:
    - A single value (e.g. "3" or 3) means *exactly* IROOT=3
    - A comma/whitespace separated list (e.g. "1,2,3") means those IROOTs
    - A simple range "1-3" expands to 1,2,3
    """
    if raw is None:
        return list(default)
    if isinstance(raw, (int, float)):
        val = int(raw)
        return [val] if val > 0 else list(default)
    text = str(raw).strip()
    if not text:
        return list(default)

    # Range "a-b"
    import re

    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", text)
    if m:
        start = int(m.group(1))
        end = int(m.group(2))
        if start <= 0 or end <= 0:
            return list(default)
        if start <= end:
            return list(range(start, end + 1))
        return list(range(start, end - 1, -1))

    # Single integer token -> single IROOT
    if text.isdigit() or (text.startswith("-") and text[1:].isdigit()):
        val = int(text)
        return [val] if val > 0 else list(default)

    tokens = text.replace(";", ",").replace(" ", ",").split(",")
    values: List[int] = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        try:
            val = int(tok)
        except Exception:  # noqa: BLE001
            continue
        if val > 0 and val not in values:
            values.append(val)
    return values or list(default)


def _resolve_temperature_K(config: Dict[str, Any], default: float = 298.15) -> float:
    """Resolve temperature in K from CONTROL config."""
    for key in ("temperature", "temperature_K", "TEMP"):
        val = config.get(key)
        if val is None:
            continue
        try:
            return float(val)
        except Exception:  # noqa: BLE001
            continue
    return float(default)


def _apply_esd_newgto(
    coord_lines: List[str],
    *,
    found_metals: List[str],
    metal_basisset: Optional[str],
    config: Dict[str, Any],
) -> List[str]:
    """Apply per-atom NewGTO tagging to an XYZ coordinate block for ESD inputs.

    Mirrors the behavior used for classic initial inputs:
    - Always tags detected metals with NewGTO "metal_basisset" end
    - Optionally tags the first coordination sphere when enabled in CONTROL
    """
    if not metal_basisset:
        return coord_lines

    enable_first = str(config.get("first_coordination_sphere_metal_basisset", "no")).lower() in (
        "yes",
        "true",
        "1",
        "on",
    )
    sphere_scale_raw = str(config.get("first_coordination_sphere_scale", "")).strip()
    load_radii = enable_first and not sphere_scale_raw

    radii_all = None
    if load_radii:
        try:
            from delfin.xyz_io import _load_covalent_radii
            radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009"))
        except Exception:
            radii_all = None

    try:
        from delfin.xyz_io import _apply_per_atom_newgto
        return _apply_per_atom_newgto(coord_lines, found_metals, metal_basisset, config, radii_all)
    except Exception:
        # Fail safe: never break input generation because of NewGTO tagging.
        return coord_lines


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
    """Generate ORCA input for a state, respecting ESD_modus (deltaSCF|TDDFT|hybrid1)."""
    mode = str(config.get("ESD_modus", "TDDFT")).strip().lower()
    # If pipe-separated options (e.g., "TDDFT|deltaSCF|hybrid1"), take first as default
    if "|" in mode:
        mode = mode.split("|")[0].strip()
    if mode == "tddft":
        input_file = _create_state_input_tddft(
            state=state,
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
        )
    elif mode == "hybrid1":
        input_file = _create_state_input_hybrid1(
            state=state,
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
        )
    else:
        input_file = _create_state_input_delta_scf(
            state=state,
            esd_dir=esd_dir,
            charge=charge,
            solvent=solvent,
            metals=metals,
            main_basisset=main_basisset,
            metal_basisset=metal_basisset,
            config=config,
        )

    # Add properties_of_interest jobs (IP, EA) only for S0 when method=classic and calc_prop_of_interest=yes
    if state.strip().upper() == "S0":
        method = str(config.get('method', '')).strip().lower()
        calc_prop = str(config.get('calc_prop_of_interest', 'no')).strip().lower()
        properties = config.get('properties_of_interest', '')
        if calc_prop in ('yes', 'true', '1', 'on') and properties and method == 'classic':
            # Determine which properties to keep based on requested red/ox steps
            ox_steps = _parse_step_set(config.get('oxidation_steps'))
            red_steps = _parse_step_set(config.get('reduction_steps'))
            prop_tokens = str(properties).strip()
            prop_tokens = prop_tokens.strip('[]').replace("'", "").replace('"', '')
            prop_set = {p.strip().upper() for p in prop_tokens.split(',') if p.strip()}

            filtered_props = []
            if 'IP' in prop_set and (not ox_steps or 1 in ox_steps):
                filtered_props.append('IP')
            if 'EA' in prop_set and (not red_steps or 1 in red_steps):
                filtered_props.append('EA')

            if filtered_props:
                # Get multiplicity for S0
                multiplicity = config.get('multiplicity_0') or config.get('multiplicity')
                if multiplicity is None:
                    # Calculate from charge and total electrons (simplified)
                    multiplicity = 1  # Default singlet for S0

                xyz_file = str(esd_dir / "S0.xyz")
                append_properties_of_interest_jobs(
                    inp_file=input_file,
                    xyz_file=xyz_file,
                    base_charge=charge,
                    base_multiplicity=multiplicity,
                    properties=",".join(filtered_props),
                    config=config,
                    solvent=solvent,
                    metals=metals,
                    main_basisset=main_basisset,
                    metal_basisset=metal_basisset,
                )

    return input_file


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
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "S3":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "S4":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "S5":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "S6":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "T1":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = False
    elif state_upper == "T2":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    elif state_upper == "T3":
        xyz_file = "S0.xyz"
        moinp_gbw = "S0.gbw"
        use_deltascf = True
    else:
        raise ValueError(f"Unknown state: {state}")

    # Functional and basis set
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')

    # Solvation
    implicit_solvation = config.get('implicit_solvation_model', '')

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
    # For deltaSCF calculations: use RKS for S0 (closed-shell), UKS for excited states
    # This applies to deltaSCF mode and hybrid1 mode (T1, second steps of S1/T2/etc.)
    scf_type = "RKS" if state_upper == "S0" else "UKS"

    keywords = [
        functional,
        scf_type,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
    ]

    # Add solvation keyword only if model is set
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)
    if solvation_kw:
        keywords.append(solvation_kw)

    if geom_token:
        keywords.append(geom_token)

    # Only add FREQ/numFREQ if frequency calculations are enabled
    # S0 always uses FREQ, all other states use freq_type (typically numFREQ)
    if esd_frequency_enabled:
        if state_upper == "S0":
            keywords.append("FREQ")
        else:
            keywords.append(freq_type)

    if use_deltascf:
        keywords.append("deltaSCF")

    if moinp_gbw:
        keywords.append("MOREAD")

    # Only add initial_guess if not PModel (MOREAD replaces PModel when reading orbitals)
    if initial_guess and not (moinp_gbw and initial_guess == "PModel"):
        keywords.append(initial_guess)

    # Add deltaSCF-specific keywords
    if use_deltascf:
        deltascf_keywords = config.get('deltaSCF_keywords', '')
        if deltascf_keywords:
            keywords.extend(deltascf_keywords.split())

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
        maxiter = config.get('deltaSCF_maxiter', 1000)
        soscf_convfactor = config.get('deltaSCF_SOSCFConvFactor', 500)
        soscf_maxstep = config.get('deltaSCF_SOSCFMaxStep', 0.1)

        scf_block = [
            "%scf",
            f"  DOMOM {domom}",
            f"  pmom {pmom}",
            f"  keepinitialref {keepinitialref}",
        ]

        # State-specific orbital configurations - derived from TD-DFT if available
        alphaconf = None
        betaconf = None

        # Try to derive ALPHACONF/BETACONF from S0.out TD-DFT results
        s0_out_path = esd_dir / "S0.out"
        if s0_out_path.exists():
            result = _parse_tddft_and_derive_deltascf(s0_out_path, state_upper)
            if result:
                alphaconf, betaconf = result

        # Fallback to hardcoded values if parsing failed
        if alphaconf is None:
            if state_upper == "S1":
                alphaconf = "0,1"
                betaconf = "0"
            elif state_upper == "T2":
                alphaconf = "0,1"
                betaconf = "0,1"  # Triplet spin-flip
            else:
                # Default for other states
                alphaconf = "0,1"
                betaconf = "0"

        scf_block.extend([
            f"  alphaconf {alphaconf}",
            f"  betaconf {betaconf}",
        ])

        scf_block.append(f"  SOSCFHESSUP {soscfhessup}")
        scf_block.append(f"  maxiter {maxiter}")
        scf_block.append(f"  SOSCFConvFactor {soscf_convfactor}")
        scf_block.append(f"  SOSCFMaxStep {soscf_maxstep}")
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
    coord_lines = _apply_esd_newgto(
        coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

    # Write input file
    with open(input_file, 'w', encoding='utf-8') as f:
        f.write(simple_line + "\n")
        for block in blocks:
            f.write(block + "\n")

        # Add custom additions for S0 state if specified in CONTROL
        if state_upper == "S0":
            # Generate %elprop block from elprop_properties config
            elprop_block = generate_elprop_block(config)
            if elprop_block:
                f.write(f"{elprop_block}\n")

            # Legacy addition_S0 support (can be used together with elprop_properties)
            addition_s0 = config.get('addition_S0', '').strip()
            if addition_s0:
                f.write(f"{addition_s0}\n")

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
            ]
            # Add solvation keyword only if model is set
            if solvation_kw:
                tddft_keywords.append(solvation_kw)
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
            # Read ESD_SOC setting
            dosoc_flag = str(config.get('ESD_SOC', 'false')).strip().lower()
            dosoc_value = "true" if dosoc_flag in ('yes', 'true', '1', 'on') else "false"
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            if tddft_maxiter is not None:
                f.write(f"  maxiter {tddft_maxiter}\n")
            f.write("  triplets true\n")
            f.write(f"  dosoc {dosoc_value}\n")
            f.write("end\n")

            # Geometry reference - use optimized S0.xyz
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 S0.xyz\n")

            logger.info(f"Added TDDFT check job to S0 input for state identification")

        # Add TDDFT check job for T1 only
        elif state_upper == "T1":
            # Determine XYZ file for this state (S1.xyz, T1.xyz, etc.)
            state_xyz_file = f"{state_upper}.xyz"

            f.write("\n")
            f.write("#==========================================\n")
            f.write(f"# TDDFT Check: Transitions from {state_upper}\n")
            f.write("#==========================================\n")
            f.write("\n")
            f.write("$new_job\n")

            # TDDFT keyword line - always RKS (closed shell) for TDDFT check jobs
            tddft_keywords = [
                functional,
                "RKS",
                main_basisset,
                disp_corr,
                ri_jkx,
                aux_jk,
            ]
            # Add solvation keyword only if model is set
            if solvation_kw:
                tddft_keywords.append(solvation_kw)
            f.write("! " + " ".join(tddft_keywords) + "\n")

            # Base block for TDDFT check
            f.write(f'%base "{state_upper}_TDDFT"\n')

            # PAL and maxcore
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # TDDFT block - NO followiroot for excited state checks
            nroots = config.get('ESD_nroots', 15)
            tda_flag = str(config.get('TDA', 'FALSE')).upper()
            esd_maxdim = config.get('ESD_maxdim', None)
            maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
            dosoc_flag = str(config.get('ESD_SOC', 'false')).strip().lower()
            dosoc_value = "true" if dosoc_flag in ('yes', 'true', '1', 'on') else "false"
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            if tddft_maxiter is not None:
                f.write(f"  maxiter {tddft_maxiter}\n")
            f.write("  triplets true\n")
            f.write(f"  dosoc {dosoc_value}\n")
            f.write("end\n")

            # Geometry reference - always multiplicity 1 (closed shell) for TDDFT
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 {state_xyz_file}\n")

            logger.info(f"Added TDDFT check job to {state_upper} input for transition analysis")

    def _write_tddft_check_input(state_label: str) -> None:
        if not (state_label.startswith("S") and state_label != "S0"):
            return

        tddft_input = esd_dir / f"{state_label}_TDDFT.inp"
        state_xyz_file = f"{state_label}.xyz"

        with open(tddft_input, "w", encoding="utf-8") as f:
            f.write("#==========================================\n")
            f.write(f"# TDDFT Check: Transitions from {state_label}\n")
            f.write("#==========================================\n")
            f.write("\n")

            # TDDFT keyword line - always RKS (closed shell) for TDDFT check jobs
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
            f.write("! " + " ".join(tddft_keywords) + "\n")

            # Base block for TDDFT check
            f.write(f'%base "{state_label}_TDDFT"\n')

            # PAL and maxcore
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # TDDFT block - NO followiroot for excited state checks
            nroots = config.get('ESD_nroots', 15)
            tda_flag = str(config.get('TDA', 'FALSE')).upper()
            esd_maxdim = config.get('ESD_maxdim', None)
            maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
            dosoc_flag = str(config.get('ESD_SOC', 'false')).strip().lower()
            dosoc_value = "true" if dosoc_flag in ('yes', 'true', '1', 'on') else "false"
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            if tddft_maxiter is not None:
                f.write(f"  maxiter {tddft_maxiter}\n")
            f.write("  triplets true\n")
            f.write(f"  dosoc {dosoc_value}\n")
            f.write("end\n")

            # Geometry reference - always multiplicity 1 (closed shell) for TDDFT
            f.write("\n")
            f.write(f"* xyzfile {charge} 1 {state_xyz_file}\n")

        logger.info(f"Created TDDFT check input: {tddft_input}")

    _write_tddft_check_input(state_upper)

    logger.info(f"Created ESD state input: {input_file}")
    return str(input_file)


def _create_state_input_hybrid1(
    state: str,
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
) -> str:
    """Generate hybrid TDDFT→deltaSCF input files for excited states.

    For S0: Creates standard TDDFT input
    For T1: Creates simple UKS mult=3 optimization (lowest triplet, stable without deltaSCF)
    For excited states (S1, S2, T2, T3, etc.):
        1. {state}_first_TDDFT.inp: TDDFT OPT without FREQ
        2. {state}_second.inp: deltaSCF OPT with FREQ, reads from first step

    This two-step approach prevents deltaSCF collapse for higher excited states.
    """
    state_upper = state.upper()

    # For S0, just create a standard TDDFT input
    if state_upper == "S0":
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

    # For T1 (lowest triplet), use simple UKS mult=3 like S0 (no two-step needed)
    if state_upper == "T1":
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

    # For excited states (S1+, T2+): create two-step inputs
    # Step 1: TDDFT OPT (no FREQ)
    config_step1 = config.copy()
    config_step1['ESD_frequency'] = 'no'  # Disable FREQ for first step

    # Temporarily modify state to create first input
    first_input_base = f"{state_upper}_first_TDDFT"

    # Create TDDFT input for first optimization
    # We need to manually build this since we need custom filename
    from pathlib import Path
    input_file_first = esd_dir / f"{first_input_base}.inp"

    # Build first input using TDDFT logic (simplified)
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')
    implicit_solvation = config.get('implicit_solvation_model', '')
    geom_token = str(config.get('geom_opt', 'OPT')).strip() or "OPT"
    pal = config.get('PAL', 12)
    maxcore = config.get('maxcore', 6000)
    nroots = _get_tddft_param(config, 'nroots', 15)
    tda_flag = str(_get_tddft_param(config, 'TDA', 'FALSE')).upper()
    esd_maxdim = _get_tddft_param(config, 'maxdim', None)
    maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
    tddft_maxiter = _resolve_tddft_maxiter(config)
    followiroot = str(_get_tddft_param(config, 'followiroot', 'true')).lower() in ('true', 'yes', '1', 'on')

    # Determine multiplicity and iroot
    # For TDDFT (first step): always use multiplicity 1 (like pure TDDFT mode)
    # TDDFT computes triplet excitations via irootmult=triplet, not via multiplicity 3
    multiplicity = 1
    # Extract root number: T2→2, S1→1, etc.
    root_num = int(state_upper[1:]) if len(state_upper) > 1 else 1

    # Build solvation keyword
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)

    # Build keyword line for first step (TDDFT OPT, no FREQ)
    # Use RKS for TDDFT (required for dosoc, identical to pure TDDFT mode)
    keywords_first = [functional, "RKS", main_basisset, disp_corr, ri_jkx, aux_jk]
    if solvation_kw:
        keywords_first.append(solvation_kw)
    keywords_first.append(geom_token)
    keywords_first.append("MOREAD")  # Read orbitals from S0.gbw for better initial guess
    # NO FREQ/numFREQ for first step

    # Read coordinates from S0.xyz
    xyz_path = esd_dir / "S0.xyz"
    try:
        with open(xyz_path, 'r', encoding='utf-8') as f:
            all_lines = f.readlines()
            coord_lines = all_lines[2:]  # Skip atom count and comment
    except FileNotFoundError:
        logger.error(f"Coordinate file not found: {xyz_path}")
        raise

    coord_lines = _apply_esd_newgto(
        coord_lines, found_metals=metals, metal_basisset=metal_basisset, config=config
    )

    # Write first TDDFT input
    with open(input_file_first, 'w', encoding='utf-8') as f:
        f.write("! " + " ".join(keywords_first) + "\n")
        f.write(f'%base "{first_input_base}"\n')
        f.write('%moinp "S0.gbw"\n')
        f.write(f"%pal nprocs {pal} end\n")
        f.write(f"%maxcore {maxcore}\n")

        # Add %geom block for first TDDFT step (configurable via hybrid1_geom_MaxIter)
        geom_maxiter = config.get('hybrid1_geom_MaxIter', 2)
        f.write("%geom\n")
        f.write(f"  MaxIter {geom_maxiter}\n")
        f.write("end\n")

        # Output blocks
        for block in collect_output_blocks(config, allow=True):
            f.write(block + "\n")

        # TDDFT block (identical to TDDFT mode)
        # Read TDDFT_SOC setting (same as in TDDFT mode)
        dosoc_flag = str(_get_tddft_param(config, 'SOC', 'false')).strip().lower()
        dosoc_value = "true" if dosoc_flag in ('yes', 'true', '1', 'on') else "false"

        # Determine irootmult based on state type
        if state_upper.startswith('T'):
            irootmult = "triplet"
        else:
            irootmult = "singlet"

        f.write("\n%tddft\n")
        f.write(f"  nroots {nroots}\n")
        f.write(f"  maxdim {maxdim}\n")
        f.write(f"  tda {tda_flag}\n")
        if tddft_maxiter is not None:
            f.write(f"  maxiter {tddft_maxiter}\n")
        if state_upper.startswith('T'):
            f.write("  triplets true\n")
        f.write(f"  iroot {root_num}\n")
        f.write(f"  irootmult {irootmult}\n")
        if followiroot:
            f.write("  followiroot true\n")
        f.write(f"  dosoc {dosoc_value}\n")
        f.write("end\n")

        # Coordinates
        f.write(f"\n* xyz {charge} {multiplicity}\n")
        for line in coord_lines:
            f.write(line)
        f.write("*\n")

    logger.info(f"Created hybrid step 1 (TDDFT): {input_file_first}")

    # Step 2: deltaSCF OPT with FREQ
    input_file_second = esd_dir / f"{state_upper}_second.inp"

    # Now create deltaSCF input that reads from first step
    # Determine if this state uses deltaSCF
    use_deltascf = True  # All excited states in hybrid mode use deltaSCF for second step

    # Get frequency settings
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')
    freq_type = str(config.get('freq_type', 'FREQ')).strip().upper()
    if freq_type not in ('FREQ', 'NUMFREQ'):
        freq_type = 'FREQ'

    # Build keywords for second step (deltaSCF)
    keywords_second = [functional, "UKS", main_basisset, disp_corr, ri_jkx, aux_jk]
    if solvation_kw:
        keywords_second.append(solvation_kw)
    keywords_second.append(geom_token)

    # Add FREQ for second step
    if esd_frequency_enabled:
        keywords_second.append(freq_type)

    keywords_second.append("deltaSCF")
    keywords_second.append("MOREAD")

    # Add deltaSCF-specific keywords
    deltascf_keywords = config.get('deltaSCF_keywords', '')
    if deltascf_keywords:
        keywords_second.extend(deltascf_keywords.split())

    # Get deltaSCF SCF parameters
    domom = str(config.get('deltaSCF_DOMOM', 'true')).lower()
    pmom = str(config.get('deltaSCF_PMOM', 'true')).lower()
    keepinitialref = str(config.get('deltaSCF_keepinitialref', 'true')).lower()
    soscfhessup = config.get('deltaSCF_SOSCFHESSUP', 'LSR1')
    maxiter_scf = config.get('deltaSCF_maxiter', 1000)
    soscf_convfactor = config.get('deltaSCF_SOSCFConvFactor', 500)
    soscf_maxstep = config.get('deltaSCF_SOSCFMaxStep', 0.1)

    # Try to derive alphaconf/betaconf from {state}_first_TDDFT.out (not S0.out!)
    alphaconf = None
    betaconf = None
    first_out_path = esd_dir / f"{first_input_base}.out"
    if first_out_path.exists():
        result = _parse_tddft_and_derive_deltascf(first_out_path, state_upper)
        if result:
            alphaconf, betaconf = result

    # Fallback configs - should rarely be used if TDDFT parsing works
    if alphaconf is None:
        logger.warning(f"Could not parse ALPHACONF/BETACONF from {first_out_path}, using fallback")
        if state_upper.startswith('S'):
            alphaconf = "0,1"
            betaconf = "0"
        elif state_upper.startswith('T'):
            # Fallback for triplets: assume HOMO→LUMO (like T1)
            alphaconf = "1,1"
            betaconf = "0"
        else:
            alphaconf = "0,1"
            betaconf = "0"

    # Determine multiplicity for deltaSCF step
    if state_upper.startswith('T'):
        multiplicity_deltascf = 3
    else:
        multiplicity_deltascf = 1

    # Write second deltaSCF input
    with open(input_file_second, 'w', encoding='utf-8') as f:
        f.write("! " + " ".join(keywords_second) + "\n")
        f.write(f'%base "{state_upper}_second_deltaSCF"\n')
        f.write(f'%moinp "{first_input_base}.gbw"\n')
        f.write(f"%pal nprocs {pal} end\n")
        f.write(f"%maxcore {maxcore}\n")

        # Output blocks
        for block in collect_output_blocks(config, allow=True):
            f.write(block + "\n")

        # SCF block for deltaSCF
        f.write("\n%scf\n")
        f.write(f"  DOMOM {domom}\n")
        f.write(f"  pmom {pmom}\n")
        f.write(f"  keepinitialref {keepinitialref}\n")
        f.write(f"  alphaconf {alphaconf}\n")
        f.write(f"  betaconf {betaconf}\n")
        f.write(f"  SOSCFHESSUP {soscfhessup}\n")
        f.write(f"  maxiter {maxiter_scf}\n")
        f.write(f"  SOSCFConvFactor {soscf_convfactor}\n")
        f.write(f"  SOSCFMaxStep {soscf_maxstep}\n")
        f.write("end\n")

        # Coordinates - read from first step output
        # For now, we'll use the same coords as first step; ORCA will update from .gbw
        f.write(f"\n* xyzfile {charge} {multiplicity_deltascf} {first_input_base}.xyz\n")

    logger.info(f"Created hybrid step 2 (deltaSCF): {input_file_second}")

    # Return the second input file as the "main" one
    return str(input_file_second)


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
    implicit_solvation = config.get("implicit_solvation_model", "")
    geom_token_raw = config.get("geom_opt", "OPT")
    geom_token = str(geom_token_raw).strip() or "OPT"
    pal = config.get("PAL", 12)
    maxcore = config.get("maxcore", 6000)
    nroots = _get_tddft_param(config, "nroots", 15)
    tda_flag = str(_get_tddft_param(config, "TDA", config.get("TDA", "FALSE"))).upper()
    # Use TDDFT_maxdim if set, otherwise default to nroots/2 (min 5)
    esd_maxdim = _get_tddft_param(config, "maxdim", None)
    maxdim = esd_maxdim if esd_maxdim is not None else max(5, int(nroots / 2))
    tddft_maxiter = _resolve_tddft_maxiter(config)
    followiroot = str(_get_tddft_param(config, "followiroot", "true")).lower() in ("true", "yes", "1", "on")
    esd_frequency_enabled = str(config.get('ESD_frequency', 'yes')).strip().lower() in ('yes', 'true', '1', 'on')
    output_blocks = collect_output_blocks(config, allow=True)

    # Build solvation keyword once
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)

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
        ]
        if solvation_kw:
            kw.append(solvation_kw)
        kw.append(geom_token)
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
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "S3":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "S4":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "S5":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "S6":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T1":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T2":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T3":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T4":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T5":
        xyz_path = esd_dir / "S0.xyz"
        skip_lines = 2
    elif state_upper == "T6":
        xyz_path = esd_dir / "S0.xyz"
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

    coord_lines = _apply_esd_newgto(
        coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

    def _write_tddft_block(
        fh,
        iroot: Optional[int] = None,
        irootmult: Optional[str] = None,
        *,
        triplets: bool = False,
    ) -> None:
        # Read TDDFT_SOC setting (legacy: ESD_SOC)
        dosoc_flag = str(_get_tddft_param(config, 'SOC', 'false')).strip().lower()
        dosoc_value = "true" if dosoc_flag in ('yes', 'true', '1', 'on') else "false"

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
        fh.write(f"  dosoc {dosoc_value}\n")
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
            ]
            if solvation_kw:
                keywords.append(solvation_kw)
            if geom_token:
                keywords.append(geom_token)
            if esd_frequency_enabled:
                keywords.append("FREQ")
            f.write("! " + " ".join(keywords) + "\n")
            f.write('%base "S0"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)

            # Add custom additions for S0 state if specified in CONTROL
            # Generate %elprop block from elprop_properties config
            elprop_block = generate_elprop_block(config)
            if elprop_block:
                f.write(f"{elprop_block}\n")

            # Legacy addition_S0 support (can be used together with elprop_properties)
            addition_s0 = config.get('addition_S0', '').strip()
            if addition_s0:
                f.write(f"{addition_s0}\n")

            f.write(f"\n* xyz {charge} 1\n")
            for line in coord_lines:
                f.write(line)
            f.write("*\n\n")

            f.write("$new_job\n")
            # Build TDDFT keywords without OPT and numFREQ
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
            f.write("! " + " ".join(tddft_keywords) + "\n")
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
            # T1 (lowest triplet) uses simple UKS + mult 3 optimization (like S0)
            # NO deltaSCF - just a ground state triplet optimization
            f.write("! " + _join_keywords(_build_keywords("UKS")) + " MOREAD\n")
            f.write('%base "T1"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            f.write(f"\n* xyz {charge} 3\n")  # Multiplicity 3 for lowest triplet
            for line in coord_lines:
                f.write(line)
            f.write("*\n")

            # Add TDDFT check job for T1 (like in hybrid1/deltaSCF mode)
            f.write("\n")
            f.write("#==========================================\n")
            f.write("# TDDFT Check: Transitions from T1\n")
            f.write("#==========================================\n")
            f.write("\n")
            f.write("$new_job\n")

            # TDDFT keyword line - RKS for TDDFT check
            tddft_keywords_check = [
                functional,
                "RKS",
                main_basisset,
                disp_corr,
                ri_jkx,
                aux_jk,
            ]
            if solvation_kw:
                tddft_keywords_check.append(solvation_kw)
            f.write("! " + " ".join(tddft_keywords_check) + "\n")

            f.write('%base "T1_TDDFT"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # TDDFT block for check job
            f.write("\n%tddft\n")
            f.write(f"  nroots {nroots}\n")
            f.write(f"  maxdim {maxdim}\n")
            f.write(f"  tda {tda_flag}\n")
            if tddft_maxiter is not None:
                f.write(f"  maxiter {tddft_maxiter}\n")
            f.write("  triplets true\n")
            f.write("  dosoc false\n")
            f.write("end\n")

            f.write(f"\n* xyzfile {charge} 1 T1.xyz\n")
        elif state_upper == "T2":
            f.write("! " + _join_keywords(_build_keywords("RKS")) + " MOREAD\n")
            f.write('%base "T2"\n')
            f.write('%moinp "S0.gbw"\n')
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")
            _write_output_blocks(f)
            _write_tddft_block(f, iroot=2, irootmult="triplet", triplets=True)
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
            _write_tddft_block(f, iroot=3, irootmult="triplet", triplets=True)
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
            _write_tddft_block(f, iroot=4, irootmult="triplet", triplets=True)
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
            _write_tddft_block(f, iroot=5, irootmult="triplet", triplets=True)
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
            _write_tddft_block(f, iroot=6, irootmult="triplet", triplets=True)
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

    # Get ESD mode to resolve correct file names
    esd_mode = str(config.get('ESD_modus', 'tddft')).strip().lower()
    if "|" in esd_mode:
        esd_mode = esd_mode.split("|")[0].strip()

    # Determine source geometry (use optimized geometry of FINAL state, per ORCA manual)
    # In hybrid1 mode: S1_second.xyz, T2_second.xyz, etc.
    xyz_file = _resolve_state_filename(final_state, 'xyz', esd_mode)

    # Use restricted (closed-shell) reference for SOC; keep multiplicity 1 to avoid UKS
    final_multiplicity = 1

    # Calculate adiabatic energy difference (DELE) for ISC
    # DELE = E(initial) - E(final) in cm^-1
    # In hybrid1 mode: use S1_second.out, T2_second.out, etc.
    initial_out = _resolve_state_filename(initial_state, 'out', esd_mode)
    final_out = _resolve_state_filename(final_state, 'out', esd_mode)
    dele = calculate_dele_cm1(
        str(esd_dir / initial_out),
        str(esd_dir / final_out),
    )

    # Build input
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')
    implicit_solvation = config.get('implicit_solvation_model', '')

    # Simple keyword line (restricted reference)
    keywords = [
        "RKS",
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
    ]

    # Add solvation keyword only if model is set
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)
    if solvation_kw:
        keywords.append(solvation_kw)

    keywords.append("ESD(ISC)")

    simple_line = "! " + " ".join(keywords)

    # Blocks
    blocks = []

    # Base
    blocks.append(f'%base "{job_name}"')

    # TDDFT block (aligned with reference layout)
    nroots = config.get('ESD_ISC_NROOTS', config.get('ESD_nroots', 10))  # Use ESD_nroots as fallback

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
    temperature = _resolve_temperature_K(config, default=298.15)
    doht_flag = str(config.get('DOHT', 'TRUE')).upper()
    lines = str(config.get("ESD_LINES", "LORENTZ")).strip().upper() or "LORENTZ"
    linew = str(config.get("ESD_LINEW", 50)).strip()
    inlinew = str(config.get("ESD_INLINEW", 250)).strip()
    npoints = str(config.get("ESD_NPOINTS", 131072)).strip()
    maxtime = str(config.get("ESD_MAXTIME", 12000)).strip()

    # Resolve Hessian file names for hybrid1 mode
    initial_hess = _resolve_state_filename(initial_state, 'hess', esd_mode)
    final_hess = _resolve_state_filename(final_state, 'hess', esd_mode)

    esd_block = [
        "%ESD",
        f'  ISCISHESS       "{initial_hess}"',
        f'  ISCFSHESS       "{final_hess}"',
        "  USEJ            TRUE",
        f"  DOHT            {doht_flag}",
        f"  LINES           {lines}",
        f"  LINEW           {linew}",
        f"  INLINEW         {inlinew}",
        f"  NPOINTS         {npoints}",
        f"  MAXTIME         {maxtime}",
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
    coord_lines = _apply_esd_newgto(
        coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

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

    # Get ESD mode to resolve correct file names
    esd_mode = str(config.get('ESD_modus', 'tddft')).strip().lower()
    if "|" in esd_mode:
        esd_mode = esd_mode.split("|")[0].strip()

    # Determine source geometry (use lower-state geometry for IC, per ORCA manual)
    # For S1>S0: use S0.xyz; for Tn>T1: use T1.xyz
    # In hybrid1 mode: S1_second.xyz, T2_second.xyz, etc.
    xyz_file = _resolve_state_filename(final_state, 'xyz', esd_mode)

    # Multiplicity follows final state (triplet -> 3, singlet -> 1)
    final_type, _ = _parse_state_root(final_state)
    final_multiplicity = 3 if final_type == "T" else 1

    # Build input (same as ISC but labeled as IC)
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    aux_jk = config.get('aux_jk', 'def2/J')
    implicit_solvation = config.get('implicit_solvation_model', '')

    # Simple keyword line (no RKS/UKS flag - let ORCA decide based on multiplicity)
    keywords = [
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
    ]

    # Add solvation keyword only if model is set
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)
    if solvation_kw:
        keywords.append(solvation_kw)

    keywords.append("ESD(IC)")

    simple_line = "! " + " ".join(keywords)

    # Blocks
    blocks = []

    # Base
    blocks.append(f'%base "{job_name}"')

    # TDDFT block tailored for IC calculations
    nroots = config.get('ESD_IC_NROOTS', config.get('ESD_nroots', 10))  # Use ESD_nroots as fallback

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
    temperature = _resolve_temperature_K(config, default=298.15)
    lines = str(config.get("ESD_LINES", "LORENTZ")).strip().upper() or "LORENTZ"
    linew = str(config.get("ESD_LINEW", 50)).strip()
    inlinew = str(config.get("ESD_INLINEW", 250)).strip()
    npoints = str(config.get("ESD_NPOINTS", 131072)).strip()
    maxtime = str(config.get("ESD_MAXTIME", 12000)).strip()

    # Resolve Hessian file names for hybrid1 mode
    final_hess = _resolve_state_filename(final_state, 'hess', esd_mode)
    initial_hess = _resolve_state_filename(initial_state, 'hess', esd_mode)

    esd_block = [
        "%ESD",
        f'  GSHESSIAN       "{final_hess}"',
        f'  ESHESSIAN       "{initial_hess}"',
        "  USEJ            TRUE",
        f"  LINES           {lines}",
        f"  LINEW           {linew}",
        f"  INLINEW         {inlinew}",
        f"  NPOINTS         {npoints}",
        f"  MAXTIME         {maxtime}",
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
    coord_lines = _apply_esd_newgto(
        coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

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


def create_fluor_input(
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
    *,
    initial_state: str = "S1",
    final_state: str = "S0",
) -> str:
    """Generate ORCA input file for fluorescence rate calculation (ESD(FLUOR)).

    This is primarily used when CONTROL sets emission_rates=f.

    Args:
        esd_dir: ESD working directory
        charge: Molecular charge
        solvent: Solvent name
        metals: List of metal atoms (kept for signature consistency)
        main_basisset: Main basis set
        metal_basisset: Metal basis set (kept for signature consistency)
        config: Configuration dictionary
        initial_state: Excited singlet state (default: S1)
        final_state: Ground singlet state (default: S0)

    Returns:
        Path to generated input file
    """
    initial_state = initial_state.strip().upper()
    final_state = final_state.strip().upper()
    init_type, init_root = _parse_state_root(initial_state)
    final_type, _ = _parse_state_root(final_state)

    if init_type != "S" or final_state != "S0" or final_type != "S":
        raise ValueError(f"Fluorescence only supported for Sn→S0 (got {initial_state}→{final_state})")

    job_name = f"{initial_state}_{final_state}_FLUOR"
    input_file = esd_dir / f"{job_name}.inp"

    functional = config.get("functional", "PBE0")
    disp_corr = config.get("disp_corr", "D4")
    ri_jkx = config.get("ri_jkx", "RIJCOSX")
    aux_jk = config.get("aux_jk", "def2/J")
    implicit_solvation = config.get("implicit_solvation_model", "")

    keywords = [
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
    ]
    fluor_keywords = str(config.get("fluor_keywords", "")).strip()
    if fluor_keywords:
        keywords.append(fluor_keywords)

    keywords = [k for k in keywords if str(k).strip()]

    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)
    if solvation_kw:
        keywords.append(solvation_kw)

    keywords.append("ESD(FLUOR)")
    simple_line = "! " + " ".join(keywords)

    blocks: list[str] = []
    blocks.append(f'%base "{job_name}"')

    # TDDFT block (minimal, with configurable nroots)
    nroots = int(config.get("ESD_FLUOR_NROOTS", config.get("ESD_nroots", 15)))
    iroot = int(config.get("ESD_FLUOR_IROOT", init_root))
    tddft_block = [
        "%TDDFT",
        f"  NROOTS     {nroots}",
        f"  IROOT      {iroot}",
        "END",
    ]
    tddft_maxiter = _resolve_tddft_maxiter(config)
    if tddft_maxiter is not None:
        tddft_block.insert(-1, f"  maxiter    {tddft_maxiter}")
    blocks.append("\n".join(tddft_block))

    # Get ESD mode to resolve correct file names for hybrid1
    esd_mode = str(config.get('ESD_modus', 'tddft')).strip().lower()
    if "|" in esd_mode:
        esd_mode = esd_mode.split("|")[0].strip()

    # Resolve Hessian file names for hybrid1 mode
    final_hess = _resolve_state_filename(final_state, 'hess', esd_mode)
    initial_hess = _resolve_state_filename(initial_state, 'hess', esd_mode)

    # ESD block
    doht_flag = str(config.get("DOHT", "TRUE")).upper()
    lines = str(config.get("ESD_LINES", "LORENTZ")).strip().upper() or "LORENTZ"
    linew = str(config.get("ESD_LINEW", 50)).strip()
    inlinew = str(config.get("ESD_INLINEW", 250)).strip()
    temperature = _resolve_temperature_K(config, default=298.15)
    esd_block = [
        "%ESD",
        f'  GSHESSIAN  "{final_hess}"',
        f'  ESHESSIAN  "{initial_hess}"',
        f"  DOHT       {doht_flag}",
        f"  LINES      {lines}",
        f"  LINEW      {linew}",
        f"  INLINEW    {inlinew}",
        f"  TEMP       {temperature}",
        "END",
    ]
    blocks.append("\n".join(esd_block))

    # PAL and maxcore
    pal = config.get("PAL", 12)
    maxcore = config.get("maxcore", 6000)
    blocks.append(f"%pal nprocs {pal} end")
    blocks.append(f"%maxcore {maxcore}")

    # Geometry uses the optimized ground-state geometry from the ESD directory.
    # Use inline coordinates so we can attach per-atom NewGTO for metals (like classic inputs).
    # In hybrid1 mode: use S0.xyz (S0 is always simple naming)
    final_xyz = _resolve_state_filename(final_state, 'xyz', esd_mode)
    xyz_path = esd_dir / final_xyz
    with open(xyz_path, "r", encoding="utf-8") as f:
        all_lines = f.readlines()
    coord_lines = all_lines[2:]  # skip XYZ header
    coord_lines = _apply_esd_newgto(
        coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

    with open(input_file, "w", encoding="utf-8") as f:
        f.write(simple_line + "\n")
        for block in blocks:
            f.write(block + "\n")
        f.write("\n")
        f.write(f"* xyz {charge} 1\n")
        for line in coord_lines:
            f.write(line if line.endswith("\n") else line + "\n")
        f.write("*\n")

    logger.info(f"Created FLUOR input: {input_file}")
    return str(input_file)


def create_phosp_input(
    esd_dir: Path,
    charge: int,
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
    config: Dict[str, Any],
    *,
    initial_state: str = "T1",
    final_state: str = "S0",
) -> str:
    """Generate ORCA input file for phosphorescence rate/spectrum calculation (ESD(PHOSP)).

    Notes (ORCA manual):
    - Requires SOC and the ESD(PHOSP) module.
    - Needs GSHESSIAN (ground singlet) and TSHESSIAN (triplet Hessian).
    - DELE should be the adiabatic energy difference between S0 and T1 at their respective geometries
      (electronic energies, no ZPE correction), given in cm^-1.
    - For triplet sub-levels after SOC splitting, compute IROOT=1,2,3 (often via $new_job).
    """
    initial_state = initial_state.strip().upper()
    final_state = final_state.strip().upper()
    init_type, init_root = _parse_state_root(initial_state)
    final_type, _ = _parse_state_root(final_state)

    if init_type != "T" or initial_state != "T1" or final_state != "S0" or final_type != "S":
        raise ValueError(f"Phosphorescence only supported for T1→S0 (got {initial_state}→{final_state})")

    job_name = f"{initial_state}_{final_state}_PHOSP"
    input_file = esd_dir / f"{job_name}.inp"

    functional = config.get("functional", "PBE0")
    disp_corr = config.get("disp_corr", "D4")
    ri_jkx = config.get("ri_jkx", "RIJCOSX")
    aux_jk = config.get("aux_jk", "def2/J")
    implicit_solvation = config.get("implicit_solvation_model", "")

    keywords = [
        functional,
        main_basisset,
        disp_corr,
        ri_jkx,
        aux_jk,
        "ESD(PHOSP)",
    ]
    # Optional extra keywords for the main keyword line (e.g. "RI-SOMF(1X)")
    phosp_keywords = str(config.get("phosp_keywords", "")).strip()
    if phosp_keywords:
        keywords.append(phosp_keywords)

    keywords = [k for k in keywords if str(k).strip()]

    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)
    if solvation_kw:
        # Keep solvation close to the method tokens (before ESD(PHOSP))
        try:
            esd_idx = keywords.index("ESD(PHOSP)")
        except ValueError:
            esd_idx = len(keywords)
        keywords.insert(esd_idx, solvation_kw)

    simple_line = "! " + " ".join(keywords)

    # Compute DELE (cm^-1) from electronic energies at optimized geometries (no ZPE)
    # Use T1.out and S0.out produced by state jobs.
    t1_out = esd_dir / "T1.out"
    s0_out = esd_dir / "S0.out"
    dele_cm1 = calculate_dele_cm1(str(t1_out), str(s0_out))
    dele_int = int(round(dele_cm1)) if dele_cm1 is not None else None

    # Shared settings
    doht_flag = str(config.get("DOHT", "TRUE")).upper()
    temperature = _resolve_temperature_K(config, default=298.15)
    tda_flag = str(config.get("ESD_TDA", config.get("TDA", "FALSE"))).upper()
    # Use the general ESD_nroots by default (CONTROL), allow PHOSP override if desired.
    nroots = int(config.get("ESD_PHOSP_NROOTS", config.get("ESD_nroots", 15)))
    tddft_maxiter = _resolve_tddft_maxiter(config)

    # ORCA recommends DOSOC TRUE for phosphorescence
    dosoc_flag = str(config.get("ESD_PHOSP_DOSOC", "TRUE")).upper()
    lines = str(config.get("ESD_LINES", "LORENTZ")).strip().upper() or "LORENTZ"
    linew = str(config.get("ESD_LINEW", 50)).strip()
    inlinew = str(config.get("ESD_INLINEW", 250)).strip()

    # Which IROOT subjobs to run (triplet SOC-split components).
    # Default is 1,2,3, but allow overriding via CONTROL.
    iroots = _parse_iroot_spec(config.get("phosp_IROOT", None), default=[1, 2, 3])

    # Get ESD mode to resolve correct file names for hybrid1
    esd_mode = str(config.get('ESD_modus', 'tddft')).strip().lower()
    if "|" in esd_mode:
        esd_mode = esd_mode.split("|")[0].strip()

    # Geometry uses the ground-state geometry (S0.xyz), multiplicity 1.
    # Use inline coordinates so we can attach per-atom NewGTO for metals (like classic inputs).
    # Note: For PHOSP (T1→S0), both states use simple naming even in hybrid1
    final_xyz = _resolve_state_filename(final_state, 'xyz', esd_mode)
    xyz_path = esd_dir / final_xyz
    with open(xyz_path, "r", encoding="utf-8") as f:
        all_lines = f.readlines()
    geom_coord_lines = all_lines[2:]  # skip XYZ header
    geom_coord_lines = _apply_esd_newgto(
        geom_coord_lines,
        found_metals=metals,
        metal_basisset=metal_basisset,
        config=config,
    )

    # Resolve Hessian file names for hybrid1 mode
    # Note: For PHOSP (T1→S0), both states use simple naming even in hybrid1
    final_hess = _resolve_state_filename(final_state, 'hess', esd_mode)
    initial_hess = _resolve_state_filename(initial_state, 'hess', esd_mode)

    def _job_block(iroot: int) -> str:
        blocks: list[str] = []
        blocks.append(simple_line)
        blocks.append(f'%base "{job_name}_iroot{iroot}"')
        tddft = [
            "%TDDFT",
            f"  NROOTS  {nroots}",
            f"  DOSOC   {dosoc_flag}",
            f"  TDA     {tda_flag}",
            f"  IROOT   {iroot}",
        ]
        if tddft_maxiter is not None:
            tddft.append(f"  maxiter {tddft_maxiter}")
        tddft.append("END")
        blocks.append("\n".join(tddft))

        esd = [
            "%ESD",
            f'  GSHESSIAN       "{final_hess}"',
            f'  TSHESSIAN       "{initial_hess}"',
            f"  DOHT            {doht_flag}",
            f"  LINES           {lines}",
            f"  LINEW           {linew}",
            f"  INLINEW         {inlinew}",
            f"  TEMP            {temperature}",
        ]
        if dele_int is not None:
            esd.append(f"  DELE            {dele_int}")
        esd.append("END")
        blocks.append("\n".join(esd))

        blocks.append(f"* xyz {charge} 1")
        blocks.extend([ln.rstrip("\n") for ln in geom_coord_lines])
        blocks.append("*")
        return "\n".join(blocks)

    pal = config.get("PAL", 12)
    maxcore = config.get("maxcore", 6000)

    with open(input_file, "w", encoding="utf-8") as f:
        first = True
        for iroot in iroots:
            if not first:
                f.write("\n$new_job\n\n")
            first = False
            f.write(_job_block(iroot) + "\n")
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

    logger.info(f"Created PHOSP input: {input_file}")
    return str(input_file)


def append_properties_of_interest_jobs(
    inp_file: str,
    xyz_file: str,
    base_charge: int,
    base_multiplicity: int,
    properties: str,
    config: Dict[str, Any],
    solvent: str,
    metals: List[str],
    main_basisset: str,
    metal_basisset: str,
) -> None:
    """Append $new_job blocks for IP and EA calculations to an existing input file.

    Args:
        inp_file: Path to the input file to modify
        xyz_file: Path to the xyz geometry file (e.g., 'S0.xyz')
        base_charge: Base charge of the system
        base_multiplicity: Base multiplicity of the system
        properties: Comma-separated list of properties (e.g., 'IP,EA')
        config: Configuration dictionary
        solvent: Solvent name
        metals: List of metal atoms
        main_basisset: Main basis set
        metal_basisset: Metal basis set
    """
    if not properties:
        return

    # Parse properties list - handle both string and list input
    if isinstance(properties, (list, tuple)):
        prop_list = [p.strip().upper() for p in properties]
    else:
        # Handle string that might look like "['IP', 'EA']" or "IP,EA"
        prop_str = str(properties).strip()
        # Remove list brackets and quotes if present
        prop_str = prop_str.strip('[]').replace("'", "").replace('"', '')
        prop_list = [p.strip().upper() for p in prop_str.split(',')]

    if not any(p in ['IP', 'EA'] for p in prop_list):
        return

    logger.info(f"Adding properties_of_interest jobs to {inp_file}: {prop_list}")

    # Get settings from config
    functional = config.get('functional', 'PBE0')
    disp_corr = config.get('disp_corr', 'D4')
    ri_jkx = config.get('ri_jkx', 'RIJCOSX')
    relativity = str(config.get('relativity', 'none')).strip().lower()
    implicit_solvation = config.get('implicit_solvation_model', '')
    pal = config.get('PAL', 6)
    maxcore = config.get('maxcore', 6000)
    maxiter = config.get('maxiter', 125)

    # Build solvation keyword
    solvation_kw = _build_solvation_keyword(implicit_solvation, solvent)

    # Determine basis sets and aux basis based on relativity
    if relativity in ['zora', 'dkh', 'dkh2']:
        # Use relativistic basis sets
        if relativity == 'zora':
            actual_main_basis = config.get('main_basisset_rel', f'ZORA-{main_basisset}')
            aux_jk = config.get('aux_jk_rel', 'SARC/J')
        else:  # DKH
            actual_main_basis = config.get('main_basisset_rel', f'DKH-{main_basisset}')
            aux_jk = config.get('aux_jk', 'def2/J')
        relativity_keyword = relativity.upper()
    else:
        actual_main_basis = main_basisset
        aux_jk = config.get('aux_jk', 'def2/J')
        relativity_keyword = None

    # Calculate total electrons for multiplicity determination
    # For now, we use simple rules: odd electrons -> mult=2, even -> mult=1
    # This is a simplification and may need refinement

    jobs_to_add = []

    if 'IP' in prop_list:
        # Ionization Potential: remove one electron (charge +1)
        ip_charge = base_charge + 1
        # Always use doublet (mult=2) for IP/EA
        ip_mult = 2
        jobs_to_add.append(('IP', ip_charge, ip_mult))

    if 'EA' in prop_list:
        # Electron Affinity: add one electron (charge -1)
        ea_charge = base_charge - 1
        # Always use doublet (mult=2) for IP/EA
        ea_mult = 2
        jobs_to_add.append(('EA', ea_charge, ea_mult))

    # Build output blocks (same as main job)
    output_blocks = collect_output_blocks(config)

    # Append to input file
    with open(inp_file, 'a') as f:
        for prop_name, charge, mult in jobs_to_add:
            f.write("\n")
            f.write(f"$new_job\n")

            # Build keyword line (without OPT/FREQ)
            keywords = [functional, actual_main_basis, disp_corr, ri_jkx, aux_jk]
            if relativity_keyword:
                keywords.append(relativity_keyword)
            # Add solvation keyword only if model is set
            if solvation_kw:
                keywords.append(solvation_kw)

            f.write("! " + " ".join(keywords) + "\n")

            # Base name for output files
            base_name = Path(inp_file).stem  # e.g., 'S0'
            f.write(f'%base "{base_name}_{prop_name}"\n')

            # PAL and maxcore
            f.write(f"%pal nprocs {pal} end\n")
            f.write(f"%maxcore {maxcore}\n")

            # SCF settings
            f.write(f"%scf maxiter {maxiter} end\n")

            # Output blocks
            for block in output_blocks:
                f.write(block + "\n")

            # Add %basis block with NewGTO for metals if present
            if metals and metal_basisset:
                f.write("\n%basis\n")
                for metal in metals:
                    f.write(f'  NewGTO {metal} "{metal_basisset}" end\n')
                f.write("end\n")

            # Geometry - use xyzfile to reference optimized geometry (filename only, no path)
            xyz_filename = Path(xyz_file).name
            f.write(f"\n* xyzfile {charge} {mult} {xyz_filename}\n")

    logger.info(f"Added {len(jobs_to_add)} properties_of_interest job(s) to {inp_file}")
