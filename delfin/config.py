import ast
from typing import Dict, Any, Optional

from delfin.common.control_validator import validate_control_config

from delfin.common.logging import get_logger

logger = get_logger(__name__)

def read_control_file(file_path: str) -> Dict[str, Any]:
    """Parse CONTROL.txt file and return configuration dictionary.

    Supports:
    - Key=value pairs with type inference
    - Multi-line lists in [...] format
    - Comma-separated values converted to lists
    - Comments starting with # or --- or ***

    Args:
        file_path: Path to CONTROL.txt file

    Returns:
        Dictionary containing parsed configuration parameters
    """
    config = {}
    multi_key = None
    multi_val = ""

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            # Skip comments / blank lines
            if not line or line.startswith('#') or line.startswith('---') or line.startswith('***'):
                continue

            # Continuation of multi-line list
            if multi_key:
                multi_val += line + '\n'
                if line.endswith(']'):
                    try:
                        config[multi_key] = ast.literal_eval(multi_val)
                    except Exception:
                        config[multi_key] = []
                    multi_key, multi_val = None, ""
                continue

            # Normal key=value lines
            if '=' in line:
                key, value = [x.strip() for x in line.split('=', 1)]

                # ---------- NEW: Ox/Red‑Steps always as string -----------------
                if key in ('oxidation_steps', 'reduction_steps'):
                    config[key] = value                # no type conversion
                    continue
                # ----------------------------------------------------------------

                # Start of a multi-line list
                if value.startswith('[') and not value.endswith(']'):
                    multi_key, multi_val = key, value + '\n'
                    continue

                # Comma separated values → List of strings
                if ',' in value and not value.startswith('{') and not value.startswith('['):
                    config[key] = [v.strip() for v in value.split(',') if v.strip()]
                    continue

                # Everything else: try to parse (int, float, dict …)
                try:
                    config[key] = ast.literal_eval(value)
                except Exception:
                    config[key] = value
                continue

            # Ignore section headings (with colon)
            elif ':' in line:
                continue

    validated = validate_control_config(config)
    return validated

def OCCUPIER_parser(path: str) -> Dict[str, Any]:
    """Parse OCCUPIER-specific configuration file.

    Similar to read_control_file but with specialized handling for OCCUPIER workflow.

    Args:
        path: Path to configuration file

    Returns:
        Dictionary containing parsed OCCUPIER configuration
    """
    config = {}
    multi_key = None
    multi_val = ""

    with open(path, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip comments, separators, or empty lines
            if not line or line.startswith('#') or line.startswith('---') or line.startswith('***'):
                continue

            # Handle continuation of a multi-line list
            if multi_key:
                multi_val += line + '\n'
                if line.endswith(']'):
                    try:
                        parsed = ast.literal_eval(multi_val)
                        config[multi_key] = parsed
                    except Exception as e:
                        logger.error(f"Could not parse list for {multi_key}: {e}")
                        config[multi_key] = []
                    multi_key = None
                    multi_val = ""
                continue

            # Normal key=value line
            if '=' in line:
                key, value = line.split('=', 1)
                key: str = key.strip()
                value: str = value.strip()

                # Start of a multiline list
                if value.startswith('[') and not value.endswith(']'):
                    multi_key = key
                    multi_val = value + '\n'
                    continue

                # Convert comma-separated values to list of strings
                if ',' in value and not value.startswith('{') and not value.startswith('['):
                    config[key] = [v.strip() for v in value.split(',') if v.strip()]
                else:
                    try:
                        config[key] = ast.literal_eval(value)
                    except Exception:
                        config[key] = value

            # Optional: Skip section headers like "odd electron number:"
            elif ':' in line:
                continue

    validated = validate_control_config(config)
    return validated


def _coerce_float(val: Any) -> Optional[float]:
    """Convert various types to float with robust error handling.

    Handles:
    - Integers and floats
    - String representations (including comma as decimal separator)
    - Boolean values (returns None)
    - Infinity and NaN checks

    Args:
        val: Value to convert to float

    Returns:
        Float value or None if conversion fails
    """
    if val is None:
        return None
    if isinstance(val, bool):
        return None
    if isinstance(val, (int, float)):
        try:
            from math import isfinite
            f = float(val)
            return f if isfinite(f) else None
        except Exception:
            return None
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return None
        s = s.replace(",", ".")
        try:
            return float(s)
        except ValueError:
            return None
    return None


def _is_oniom_calculation(config: Dict[str, Any]) -> bool:
    """Check if this is an ONIOM (QM/QM2) calculation.

    Detects ONIOM by checking:
    1. If input_file exists and contains ONIOM markers
    2. If any generated ORCA input files contain QM/QM2 keywords

    Args:
        config: Configuration dictionary

    Returns:
        True if ONIOM calculation is detected, False otherwise
    """
    import os
    from pathlib import Path

    # QM/MM method patterns to detect
    qmmm_patterns = ['QM/XTB', 'QM/MM', 'QM/QM2', 'QM/PBEH-3C', 'QM/HF-3C', 'QM/r2SCAN-3C']

    # Check input_file specified in config
    input_file = config.get('input_file', 'input.txt')
    if os.path.exists(input_file):
        try:
            with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read(2000)  # Read first 2000 chars
                # Check for ONIOM markers
                if any(pattern in content for pattern in qmmm_patterns):
                    return True
        except Exception:
            pass

    # Check generated ORCA input files (initial.inp, etc.)
    input_files_to_check = [
        'initial.inp',
        'ox_step_1.inp',
        'red_step_1.inp',
    ]

    for fname in input_files_to_check:
        if os.path.exists(fname):
            try:
                with open(fname, 'r', encoding='utf-8', errors='ignore') as f:
                    first_line = f.readline()
                    # ORCA input files start with "! keywords"
                    if any(pattern in first_line for pattern in qmmm_patterns):
                        return True
            except Exception:
                pass

    return False


def get_E_ref(config: Dict[str, Any]) -> float:
    """Get reference electrode potential for redox calculations.

    Returns user-specified E_ref if available, otherwise looks up
    solvent-specific reference potentials vs. SHE.

    For ONIOM calculations (QM/QM2), automatically uses adjusted E_ref values
    to account for systematic DFT errors in electrostatic interactions.

    Args:
        config: Configuration dictionary containing 'E_ref' and 'solvent'

    Returns:
        Reference electrode potential in V vs. SHE (default: 4.345 V)
    """
    e_ref = _coerce_float(config.get('E_ref', None))
    if e_ref is not None:
        return e_ref

    solvent_raw = config.get('solvent', '')
    solvent_key = solvent_raw.strip().lower() if isinstance(solvent_raw, str) else ''

    # Check if this is an ONIOM calculation
    is_oniom = _is_oniom_calculation(config)

    if is_oniom:
        # ONIOM-specific E_ref values (can be customized per solvent)
        solvent_E_ref_oniom = {
            "dmf": -3.31, "n,n-dimethylformamide": -3.31,
            "dcm": -3.31, "ch2cl2": -3.31, "dichloromethane": -3.31,
            "acetonitrile": -3.31, "mecn": -3.31,
            "thf": -3.31, "tetrahydrofuran": -3.31,
            "dmso": -3.31, "dimethylsulfoxide": -3.31,
            "dme": -3.31, "dimethoxyethane": -3.31,
            "acetone": -3.31, "propanone": -3.31,
        }

        e_ref_value = solvent_E_ref_oniom.get(solvent_key, -3.31)

        # Log the automatic adjustment
        logger.info(f"ONIOM calculation detected: Using E_ref = {e_ref_value:.3f} V for {solvent_raw or 'default solvent'}")
        logger.info(f"  → This accounts for systematic DFT errors in QM/QM2 electrostatic interactions")
        logger.info(f"  → To override, set E_ref manually in CONTROL.txt")

        return e_ref_value
    else:
        # Standard E_ref values for non-ONIOM calculations
        solvent_E_ref = {
            "dmf": 4.795, "n,n-dimethylformamide": 4.795,
            "dcm": 4.805, "ch2cl2": 4.805, "dichloromethane": 4.805,
            "acetonitrile": 4.745, "mecn": 4.745,
            "thf": 4.905, "tetrahydrofuran": 4.905,
            "dmso": 4.780, "dimethylsulfoxide": 4.780,
            "dme": 4.855, "dimethoxyethane": 4.855,
            "acetone": 4.825, "propanone": 4.825,
        }

        return solvent_E_ref.get(solvent_key, 4.345)


