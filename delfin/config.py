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


def get_E_ref(config: Dict[str, Any]) -> float:
    """Get reference electrode potential for redox calculations.

    Returns user-specified E_ref if available, otherwise looks up
    solvent-specific reference potentials vs. SHE.

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


