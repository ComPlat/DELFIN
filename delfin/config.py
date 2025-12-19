import ast
import re
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from delfin.common.control_validator import validate_control_config
from delfin.occupier_sequences import remove_existing_sequence_blocks
from delfin.define import TEMPLATE as CONTROL_TEMPLATE

from delfin.common.logging import get_logger

logger = get_logger(__name__)


_SEQUENCE_BLOCK_HEADER = re.compile(r"^\s*([+\-0-9,\s]+)=\[\s*$")
_TEMPLATE_DEFAULTS_CACHE: Optional[Dict[str, Any]] = None
_TEMPLATE_REQUIRED_KEYS: Set[str] = set()
_PLACEHOLDER_VALIDATION_VALUES: Dict[str, Any] = {
    "charge": 0,
    "solvent": "DMF",
}


def _is_placeholder_value(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    stripped = value.strip()
    return bool(stripped.startswith("[") and stripped.endswith("]") and len(stripped) > 2)


def _collect_literal_list(lines: List[str], start_idx: int, initial_value: str) -> Tuple[List[Any], int]:
    """Collect a multi-line Python literal (list/dict) starting at start_idx."""
    buffer = initial_value + "\n"
    depth = initial_value.count('[') - initial_value.count(']')
    idx = start_idx + 1

    while idx < len(lines) and depth > 0:
        line = lines[idx]
        buffer += line
        depth += line.count('[') - line.count(']')
        idx += 1

    try:
        parsed = ast.literal_eval(buffer)
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Failed to parse list literal near line {start_idx + 1}: {exc}")
        parsed = []

    return parsed if isinstance(parsed, list) else [], idx


def _parse_delta_tokens(raw: str) -> List[int]:
    deltas: List[int] = []
    seen: set[int] = set()
    for token in raw.split(','):
        text = token.strip()
        if not text:
            continue
        try:
            value = int(text)
        except ValueError:
            continue
        if value in seen:
            continue
        seen.add(value)
        deltas.append(value)
    return deltas


def _parse_sequence_block(lines: List[str], start_idx: int) -> Tuple[Dict[str, Any], int]:
    block: Dict[str, Any] = {"even_seq": [], "odd_seq": []}
    idx = start_idx

    while idx < len(lines):
        line = lines[idx]
        stripped = line.strip()

        if not stripped or stripped.startswith("#") or stripped.startswith("---") or stripped.startswith("***"):
            idx += 1
            continue
        if stripped.endswith(":"):
            idx += 1
            continue
        if stripped == "]":
            return block, idx + 1

        if stripped.startswith("even_seq"):
            parts = line.split("=", 1)
            if len(parts) != 2:
                idx += 1
                continue
            value = parts[1].strip()
            if value.startswith('[') and not value.endswith(']'):
                parsed, idx = _collect_literal_list(lines, idx, value)
            else:
                try:
                    parsed = ast.literal_eval(value)
                    idx += 1
                except Exception as exc:  # noqa: BLE001
                    logger.error(f"Failed to parse even_seq block near line {idx + 1}: {exc}")
                    parsed = []
                    idx += 1
            block["even_seq"] = parsed if isinstance(parsed, list) else []
            continue

        if stripped.startswith("odd_seq"):
            parts = line.split("=", 1)
            if len(parts) != 2:
                idx += 1
                continue
            value = parts[1].strip()
            if value.startswith('[') and not value.endswith(']'):
                parsed, idx = _collect_literal_list(lines, idx, value)
            else:
                try:
                    parsed = ast.literal_eval(value)
                    idx += 1
                except Exception as exc:  # noqa: BLE001
                    logger.error(f"Failed to parse odd_seq block near line {idx + 1}: {exc}")
                    parsed = []
                    idx += 1
            block["odd_seq"] = parsed if isinstance(parsed, list) else []
            continue

        idx += 1

    logger.warning("Unterminated OCCUPIER sequence block; reached EOF without closing ']'.")
    return block, idx


def _load_template_defaults() -> Dict[str, Any]:
    """Parse define.py TEMPLATE and cache defaults for missing CONTROL keys."""
    global _TEMPLATE_DEFAULTS_CACHE
    if _TEMPLATE_DEFAULTS_CACHE is None:
        parsed = _parse_control_file("<template>", keep_steps_literal=True, content=CONTROL_TEMPLATE)
        sanitized: Dict[str, Any] = {}
        for key, value in parsed.items():
            if isinstance(value, str) and "|" in value:
                value = value.split("|", 1)[0].strip()

            if _is_placeholder_value(value):
                _TEMPLATE_REQUIRED_KEYS.add(key)
                continue

            sanitized[key] = value

        validation_seed = dict(sanitized)
        for key in _TEMPLATE_REQUIRED_KEYS:
            if key in validation_seed:
                continue
            if key not in _PLACEHOLDER_VALIDATION_VALUES:
                raise ValueError(f"No validation fallback configured for placeholder key '{key}'")
            validation_seed[key] = _PLACEHOLDER_VALIDATION_VALUES[key]

        defaults = validate_control_config(validation_seed)
        if "_occupier_sequence_blocks" in parsed:
            defaults["_occupier_sequence_blocks"] = parsed["_occupier_sequence_blocks"]

        # Override ESD-specific defaults to empty (opt-in only)
        # Users must explicitly set these in CONTROL.txt to enable ESD calculations
        esd_opt_in_keys = ['states', 'ISCs', 'ICs', 'emission_rates']
        for key in esd_opt_in_keys:
            if key in defaults:
                defaults[key] = ''

        _TEMPLATE_DEFAULTS_CACHE = defaults

    return deepcopy(_TEMPLATE_DEFAULTS_CACHE)


def _collect_missing_required_keys(user_keys: Set[str], placeholder_keys: Set[str]) -> List[str]:
    required_missing = {key for key in _TEMPLATE_REQUIRED_KEYS if key not in user_keys}
    placeholder_required = _TEMPLATE_REQUIRED_KEYS & placeholder_keys
    return sorted(required_missing | placeholder_required)


def _raise_if_missing_required(user_keys: Set[str], placeholder_keys: Set[str]) -> None:
    missing_all = _collect_missing_required_keys(user_keys, placeholder_keys)
    if missing_all:
        joined = ", ".join(missing_all)
        raise ValueError(
            f"Missing required CONTROL values for: {joined}. "
            "Replace template placeholders (e.g. [CHARGE], [SOLVENT]) with actual values."
        )


def _parse_control_file(file_path: str, *, keep_steps_literal: bool, content: Optional[str] = None) -> Dict[str, Any]:
    config: Dict[str, Any] = {}
    sequence_blocks: List[Dict[str, Any]] = []

    if content is None:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    else:
        lines = content.splitlines(keepends=True)

    idx = 0
    total_lines = len(lines)

    while idx < total_lines:
        line = lines[idx]
        stripped = line.strip()

        upper = stripped.upper()
        # Stop parsing before the documentation block at the end of CONTROL.txt.
        # Require the exact "ESD MODULE:" heading so we don't accidentally stop
        # at the real ESD configuration section ("ESD module (excited state dynamics):").
        if upper.startswith("INFOS:") or upper.startswith("ESD MODULE:"):
            # Stop parsing before template/documentation blocks to avoid overriding user values.
            break

        # Skip comments / separators / blanks
        if not stripped or stripped.startswith('#') or stripped.startswith('---') or stripped.startswith('***'):
            idx += 1
            continue

        # Sequence block definition (e.g. "-1,0,+1=[")
        block_match = _SEQUENCE_BLOCK_HEADER.match(stripped)
        if block_match:
            deltas = _parse_delta_tokens(block_match.group(1))
            block, next_idx = _parse_sequence_block(lines, idx + 1)
            if deltas and (block.get("even_seq") or block.get("odd_seq")):
                block["deltas"] = deltas
                sequence_blocks.append(block)
            idx = next_idx
            continue

        # Ignore headings like "odd electron number:"
        if ':' in stripped and not '=' in stripped:
            idx += 1
            continue

        if '=' not in line:
            idx += 1
            continue

        key_raw, value_raw = line.split('=', 1)
        key = key_raw.strip()
        value = value_raw.strip()

        if keep_steps_literal and key in ('oxidation_steps', 'reduction_steps'):
            config[key] = value
            idx += 1
            continue

        if value.startswith('[') and not value.endswith(']'):
            parsed, next_idx = _collect_literal_list(lines, idx, value)
            config[key] = parsed
            idx = next_idx
            continue

        if ',' in value and not value.startswith('{') and not value.startswith('['):
            config[key] = [v.strip() for v in value.split(',') if v.strip()]
            idx += 1
            continue

        try:
            config[key] = ast.literal_eval(value)
        except Exception:
            config[key] = value
        idx += 1

    if sequence_blocks:
        config["_occupier_sequence_blocks"] = sequence_blocks
    return config


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
    sanitized = remove_existing_sequence_blocks(Path(file_path), persist=False)
    text = sanitized or Path(file_path).read_text(encoding="utf-8")
    config = _parse_control_file(file_path, keep_steps_literal=True, content=text)
    user_keys = set(config.keys())
    placeholder_keys = {key for key, value in config.items() if _is_placeholder_value(value)}
    defaults = _load_template_defaults()
    _raise_if_missing_required(user_keys, placeholder_keys)
    validated = validate_control_config(config)
    if "_occupier_sequence_blocks" in config:
        validated["_occupier_sequence_blocks"] = config["_occupier_sequence_blocks"]
        user_keys.add("_occupier_sequence_blocks")

    merged = dict(validated)
    for key, value in defaults.items():
        if key not in user_keys:
            merged[key] = value

    return merged

def OCCUPIER_parser(path: str) -> Dict[str, Any]:
    """Parse OCCUPIER-specific configuration file.

    Similar to read_control_file but with specialized handling for OCCUPIER workflow.

    Args:
        path: Path to configuration file

    Returns:
        Dictionary containing parsed OCCUPIER configuration
    """
    sanitized = remove_existing_sequence_blocks(Path(path), persist=False)
    text = sanitized or Path(path).read_text(encoding="utf-8")
    config = _parse_control_file(path, keep_steps_literal=False, content=text)
    user_keys = set(config.keys())
    placeholder_keys = {key for key, value in config.items() if _is_placeholder_value(value)}
    defaults = _load_template_defaults()
    _raise_if_missing_required(user_keys, placeholder_keys)
    validated = validate_control_config(config)
    if "_occupier_sequence_blocks" in config:
        validated["_occupier_sequence_blocks"] = config["_occupier_sequence_blocks"]
        user_keys.add("_occupier_sequence_blocks")

    merged = dict(validated)
    for key, value in defaults.items():
        if key not in user_keys:
            merged[key] = value

    return merged


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
