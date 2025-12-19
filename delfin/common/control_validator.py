"""Schema validation for CONTROL.txt configurations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Mapping, MutableMapping


@dataclass(frozen=True)
class FieldSpec:
    name: str
    coerce: Callable[[Any], Any]
    required: bool = False
    default: Any = None
    allow_none: bool = False


def _as_int(value: Any) -> int:
    if value is None or value == "":
        raise ValueError("must be an integer")
    return int(value)


def _as_float(value: Any) -> float:
    if value is None or value == "":
        raise ValueError("must be a float")
    return float(value)


def _as_str(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _as_yes_no(value: Any) -> str:
    text = str(value or "no").strip().lower()
    return "yes" if text in {"yes", "true", "1", "on"} else "no"


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        return [item.strip() for item in value.split(',') if item.strip()]
    raise ValueError("must be a list or comma-separated string")

def _as_parallel_strategy(value: Any) -> str:
    """Coerce user value into a known ORCA parallel strategy token."""
    text = str(value or "auto").strip().lower()
    if text in {"threads", "serial", "auto"}:
        return text
    if text in {"mpi", "default"}:
        return "auto"
    raise ValueError("must be one of: auto, threads, serial")

def _as_imag_scope(value: Any) -> str:
    """Coerce user value into a known IMAG scope."""
    text = str(value or "initial").strip().lower()
    if text in {"initial", "all"}:
        return text
    raise ValueError("must be one of: initial, all")

def _as_imag_option(value: Any) -> int:
    """Coerce IMAG scheduler behaviour selector."""
    if value is None or value == "":
        return 2
    try:
        option = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("must be 1 or 2") from exc
    if option not in (1, 2):
        raise ValueError("must be 1 or 2")
    return option


def _as_occupier_method(value: Any) -> str:
    text = str(value or "auto").strip().lower()
    if text in {"manual", "manually"}:
        return "manually"
    if text == "auto":
        return "auto"
    # Default to "auto" for invalid values instead of raising error
    return "auto"


def _as_occupier_tree(value: Any) -> str:
    """Coerce user value into a known tree token (flat, deep2, deep3, deep, own)."""
    # If empty/None, use default "deep"
    if not value or str(value).strip() == "":
        return "deep"

    text = str(value).strip().lower()
    if text in {"deep", "tree"}:
        return "deep"
    if text in {"flat", "flatt", "legacy"}:
        return "flat"
    if text == "deep2":
        return "deep2"
    if text == "deep3":
        return "deep3"
    if text in {"deep4", "dee4"}:
        return "own"
    if text in {"own", "custom"}:
        return "own"
    # Legacy aliases for backwards compatibility - map to "deep"
    if text in {"deep5", "dee5", "deep6", "dee6"}:
        return "deep"
    # Default to "own" for invalid/unrecognized values (including multi-value syntax like "deep|flat|own")
    return "own"


def _as_ap_method(value: Any) -> int | None:
    """Coerce APMethod tokens into ORCA-compatible integers."""
    if value is None:
        return None

    text = str(value).strip()
    lowered = text.lower()

    if not text:
        return None
    if lowered in {"none", "off", "disable", "disabled", "false"}:
        return None
    if lowered == "apbs":
        return 2
    candidate = text
    if lowered.startswith("ap") and lowered[2:].isdigit():
        candidate = lowered[2:]

    try:
        parsed = int(candidate)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("must be an integer (e.g. 1-4) or 'none'") from exc

    if parsed <= 0:
        return None
    return parsed


CONTROL_FIELD_SPECS: Iterable[FieldSpec] = (
    FieldSpec("NAME", _as_str, default=""),
    FieldSpec("SMILES", _as_str, default=""),
    FieldSpec("charge", _as_int, required=True),
    FieldSpec("multiplicity_global_opt", _as_int, allow_none=True),
    FieldSpec("PAL", _as_int, default=6),
    FieldSpec("number_explicit_solv_molecules", _as_int, default=0),
    FieldSpec("method", _as_str, default="", allow_none=True),
    FieldSpec("frequency_calculation", _as_yes_no, default="no"),
    FieldSpec("frequency_calculation_OCCUPIER", _as_yes_no, default="no"),
    FieldSpec("xTB_method", _as_str, default="GFN2xTB"),
    FieldSpec("functional", _as_str, default="PBE0"),
    FieldSpec("main_basisset", _as_str, default="def2-SVP"),
    FieldSpec("metal_basisset", _as_str, default=""),
    FieldSpec("initial_guess", _as_str, default="PModel"),
    FieldSpec("relativity", _as_str, default="none"),
    FieldSpec("geom_opt", _as_str, default="OPT"),
    FieldSpec("freq_type", _as_str, default="FREQ"),
    FieldSpec("orca_parallel_strategy", _as_parallel_strategy, default="auto"),
    FieldSpec("IMAG_scope", _as_imag_scope, default="initial"),
    FieldSpec("IMAG_option", _as_imag_option, default=2),
    FieldSpec("OCCUPIER_method", _as_occupier_method, default="auto"),
    FieldSpec("OCCUPIER_tree", _as_occupier_tree, default="deep"),
    FieldSpec("OWN_progressive_from", _as_yes_no, default="no"),
    FieldSpec("OWN_TREE_PURE_WINDOW", _as_int, default=None, allow_none=True),
    FieldSpec("approximate_spin_projection_APMethod", _as_ap_method, default=2),
    FieldSpec("ESD_nroots", _as_int, default=15),
    FieldSpec("ESD_maxdim", _as_int, default=None, allow_none=True),
    FieldSpec("ESD_SOC", _as_yes_no, default="false"),
    FieldSpec("properties_of_interest", _as_str, default=""),
)


def validate_control_config(config: MutableMapping[str, Any]) -> dict[str, Any]:
    """Validate and coerce CONTROL configuration values."""
    errors: list[str] = []
    validated: dict[str, Any] = dict(config)

    # First pass: validate OCCUPIER_method to determine if tree validation is needed
    occupier_method_raw = config.get("OCCUPIER_method", None)
    occupier_method = "manually"  # default
    if occupier_method_raw is not None and occupier_method_raw != "":
        try:
            occupier_method = _as_occupier_method(occupier_method_raw)
        except Exception:  # noqa: BLE001
            occupier_method = "manually"

    for spec in CONTROL_FIELD_SPECS:
        # Skip OCCUPIER_tree validation if OCCUPIER_method is 'manually'
        if spec.name == "OCCUPIER_tree" and occupier_method == "manually":
            # Just set default without validation when manually mode
            validated[spec.name] = "deep"
            continue

        raw = config.get(spec.name, None)
        if raw is None or raw == "":
            if spec.required and spec.default is None:
                errors.append(f"Missing required key: {spec.name}")
                continue
            if spec.default is not None:
                validated[spec.name] = spec.default
                continue
            if spec.allow_none:
                validated[spec.name] = None
                continue
        try:
            validated[spec.name] = spec.coerce(raw)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"Invalid value for {spec.name}: {exc}")

    def _validate_sequence_list(seq_value: Any, label: str) -> None:
        if not isinstance(seq_value, list):
            errors.append(f"{label} must be a list of mappings")
            return
        for idx, item in enumerate(seq_value, start=1):
            if not isinstance(item, Mapping):
                errors.append(f"{label}[{idx}] must be a mapping")
                continue
            if "index" not in item or "m" not in item:
                errors.append(f"{label}[{idx}] must define 'index' and 'm'")
                continue
            try:
                int(item["index"])
                int(item["m"])
            except Exception:  # noqa: BLE001
                errors.append(f"{label}[{idx}] has non-integer 'index' or 'm'")

    # ensure electron sequences have expected structure if present
    for seq_key in ("even_seq", "odd_seq"):
        if seq_key in config:
            _validate_sequence_list(config[seq_key], seq_key)

    blocks = config.get("_occupier_sequence_blocks")
    if blocks is not None:
        if not isinstance(blocks, list):
            errors.append("_occupier_sequence_blocks must be a list")
        else:
            for block_idx, block in enumerate(blocks, start=1):
                if not isinstance(block, Mapping):
                    errors.append(f"sequence block #{block_idx} must be a mapping")
                    continue
                deltas = block.get("deltas")
                if not isinstance(deltas, list) or not all(isinstance(d, int) for d in deltas):
                    errors.append(f"sequence block #{block_idx} has invalid 'deltas'")
                if "even_seq" in block:
                    _validate_sequence_list(block["even_seq"], f"sequence block #{block_idx} even_seq")
                if "odd_seq" in block:
                    _validate_sequence_list(block["odd_seq"], f"sequence block #{block_idx} odd_seq")

    if errors:
        raise ValueError("; ".join(errors))

    return validated
