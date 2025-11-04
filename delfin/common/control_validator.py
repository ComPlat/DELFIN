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
    FieldSpec("orca_parallel_strategy", _as_parallel_strategy, default="auto"),
    FieldSpec("IMAG_scope", _as_imag_scope, default="initial"),
    FieldSpec("IMAG_option", _as_imag_option, default=2),
)


def validate_control_config(config: MutableMapping[str, Any]) -> dict[str, Any]:
    """Validate and coerce CONTROL configuration values."""
    errors: list[str] = []
    validated: dict[str, Any] = dict(config)

    for spec in CONTROL_FIELD_SPECS:
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

    # ensure electron sequences have expected structure if present
    for seq_key in ("even_seq", "odd_seq"):
        if seq_key in config:
            seq_value = config[seq_key]
            if not isinstance(seq_value, list):
                errors.append(f"{seq_key} must be a list of mappings")
                continue
            for idx, item in enumerate(seq_value, start=1):
                if not isinstance(item, Mapping):
                    errors.append(f"{seq_key}[{idx}] must be a mapping")
                    continue
                if "index" not in item or "m" not in item:
                    errors.append(f"{seq_key}[{idx}] must define 'index' and 'm'")
                    continue
                try:
                    int(item["index"])
                    int(item["m"])
                except Exception:  # noqa: BLE001
                    errors.append(f"{seq_key}[{idx}] has non-integer 'index' or 'm'")

    if errors:
        raise ValueError("; ".join(errors))

    return validated
