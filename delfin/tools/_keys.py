"""Central key vocabulary — define well-known parameters once, reuse everywhere.

A *key* is a named, well-known parameter (``functional``, ``basis``, ``solvent``,
``dispersion``, ``relativity`` …) with its type, allowed values (enum), default
and description.  Defining them in one place means adapters, applications and the
exported manifest all agree on what ``functional`` means and which values are
valid — instead of each adapter re-declaring them ad hoc.

The allowed-value lists are **referenced from the single source of truth**
(:mod:`delfin.common.control_validator`: ``ORCA_FUNCTIONALS``, ``ORCA_SOLVENTS``,
``ORCA_BASIS_SETS`` …), never duplicated, so the tool framework and CONTROL.txt
validation can never drift apart.

    from delfin.tools._keys import key

    params = (
        key("charge", required=True),
        key("functional"),            # enum = the canonical ORCA functional list
        key("basis"),
        key("solvent"),
    )
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

from delfin.tools._spec import ParamSpec


@dataclass(frozen=True)
class KeySpec:
    """The vocabulary entry for one well-known key (usage-independent)."""

    name: str
    type: str = "str"
    description: str = ""
    default: Any = None
    enum: Optional[Tuple[Any, ...]] = None
    unit: str = ""
    enum_source: str = ""           # provenance of the allowed-value list


_UNSET = object()
_REGISTRY: Dict[str, KeySpec] = {}
_BUILT = False


def _load_allowed_values() -> Dict[str, Tuple[Any, ...]]:
    """Pull the canonical allowed-value lists from control_validator (best-effort)."""
    out: Dict[str, Tuple[Any, ...]] = {}
    try:
        from delfin.common import control_validator as cv
        out["functional"] = tuple(cv.ORCA_FUNCTIONALS)
        out["basis"] = tuple(cv.ORCA_BASIS_SETS)
        out["solvent"] = tuple(cv.ORCA_SOLVENTS)
        out["relativity"] = tuple(sorted(cv._RELATIVITY_VALUES) + ["NONE"])
        out["dispersion"] = tuple(sorted(v for v in cv.DISP_CORR_VALUES if v))
    except Exception:
        pass
    return out


def _build() -> None:
    global _BUILT
    _BUILT = True
    allowed = _load_allowed_values()

    def add(spec: KeySpec) -> None:
        _REGISTRY[spec.name] = spec

    # --- system / electronic keys ---
    add(KeySpec("charge", "int", "Molecular charge"))
    add(KeySpec("mult", "int", "Spin multiplicity (2S+1)", default=1))

    # --- level of theory ---
    add(KeySpec("functional", "str", "DFT functional", default="PBE0",
                enum=allowed.get("functional"),
                enum_source="control_validator.ORCA_FUNCTIONALS"))
    # `method` is the adapter-side name for the functional (ORCA's bang line).
    add(KeySpec("method", "str", "DFT functional / method", default="B3LYP",
                enum=allowed.get("functional"),
                enum_source="control_validator.ORCA_FUNCTIONALS"))
    add(KeySpec("basis", "str", "Main basis set", default="def2-SVP",
                enum=allowed.get("basis"),
                enum_source="control_validator.ORCA_BASIS_SETS"))
    add(KeySpec("metal_basisset", "str", "Basis set for metal centres",
                default="def2-TZVP", enum=allowed.get("basis"),
                enum_source="control_validator.ORCA_BASIS_SETS"))
    add(KeySpec("dispersion", "str", "Dispersion correction",
                enum=allowed.get("dispersion"),
                enum_source="control_validator.DISP_CORR_VALUES"))
    add(KeySpec("relativity", "str", "Relativistic treatment",
                enum=allowed.get("relativity"),
                enum_source="control_validator._RELATIVITY_VALUES"))

    # --- solvation ---
    add(KeySpec("solvent", "str", "Implicit solvent",
                enum=allowed.get("solvent"),
                enum_source="control_validator.ORCA_SOLVENTS"))
    add(KeySpec("solvent_model", "str", "Implicit solvation model",
                default="CPCM", enum=("CPCM", "SMD")))


def _ensure() -> None:
    if not _BUILT:
        _build()


def list_keys() -> Dict[str, KeySpec]:
    """All registered well-known keys."""
    _ensure()
    return dict(_REGISTRY)


def get_key(name: str) -> Optional[KeySpec]:
    """The :class:`KeySpec` for *name*, or ``None`` if it is not a known key."""
    _ensure()
    return _REGISTRY.get(name)


def register_key(spec: KeySpec) -> KeySpec:
    """Add or override a well-known key at runtime."""
    _ensure()
    _REGISTRY[spec.name] = spec
    return spec


def key(
    name: str,
    *,
    required: bool = False,
    default: Any = _UNSET,
    description: Optional[str] = None,
) -> ParamSpec:
    """Return a :class:`ParamSpec` for a well-known key, with optional overrides.

    Adapters and applications use this so a parameter named ``functional`` always
    carries the same type, allowed values and description.  Unknown names yield a
    plain string ParamSpec so callers are never blocked.
    """
    ks = get_key(name)
    if ks is None:
        return ParamSpec(
            name=name, type="str", required=required,
            default=(None if default is _UNSET else default),
            description=description or "",
        )
    return ParamSpec(
        name=ks.name,
        type=ks.type,
        required=required,
        default=(ks.default if default is _UNSET else default),
        enum=ks.enum,
        description=description or ks.description,
        unit=ks.unit,
    )


__all__ = ["KeySpec", "list_keys", "get_key", "register_key", "key"]
