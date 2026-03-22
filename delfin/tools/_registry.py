"""Step adapter registry with lazy discovery."""

from __future__ import annotations

from typing import Dict, Optional

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter

logger = get_logger(__name__)

_REGISTRY: Dict[str, StepAdapter] = {}
_DISCOVERED = False


def register(adapter: StepAdapter) -> StepAdapter:
    """Register a step adapter instance."""
    if not adapter.name:
        raise ValueError(f"Adapter {adapter!r} has no 'name'")
    _REGISTRY[adapter.name] = adapter
    logger.debug("Registered step adapter: %s", adapter.name)
    return adapter


def get(name: str) -> Optional[StepAdapter]:
    """Look up a step adapter by name.  Triggers lazy discovery."""
    if not _DISCOVERED:
        _discover_builtin_adapters()
    return _REGISTRY.get(name)


def list_steps() -> Dict[str, StepAdapter]:
    """Return all registered step adapters."""
    if not _DISCOVERED:
        _discover_builtin_adapters()
    return dict(_REGISTRY)


def _discover_builtin_adapters() -> None:
    """Lazy-import built-in adapter modules to populate the registry."""
    global _DISCOVERED
    _DISCOVERED = True
    for mod in (
        "delfin.tools.adapters.smiles",
        "delfin.tools.adapters.morfeus",
        "delfin.tools.adapters.cclib",
        "delfin.tools.adapters.xtb",
        "delfin.tools.adapters.orca",
        "delfin.tools.adapters.crest",
        "delfin.tools.adapters.censo",
        "delfin.tools.adapters.packmol",
        "delfin.tools.adapters.multiwfn",
        "delfin.tools.adapters.imag",
        "delfin.tools.adapters.guppy",
        "delfin.tools.adapters.mlp",
        "delfin.tools.adapters.xtb_solvator",
        "delfin.tools.adapters.genarris",
        "delfin.tools.adapters.uv_vis",
        "delfin.tools.adapters.afp",
        "delfin.tools.adapters.reporting",
        "delfin.tools.adapters.python_func",
        "delfin.tools.adapters.ase_adapter",
        "delfin.tools.adapters.turbomole",
        "delfin.tools.adapters.openbabel",
    ):
        _try_import(mod)


def _try_import(module_path: str) -> None:
    try:
        __import__(module_path)
    except ImportError:
        logger.debug("Optional adapter module not available: %s", module_path)
