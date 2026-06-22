"""Central workflow registry for DELFIN.

Provides :func:`register`, :func:`get`, and :func:`list_workflows` so that
CLI, Python API, and dashboard share a single discovery mechanism.
"""

from __future__ import annotations

from typing import Dict, Optional

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# Internal registry: name -> Workflow instance
_REGISTRY: Dict[str, object] = {}
_DISCOVERED = False


def register(workflow: object) -> object:
    """Register a workflow instance (or use as a decorator on a class).

    The object must have a ``name`` attribute.
    """
    name = getattr(workflow, "name", None)
    if name is None:
        raise ValueError(f"Workflow {workflow!r} has no 'name' attribute")
    _REGISTRY[name] = workflow
    logger.debug("Registered workflow: %s", name)
    return workflow


def get(name: str) -> Optional[object]:
    """Look up a workflow by name.  Triggers lazy discovery on first call."""
    if not _DISCOVERED:
        _discover_builtin_workflows()
    return _REGISTRY.get(name)


def list_workflows() -> Dict[str, object]:
    """Return all registered workflows.  Triggers lazy discovery."""
    if not _DISCOVERED:
        _discover_builtin_workflows()
    return dict(_REGISTRY)


def _discover_builtin_workflows() -> None:
    """Lazy-import built-in workflow wrappers to populate the registry."""
    global _DISCOVERED
    _DISCOVERED = True

    _try_import("delfin.workflows.contrib.tadf_xtb_workflow")
    _try_import("delfin.workflows.contrib.hyperpol_workflow")
    _try_import("delfin.workflows.contrib.occupier_workflow")
    _try_import("delfin.workflows.contrib.esd_workflow")
    _try_import("delfin.workflows.contrib.imag_workflow")
    _try_import("delfin.workflows.contrib.co2_workflow")


def _try_import(module_path: str) -> None:
    try:
        __import__(module_path)
    except ImportError:
        logger.debug("Optional workflow module not found: %s", module_path)
