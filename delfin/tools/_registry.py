"""Step adapter registry with lazy discovery."""

from __future__ import annotations

from typing import Dict, Optional

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter

logger = get_logger(__name__)

_REGISTRY: Dict[str, StepAdapter] = {}
_DISCOVERED = False

# Modules kept for backwards compatibility that the directory walk cannot
# reach (currently none — every adapter lives as <stem>.py under adapters/).
# Kept as a hook so a renamed/legacy module still loads after discovery.
_LEGACY_ADAPTER_MODULES: tuple = ()


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
    """Lazy-import built-in adapter modules to populate the registry.

    Every module file under ``delfin/tools/adapters/`` is imported, so a new
    adapter that registers itself on import is reachable without editing
    this list (the wave-20 OR1 defect: three adapters registered on import
    but their module names were missing from a hardcoded tuple and nothing
    imported them).  The explicit tuple below is kept ONLY for modules whose
    names no longer match their file (legacy); prefer the directory walk.
    """
    global _DISCOVERED
    _DISCOVERED = True
    from pathlib import Path

    adapters_dir = Path(__file__).parent / "adapters"
    if adapters_dir.is_dir():
        for path in sorted(adapters_dir.glob("*.py")):
            if path.stem == "__init__":
                continue
            mod_name = f"{__name__.rsplit('.', 1)[0]}.adapters.{path.stem}"
            _reexec_module(mod_name, path)
    for mod in _LEGACY_ADAPTER_MODULES:
        _try_import(mod)
    _load_user_adapters()


def _reexec_module(mod_name: str, path) -> None:
    """(Re-)execute an adapter module so registration always happens.

    Import from sys.modules cache re-runs no module body, so a SECOND
    discovery in one process (e.g. after a registry reset) registered
    nothing and the registry stayed empty.  Re-executing the file
    guarantees every module-level ``register(...)`` call runs again.
    User adapters never enter sys.modules (spec_from_file_location),
    so only the builtin walk needs this.
    """
    import sys

    try:
        module = sys.modules.pop(mod_name, None)
        if module is None:
            __import__(mod_name)
            return
        # The module must be in sys.modules while its body runs: @dataclass
        # (and typing.get_type_hints) look the class's module up there.
        sys.modules[mod_name] = module
        module.__spec__.loader.exec_module(module)
    except ImportError:
        logger.debug("Optional adapter module not available: %s", mod_name)


def _user_adapter_dirs():
    import os
    from pathlib import Path

    dirs = []
    env = os.environ.get("DELFIN_ADAPTERS_DIR")
    if env:
        dirs.append(Path(env))
    dirs.append(Path.home() / ".delfin" / "adapters")
    return dirs


def load_user_adapter_file(path) -> None:
    """Execute a single user adapter file (it calls ``register(...)`` itself).

    Raises on any import/exec error so callers that want feedback (e.g. the
    agent integrating a freshly built module) can surface it.
    """
    import importlib.util
    from pathlib import Path

    path = Path(path)
    spec = importlib.util.spec_from_file_location(f"delfin_user_adapter_{path.stem}", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load adapter file {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


def _load_user_adapters() -> None:
    """Load user-defined building blocks from ~/.delfin/adapters/*.py.

    Each file is a normal adapter module that calls ``register(...)``; drop one
    there and it becomes a registered capability without editing the package.
    Failures per file are skipped so one bad file never blocks discovery.
    """
    for directory in _user_adapter_dirs():
        try:
            if not directory.is_dir():
                continue
            for f in sorted(directory.glob("*.py")):
                try:
                    load_user_adapter_file(f)
                except Exception:
                    logger.debug("Failed to load user adapter %s", f, exc_info=True)
        except Exception:
            continue


def _try_import(module_path: str) -> None:
    try:
        __import__(module_path)
    except ImportError:
        logger.debug("Optional adapter module not available: %s", module_path)
