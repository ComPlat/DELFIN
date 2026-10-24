"""WB3 phase 1 / F1: registry discovery must not rely on a hardcoded list.

Red test for the wave-20 OR1 defect class: adapters that register
themselves on import are unreachable at runtime when their module name
is missing from the hardcoded tuple in
``delfin.tools._registry._discover_builtin_adapters``.  This test walks
the real adapters directory and checks, per module, that the adapter
instances it registers on import appear in ``list_steps()``.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import delfin.tools._registry as reg
from delfin.tools._registry import list_steps
from delfin.tools._base import StepAdapter


def test_discovery_is_rerunnable_after_the_registry_is_cleared():
    """A consumer that resets (_DISCOVERED, _REGISTRY) gets a full registry.

    Red control for the operator's batch-order leak: discovery imported
    adapter modules through the sys.modules cache, so a SECOND discovery
    in one process re-executed no module and the registry stayed empty
    ("Unknown step 'esd'. Available:" blank).  Discovery must bypass the
    module cache and re-execute every adapters/*.py.
    """
    saved = (reg._DISCOVERED, dict(reg._REGISTRY))
    try:
        reg._DISCOVERED = False
        reg._REGISTRY.clear()
        # trigger discovery via the public path (cache is already warm
        # from earlier imports in this process — that is the defect case)
        names = set(list_steps())
        assert names, (
            "re-discovery after a registry reset must repopulate the "
            "registry; sys.modules cache made it a no-op"
        )
    finally:
        reg._DISCOVERED, snapshot = saved
        reg._REGISTRY.clear()
        reg._REGISTRY.update(snapshot)

    """Every adapters/ module that registers an instance must be reachable.

    For each module file in delfin/tools/adapters/ (excluding __init__),
    import it (exactly as discovery does), collect its ``StepAdapter``
    instances (i.e. module-level ``register(...)`` calls), and require
    every such name to appear in ``list_steps()``.  The registry is
    snapshotted and restored so this test leaves no state behind in a
    full-suite run.
    """
    saved = (reg._DISCOVERED, dict(reg._REGISTRY))
    try:
        reg._DISCOVERED = False
        reg._REGISTRY.clear()
        names = set(list_steps())

        adapters_dir = Path(reg.__file__).parent / "adapters"
        missing = []
        for path in sorted(adapters_dir.glob("*.py")):
            if path.stem == "__init__":
                continue
            mod_name = f"delfin.tools.adapters.{path.stem}"
            try:
                module = importlib.import_module(mod_name)
            except ImportError:
                # optional-dependency module that cannot even import here
                continue
            instance_names = [
                getattr(obj, "name")
                for obj in vars(module).values()
                if isinstance(obj, StepAdapter) and getattr(obj, "name", "")
            ]
            for n in instance_names:
                if n not in names:
                    missing.append(f"{mod_name}: {n}")
        assert missing == [], (
            "modules under delfin/tools/adapters/ register adapters that the "
            f"registry discovery never reaches: {missing}"
        )
    finally:
        reg._DISCOVERED, snapshot = saved
        reg._REGISTRY.clear()
        reg._REGISTRY.update(snapshot)
