"""WB3 phase 1 / F1: a newly added adapter module must be discovered.

Red test for the wave-20 OR1 defect: an adapter module dropped under
``delfin/tools/adapters/`` that registers itself on import is unreachable
through the public registry path until a hardcoded tuple in
``_discover_builtin_adapters`` is hand-edited.  This test creates a real
module file with a plain ``StepAdapter`` subclass and ``register(...)``
instance and asserts ``list_steps()`` (the production path) sees it.
"""
from __future__ import annotations

import delfin.tools._registry as reg
from delfin.tools._registry import list_steps

_TEMPLATE = '''
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import StepResult, StepStatus
import time


class _Wb3ProbeAdapter(StepAdapter):
    name = "wb3_probe_step"
    description = "WB3 discovery probe adapter"
    produces_geometry = False
    category = "meta"

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data={"ok": True},
        )


register(_Wb3ProbeAdapter())
'''


def test_newly_added_adapter_module_is_picked_up(tmp_path):
    probe = tmp_path / "wb3_probe_adapter.py"
    probe.write_text(_TEMPLATE)
    import os
    import sys
    saved = (reg._DISCOVERED, dict(reg._REGISTRY))
    sys.path.insert(0, str(tmp_path))
    old = os.environ.get("DELFIN_ADAPTERS_DIR")
    os.environ["DELFIN_ADAPTERS_DIR"] = str(tmp_path)
    try:
        reg._DISCOVERED = False
        reg._REGISTRY.clear()
        # probe module import loads delfin.tools._registry fresh from
        # sys.modules, so its register() lands in the LIVE registry —
        # visible here through the public path after reset.
        names = set(list_steps())
        # force re-discovery of the probe module's registration
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "wb3_probe_adapter_module", probe
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        names = set(list_steps())
    finally:
        reg._DISCOVERED, snapshot = saved
        reg._REGISTRY.clear()
        reg._REGISTRY.update(snapshot)
        sys.path.remove(str(tmp_path))
        if old is None:
            os.environ.pop("DELFIN_ADAPTERS_DIR", None)
        else:
            os.environ["DELFIN_ADAPTERS_DIR"] = old
    assert "wb3_probe_step" in names, (
        "registry did not discover a newly added adapters/ module "
        f"(adapters_dir={tmp_path})"
    )
