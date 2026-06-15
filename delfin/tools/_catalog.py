"""Discovery / introspection over the registered building blocks.

This is the read-side counterpart to the serialization layer: it lets a UI,
an agent, or a human *browse* the available Bausteine and ask "what can I plug
in after this one?" without reading adapter source code.  Everything here is
derived from the optional :class:`~delfin.tools._spec.StepContract` each adapter
exposes via :meth:`StepAdapter.contract`, so undeclared adapters simply show up
with empty ports.

    from delfin.tools import describe, catalog, compatible_successors

    describe("orca_sp")                 # full contract for one block
    catalog(by="category")              # {"dft": [...], "structure": [...], ...}
    compatible_successors("orca_freq")  # blocks whose inputs orca_freq satisfies
"""

from __future__ import annotations

from typing import Dict, List, Optional

from delfin.tools._registry import get, list_steps
from delfin.tools._spec import StepContract


def describe(step_name: str) -> Optional[StepContract]:
    """Return the :class:`StepContract` for *step_name*, or ``None`` if unknown."""
    adapter = get(step_name)
    return adapter.contract() if adapter is not None else None


def catalog(*, by: str = "category") -> Dict[str, List[StepContract]]:
    """Return all registered building blocks, grouped.

    Parameters
    ----------
    by:
        ``"category"`` — group by the adapter's declared ``category``
        (``"uncategorized"`` for blocks that declare none).
        ``"produces"`` — index by each produced capability tag.
        ``"consumes"`` — index by each consumed capability tag (the
        "what can follow a block that yields X" view).
    """
    out: Dict[str, List[StepContract]] = {}
    for name, adapter in sorted(list_steps().items()):
        c = adapter.contract()
        if by == "produces":
            keys = sorted(c.produces) or ["(none)"]
        elif by == "consumes":
            keys = sorted(c.consumes) or ["(none)"]
        else:  # "category"
            keys = [c.category or "uncategorized"]
        for k in keys:
            out.setdefault(k, []).append(c)
    return out


def compatible_successors(step_name: str) -> List[str]:
    """Names of blocks whose required inputs *step_name* can satisfy.

    A successor is compatible when everything it ``consumes`` is contained in
    what *step_name* ``produces``.  This is the core "plug something in after
    this" query a builder needs.  Blocks that consume nothing are always
    compatible (they need no upstream port).
    """
    src = describe(step_name)
    if src is None:
        return []
    available = set(src.produces)
    out: List[str] = []
    for name, adapter in list_steps().items():
        if adapter.contract().consumes <= available:
            out.append(name)
    return sorted(out)


__all__ = ["describe", "catalog", "compatible_successors"]
