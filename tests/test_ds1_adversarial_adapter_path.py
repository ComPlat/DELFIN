"""DS1 adversarial: the ADAPTER path through the real registry entry point.

History: first written against the adapter CLASS imported directly, that
process never ran the registry's lazy discovery (delfin/tools/_registry.py:44
walks adapters/*.py and imports morfeus_realism.py, whose install() attaches
the physics gate), so the gate was absent and the test measured an artifact
of my own entry point. Corrected to the real platform path (registry get()),
which is what a workflow step and platform.resolve use.

The red control for the underlying defect is tests/test_ds1_index_validation.py
(builder, red on a19f5a8a). This file pins the fixed behavior through the real
caller path: a negative index must raise a ValueError naming the physics
BEFORE morfeus is invoked, with no compute spent.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


def _small_xyz(directory: Path):
    xyz = directory / "c.xyz"
    xyz.write_text("3\nnone\nPt 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n")
    return xyz


def test_adapter_buried_volume_negative_metal_index_is_a_clear_error():
    """Negative index through the registry path must be refused pre-morfeus.

    The wrapper adds +1 (0-based -> 1-based, morfeus_wrapper.py:85), so -1
    would reach morfeus as 0 and die with the misleading 'indices should
    not be 0 (1-indexed)' IndexError. The physics gate must refuse it first.
    """
    pytest.importorskip("morfeus")
    from delfin.tools._registry import get

    with tempfile.TemporaryDirectory() as d:
        xyz = _small_xyz(Path(d))
        adapter = get("morfeus_buried_volume")
        assert adapter is not None
        with pytest.raises(ValueError, match="metal_index"):
            adapter.execute(Path(d), geometry=xyz, metal_index=-1)


def test_adapter_indexerror_from_morfeus_becomes_named_physics_error():
    """An out-of-range index that slips past the value gate must not crash
    the step with a bare IndexError: the wrapped execute re-raises it as a
    ValueError naming the physics (morfeus_realism.py:162)."""
    pytest.importorskip("morfeus")
    import delfin.analysis_tools.morfeus_wrapper as mw
    from delfin.tools._registry import get

    def _boom(*args, **kwargs):
        raise IndexError("Atom indices should not be 0 (1-indexed).")

    original = mw.buried_volume
    mw.buried_volume = _boom
    try:
        with tempfile.TemporaryDirectory() as d:
            xyz = _small_xyz(Path(d))
            adapter = get("morfeus_buried_volume")
            with pytest.raises(ValueError, match="index rejected by morfeus"):
                adapter.execute(Path(d), geometry=xyz, metal_index=0)
    finally:
        mw.buried_volume = original
