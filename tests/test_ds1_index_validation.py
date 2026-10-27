"""DS1 R2 (reviewer nacht-s4): the negative metal index must be refused.

Red on a19f5a8a: a negative metal_index slips through validate_physics
(it is present) and the wrapper's +1 conversion turns -1 into index 0,
which morfeus rejects with a misleading "indices should not be 0
(1-indexed)" IndexError — the user's error is not named at the layer that
could name it. The adapter's execute-time physics gate must refuse a
negative index before the tool is invoked (the wrapper is stubbed, so
only the gate can fail).

Also: a positive index must still pass (the "does not fire" direction).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.tools._registry import get


@pytest.fixture(autouse=True)
def stub_wrappers(monkeypatch):
    import delfin.analysis_tools.morfeus_wrapper as mw

    def fake(*args, **kwargs):
        return {"ok": True}

    monkeypatch.setattr(mw, "buried_volume", fake)
    monkeypatch.setattr(mw, "cone_angle", fake)
    monkeypatch.setattr(mw, "sterimol", fake)
    monkeypatch.setattr(mw, "full_steric_report", fake)


def _run(name, **kwargs):
    return get(name).execute(Path("."), **kwargs)


@pytest.mark.parametrize("name", [
    "morfeus_buried_volume", "morfeus_cone_angle", "morfeus_full_report",
])
def test_negative_metal_index_is_refused(name):
    with pytest.raises(ValueError, match="metal_index"):
        _run(name, metal_index=-1)


def test_negative_sterimol_index_is_refused():
    with pytest.raises(ValueError, match="atom1_index"):
        _run("morfeus_sterimol", atom1_index=-2, atom2_index=3)


def test_index_zero_is_still_accepted():
    """0 is a valid 0-based index — the old refusal must not over-fire."""
    _run("morfeus_buried_volume", metal_index=0)


def test_out_of_range_indexerror_becomes_named_physics_error(monkeypatch):
    """If a bad index still reaches morfeus, its bare IndexError is
    re-raised as a ValueError naming the physics, not crashing."""
    import delfin.analysis_tools.morfeus_wrapper as mw

    def fake(*args, **kwargs):
        raise IndexError("indices should not be 0 (1-indexed)")

    monkeypatch.setattr(mw, "buried_volume", fake)
    with pytest.raises(ValueError, match="index rejected by morfeus"):
        _run("morfeus_buried_volume", metal_index=0)
