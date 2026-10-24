"""DS1 adversarial: the ADAPTER path (the real caller's path).

The builder's physics gate lives in the wrapped execute() of the adapters,
not in the delfin/analysis_tools wrappers. These probes walk the adapter
surface (the path a workflow and a user step take) with hostile params:
negative and out-of-range indices, empty donor lists, huge radius.
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
    """Negative index through the ADAPTER must be refused before morfeus.

    The wrapper adds +1 (0-based -> 1-based, morfeus_wrapper.py:85), so -1
    reaches morfeus as 0 and dies with the misleading 'indices should not
    be 0 (1-indexed)' IndexError. The phase-1 physics gate checks presence
    and radius but not the index value.
    """
    pytest.importorskip("morfeus")
    from delfin.tools.adapters.morfeus import MorfeusBuriedVolumeAdapter

    with tempfile.TemporaryDirectory() as d:
        xyz = _small_xyz(Path(d))
        result = MorfeusBuriedVolumeAdapter().execute(
            Path(d), geometry=xyz, metal_index=-1
        )
        assert result.status.name in {"FAILED", "ERROR"}, result.status
        assert "metal_index" in (result.error or ""), result.error
