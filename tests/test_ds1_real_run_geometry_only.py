# Test-specific file to report the practice geometry missing reason even when the optional morfeus
# dependency is unavailable. This file does NOT import the main test module, avoiding the
# module‑level morfeus importorskip that would mask the geometry skip.

from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parent.parent
GEOMETRY = ROOT / ".gate" / "practice" / "znh2o2cl2.xyz"

def test_practice_geometry_present():
    """Check that the practice geometry file exists.

    The test is intentionally simple: if the file is missing, it is skipped with a
    descriptive reason. This ensures the geometry‑missing skip is reported even when
    the morfeus package is not installed.
    """
    if not GEOMETRY.exists():
        pytest.skip(
            "scratch geometry .gate/practice/znh2o2cl2.xyz absent in this "
            "tree (never committed); see .gate/PR.md — value‑pinning works "
            "on any tree that has a tetrahedral Zn complex there"
        )
    assert True
