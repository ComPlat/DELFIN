"""DS1 adversarial tests (reviewer s4).

Probes the descriptor wrappers the way a user with a typo would, not the
way the contract tests do.  Targets (verified on fcf63750):

  M1: multiwfn_wrapper.bond_order_analysis silently falls back to Mayer
      for an unknown method (method_map.get(..., "1"), wrapper:116) --
      a wrong quantity that looks fine.
  M2: morfeus_wrapper.buried_volume does not validate index/radius before
      morfeus raises an opaque IndexError deep inside numpy.
  C1: cclib_wrapper.parse_output uses bare cclib.io.ccread, although
      delfin/agent/energy_profile.py:264 documents that cclib 1.8.1
      raises IndexError on ORCA 6 outputs and carries its own fallback.
"""

from __future__ import annotations

import gzip  # noqa: F401  (fixture loader twin of test_energy_profile_layout.py)
import textwrap  # noqa: F401  (kept for upcoming C1 fixture shaping)
from pathlib import Path

import pytest

_FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------- M1: Multiwfn


def test_m1_unknown_bond_order_method_is_refused_not_silently_swapped():
    """An unknown method must raise, not silently compute Mayer instead.

    The adapter layer enforces the enum (delfin/tools/_validate.py:204), but
    the wrapper is also called directly (docstring usage examples) and by
    future steps; there the silent swap returns a plausible wrong number.
    """
    from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis
    import delfin.analysis_tools.multiwfn_wrapper as mw

    def _fake_run(input_file, commands, *, timeout=300):
        # record the menu commands the wrapper chose
        _fake_run.seen = commands
        return "1(C ) -- 2(C ) 1.0000"

    _fake_run.seen = None
    mw._run_multiwfn = _fake_run

    with pytest.raises(ValueError):
        bond_order_analysis("x.molden", method="nonsense")
    # and a harmless method must still work
    out = bond_order_analysis("x.molden", method="mayer")
    assert out["method"] == "mayer"
    assert _fake_run.seen.startswith("9\n")


def test_m1_none_method_is_refused():
    from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis
    import delfin.analysis_tools.multiwfn_wrapper as mw

    mw._run_multiwfn = lambda *a, **k: "1(C ) -- 2(C ) 1.0000"
    with pytest.raises(ValueError):
        bond_order_analysis("x.molden", method=None)


def test_m1_unicode_lookalike_method_is_refused():
    from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis
    import delfin.analysis_tools.multiwfn_wrapper as mw

    mw._run_multiwfn = lambda *a, **k: "1(C ) -- 2(C ) 1.0000"
    # Cyrillic p (U+043F) instead of Latin p in "mayer"
    with pytest.raises(ValueError):
        bond_order_analysis("x.molden", method="may\u043fer")


# ------------------------------------------------------------ M2: morfeus index


def _small_xyz(directory: Path):
    """3-atom dummy complex: M at 0, two ligand atoms."""
    xyz = directory / "c.xyz"
    xyz.write_text(
        "3\nnone\nPt 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n"
    )
    return xyz


def test_m2_metal_index_out_of_range_is_a_clear_error():
    from delfin.analysis_tools.morfeus_wrapper import buried_volume

    pytest.importorskip("morfeus")
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        xyz = _small_xyz(Path(d))
        with pytest.raises((ValueError, IndexError)) as exc_info:
            buried_volume(xyz, metal_index=99)
        # the message must name the problem, not be a bare numpy IndexError
        assert "99" in str(exc_info.value)


def test_m2_negative_metal_index_is_refused():
    from delfin.analysis_tools.morfeus_wrapper import buried_volume

    pytest.importorskip("morfeus")
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        xyz = _small_xyz(Path(d))
        with pytest.raises(ValueError):
            buried_volume(xyz, metal_index=-1)


def test_m2_donor_count_for_bite_angle_is_two():
    from delfin.analysis_tools.morfeus_wrapper import bite_angle

    pytest.importorskip("morfeus")
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        xyz = _small_xyz(Path(d))
        with pytest.raises(ValueError):
            bite_angle(xyz, 0, [1])
