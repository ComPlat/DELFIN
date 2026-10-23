"""g-xTB has to be the program that answers, not the one that agrees to.

An ordinary xtb **accepts** ``--gxtb``.  It prints no warning, exits zero, and
returns the GFN2 energy.  Measured on this machine with a two-atom probe:

    ordinary xtb   --gxtb  -0.981983694723   --gfn 2  -0.981983694723
    xtb-gxtb       --gxtb  -1.163058822224   --gfn 2  -0.981983694723

So a ranking asked for in g-xTB and served by the xtb next to it is a GFN2
ranking wearing a g-xTB label, and nothing anywhere says so.  The only way to
tell is to ask the same binary both questions and see whether it gives two
answers.

The test is skipped where no xtb is installed, and the discriminating half is
skipped where no g-xTB is -- but where both are present it is the check that
the guard actually guards.
"""

from __future__ import annotations

import pytest

from delfin.manta import _gfnff_rank as ranking


def test_gxtb_is_a_method_the_ranker_knows():
    assert "gxtb" in ranking._METHOD_FLAGS
    assert ranking._METHOD_FLAGS["gxtb"] == ["--gxtb"]


def test_gxtb_needs_its_own_binary():
    """It is not reachable by handing a flag to the xtb beside it."""
    assert "gxtb" in ranking._NEEDS_GXTB


@pytest.mark.skipif(ranking._XTB is None, reason="no xtb installed")
def test_an_ordinary_xtb_is_refused_as_gxtb():
    """The whole point: this binary answers ``--gxtb``, and must not be believed."""
    ranking._GXTB_VERIFIED.clear()
    try:
        assert ranking._binary_is_really_gxtb(ranking._XTB) is False
    finally:
        ranking._GXTB_VERIFIED.clear()


@pytest.mark.skipif(ranking._GXTB is None, reason="no g-xTB build installed")
def test_a_real_gxtb_build_is_accepted():
    ranking._GXTB_VERIFIED.clear()
    try:
        assert ranking._binary_is_really_gxtb(ranking._GXTB) is True
    finally:
        ranking._GXTB_VERIFIED.clear()


@pytest.mark.skipif(ranking._XTB is None, reason="no xtb installed")
def test_the_two_hamiltonians_disagree_on_the_probe():
    """If they ever agreed, the probe could not tell the programs apart and the
    guard would be a coin toss.  Two atoms are enough for them to disagree."""
    with_gxtb = ranking._probe_energy(
        ranking._XTB, ranking._METHOD_FLAGS["gfnff"], 60.0)
    with_gfn2 = ranking._probe_energy(
        ranking._XTB, ranking._METHOD_FLAGS["gfn2"], 60.0)
    assert with_gxtb is not None and with_gfn2 is not None
    assert abs(with_gxtb - with_gfn2) > 1e-8


def test_a_refused_gxtb_yields_no_binary_rather_than_a_substitute():
    """``binary_for`` returns None instead of quietly handing back the ordinary
    xtb.  The caller then has to fall back and say so."""
    saved_binary, saved_cache = ranking._GXTB, dict(ranking._GXTB_VERIFIED)
    try:
        ranking._GXTB = None
        assert ranking.binary_for("gxtb") is None
        assert ranking.available("gxtb") is False
        # every other method is unaffected
        if saved_binary is not None or ranking._XTB is not None:
            assert ranking.binary_for("gfn2") == ranking._XTB
    finally:
        ranking._GXTB = saved_binary
        ranking._GXTB_VERIFIED.clear()
        ranking._GXTB_VERIFIED.update(saved_cache)


def test_the_control_key_accepts_gxtb_and_refuses_nonsense():
    from delfin.common.control_validator import _as_manta_rank

    assert _as_manta_rank("gxtb") == "gxtb"
    assert _as_manta_rank("") == "gfn2"
    with pytest.raises(ValueError):
        _as_manta_rank("dft")
