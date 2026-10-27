"""XT1: build_xtb_cmd must reject unbuildable xtb input, not guess it.

Reviewer findings F1-F4 (nacht-s17): an unknown GFN level silently fell
back to --gfn 2, uhf was built from an unvalidated multiplicity,
solvent_model passed straight onto the command line, and extra_args could
both arrive as a string and override the validated charge/spin flags.
A wrong level or spin produces plausible-but-wrong numbers, so the
command builder raises instead of guessing.
"""

import pytest

from delfin.tools.adapters.xtb_native import build_xtb_cmd


def test_known_method_still_builds():
    cmd = build_xtb_cmd("in.xyz", method="gfn2", charge=0, uhf=0)
    assert cmd[:4] == ["xtb", "in.xyz", "--gfn", "2"]


def test_unknown_method_is_rejected_not_silently_gfn2():
    with pytest.raises(ValueError, match="gfn3"):
        build_xtb_cmd("in.xyz", method="gfn3", charge=0, uhf=0)


def test_method_none_is_rejected():
    with pytest.raises(ValueError, match="method"):
        build_xtb_cmd("in.xyz", method=None, charge=0, uhf=0)


def test_gfnff_is_a_known_method():
    cmd = build_xtb_cmd("in.xyz", method="gfnff", charge=0, uhf=0)
    assert "--gfnff" in cmd


def test_uhf_negative_is_rejected():
    with pytest.raises(ValueError, match="uhf"):
        build_xtb_cmd("in.xyz", charge=0, uhf=-1)


def test_uhf_absurdly_large_is_rejected():
    with pytest.raises(ValueError, match="uhf"):
        build_xtb_cmd("in.xyz", charge=0, uhf=999)


def test_charge_non_integer_is_rejected():
    with pytest.raises((ValueError, TypeError)):
        build_xtb_cmd("in.xyz", charge="a", uhf=0)


def test_unknown_solvent_model_is_rejected():
    with pytest.raises(ValueError, match="solvent_model"):
        build_xtb_cmd("in.xyz", charge=0, uhf=0,
                      solvent="water", solvent_model="cpcm")


def test_alpb_solvent_still_builds():
    cmd = build_xtb_cmd("in.xyz", charge=0, uhf=0,
                        solvent="water", solvent_model="alpb")
    assert "--alpb" in cmd and "water" in cmd


def test_extra_args_string_is_rejected_not_silently_split():
    with pytest.raises((ValueError, TypeError)):
        build_xtb_cmd("in.xyz", charge=0, uhf=0, extra_args="--gfn 1")


def test_extra_args_cannot_override_chrg():
    with pytest.raises(ValueError, match="chrg"):
        build_xtb_cmd("in.xyz", charge=1, uhf=0,
                      extra_args=["--chrg", "-2"])


def test_extra_args_cannot_override_uhf():
    with pytest.raises(ValueError, match="uhf"):
        build_xtb_cmd("in.xyz", charge=0, uhf=0,
                      extra_args=["--uhf", "3"])


def test_extra_args_cannot_override_method():
    with pytest.raises(ValueError, match="gfn"):
        build_xtb_cmd("in.xyz", method="gfn1", charge=0, uhf=0,
                      extra_args=["--gfn", "2"])


def test_benign_extra_args_pass_through():
    cmd = build_xtb_cmd("in.xyz", charge=0, uhf=0,
                        extra_args=["--nofinalnames"])
    assert "--nofinalnames" in cmd
