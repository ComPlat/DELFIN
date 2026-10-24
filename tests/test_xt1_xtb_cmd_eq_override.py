"""XT1 F6: extra_args cannot shadow validated flags via the ``=`` form.

xTB parses ``--chrg=-2`` the same as ``--chrg -2``, and F4's set check
(``a in protected``) only matches the exact bare token — so the
validated charge, unpaired electrons, GFN level and solvation model
were still overridable by one extra token. The builder must reject the
``=``-attached form of every protected flag as well.
"""

import pytest

from delfin.tools.adapters.xtb_native import build_xtb_cmd


def test_extra_args_cannot_override_chrg():
    with pytest.raises(ValueError):
        build_xtb_cmd("mol.xyz", charge=0, uhf=0, extra_args=("--chrg=-2",))
    with pytest.raises(ValueError):
        build_xtb_cmd("mol.xyz", charge=0, uhf=0, extra_args=("--chrg=2",))


def test_extra_args_cannot_override_uhf_gfn_solvent():
    for token in ("--uhf=2", "--gfn=1", "--gfnff=1", "--alpb=water",
                  "--gbsa=water"):
        with pytest.raises(ValueError):
            build_xtb_cmd("mol.xyz", charge=0, uhf=0, extra_args=(token,))


def test_harmless_flags_with_values_are_still_allowed():
    cmd = build_xtb_cmd("mol.xyz", charge=0, uhf=0, extra_args=("--etemp=300",))
    assert "--etemp=300" in cmd


def test_bare_protected_tokens_still_refused():
    with pytest.raises(ValueError):
        build_xtb_cmd("mol.xyz", charge=0, uhf=0,
                      extra_args=("--chrg", "-2"))
