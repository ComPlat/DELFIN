"""XT1 F5: the ORCA-xTB input builder honours its declared enum.

``_build_orca_xtb_input`` feeds its ``method`` straight into the ORCA
keyword line while the ParamSpec at delfin/tools/adapters/xtb.py:18
promises XTB1/XTB2 only. "XTB3" (or any other string) built an input
ORCA would refuse or misread. The raw builder now enforces the enum.
"""

import pytest

from delfin.tools.adapters.xtb import _build_orca_xtb_input


def test_xtb2_still_builds():
    inp = _build_orca_xtb_input("XTB2", "OPT", 0, 1, 1, 1000, "C 0 0 0")
    assert "XTB2" in inp


def test_xtb1_still_builds():
    inp = _build_orca_xtb_input("XTB1", "OPT", 0, 1, 1, 1000, "C 0 0 0")
    assert "XTB1" in inp


def test_xtb3_is_rejected():
    with pytest.raises(ValueError, match="XTB3"):
        _build_orca_xtb_input("XTB3", "OPT", 0, 1, 1, 1000, "C 0 0 0")


def test_none_method_is_rejected():
    with pytest.raises(ValueError, match="method"):
        _build_orca_xtb_input(None, "OPT", 0, 1, 1, 1000, "C 0 0 0")


def test_lowercase_xtb2_is_not_accepted_by_the_orca_spelling():
    with pytest.raises(ValueError, match="XTB2|method"):
        _build_orca_xtb_input("xtb2", "OPT", 0, 1, 1, 1000, "C 0 0 0")


def test_goat_keyword_is_untouched():
    inp = _build_orca_xtb_input("XTB2", "GOAT", 0, 1, 1, 1000, "C 0 0 0")
    assert "GOAT" in inp
