"""Adversarial tests for the ORCA-driven xTB adapters (XT1 review).

Reviewer nacht-s17, package XT1. Attacks the ORCA input builders in
``delfin/tools/adapters/xtb.py`` (xtb_opt, xtb_goat) and
``xtb_solvator.py`` — the text DELFIN writes for ORCA's built-in xTB.
"""

import pytest

from delfin.tools.adapters.xtb import _build_orca_xtb_input
from delfin.tools.adapters.xtb import _read_xyz_coords


class TestOrcaXtbInput:
    def test_input_shape_opt(self):
        inp = _build_orca_xtb_input("XTB2", "OPT", charge=1, mult=2, cores=4,
                                    maxcore=1000, coords="O 0 0 0\nH 1 0 0\nH 0 1 0")
        assert inp.startswith("!XTB2 OPT\n")
        assert "*xyz 1 2\n" in inp
        assert "%pal nprocs 4 end\n" in inp
        assert inp.rstrip().endswith("*")

    def test_only_xtb1_or_xtb2_accepted(self):
        # The ParamSpec enum (xtb.py:18) promises XTB1/XTB2 only; the raw
        # builder must not silently accept e.g. "XTB3" (no such ORCA keyword).
        with pytest.raises((ValueError, AssertionError)):
            _build_orca_xtb_input("XTB3", "OPT", 0, 1, 1, 1000, "H 0 0 0")


class TestReadXyzCoords:
    def test_standard_xyz_header_stripped(self, tmp_path):
        xyz = tmp_path / "water.xyz"
        xyz.write_text("3\nwater\nO 0 0 0\nH 1 0 0\nH 0 1 0\n")
        coords = _read_xyz_coords(xyz)
        assert "O 0 0 0" in coords
        assert "water" not in coords

    def test_empty_file(self, tmp_path):
        # Empty geometry: current behavior passes the empty string through
        # (ORCA would reject the input). Pinned as characterization; a hard
        # validation here would be better but is not a correctness defect.
        xyz = tmp_path / "empty.xyz"
        xyz.write_text("")
        assert _read_xyz_coords(xyz).strip() == ""

    def test_headerless_coords_pass_through(self, tmp_path):
        xyz = tmp_path / "raw.xyz"
        xyz.write_text("O 0 0 0\nH 1 0 0\n")
        assert "O 0 0 0" in _read_xyz_coords(xyz)
