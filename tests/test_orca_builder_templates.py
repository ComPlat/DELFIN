"""Unit tests for the ORCA Builder template helpers (coordinate stripping)."""

from delfin.dashboard.tab_orca_builder import strip_coord_block


def test_strip_coord_block_inline_xyz():
    text = (
        "! PBE0 OPT def2-SVP\n\n"
        "manual line\n\n"
        "%pal\n  nprocs 12\nend\n\n"
        "%maxcore 6000\n\n"
        "* xyz 0 1\nC 0.0 0.0 0.0\nO 0.0 0.0 1.2\n*\n"
    )
    body = strip_coord_block(text)
    assert body == (
        "! PBE0 OPT def2-SVP\n\n"
        "manual line\n\n"
        "%pal\n  nprocs 12\nend\n\n"
        "%maxcore 6000"
    )
    assert "* xyz" not in body
    assert "C 0.0 0.0 0.0" not in body


def test_strip_coord_block_xyzfile():
    assert strip_coord_block("! PBE0\n\n%maxcore 6000\n\n* xyzfile 0 1 mol.xyz\n") == (
        "! PBE0\n\n%maxcore 6000"
    )


def test_strip_coord_block_preserves_manual_free_lines():
    text = "! PBE0\n\nMANUAL_A\nMANUAL_B\n\n* xyz 0 1\nH 0 0 0\n*"
    body = strip_coord_block(text)
    assert "MANUAL_A" in body and "MANUAL_B" in body
    assert "H 0 0 0" not in body


def test_strip_coord_block_no_coords_and_empty():
    assert strip_coord_block("! PBE0\n%maxcore 6000") == "! PBE0\n%maxcore 6000"
    assert strip_coord_block("") == ""
    assert strip_coord_block(None) == ""
