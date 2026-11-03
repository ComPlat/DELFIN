"""The rotation-movie glue: ORCA .out -> multi-frame xyz.

A dihedral scan, an optimisation or an IRC each print one
``CARTESIAN COORDINATES (ANGSTROEM)`` block per geometry. The glue
collects every block into one multi-frame xyz so the agent's native
``show_molecule`` renders the rotation movie. No new renderer.

No OPI dependency; runs under the py3.11 CI interpreter.
"""
import pytest

from delfin.tools.geometry_edit import (
    GeometryError,
    parse_xyz,
    trajectory_xyz,
)

SYN = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"

# two-frame fake output: anti (phi=180) then gauche (phi=60). Hand-built
# arithmetic, no CSD data.
_AT = ("C   0.000000000    0.000000000    0.000000000\n"
       "C   0.000000000    0.000000000    1.526000000\n"
       "C   1.415000000    0.000000000    2.034000000\n"
       "C   2.830000000    0.000000000    1.012000000\n")
_GA = ("C   0.000000000    0.000000000    0.000000000\n"
       "C   0.000000000    0.000000000    1.526000000\n"
       "C   1.415000000    0.700000000    2.034000000\n"
       "C   2.830000000    1.400000000    1.012000000\n")


def _out(n_blocks: int = 2) -> str:
    parts = [SYN, "! r2SCAN-3c TightOpt TightSCF\n"]
    for i, body in enumerate([_AT, _GA][:n_blocks]):
        parts.append(f"CARTESIAN COORDINATES (ANGSTROEM)\n"
                     f"  0  1\n  0  1\n{body}")
    return "".join(parts)


def test_collects_every_cartesian_block_as_frames(tmp_path):
    out = tmp_path / "scan.out"
    out.write_text(_out())
    xyz = trajectory_xyz(str(out))
    # multi-frame: the file holds n_blocks count+comment headers
    assert xyz.count("frame") == 2
    g = parse_xyz(xyz.split("\n\n")[0])
    assert g.symbols == ["C", "C", "C", "C"]


def test_returns_xyz_string_that_show_molecule_can_open(tmp_path):
    out = tmp_path / "scan.out"
    out.write_text(_out())
    xyz = trajectory_xyz(str(out))
    # multi-frame xyz: alternating count line and comment lines
    lines = xyz.splitlines()
    assert lines[0] == "4"
    assert lines[1].startswith("frame")
    assert lines[5] == "4"
    # second frame holds the gauche coords
    assert lines[9].split()[2].startswith("1.400")


def test_refuses_output_without_geometry_blocks(tmp_path):
    out = tmp_path / "sp.out"
    out.write_text(SYN + "FINAL SINGLE POINT ENERGY  -157.41\n")
    with pytest.raises(GeometryError):
        trajectory_xyz(str(out))
