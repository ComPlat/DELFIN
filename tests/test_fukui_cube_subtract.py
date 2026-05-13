"""Unit tests for delfin.fukui cube subtraction."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delfin import fukui


def _write_cube(path: Path, data_values: list[float], *, n_atoms: int = 1) -> None:
    """Write a tiny 2×2×2 cube with the given voxel data."""
    header = textwrap.dedent(f"""\
        Test cube
        comment line
        {n_atoms:5d}    0.000000    0.000000    0.000000
            2    1.000000    0.000000    0.000000
            2    0.000000    1.000000    0.000000
            2    0.000000    0.000000    1.000000
            6    6.000000    0.500000    0.500000    0.500000
        """).rstrip() + "\n"
    # 8 voxels = 2 lines of 6 + 2
    data_str = " ".join(f"{v: .5E}" for v in data_values)
    path.write_text(header + data_str + "\n", encoding="utf-8")


def test_read_cube_round_trip(tmp_path):
    voxels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    p = tmp_path / "test.cube"
    _write_cube(p, voxels)

    cube = fukui.read_cube(p)
    assert cube.n_atoms == 1
    assert cube.origin == (0.0, 0.0, 0.0)
    assert cube.axes == [(2, 1.0, 0.0, 0.0), (2, 0.0, 1.0, 0.0), (2, 0.0, 0.0, 1.0)]
    assert cube.data == pytest.approx(voxels)


def test_subtract_cubes_basic(tmp_path):
    cube_a = tmp_path / "a.cube"
    cube_b = tmp_path / "b.cube"
    out = tmp_path / "diff.cube"

    _write_cube(cube_a, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    _write_cube(cube_b, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

    fukui.subtract_cubes(cube_a, cube_b, out)

    diff = fukui.read_cube(out)
    assert diff.data == pytest.approx([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])


def test_subtract_cubes_with_scale_half(tmp_path):
    """f_zero = (rho(anion) - rho(cation)) / 2 uses scale=0.5."""
    cube_a = tmp_path / "a.cube"
    cube_b = tmp_path / "b.cube"
    out = tmp_path / "diff.cube"

    _write_cube(cube_a, [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
    _write_cube(cube_b, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    fukui.subtract_cubes(cube_a, cube_b, out, scale=0.5)

    diff = fukui.read_cube(out)
    assert diff.data == pytest.approx([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])


def test_subtract_cubes_grid_mismatch_raises(tmp_path):
    cube_a = tmp_path / "a.cube"
    cube_b = tmp_path / "b.cube"
    out = tmp_path / "diff.cube"
    _write_cube(cube_a, [1.0] * 8)

    # Cube B with a different origin
    bad = textwrap.dedent("""\
            1    0.500000    0.000000    0.000000
                2    1.000000    0.000000    0.000000
                2    0.000000    1.000000    0.000000
                2    0.000000    0.000000    1.000000
                6    6.000000    0.500000    0.500000    0.500000
        """)
    cube_b.write_text("Test cube\ncomment\n" + bad + "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")

    with pytest.raises(ValueError, match="grid mismatch"):
        fukui.subtract_cubes(cube_a, cube_b, out)


def test_subtract_cubes_voxel_count_mismatch_raises(tmp_path):
    """Two cubes with same grid axes but different data lengths."""
    cube_a = tmp_path / "a.cube"
    cube_b = tmp_path / "b.cube"
    out = tmp_path / "diff.cube"
    _write_cube(cube_a, [1.0] * 8)
    _write_cube(cube_b, [0.5] * 4)

    with pytest.raises(ValueError, match="voxel count mismatch"):
        fukui.subtract_cubes(cube_a, cube_b, out)
