"""Recalc geometry recovery in delfin.define.

Reproduces the MarcDy3 recalc failure: CONTROL.txt names ``input_file=input.xyz``
but only the bare coordinate block (``input.txt`` / ``start.txt``) was staged to
the run directory. The geometry resolver must rebuild the missing ``input.xyz``
from that block (instead of clobbering it with an empty file), so the downstream
electron-count / OCCUPIER steps can run instead of aborting on empty geometry.
"""

from delfin.define import convert_input_txt_to_xyz, convert_xyz_to_input_txt
from delfin.utils import calculate_total_electrons_txt

# Water-like 3-atom block: 8 + 1 + 1 = 10 electrons (even -> multiplicity 1).
COORDS = "O 0.000000 0.000000 0.000000\nH 0.000000 0.000000 0.957000\nH 0.926000 0.000000 -0.239000\n"


def test_recovery_rebuilds_missing_xyz_from_input_txt(tmp_path):
    control = tmp_path / "CONTROL.txt"
    input_txt = tmp_path / "input.txt"
    input_xyz = tmp_path / "input.xyz"

    control.write_text("input_file=input.xyz\ncharge=0\n", encoding="utf-8")
    input_txt.write_text(COORDS, encoding="utf-8")
    assert not input_xyz.exists()

    # Mirrors pipeline geometry-prep: convert the CONTROL-named .xyz -> input.txt.
    convert_xyz_to_input_txt(str(input_xyz), str(input_txt))

    # The missing .xyz is rebuilt with a proper 2-line XYZ header.
    assert input_xyz.exists()
    lines = input_xyz.read_text(encoding="utf-8").splitlines()
    assert lines[0].strip() == "3"
    assert sum(1 for ln in lines[2:] if ln.strip()) == 3

    # And electron counting (the step that aborted before) now succeeds.
    result = calculate_total_electrons_txt(str(control))
    assert result is not None
    total_electrons, multiplicity = result
    assert total_electrons == 10
    assert multiplicity == 1


def test_recovery_rebuilds_missing_xyz_from_start_txt(tmp_path):
    input_txt = tmp_path / "input.txt"
    start_txt = tmp_path / "start.txt"
    input_xyz = tmp_path / "input.xyz"

    # input.txt is empty; only start.txt carries the surviving geometry.
    input_txt.write_text("", encoding="utf-8")
    start_txt.write_text(COORDS, encoding="utf-8")

    convert_xyz_to_input_txt(str(input_xyz), str(input_txt))

    assert input_xyz.read_text(encoding="utf-8").splitlines()[0].strip() == "3"
    # input.txt is backfilled with the bare coordinate block too.
    assert sum(1 for ln in input_txt.read_text(encoding="utf-8").splitlines() if ln.strip()) == 3


def test_no_recovery_keeps_empty_fallback(tmp_path):
    """Neither input.txt nor start.txt has coordinates -> old behaviour:
    empty input.txt is created, no input.xyz, and no exception is raised."""
    input_txt = tmp_path / "input.txt"
    input_xyz = tmp_path / "input.xyz"
    input_txt.write_text("", encoding="utf-8")

    convert_xyz_to_input_txt(str(input_xyz), str(input_txt))

    assert input_txt.exists()
    assert not input_xyz.exists()


def test_convert_input_txt_to_xyz_is_inverse_of_convert(tmp_path):
    src_txt = tmp_path / "input.txt"
    dst_xyz = tmp_path / "input.xyz"
    src_txt.write_text(COORDS, encoding="utf-8")

    convert_input_txt_to_xyz(str(src_txt), str(dst_xyz))
    xyz_lines = dst_xyz.read_text(encoding="utf-8").splitlines()
    assert xyz_lines[0].strip() == "3"

    # Round-trip back to a coordinate block reproduces the original atoms.
    back_txt = tmp_path / "back.txt"
    convert_xyz_to_input_txt(str(dst_xyz), str(back_txt))
    original_atoms = [ln.split()[0] for ln in COORDS.splitlines() if ln.strip()]
    back_atoms = [ln.split()[0] for ln in back_txt.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert back_atoms == original_atoms


def test_convert_input_txt_to_xyz_accepts_full_xyz_input(tmp_path):
    """The builder also tolerates an already-headered XYZ as source (header
    lines are not coordinate rows and are dropped, not double-counted)."""
    src_xyz = tmp_path / "geom.xyz"
    dst_xyz = tmp_path / "out.xyz"
    src_xyz.write_text(f"3\nsome comment\n{COORDS}", encoding="utf-8")

    convert_input_txt_to_xyz(str(src_xyz), str(dst_xyz))
    assert dst_xyz.read_text(encoding="utf-8").splitlines()[0].strip() == "3"
