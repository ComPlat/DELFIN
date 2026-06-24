from delfin.utils import calculate_total_electrons_txt


def test_calculate_total_electrons_txt_uses_existing_start_geometry_for_smiles(tmp_path):
    control = tmp_path / "CONTROL.txt"
    input_file = tmp_path / "input.txt"
    start_file = tmp_path / "start.txt"

    control.write_text("input_file=input.txt\ncharge=0\n", encoding="utf-8")
    input_file.write_text("C\n", encoding="utf-8")
    start_file.write_text("1\ncomment\nC 0.0 0.0 0.0\n", encoding="utf-8")

    total_electrons, multiplicity = calculate_total_electrons_txt(str(control))

    assert total_electrons == 6
    assert multiplicity == 1


def test_calculate_total_electrons_txt_counts_smiles_without_xyz_conversion(tmp_path):
    control = tmp_path / "CONTROL.txt"
    input_file = tmp_path / "input.txt"

    control.write_text("input_file=input.txt\ncharge=0\n", encoding="utf-8")
    input_file.write_text("[Ag+]\n", encoding="utf-8")

    total_electrons, multiplicity = calculate_total_electrons_txt(str(control))

    assert total_electrons == 47
    assert multiplicity == 2


def test_calculate_total_electrons_txt_returns_none_when_input_missing(tmp_path):
    """Reproduces the smart-recalc crash: a missing input geometry must yield
    None (a clean signal the caller can guard) — never a partial/garbage tuple.

    The run path in cli.py relies on this: None -> log + _finalize(1), instead
    of unpacking None into a TypeError or computing on an empty geometry.
    """
    control = tmp_path / "CONTROL.txt"
    # No input file is written: CONTROL points at a geometry that isn't there.
    control.write_text("input_file=input.xyz\ncharge=0\n", encoding="utf-8")

    result = calculate_total_electrons_txt(str(control))

    assert result is None


def test_calculate_total_electrons_txt_falls_back_to_input_txt(tmp_path):
    """Recalc safety net (MarcDy3): CONTROL names ``input.xyz`` (a derived file
    that was not staged), but the real coordinates survive as ``input.txt``.
    The electron count must use that sibling instead of aborting with None.
    """
    control = tmp_path / "CONTROL.txt"
    control.write_text("input_file=input.xyz\ncharge=0\n", encoding="utf-8")
    # input.xyz is absent; geometry survives as a bare coordinate block.
    (tmp_path / "input.txt").write_text(
        "Dy 0.0 0.0 0.0\nP 0.0 0.0 2.3\n", encoding="utf-8"
    )

    result = calculate_total_electrons_txt(str(control))

    assert result is not None
    total_electrons, multiplicity = result
    # Dy(66) + P(15) = 81 electrons (odd -> multiplicity 2).
    assert total_electrons == 81
    assert multiplicity == 2


def test_calculate_total_electrons_txt_falls_back_to_start_txt(tmp_path):
    """If neither ``input.xyz`` nor a coordinate-bearing ``input.txt`` exist,
    the resolver falls back to ``start.txt`` before giving up.
    """
    control = tmp_path / "CONTROL.txt"
    control.write_text("input_file=input.xyz\ncharge=0\n", encoding="utf-8")
    (tmp_path / "start.txt").write_text("C 0.0 0.0 0.0\n", encoding="utf-8")

    result = calculate_total_electrons_txt(str(control))

    assert result == (6, 1)
