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
