from pathlib import Path

from delfin.scheduler_profiles import (
    get_molecule_features,
    infer_stage_key,
    parse_molecule_text,
    predict_stage_profile,
)


def test_infer_stage_key_handles_esd_jobs():
    assert infer_stage_key("esd_S0", "ESD S0 optimization") == "esd_S0"
    assert infer_stage_key("esd_T1", "ESD T1 optimization") == "esd_T1"
    assert infer_stage_key("job", "Fluorescence S1→S0") == "fluorescence_S1_S0"


def test_parse_molecule_text_counts_xyz_atoms():
    parsed = parse_molecule_text("3\ncomment\nC 0 0 0\nH 0 0 1\nH 0 1 0\n")

    assert parsed is not None
    assert parsed["atom_count"] == 3
    assert parsed["elements"] == ["C", "H", "H"]


def test_parse_molecule_text_counts_smiles_atoms():
    parsed = parse_molecule_text("C1=CC=CC=C1")

    assert parsed is not None
    assert parsed["atom_count"] >= 6


def test_get_molecule_features_prefers_input_txt(tmp_path: Path):
    workdir = tmp_path / "case" / "ESD"
    workdir.mkdir(parents=True)
    (workdir.parent / "input.txt").write_text("2\n\nFe 0 0 0\nCl 0 0 2\n", encoding="utf-8")

    features = get_molecule_features(workdir)

    assert features["atom_count"] == 2
    assert features["has_transition_metal"] is True


def test_predict_stage_profile_scales_with_atom_count():
    small = predict_stage_profile("ox_step_1", atom_count=30)
    large = predict_stage_profile("ox_step_1", atom_count=90)

    assert small["duration_s"] is not None
    assert large["duration_s"] is not None
    assert large["duration_s"] > small["duration_s"]
    assert large["observed_avg_cores"] >= small["observed_avg_cores"]
