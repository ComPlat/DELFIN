import json

import pytest

from delfin.tools import run_step
from delfin.tools.bio import BindingBox, StructureFormatError, read_structure
from delfin.tools.bio._confidence import parse_pae_json
from test_how_sure_a_predicted_model_is_is_read_from_the_file import PREDICTED_CIF, PREDICTED_PDB


def test_uniform_b_factors_are_not_automatically_plddt(tmp_path):
    path = tmp_path / "unknown.pdb"
    path.write_text(PREDICTED_PDB)
    assert read_structure(path).plddt_per_residue() == ()
    result = run_step("structure_confidence", structure_file=str(path), work_dir=tmp_path / "run")
    assert result.ok
    assert result.data["has_confidence"] is False
    assert "origin is not established" in result.data["statement"]


def test_explicit_b_factor_interpretation_is_recorded(tmp_path):
    path = tmp_path / "predicted.pdb"
    path.write_text(PREDICTED_PDB)
    result = run_step("structure_confidence", structure_file=str(path),
                      b_factor_is_plddt=True, work_dir=tmp_path / "run")
    assert result.ok
    assert result.data["plddt_source"] == "pdb_b_factor"
    assert result.data["confidence_interpretation"] == "explicit_b_factor"
    assert result.data["model_id"] == 1


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, True])
def test_pae_entries_must_be_finite_nonnegative_numbers(value):
    with pytest.raises(StructureFormatError):
        parse_pae_json({"predicted_aligned_error": [[value]], "max_predicted_aligned_error": 5})


def test_pae_dimensions_must_match_the_selected_structure(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text(PREDICTED_PDB)
    pae = tmp_path / "pae.json"
    pae.write_text(json.dumps({"predicted_aligned_error": [[0]], "max_predicted_aligned_error": 5}))
    result = run_step("structure_confidence", structure_file=str(path),
                      pae_file=str(pae), work_dir=tmp_path / "run")
    assert not result.ok
    assert "residue" in result.error.lower()
    assert not (tmp_path / "run" / "structure_confidence.json").exists()


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_binding_boxes_refuse_nonfinite_coordinates_and_sizes(value):
    with pytest.raises(ValueError):
        BindingBox(center=(value, 0, 0), size=(1, 1, 1))
    with pytest.raises(ValueError):
        BindingBox(center=(0, 0, 0), size=(value, 1, 1))


def test_binding_box_units_are_not_silently_reinterpreted():
    with pytest.raises(ValueError, match="unit"):
        BindingBox.from_dict({"center": [0, 0, 0], "size": [1, 1, 1], "unit": "nm"})


def test_canonical_metric_fields_follow_author_residue_numbering(tmp_path):
    text = PREDICTED_CIF.replace("ma_qa_metric_id", "metric_id").replace(
        "ma_qa_metric.value", "metric_value",
    )
    text = text.replace("MET A A 1 1 ?", "MET A Z 1 41 A")
    text = text.replace("ALA A A 2 2 ?", "ALA A Z 2 42 ?")
    text = text.replace("GLY A A 3 3 ?", "GLY A Z 3 43 ?")
    path = tmp_path / "renumbered.cif"
    path.write_text(text)
    values = read_structure(path).plddt_per_residue()
    assert [(residue.chain_id, residue.seq_id, residue.insertion_code, score)
            for residue, score in values] == [("Z", 41, "A", 90.5), ("Z", 42, "", 40), ("Z", 43, "", 20)]


def test_model_specific_metrics_are_not_mixed(tmp_path):
    text = """data_models
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N MET A X 1 10 0 0 0 1
ATOM 2 N N MET A X 1 10 1 0 0 2
#
loop_
_ma_qa_metric.id
_ma_qa_metric.type
1 pLDDT
#
loop_
_ma_qa_metric_local.model_id
_ma_qa_metric_local.metric_id
_ma_qa_metric_local.label_asym_id
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_value
1 1 A 1 90
2 1 A 1 20
#
"""
    path = tmp_path / "models.cif"
    path.write_text(text)
    result = run_step("structure_confidence", structure_file=str(path),
                      model=2, work_dir=tmp_path / "run")
    assert result.ok, result.error
    assert result.data["model_id"] == 2
    assert result.data["mean_plddt"] == 20
    assert result.data["plddt_source"] == "mmcif_ma_qa_metric_local"
