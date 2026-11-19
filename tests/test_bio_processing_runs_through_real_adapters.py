import json
from pathlib import Path
import shutil

import pytest

from delfin.tools import run_step
from delfin.tools.bio.artifacts import validate_artifact


@pytest.fixture
def files(tmp_path):
    root = Path(__file__).resolve().parents[1] / "examples" / "biomolecular_foundation"
    for name in ("receptor.pdb", "reference.sdf", "poses.sdf"):
        shutil.copyfile(root / name, tmp_path / name)
    return tmp_path


def test_properties_and_depiction_preserve_the_source(files):
    original = (files / "poses.sdf").read_bytes()
    result = run_step("bio_ligand_properties", ligands_file=files / "poses.sdf", work_dir=files / "properties")
    assert result.ok, result.error
    rows = result.data["records"]
    assert len(rows) == 3
    assert len({row["chemical_state_id"] for row in rows}) == 1
    assert rows[0]["molecular_weight"] > 0
    depiction = run_step("bio_ligand_depict", ligands_file=files / "poses.sdf", work_dir=files / "depiction")
    assert depiction.ok, depiction.error
    assert "<svg" in depiction.artifacts["svg"].read_text()
    assert (files / "poses.sdf").read_bytes() == original


def test_tautomers_have_parent_identity_and_explicit_no_ph_claim(files):
    source = files / "ligand.smi"
    source.write_text("CC(=O)CC(=O)C\n")
    result = run_step("bio_ligand_states", ligand_file=source, work_dir=files / "states")
    assert result.ok, result.error
    assert len(result.data["states"]) > 1
    assert result.data["protonation"] == "no pH model"
    validate_artifact(result.artifacts["ligand_state_set"], "ligand_state_set")


def test_structure_selection_is_explicit_and_not_receptor_preparation(files):
    result = run_step("bio_structure_select", structure_file=files / "receptor.pdb", chains=["A"],
                      waters="drop", work_dir=files / "select")
    assert result.ok, result.error
    assert not result.data["receptor_prepared"]
    assert len(result.data["selected"]) == 5
    validate_artifact(result.artifacts["macromolecule"], "macromolecule")
    missing = run_step("bio_structure_select", structure_file=files / "receptor.pdb", chains=["Z"], work_dir=files / "missing")
    assert not missing.ok


def test_prediction_import_and_selection_verify_the_actual_files(files):
    imported = run_step("bio_prediction_import", source_dir=files,
                        models=[{"id": "model_7", "path": "receptor.pdb"}], engine="fixture", engine_version="1",
                        work_dir=files / "import")
    assert imported.ok, imported.error
    bundle = imported.artifacts["prediction_bundle"]
    selected = run_step("bio_model_select", bundle_file=bundle, model_id="model_7", work_dir=files / "model")
    assert selected.ok, selected.error
    assert selected.artifacts["macromolecule"].read_bytes() == (files / "receptor.pdb").read_bytes()
    (bundle.parent / "model_1.pdb").write_text((files / "receptor.pdb").read_text() + "REMARK modified\n")
    refused = run_step("bio_model_select", bundle_file=bundle, model_id="model_7", work_dir=files / "tampered")
    assert not refused.ok and "checksum" in refused.error


def test_prediction_import_cannot_escape_its_source_directory(files):
    result = run_step("bio_prediction_import", source_dir=files / "empty",
                      models=[{"id": "x", "path": "../receptor.pdb"}], engine="fixture", engine_version="1",
                      work_dir=files / "bad_import")
    assert not result.ok and "escapes" in result.error


def test_visualization_is_self_contained_and_exports_both_desktop_viewers(files):
    result = run_step("bio_visualize", structure_file=files / "receptor.pdb", poses_file=files / "poses.sdf",
                      work_dir=files / "view")
    assert result.ok, result.error
    root = result.artifacts["html_view"].parent
    assert (root / "3Dmol-min.js").is_file()
    for filename in ("view_pymol.py", "view_chimerax.py"):
        compile((root / filename).read_text(), filename, "exec")
    page = (root / "view.html").read_text()
    assert 'src="3Dmol-min.js"' in page
    assert 'src="https:' not in page
    data = json.loads(page.split('<script id="data" type="application/json">')[1].split('</script>')[0])
    assert len(data["poses"]) == 3
    assert result.data["external_viewers_launched"] is False


def test_confidence_integrity_is_checked_on_model_selection(files):
    (files / 'confidence.json').write_text('{"plddt": [70]}')
    imported = run_step('bio_prediction_import', source_dir=files,
                        models=[{'id': 'one', 'path': 'receptor.pdb', 'confidence': 'confidence.json'}],
                        engine='fixture', engine_version='1', work_dir=files / 'import')
    assert imported.ok, imported.error
    bundle = imported.artifacts['prediction_bundle']
    (bundle.parent / 'confidence_1.json').write_text('{"plddt": [99]}')
    with pytest.raises(ValueError, match='confidence checksum'):
        validate_artifact(bundle, 'prediction_bundle')


@pytest.mark.parametrize('number', ['NaN', 'Infinity', '-Infinity', '1e999'])
def test_nonfinite_json_is_rejected(files, number):
    from delfin.tools.bio.artifacts import read_json
    source = files / 'invalid.json'
    source.write_text('{"value": ' + number + '}')
    with pytest.raises(ValueError, match='non-finite'):
        read_json(source)


def test_pose_selection_and_box_preserve_the_selected_coordinate_frame(files):
    from delfin.tools.bio.processing import load_ligands
    selected = run_step('bio_ligand_select', ligands_file=files / 'poses.sdf', record_id=2, work_dir=files / 'selected')
    assert selected.ok, selected.error
    before = load_ligands(files / 'poses.sdf')[1].GetConformer().GetPositions()
    after = load_ligands(selected.artifacts['ligand'])[0].GetConformer().GetPositions()
    assert (before == after).all()
    boxed = run_step('bio_box_from_ligand', ligand_file=selected.artifacts['ligand'], padding=3.0, work_dir=files / 'box')
    assert boxed.ok, boxed.error
    box = boxed.data['box']
    for point in after:
        assert all(abs(point[i] - box['center'][i]) <= box['size'][i] / 2 for i in range(3))
    invalid = run_step('bio_ligand_select', ligands_file=files / 'poses.sdf', record_id=4, work_dir=files / 'invalid')
    assert not invalid.ok


def test_box_refuses_unembedded_smiles(files):
    source = files / 'ligand.smi'
    source.write_text('CCO\n')
    result = run_step('bio_box_from_ligand', ligand_file=source, work_dir=files / 'no_coords')
    assert not result.ok and '3D' in result.error
