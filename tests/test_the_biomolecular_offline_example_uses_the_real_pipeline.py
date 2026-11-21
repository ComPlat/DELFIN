import json
import runpy
from pathlib import Path

import pytest


def test_the_example_preserves_inputs_and_wires_a_ligand_reference(tmp_path):
    example = Path(__file__).resolve().parents[1] / "examples" / "biomolecular_foundation"
    inputs = {path: path.read_bytes() for path in example.glob("*.sdf")}
    inputs[example / "receptor.pdb"] = (example / "receptor.pdb").read_bytes()
    run_example = runpy.run_path(str(example / "run_example.py"))["run_example"]
    report = run_example(tmp_path / "run")
    assert [step["name"] for step in report["steps"]] == [
        "biostructure_summary", "ligand_prepare", "pose_rmsd",
    ]
    assert report["live_docking_performed"] is False
    assert report["observed_n_under_cutoff"] == 2
    artifact = tmp_path / "run" / report["steps"][-1]["artifacts"]["json_report"]
    poses = json.loads(artifact.read_text())["poses"]
    assert [pose["rmsd_A"] for pose in poses] == pytest.approx([0, 1, 3], abs=1e-4)
    assert report["steps"][1]["data"]["method"] == "coordinates-kept"
    assert all(path.read_bytes() == original for path, original in inputs.items())
    with pytest.raises(ValueError, match="empty"):
        run_example(tmp_path / "run")
