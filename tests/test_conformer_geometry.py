# SYNTHETIC test fixture written for DELFIN -- not ORCA output
"""Independent synthetic geometry and ensemble regression tests."""
import json
from pathlib import Path

import numpy as np
import pytest

from delfin.analysis_tools.conformer_shape_classification import (
    analyze_ensemble,
    ez_cip_assignment,
    custom_torsion,
    load_ensemble,
    shape_descriptors,
)
from delfin.analysis_tools.ensemble_profiling import (
    _select_and_rank,
    export_selection,
    profile_ensemble,
    resolve_selection,
    suggest_classifications,
)


def make_xyz(tmp_path: Path, *, goat=False) -> Path:
    name = "synthetic.finalensemble.xyz" if goat else "synthetic_crest.xyz"
    path = tmp_path / name
    # Original invented coordinates: H-O-H, two variations.
    energies = (0.0, 1.5) if goat else (-1.0, -1.001)
    positions = (
        ((0.0, 0.0, 0.0), (0.95, 0.0, 0.0), (-0.239, 0.927, 0.0)),
        ((0.0, 0.0, 0.0), (0.99, 0.0, 0.0), (-0.20, 1.01, 0.1)),
    )
    lines = []
    for e, frame in zip(energies, positions):
        lines.extend(["3", f"Energy {e}" if goat else f"energy: {e}"])
        lines.extend(f"{sym} {x:.3f} {y:.3f} {z:.3f}"
                     for sym, (x, y, z) in zip(("O", "H", "H"), frame))
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.mark.parametrize("goat", [False, True])
def test_reader_uses_correct_energy_units(tmp_path, goat):
    data = load_ensemble(make_xyz(tmp_path, goat=goat))
    delta_kcal = (data.energies[1] - data.energies[0]) * 627.5094740631
    assert delta_kcal == pytest.approx(1.5 if goat else -0.6275094740631)


def test_descriptor_invariants():
    pos = np.array([[0., 0., 0.], [1., 0., 0.],
                    [2., 0., 0.], [3., 0., 0.]])
    desc = shape_descriptors(pos, ["C"] * 4, "uniform")
    assert desc["dmax"].value == pytest.approx(3.)
    assert desc["npr1"].value == pytest.approx(0.)
    assert desc["kappa2"].value == pytest.approx(desc["asphericity"].value)


def test_profiling_and_suggestion_consistency(tmp_path):
    path = make_xyz(tmp_path)
    assert len(analyze_ensemble(path)["frames"]) == 2
    profile = profile_ensemble(path, use_cache=False)
    assert profile["n_conformers"] == 2
    assert profile["composition"] == {"O": 1, "H": 2}
    assert profile["descriptor_variation"]["dmax"]["spread"] > 0
    for suggestion in suggest_classifications(profile):
        if suggestion["kind"] != "continuous":
            continue
        selection = resolve_selection(suggestion["selection_query"])
        assert selection["descriptor"] == suggestion["criterion"]
        assert selection["selection"] == suggestion["selection"]
        assert selection["descending"] == suggestion["descending"]
    selection = resolve_selection("die zwei mit der größten räumlichen Ausdehnung")
    assert selection["n"] == 2 and selection["descriptor"] == "dmax"
    rows = _select_and_rank(path, "dmax", "whole", True, 2, use_cache=False)
    assert {r[0] for r in rows} == {0, 1}


def test_export_preserves_coordinates_and_mapping(tmp_path):
    source = make_xyz(tmp_path)
    original = source.read_text()
    paths = export_selection(source, [1, 0], "synthetic_selection")
    data = json.loads(paths["json"].read_text())
    assert [r["conformer_id"] for r in data["conformer_id_to_export_frame"]] == [1, 0]
    assert source.read_text() == original
    src_lines, out_lines = original.splitlines(), paths["xyz"].read_text().splitlines()
    assert out_lines[2:5] == src_lines[7:10]
    assert out_lines[7:10] == src_lines[2:5]
    with pytest.raises(FileExistsError):
        export_selection(source, [1, 0], "synthetic_selection")


def test_coincident_torsion_atoms_report_error_not_nan():
    # Coincident middle atoms make the torsion undefined; the engine must
    # say so explicitly instead of reporting NaN with status "ok"
    # (found by reviewer nacht-s16, adversarial probe on the base commit).
    coords = np.array([[0., 0., 0.], [1., 0., 0.], [1., 0., 0.], [1., 1., 0.]])
    status = custom_torsion(coords, 0, 1, 2, 3)
    assert status.status == "error", status
    assert status.value is None, status


def test_bond_order_failure_is_not_a_stereochemical_guess(monkeypatch):
    from rdkit.Chem import rdDetermineBonds
    def fail(*args, **kwargs):
        raise ValueError("synthetic charged metal complex perception failure")
    monkeypatch.setattr(rdDetermineBonds, "DetermineBonds", fail)
    result = ez_cip_assignment(["Fe", "O"],
                               np.array([[0., 0., 0.], [0., 0., 1.8]]),
                               charge=2)
    assert len(result) == 1
    assert result[0]["status"] == "unknown"
    assert result[0]["stereo"] is None
    assert "not independently validated" in result[0]["kind"]
