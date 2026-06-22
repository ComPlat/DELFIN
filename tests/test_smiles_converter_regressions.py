import pytest

from delfin import smiles_converter as sc
from delfin.dashboard import input_processing


def test_segment_distance_sq_exists_and_returns_expected_value():
    dist_sq = sc._segment_distance_sq(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (1.0, 1.0, 0.0),
    )
    assert dist_sq == pytest.approx(1.0)


def test_smiles_to_xyz_isomers_uses_weighted_ob_search_when_non_deterministic(monkeypatch):
    recorded = []

    class DummyMol:
        def RemoveAllConformers(self):
            return None

    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "OPENBABEL_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smiles: True)
    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", lambda _smiles: [])
    monkeypatch.setattr(sc, "_prepare_mol_for_embedding", lambda *_args, **_kwargs: DummyMol())
    monkeypatch.setattr(sc, "_donor_type_map", lambda _mol: {})
    monkeypatch.setattr(sc, "_embed_multiple_confs_robust", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(sc, "smiles_to_xyz", lambda *_args, **_kwargs: ("C 0.0 0.0 0.0\n", None))

    def _fake_ob(*_args, **kwargs):
        recorded.append(kwargs.get("deterministic"))
        return [], "no conformers"

    monkeypatch.setattr(sc, "_openbabel_generate_conformer_xyz", _fake_ob)

    results, error = sc.smiles_to_xyz_isomers(
        "[Fe](C)(C)(C)(C)(C)(C)",
        deterministic=False,
        include_binding_mode_isomers=False,
    )

    assert error is None
    assert results
    assert recorded
    assert all(flag is False for flag in recorded)


def test_dashboard_isomer_wrapper_passes_deterministic_flag(monkeypatch):
    recorded = []

    def _fake_converter(*_args, **kwargs):
        recorded.append(kwargs.get("deterministic"))
        return [("C 0.0 0.0 0.0\n", "isomer")], None

    monkeypatch.setattr(
        input_processing,
        "_delfin_smiles_to_xyz_isomers",
        _fake_converter,
    )

    results, error = input_processing.smiles_to_xyz_isomers(
        "C",
        deterministic=False,
    )

    assert error is None
    assert results == [("C 0.0 0.0 0.0\n", 1, "isomer")]
    assert recorded == [False]
