from pathlib import Path

from delfin import guppy_sampling


def test_sampling_accepts_relaxed_fragment_fallback(monkeypatch, tmp_path):
    start_coords = [
        "H 0.0 0.0 0.0",
        "H 0.0 0.0 0.7",
    ]

    def fake_run_orca(inp_path, out_path, working_dir=None, isolate=True):
        xyz_path = Path(working_dir) / "XTB.xyz"
        xyz_path.write_text(
            "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
            encoding="utf-8",
        )
        return True

    monkeypatch.setattr(guppy_sampling, "run_orca", fake_run_orca)
    monkeypatch.setattr(guppy_sampling, "_extract_total_energy_eh", lambda _path: -1.23)
    monkeypatch.setattr(guppy_sampling, "_fragment_topology_ok", lambda xyz, smiles: False)
    monkeypatch.setattr(
        guppy_sampling,
        "_fragment_topology_relaxed_fallback_ok",
        lambda xyz, smiles: True,
    )
    monkeypatch.setattr(guppy_sampling, "_xyz_passes_final_geometry_checks", lambda xyz, mol: True)
    monkeypatch.setattr(guppy_sampling, "_roundtrip_ring_count_ok", lambda xyz, smiles: True)
    monkeypatch.setattr(guppy_sampling, "_no_spurious_bonds", lambda xyz, smiles: True)

    ok, result, error = guppy_sampling._execute_single_sampling_run(
        run_idx=1,
        start_coords=start_coords,
        start_label="test",
        start_source="quick",
        resolved_charge=0,
        multiplicity=1,
        pal=1,
        maxcore=1000,
        method="XTB2",
        workdir=tmp_path,
        smiles="[H][H]",
        mol_template=object(),
    )

    assert ok is True
    assert error is None
    assert result is not None
    assert result[0] == -1.23


def test_sampling_rejects_relaxed_fragment_fallback_when_geometry_checks_fail(monkeypatch, tmp_path):
    start_coords = [
        "H 0.0 0.0 0.0",
        "H 0.0 0.0 0.7",
    ]

    def fake_run_orca(inp_path, out_path, working_dir=None, isolate=True):
        xyz_path = Path(working_dir) / "XTB.xyz"
        xyz_path.write_text(
            "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
            encoding="utf-8",
        )
        return True

    monkeypatch.setattr(guppy_sampling, "run_orca", fake_run_orca)
    monkeypatch.setattr(guppy_sampling, "_extract_total_energy_eh", lambda _path: -1.23)
    monkeypatch.setattr(guppy_sampling, "_fragment_topology_ok", lambda xyz, smiles: False)
    monkeypatch.setattr(
        guppy_sampling,
        "_fragment_topology_relaxed_fallback_ok",
        lambda xyz, smiles: True,
    )
    monkeypatch.setattr(guppy_sampling, "_xyz_passes_final_geometry_checks", lambda xyz, mol: False)
    monkeypatch.setattr(guppy_sampling, "_roundtrip_ring_count_ok", lambda xyz, smiles: True)
    monkeypatch.setattr(guppy_sampling, "_no_spurious_bonds", lambda xyz, smiles: True)

    ok, result, error = guppy_sampling._execute_single_sampling_run(
        run_idx=1,
        start_coords=start_coords,
        start_label="test",
        start_source="quick",
        resolved_charge=0,
        multiplicity=1,
        pal=1,
        maxcore=1000,
        method="XTB2",
        workdir=tmp_path,
        smiles="[H][H]",
        mol_template=object(),
    )

    assert ok is False
    assert result is None
    assert error == "Topology changed: final geometry checks failed after XTB"
