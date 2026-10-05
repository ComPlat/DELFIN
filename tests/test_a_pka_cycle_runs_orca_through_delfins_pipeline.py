"""Red control for the pKa ORCA execution layer (delfin/pka.py).

The layer must:
- build each species' XYZ through DELFIN's SMILES converter
  (delfin/tadf_xtb._write_smiles_inputs -> smiles_converter), never a
  second converter;
- write one ORCA Opt+Freq input per species via pka.build_opt_freq_input
  and run it through delfin.orca.run_orca;
- refine energies with a single-point input through run_orca as well;
- assemble the Gibbs energies with read_cycle_gibbs_energies (which
  reuses energies.find_gibbs_energy) and pka_from_cycle.

ORCA itself is mocked everywhere; no real ORCA runs from the tests.
"""

from pathlib import Path

import pytest

from delfin.pka import (
    KNOWN_ACIDS,
    Species,
    build_opt_freq_input,
    plan_species,
    prepare_cycle,
    run_cycle,
)

_RDKIT = pytest.importorskip("rdkit")


def test_prepare_cycle_builds_xyz_from_delfins_smiles_converter(tmp_path):
    calls = []

    def fake_converter(smiles, label, workdir):
        calls.append((smiles, label, workdir))
        xyz = workdir / f"{label}.xyz"
        xyz.write_text("3\ncomment\nC 0 0 0\nO 1 0 0\nO 2 0 0\n", encoding="utf-8")
        return workdir / "start.txt", xyz

    cycle_dir = tmp_path / "phenol"
    cycle_dir.mkdir()
    prepared = prepare_cycle(
        plan_species("phenol"),
        cycle_dir,
        write_smiles_inputs=fake_converter,
    )
    assert [p.label for p in prepared] == [
        "phenol_HA", "phenol_A", "acetic_HA", "acetic_A",
    ]
    assert len(calls) == 4
    # Each species' xyz file exists in the cycle directory.
    for species in prepared:
        assert (cycle_dir / f"{species.label}.xyz").is_file()


def test_prepare_cycle_never_produces_a_species_without_structure(tmp_path):
    cycle_dir = tmp_path / "no_structure"
    cycle_dir.mkdir()
    # A species with a SMILES but a converter that yields nothing must
    # fail loudly instead of silently leaving an ORCA run without a
    # structure.
    def empty_converter(smiles, label, workdir):
        return workdir / "start.txt", workdir / f"{label}.xyz"

    with pytest.raises(RuntimeError):
        prepare_cycle(
            [Species(label="not_an_acid_HA", charge=0, smiles="C")],
            cycle_dir,
            write_smiles_inputs=empty_converter,
        )


def test_run_cycle_executes_orca_per_species_and_anchors_to_the_reference(
    tmp_path, monkeypatch
):
    cycle_dir = tmp_path / "phenol"
    cycle_dir.mkdir()

    def fake_write_smiles_inputs(smiles, label, workdir):
        xyz = workdir / f"{label}.xyz"
        xyz.write_text("3\ncomment\nC 0 0 0\nO 1 0 0\nO 2 0 0\n", encoding="utf-8")
        return workdir / "start.txt", xyz

    def fake_run_orca(input_file_path, output_log, **kwargs):
        # Simulate an ORCA Opt+Freq run: write the thermochemistry line.
        Path(output_log).write_text(
            "Final Gibbs free energy        -307.123456 Eh\n", encoding="utf-8"
        )
        return True

    monkeypatch.setattr(
        "delfin.pka._write_smiles_inputs", fake_write_smiles_inputs
    )
    monkeypatch.setattr("delfin.pka.run_orca", fake_run_orca)

    result = run_cycle("phenol", cycle_dir)
    assert result["pka"] == pytest.approx(
        KNOWN_ACIDS["acetic"]["pka"], abs=1e-6
    )
    # Equal Gibbs energies everywhere -> the isodesmic difference is 0.
    # result["gibbs"] is keyed by cycle keys (target_HA, ...), not labels.
    for key in ("target_HA", "target_A", "reference_HA", "reference_A"):
        assert result["gibbs"][key] == pytest.approx(-307.123456)
    # One ORCA input + one ORCA output per species exist.
    for species in ("phenol_HA", "phenol_A", "acetic_HA", "acetic_A"):
        assert (cycle_dir / f"{species}.inp").is_file()
        assert (cycle_dir / f"{species}.out").is_file()


def test_run_cycle_reports_a_failed_orca_run_instead_of_guessing(
    tmp_path, monkeypatch
):
    cycle_dir = tmp_path / "phenol"
    cycle_dir.mkdir()

    def fake_write_smiles_inputs(smiles, label, workdir):
        xyz = workdir / f"{label}.xyz"
        xyz.write_text("3\ncomment\nC 0 0 0\nO 1 0 0\nO 2 0 0\n", encoding="utf-8")
        return workdir / "start.txt", xyz

    def failing_run_orca(input_file_path, output_log, **kwargs):
        Path(output_log).write_text(
            "ORCA has failed in some way\n", encoding="utf-8"
        )
        return False

    monkeypatch.setattr(
        "delfin.pka._write_smiles_inputs", fake_write_smiles_inputs
    )
    monkeypatch.setattr("delfin.pka.run_orca", failing_run_orca)

    result = run_cycle("phenol", cycle_dir)
    assert result["pka"] is None
    # Every species failed, so all four are reported - the cycle refuses
    # to guess and reports the full failure set.
    assert set(result["failures"]) == {
        "phenol_HA", "phenol_A", "acetic_HA", "acetic_A",
    }
    assert result["gibbs"] == {
        "target_HA": None, "target_A": None,
        "reference_HA": None, "reference_A": None,
    }


def test_run_cycle_never_anchors_from_a_partial_cycle(tmp_path, monkeypatch):
    cycle_dir = tmp_path / "phenol"
    cycle_dir.mkdir()

    def fake_write_smiles_inputs(smiles, label, workdir):
        xyz = workdir / f"{label}.xyz"
        xyz.write_text("3\ncomment\nC 0 0 0\nO 1 0 0\nO 2 0 0\n", encoding="utf-8")
        return workdir / "start.txt", xyz

    def partial_run_orca(input_file_path, output_log, **kwargs):
        # Only the reference acid produces thermochemistry.
        if "acetic" in Path(input_file_path).name:
            Path(output_log).write_text(
                "Final Gibbs free energy        -307.0 Eh\n", encoding="utf-8"
            )
            return True
        Path(output_log).write_text("no thermochemistry\n", encoding="utf-8")
        return True

    monkeypatch.setattr(
        "delfin.pka._write_smiles_inputs", fake_write_smiles_inputs
    )
    monkeypatch.setattr("delfin.pka.run_orca", partial_run_orca)

    result = run_cycle("phenol", cycle_dir)
    assert result["pka"] is None
    assert set(result["failures"]) == {
        "phenol_HA", "phenol_A",
    }
