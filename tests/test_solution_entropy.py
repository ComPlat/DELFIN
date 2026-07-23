from __future__ import annotations

import json
import math

from delfin.solvation_entropy import (
    MoleculeGeometry,
    calculate_reaction_entropy_correction,
    calculate_solution_entropy,
    parse_orca_thermochemistry,
    solvent_params,
    vdw_volume,
)
from delfin.tools import list_steps, run_step


def test_solvent_alias_and_custom_override():
    dmf = solvent_params("dmf")
    assert dmf.name == "dmf"
    custom = solvent_params(
        "my-solvent",
        density_g_ml=1.0,
        molar_mass_g_mol=50.0,
        solvent_vdw_volume_A3=45.0,
        permittivity=10.0,
    )
    assert custom.source == "user-provided"
    assert custom.vdw_volume_A3 == 45.0


def test_single_atom_vdw_volume_is_exact_sphere():
    geom = MoleculeGeometry(("He",), ((0.0, 0.0, 0.0),))
    assert math.isclose(vdw_volume(geom), 4.0 / 3.0 * math.pi * 1.40**3)


def test_calculate_solution_entropy_geometry_only():
    geom = MoleculeGeometry(("He",), ((0.0, 0.0, 0.0),))
    result = calculate_solution_entropy(geom, solvent="benzene", volume_samples=1000)
    assert result.solvent == "benzene"
    assert result.S_solv < 0.0
    assert result.S_soln is None
    assert result.n_imaginary_ignored == 0


def test_parse_minimal_orca_thermochemistry(tmp_path):
    out = tmp_path / "calc.out"
    out.write_text(
        """
        * O   R   C   A *
        THERMOCHEMISTRY
        Temperature         ...   298.15 K
        Point Group         C1
        Symmetry Number     1
        ENTROPY
        Rotational entropy      0.00001500 Eh/K
        freq.        120.00
        freq.        -24.00
        """,
        encoding="utf-8",
    )
    thermo = parse_orca_thermochemistry(out)
    assert thermo.program == "orca"
    assert thermo.temperature_K == 298.15
    assert thermo.symmetry_number == 1
    assert thermo.frequencies_cm1 == (120.0, -24.0)
    assert thermo.rotational_entropy_cal_mol_K is not None


def test_run_step_solution_entropy_writes_report(tmp_path):
    xyz = tmp_path / "he.xyz"
    xyz.write_text("1\nHe\nHe 0.0 0.0 0.0\n", encoding="utf-8")

    result = run_step(
        "solution_entropy",
        geometry=xyz,
        work_dir=tmp_path / "entropy",
        solvent="benzene",
        volume_samples=1000,
    )

    assert result.ok, result.error
    assert "solution_entropy_report" in result.artifacts
    assert result.data["report_path"].endswith("solution_entropy.json")
    data = json.loads(result.artifacts["solution_entropy_report"].read_text())
    assert data["S_solv"] == result.data["S_solv"]


def test_reaction_entropy_correction_from_values():
    corr = calculate_reaction_entropy_correction(
        {"A": -10.0, "B": -20.0, "C": -40.0},
        {"A": -1, "B": -1, "C": 1},
        temperature_K=300.0,
        uncorrected_delta_g_kcal_mol=5.0,
    )
    assert corr.delta_S_solv_cal_mol_K == -10.0
    assert corr.delta_G_entropy_corr_kcal_mol == 3.0
    assert corr.delta_G_corrected_kcal_mol == 8.0


def test_run_step_reaction_entropy_accepts_json_strings(tmp_path):
    result = run_step(
        "reaction_solution_entropy",
        work_dir=tmp_path / "reaction",
        species='{"A": -10.0, "B": -20.0, "C": -40.0}',
        stoichiometry='{"A": -1, "B": -1, "C": 1}',
        temperature=300.0,
    )
    assert result.ok, result.error
    assert result.data["delta_G_entropy_corr_kcal_mol"] == 3.0


def test_solution_entropy_registered_and_applications_visible():
    steps = list_steps()
    assert "solution_entropy" in steps
    assert steps["solution_entropy"].contract().category == "solvation"

    from delfin.tools import platform

    apps = set(platform.list_applications())
    assert "solution_entropy" in apps
    assert "reaction_solution_entropy" in apps
