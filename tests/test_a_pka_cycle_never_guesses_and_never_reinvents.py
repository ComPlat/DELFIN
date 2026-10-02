"""Red control for the W10/s17 pKa module.

Every case below must fail on the unmodified tree: ``delfin/pka.py`` does
not exist yet, so the import at the top makes the whole file error out.

Protocol under test: isodesmic proton transfer against a reference acid,

    pKa(target) = pKa(ref) + (dG_deprot(target) - dG_deprot(ref)) / (ln10 * R * T)

with dG_deprot = G(conjugate base) - G(acid), all species 1 M in solution
(SMD/CPCM continuum, ORCA Manual 6.1.1 section 2.13.3).  The 1 atm -> 1 M
concentration correction (Manual: +1.89 kcal/mol) cancels inside the
isodesmic difference because both acids are treated identically.
"""

import math
from pathlib import Path

import pytest

from delfin.pka import (
    KNOWN_ACIDS,
    REFERENCE_ACID,
    build_opt_freq_input,
    build_single_point_input,
    pka_from_cycle,
    plan_species,
    read_cycle_gibbs_energies,
)

KCAL_PER_HARTREE = 627.5094740631
GAS_CONSTANT_KCAL = 1.98720425864083e-3


def _pk_factor(temperature_k: float = 298.15) -> float:
    """kcal/mol per pK unit: ln(10) * R * T."""
    return math.log(10) * GAS_CONSTANT_KCAL * temperature_k


def test_the_cycle_anchors_the_target_acid_to_the_reference_pka():
    gibbs = {
        "target_HA": -700.0,
        "target_A": -699.0,
        "reference_HA": -500.0,
        "reference_A": -499.5,
    }
    pka = pka_from_cycle(gibbs, reference_pka=4.756)
    expected = 4.756 + (1.0 - 0.5) * KCAL_PER_HARTREE / _pk_factor()
    assert pka == pytest.approx(expected, abs=1e-3)


def test_the_cycle_signs_follow_the_deprotonation_energies():
    # Target deprotonates more easily than the reference -> pKa below ref.
    gibbs = {
        "target_HA": -700.0,
        "target_A": -699.6,
        "reference_HA": -500.0,
        "reference_A": -499.5,
    }
    pka = pka_from_cycle(gibbs, reference_pka=4.756)
    expected = 4.756 + (0.4 - 0.5) * KCAL_PER_HARTREE / _pk_factor()
    assert pka == pytest.approx(expected, abs=1e-3)
    assert pka < 4.756


def test_the_reference_acid_is_acetic_at_the_experimental_pka():
    assert REFERENCE_ACID == "acetic"
    assert KNOWN_ACIDS["acetic"]["pka"] == pytest.approx(4.756)
    assert KNOWN_ACIDS["acetic"]["ha_smiles"] == "CC(=O)O"
    assert KNOWN_ACIDS["acetic"]["a_smiles"] == "CC(=O)[O-]"


def test_the_plan_covers_both_species_of_target_and_reference():
    plan = plan_species("phenol")
    labels = [species.label for species in plan]
    assert labels == [
        "phenol_HA",
        "phenol_A",
        "acetic_HA",
        "acetic_A",
    ]
    for species in plan:
        assert species.multiplicity == 1
    charges = {species.label: species.charge for species in plan}
    assert charges["phenol_HA"] == 0
    assert charges["phenol_A"] == -1
    assert charges["acetic_A"] == -1


def test_the_opt_freq_input_solvates_with_smd_and_asks_for_frequencies():
    xyz_text = "2\nacetic\nC 0.0 0.0 0.0\nO 1.3 0.0 0.0\n"
    text = build_opt_freq_input(
        xyz_text, charge=0, multiplicity=1, solvent="water"
    )
    assert "B3LYP" in text
    assert "def2-SVP" in text
    assert "Opt" in text
    assert "Freq" in text
    assert "SMD(water)" in text
    assert "%maxcore" in text
    assert "%pal nprocs" in text
    assert "* xyz 0 1" in text
    assert "C 0.0 0.0 0.0" in text
    assert text.rstrip().splitlines()[-1].strip() == "*"


def test_the_single_point_input_refines_the_energy_without_reoptimizing():
    xyz_text = "2\nacetic\nC 0.0 0.0 0.0\nO 1.3 0.0 0.0\n"
    text = build_single_point_input(
        xyz_text, charge=-1, multiplicity=1, solvent="water"
    )
    assert "def2-TZVP" in text
    assert "Freq" not in text
    assert "SMD(water)" in text
    assert "* xyz -1 1" in text


def test_gibbs_energies_are_read_from_orca_outputs_not_reinvented(tmp_path):
    out = tmp_path / "acetic_HA.out"
    out.write_text(
        "Final Gibbs free energy        -700.123456 Eh\n",
        encoding="utf-8",
    )
    values = read_cycle_gibbs_energies({"acetic_HA": out})
    assert values["acetic_HA"] == pytest.approx(-700.123456)


def test_a_missing_output_is_reported_and_the_cycle_refuses_to_guess(tmp_path):
    missing = tmp_path / "does_not_exist.out"
    values = read_cycle_gibbs_energies({"phenol_HA": missing})
    assert values["phenol_HA"] is None
    gibbs = {
        "target_HA": None,
        "target_A": -699.0,
        "reference_HA": -500.0,
        "reference_A": -499.5,
    }
    with pytest.raises(ValueError):
        pka_from_cycle(gibbs, reference_pka=4.756)
