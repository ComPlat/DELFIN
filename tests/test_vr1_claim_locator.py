"""VR1 phase 1: claim_locator locates a number+unit in the real output file.

Synthetic ORCA-style fixtures only (licence rule: no real ORCA/TURBOMOLE data).
The fixture text matches DELFIN's own parser patterns (energies.py).
"""

from pathlib import Path

import pytest

from delfin.agent.claim_locator import verify_claim

# Hand-written synthetic ORCA-style output, NOT produced by ORCA.
# It reports one electronic energy, one Gibbs free energy, one ZPE,
# three vibrational frequencies and a two-atom geometry.
SYNTHETIC_OUT = """\
# Synthetic ORCA-style output -- hand-written for tests, not produced by ORCA.
FINAL SINGLE POINT ENERGY      -76.026780443768

Final Gibbs free energy                 -76.030514123   Eh

Zero point energy                      0.003733680     Eh

VIBRATIONAL FREQUENCIES (cm**(-1))
   0:     3441.3
   1:     1583.1
   2:     1021.7

CARTESIAN COORDINATES (ANGSTROEM)
    O    0.000000    0.000000    0.121000
    H    0.000000    0.000000    0.861000
"""

# A bare output that computes nothing (only a header): no source lines exist.
BARE_OUT = """\
# Synthetic ORCA-style output -- hand-written for tests, not produced by ORCA.
A job that terminated before writing any result blocks.
"""


@pytest.fixture
def calc_folder(tmp_path: Path) -> Path:
    """A calc folder holding one synthetic .out file."""
    folder = tmp_path / "water"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "water.out").write_text(SYNTHETIC_OUT, encoding="utf-8")
    return folder


@pytest.fixture
def bare_folder(tmp_path: Path) -> Path:
    """A calc folder whose .out computes no quantity."""
    folder = tmp_path / "bare"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "bare.out").write_text(BARE_OUT, encoding="utf-8")
    return folder


# --- item 1: verified, with real file:line, unit, parsed value ---------------

def test_verified_electronic_energy_Eh(calc_folder: Path) -> None:
    rec = verify_claim(-76.026780443768, "Eh", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "energy"
    assert rec.file is not None and rec.file.endswith("water.out")
    # 'FINAL SINGLE POINT ENERGY' is on line 2 (1-based) of SYNTHETIC_OUT.
    assert rec.line == 2
    assert rec.parsed_value == pytest.approx(-76.026780443768)
    assert rec.parsed_unit == "Eh"
    assert rec.tolerance is not None and rec.tolerance > 0


def test_verified_energy_in_kcal_per_mol_converts(calc_folder: Path) -> None:
    expected_kcal = -76.026780443768 * 627.509
    rec = verify_claim(expected_kcal, "kcal/mol", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "energy"
    assert rec.line == 2


def test_verified_frequency(calc_folder: Path) -> None:
    rec = verify_claim(3441.3, "cm-1", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "frequency"
    assert rec.parsed_value == pytest.approx(3441.3)
    assert rec.parsed_unit == "cm-1"


def test_true_interatomic_distance_is_verified(calc_folder: Path) -> None:
    # The true O-H separation: H z(0.861000) - O z(0.121000) = 0.740000 A.
    # A distance needs TWO atoms; this is the only interatomic distance.
    rec = verify_claim(0.740, "Angstrom", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "distance"
    assert rec.parsed_value == pytest.approx(0.740)


def test_single_coordinate_is_not_a_distance(calc_folder: Path) -> None:
    # A coordinate (z of O or H) is NOT an interatomic distance. A planted
    # wrong number must be caught, so these must be mismatch, never verified.
    for coord in (0.121, 0.861):
        rec = verify_claim(coord, "Angstrom", calc_folder)
        assert rec.verdict == "mismatch", f"{coord} A is a coordinate, not a bond"


def test_verified_distance_in_bohr_converts(calc_folder: Path) -> None:
    # The true O-H separation 0.740 Angstrom expressed in bohr.
    expected_bohr = 0.740 / 0.529177
    rec = verify_claim(expected_bohr, "bohr", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "distance"


# --- item 2: mismatch --------------------------------------------------------

def test_wrong_value_is_mismatch(calc_folder: Path) -> None:
    # -76.0 vs -76.026780443768 is far beyond the 1.0 kcal/mol energy tolerance.
    rec = verify_claim(-76.0, "Eh", calc_folder)
    assert rec.verdict == "mismatch"


def test_off_by_sign_is_mismatch(calc_folder: Path) -> None:
    rec = verify_claim(+76.026780443768, "Eh", calc_folder)
    assert rec.verdict == "mismatch"


def test_off_by_unit_is_mismatch(calc_folder: Path) -> None:
    # 0.861 bohr is 0.4556 Angstrom; the geometry holds 0.861 Angstrom.
    rec = verify_claim(0.861, "bohr", calc_folder)
    assert rec.verdict == "mismatch"


def test_frequency_off_by_sign_is_mismatch(calc_folder: Path) -> None:
    rec = verify_claim(-3441.3, "cm-1", calc_folder)
    assert rec.verdict == "mismatch"


def test_rounded_beyond_tolerance_is_mismatch(calc_folder: Path) -> None:
    # Energy tolerance is 1.0 kcal/mol; off by 5 kcal/mol is out.
    rec = verify_claim(-76.026780443768 * 627.509 + 5.0, "kcal/mol", calc_folder)
    assert rec.verdict == "mismatch"


# --- item 2: no-source -------------------------------------------------------

def test_quantity_not_computed_is_no_source(bare_folder: Path) -> None:
    # The bare job computed no electronic energy.
    rec = verify_claim(-76.0, "Eh", bare_folder)
    assert rec.verdict == "no-source"


def test_dipole_not_computed_is_no_source(calc_folder: Path) -> None:
    # The synthetic output has no DIPOLE block.
    rec = verify_claim(5.3, "Debye", calc_folder)
    assert rec.verdict == "no-source"


# --- item 3: machine-readable record ----------------------------------------

def test_record_is_machine_readable(calc_folder: Path) -> None:
    rec = verify_claim(-76.026780443768, "Eh", calc_folder)
    d = rec.to_dict()
    assert d["verdict"] == "verified"
    assert d["value"] == pytest.approx(-76.026780443768)
    assert d["unit"] == "Eh"
    assert isinstance(d["line"], int)


# --- item 4: prose numbers do not false-alarm --------------------------------

def test_count_unit_is_estimate(calc_folder: Path) -> None:
    # "3 runs" is a count, not a physical result: no quantity family.
    rec = verify_claim(3, "runs", calc_folder)
    assert rec.verdict == "estimate"


def test_year_is_estimate(calc_folder: Path) -> None:
    rec = verify_claim(2026, "a.d.", calc_folder)
    assert rec.verdict == "estimate"


def test_settings_number_is_estimate(calc_folder: Path) -> None:
    # A setting like "pal 2" / 2 cores is not an output-derived result.
    rec = verify_claim(2, "cores", calc_folder)
    assert rec.verdict == "estimate"
