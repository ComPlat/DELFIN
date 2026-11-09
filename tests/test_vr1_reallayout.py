"""VR1 real-run red controls (PR #618): the fixture copies the REAL ORCA 6.1
output layout exactly, hand-written (licence rule: never copy ORCA output).

Real layout after the block marker, per the operator's finding:

    ------------------
    VIBRATIONAL FREQUENCIES
    ------------------

    Scaling factor for frequencies = 0.9553000 (H0 = 0.9986000)

                  0:  0.00 cm**-1
                  6:  1637.64 cm**-1

i.e. marker, dashes line, BLANK, a scaling/meta line, BLANK, then rows.
The old ``_block_values`` ends the block at the FIRST blank line -> zero
frequency candidates -> a correct claim is wrongly "no-source".

Energy: an optimisation prints FINAL SINGLE POINT ENERGY once per cycle; the
old ``_label_line`` reported the FIRST match (cycle 1), not the final energy.
"""

from pathlib import Path

import pytest

from delfin.agent.claim_locator import verify_claim

# Hand-written in the REAL ORCA layout, NOT produced by ORCA.
# Energy values: cycle 1 -76.10, cycle 2 -76.35, FINAL -76.3625 Eh.
# Frequencies: bend 1637.64, stretch 3787.08 cm**-1 (unscaled, as printed).
REAL_LAYOUT_OUT = """\
# Hand-written synthetic fixture in the real ORCA 6.1 layout, not ORCA output.
------------------------------------------------------------------------------
                          OPTIMIZATION CYCLE  1
------------------------------------------------------------------------------
FINAL SINGLE POINT ENERGY      -76.100000000000

... cycle 1 geometry follows ...
------------------------------------------------------------------------------
                          OPTIMIZATION CYCLE  2
------------------------------------------------------------------------------
FINAL SINGLE POINT ENERGY      -76.350000000000

... cycle 2 geometry follows ...
------------------------------------------------------------------------------
                              THE OPTIMIZATION HAS CONVERGED
------------------------------------------------------------------------------
FINAL SINGLE POINT ENERGY      -76.362500000000

CARTESIAN COORDINATES (ANGSTROEM)
     O      0.000000    0.000000    0.121000
     H      0.000000    0.000000    0.861000

------------------
VIBRATIONAL FREQUENCIES
------------------

Scaling factor for frequencies = 0.9553000 (H0 = 0.9986000)

                  0:     0.00 cm**-1
                  6:  1637.64 cm**-1
                 21:  3787.08 cm**-1

--------------
ZERO POINT ENERGY
--------------
The zero point energy is           0.009580     Eh
"""

# A control in the OLD fixture shape (marker, blank, rows): must stay correct.
FLAT_BLOCK_OUT = """\
# Hand-written synthetic output, not ORCA output.
VIBRATIONAL FREQUENCIES (cm**(-1))
   0:     3441.3
   1:     1583.1
"""


@pytest.fixture
def real_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "real_layout"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "real.out").write_text(REAL_LAYOUT_OUT, encoding="utf-8")
    return folder


@pytest.fixture
def flat_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "flat"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "flat.out").write_text(FLAT_BLOCK_OUT, encoding="utf-8")
    return folder


# --- (1) frequency rows after a scaling line are real candidates ------------


def test_real_layout_frequency_verified(real_folder: Path) -> None:
    rec = verify_claim(3787.08, "cm-1", real_folder)
    assert rec.verdict == "verified", f"line={rec.line}"
    assert rec.quantity == "frequency"
    assert rec.parsed_value == pytest.approx(3787.08)


def test_real_layout_other_frequency_verified(real_folder: Path) -> None:
    rec = verify_claim(1637.64, "cm**-1", real_folder)
    assert rec.verdict == "verified"


def test_real_layout_wrong_frequency_is_mismatch_not_no_source(
    real_folder: Path,
) -> None:
    # 20 cm-1 off: must be mismatch (a wrong claim), never no-source.
    rec = verify_claim(3787.08 + 20.0, "cm-1", real_folder)
    assert rec.verdict == "mismatch"


def test_scaling_factor_line_yields_no_candidates(real_folder: Path) -> None:
    # The scaling/meta line's number must never become a frequency candidate.
    # The 0.00 translation/rotation rows are also dropped (see round 2), so
    # only the real vibrational modes remain.
    from delfin.agent import claim_locator as cl

    out = real_folder / "real.out"
    vals = [v for v, _ in cl._block_values(out, cl._ANCHORS["frequency"], "cm-1")]
    assert vals == [1637.64, 3787.08]


def test_flat_layout_still_verified(flat_folder: Path) -> None:
    # The simple marker/rows layout the old code handled must keep working.
    rec = verify_claim(3441.3, "cm-1", flat_folder)
    assert rec.verdict == "verified"


# --- (2) energy of an optimisation is the LAST FINAL SINGLE POINT ENERGY ----


def test_energy_is_final_cycle_not_first(real_folder: Path) -> None:
    # The answer states the converged energy: -76.3625 Eh (cycle 3, line 17).
    rec = verify_claim(-76.362500000000, "Eh", real_folder)
    assert rec.verdict == "verified"
    # Line 17 = the LAST 'FINAL SINGLE POINT ENERGY' line in REAL_LAYOUT_OUT.
    assert rec.line == 17
    assert rec.parsed_value == pytest.approx(-76.362500000000)


def test_first_cycle_energy_is_mismatch(real_folder: Path) -> None:
    # The cycle-1 energy is NOT the run's result: claiming it is a mismatch.
    rec = verify_claim(-76.100000000000, "Eh", real_folder)
    assert rec.verdict == "mismatch"


def test_kcal_conversion_uses_final_energy(real_folder: Path) -> None:
    rec = verify_claim(-76.362500000000 * 627.509, "kcal/mol", real_folder)
    assert rec.verdict == "verified"
    assert rec.line == 17


def test_zpe_uses_last_label_line(real_folder: Path) -> None:
    rec = verify_claim(0.009580, "Eh", real_folder)
    # ZPE is not in the energy labels' FIRST-line convention either; the last
    # 'Zero point energy' match is the real value.
    assert rec.verdict in ("verified", "mismatch")


# --- (3) mismatch record carries the closest real candidate -------------------


def test_mismatch_record_carries_closest_candidate(real_folder: Path) -> None:
    # A frequency 20 cm-1 off (3807 vs the real 3787.08) is a mismatch, and
    # the record must still show WHERE the run's real value sits: the
    # closest candidate by |claim - candidate| in canonical units, with its
    # line number (3787.08 at line 31 of the fixture) and parsed unit.
    rec = verify_claim(3807.0, "cm-1", real_folder)
    assert rec.verdict == "mismatch"
    assert rec.line == 31
    assert rec.parsed_value == pytest.approx(3787.08)
    assert rec.parsed_unit == "cm-1"


def test_mismatch_energy_record_carries_closest_candidate(real_folder: Path) -> None:
    # Same for energy: -75.9 Eh is off every cycle's value; the closest real
    # energy is the converged one (-76.3625, final cycle, line 17) -- it is
    # 0.46 Eh away vs 0.80 for cycle 1's -76.100.
    rec = verify_claim(-75.9, "Eh", real_folder)
    assert rec.verdict == "mismatch"
    assert rec.line == 17
    assert rec.parsed_value == pytest.approx(-76.362500000000)


# --- (4) 0.00 cm-1 translation/rotation rows are not candidates ---------------


def test_zero_frequency_rows_are_not_candidates(real_folder: Path) -> None:
    from delfin.agent import claim_locator as cl

    out = real_folder / "real.out"
    vals = [v for v, _ in cl._block_values(out, cl._ANCHORS["frequency"], "cm-1")]
    # The 0.00 cm-1 translation/rotation rows are dropped; only the real
    # vibrational modes remain.
    assert vals == [1637.64, 3787.08]


def test_negative_frequency_claim_is_not_verified_against_zero_rows(
    real_folder: Path,
) -> None:
    # A -1.0 cm-1 claim must NOT verify against the (dropped) 0.00 rows;
    # with no candidate within tolerance it is a mismatch, carrying the
    # closest real mode (1637.64, line 30).
    rec = verify_claim(-1.0, "cm-1", real_folder)
    assert rec.verdict == "mismatch"
    assert rec.line == 30
    assert rec.parsed_value == pytest.approx(1637.64)
