"""VR1-1 adversarial: a coordinate is not an interatomic distance.

Reviewer VR1-1 (nacht-s23): claim_locator treated every atomic coordinate on a
CARTESIAN COORDINATES line as a candidate "distance". A coordinate (x/y/z of
one atom) is NOT a distance -- an interatomic distance requires TWO atoms. On
the O-H fixture below the only true interatomic distance is
0.861000 - 0.121000 = 0.740000 Angstrom. Package target 4: a planted wrong
number must be caught, so a claim equal to a bare coordinate must be mismatch.
"""

from pathlib import Path

import pytest

from delfin.agent.claim_locator import verify_claim

# O at z=0.121000, H at z=0.861000 -> true O-H separation 0.740000 A.
SYNTHETIC_OUT = """\
# Synthetic ORCA-style output -- hand-written for tests, not produced by ORCA.
FINAL SINGLE POINT ENERGY      -76.026780443768

CARTESIAN COORDINATES (ANGSTROEM)
    O    0.000000    0.000000    0.121000
    H    0.000000    0.000000    0.861000
"""


@pytest.fixture
def calc_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "water"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "water.out").write_text(SYNTHETIC_OUT, encoding="utf-8")
    return folder


@pytest.mark.parametrize(
    "claimed,description",
    [
        (0.121, "O z-coordinate claimed as a distance (planted wrong number)"),
        (0.861, "H z-coordinate claimed as a distance (planted wrong number)"),
        (0.011, "nonexistent tiny 'distance' (planted wrong number)"),
    ],
)
def test_bare_coordinate_is_not_a_distance_is_mismatch(
    calc_folder: Path, claimed: float, description: str
) -> None:
    rec = verify_claim(claimed, "Angstrom", calc_folder)
    assert rec.verdict == "mismatch", description


def test_true_interatomic_distance_is_verified(calc_folder: Path) -> None:
    # The true O-H separation is the ONE real interatomic distance here.
    rec = verify_claim(0.740, "Angstrom", calc_folder)
    assert rec.verdict == "verified"
    assert rec.quantity == "distance"
    assert rec.parsed_value == pytest.approx(0.740)
