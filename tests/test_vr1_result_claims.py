"""VR1 phase 2: a reported number must be traceable to the real output FILE.

`result_claims.check` matched a number only against the evidence *text*; a
number the answer (or evidence excerpt) states but that is not actually in the
real output file still passed. Phase 2 closes that: when the evidence record's
``ref`` is a calculation folder, the number is verified with
``claim_locator.verify_claim`` against the folder's largest ``*.out``; a
mismatch or no-source number is flagged, not reported.
"""

from pathlib import Path

import pytest

from delfin.agent.claim_locator import verify_claim
from delfin.agent.result_claims import Claim, check

SYNTHETIC_OUT = (
    "# Synthetic ORCA-style output -- hand-written for tests.\n"
    "FINAL SINGLE POINT ENERGY      -76.026780443768\n"
    "\n"
    "Final Gibbs free energy                 -76.030514123   Eh\n"
    "\n"
    "Zero point energy                      0.003733680     Eh\n"
    "\n"
    "VIBRATIONAL FREQUENCIES (cm**(-1))\n"
    "   0:     3441.3\n"
    "   1:     1583.1\n"
    "   2:     1021.7\n"
    "\n"
    "CARTESIAN COORDINATES (ANGSTROEM)\n"
    "    O    0.000000    0.000000    0.121000\n"
    "    H    0.000000    0.000000    0.861000\n"
)


@pytest.fixture
def calc_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "water"
    folder.mkdir()
    (folder / "water.out").write_text(SYNTHETIC_OUT)
    return folder


def _evidence(folder: Path, text: str) -> list[dict]:
    return [{"kind": "calc", "ref": str(folder), "text": text, "flags": []}]


def _codes(claims: list[Claim]) -> list[str]:
    return [c.code for c in claims]


def test_correct_answer_passes_unchanged(calc_folder: Path) -> None:
    """A number that really is in the output file is not flagged."""
    rec = verify_claim(-76.026780443768, "Eh", calc_folder)  # sanity
    assert rec.verdict == "verified"
    text = "The total energy is -76.026780443768 Eh (B3LYP/def2-SVP)."
    claims = check(text, _evidence(calc_folder, text))
    assert claims == []


def test_planted_number_not_in_file_is_flagged(calc_folder: Path) -> None:
    """A value in the answer/evidence text but ABSENT from the real output
    file must no longer pass the text-only match."""
    # The evidence text claims a number that never appears in water.out.
    bogus = (-76.0, "Eh")
    text = f"The total energy is {bogus[0]} {bogus[1]} (B3LYP/def2-SVP)."
    claims = check(text, _evidence(calc_folder, f"energy {bogus[0]:g} {bogus[1]}"))
    assert any(c.code in ("claim-no-source", "claim-locator-mismatch",
                          "claim-locator-no-source") for c in claims), \
        f"expected the planted number to be flagged, got {_codes(claims)}"


def test_off_by_sign_is_flagged(calc_folder: Path) -> None:
    """A number in the text with the wrong SIGN is caught, not credited."""
    text = "The total energy is 76.026780443768 Eh (B3LYP/def2-SVP)."
    claims = check(text, _evidence(calc_folder, "energy -76.026780443768 Eh"))
    assert any(c.code in ("claim-no-source", "claim-locator-mismatch",
                          "claim-locator-no-source") for c in claims)


def test_off_by_unit_is_flagged(calc_folder: Path) -> None:
    """A number in the text whose UNIT differs is caught by conversion."""
    # 0.861 stated in bohr is 0.4556 Angstrom, not 0.861 Angstrom.
    text = "The O-H distance is 0.861 bohr (B3LYP/def2-SVP)."
    claims = check(text, _evidence(calc_folder, "distance 0.861 bohr"))
    assert any(c.code in ("claim-no-source", "claim-locator-mismatch",
                          "claim-locator-no-source") for c in claims)


def test_prose_numbers_do_not_false_alarm(calc_folder: Path) -> None:
    """Counts/years/settings must not be flagged (estimate-exempt)."""
    text = "the 2026 run wrote 3 files and used B3LYP/def2-SVP."
    claims = check(text, _evidence(calc_folder, text))
    assert not any("no-source" in c.code or "locator" in c.code for c in claims)
