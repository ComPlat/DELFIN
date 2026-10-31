"""OPI (orca-pi) equivalence for spin-derived reads: <S**2>, multiplicity,
charge, number of electrons.

DELFIN currently reads spin from ORCA's text output (the ``UHF SPIN
CONTAMINATION`` block via ``delfin.parser.extract_last_uhf_deviation``); charge,
multiplicity and electron count are read from the *input* (``.inp`` file header
and the CONTROL/XYZ geometry). OPI reads the same quantities from the output
(the ``<S**2>``/``S*(S+1)`` lines in the ``.out`` and the ``<base>.property.json``
next to it).

This test asserts that, on the four ORCA 6.1.1 fixture runs in
``tests/fixtures/orca611/``, DELFIN's text/input reads and OPI's reads agree
within tolerance, and that both behave sensibly on the *not-converged* run.
Where the two differ by construction, the test states which is right and why
(measure, not guess -- the numbers below were read from the fixture files).

Quantity relationships (measured on ``o2_triplet``, ORCA 6.1.1, UKS triplet
O2)::

    ORCA .out UHF SPIN CONTAMINATION block
      Expectation value of <S**2>     :     2.007206
      Ideal value S*(S+1) for S=1.0   :     2.000000
      Deviation                       :     0.007206

DELFIN's ``extract_last_uhf_deviation`` returns the ``Deviation`` line
(``<S**2> - S(S+1)``) = 0.007206. OPI's ``Output.get_s2()`` returns the tuple
``(expected, ideal)`` = (2.007206, 2.0).  Hence

    DELFIN_deviation  ==  OPI.expected - OPI.ideal

within the printed precision (the .out prints 6 decimals). ``get_charge()`` /
``get_mult()`` / ``get_nelectrons()`` read ``calculation_info`` from the
property JSON -- the launch charge/multiplicity, i.e. the same values DELFIN
reads from the ``.inp`` header.

On the closed-shell runs (h2o_optfreq, h2o_tddft) ORCA prints no UHF
contamination block: DELFIN returns ``None`` (treated as deviation 0.0 for a
singlet by OCCUPIER) and OPI's ``get_s2()`` also returns ``None``.

On the not-converged run neither parser produces values: ``get_s2()`` is None
(SCF never converged, no block) and ``get_charge()``/``get_mult()`` are None
because ORCA wrote no ``calculation_info`` to the property JSON for a failed
run. Only ``scf_converged()``/``terminated_normally()`` distinguish "no value
because closed shell" from "no value because the run failed" -- a reader must
check those before trusting a value.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.parser import extract_last_uhf_deviation

try:
    from opi.output.core import Output

    _OPI = True
except Exception:  # pragma: no cover - environment dependent
    _OPI = False


# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir as _lc1_fixture_dir  # noqa: E402, F401

#: The .out prints <S**2> to 6 decimals; anything closer is the same number.
S2_TOL = 1e-6

# (dir, basename, expected multiplicity, expected charge, expected nelectrons)
CONVERGED_CASES = [
    # h2o_optfreq :  H2O closed-shell singlet opt+freq, converged
    ("h2o_optfreq", 1, 0, 10),
    # o2_triplet  :  O2 UKS open-shell triplet, single point, converged
    ("o2_triplet", 3, 0, 16),
    # h2o_tddft   :  H2O closed-shell singlet TD-DFT, converged
    ("h2o_tddft", 1, 0, 10),
]

NOT_CONVERGED = "scf_not_converged"


@pytest.mark.skipif(not _OPI, reason="OPI (delfin[orca]) not installed")
@pytest.mark.parametrize("name,mult,charge,nelec", CONVERGED_CASES)
def test_opi_matches_on_converged_runs(name, mult, charge, nelec):
    """On converged runs DELFIN's reads equal OPI's within tolerance."""
    d = _lc1_fixture_dir(f"tests/fixtures/orca611/{name}")
    out = Output(name, working_dir=d, version_check=False)
    out.parse(read_gbw_json=False)

    assert out.scf_converged() is True
    assert out.terminated_normally() is True

    # multiplicity -- OPI from property.json, DELFIN from the .inp header
    assert out.get_mult() == mult
    assert _delfin_mult_from_inp(d / f"{name}.inp") == mult

    # charge -- OPI from property.json, DELFIN from the .inp header
    assert out.get_charge() == charge
    assert _delfin_charge_from_inp(d / f"{name}.inp") == charge

    # number of electrons -- OPI from property.json; DELFIN derives it from
    # the .inp charge + the atom list (see occupier .calculate_total_electrons)
    assert out.get_nelectrons()[0] == nelec
    assert _delfin_nelectrons_from_xyz(d / f"{name}.inp", charge) == nelec


@pytest.mark.skipif(not _OPI, reason="OPI (delfin[orca]) not installed")
def test_s2_deviation_agrees_on_open_shell():
    """DELFIN's Deviation equals OPI's (expected - ideal) on the UKS triplet."""
    d = _lc1_fixture_dir("tests/fixtures/orca611/o2_triplet")
    out = Output("o2_triplet", working_dir=d, version_check=False)
    out.parse(read_gbw_json=False)

    delfin_dev = extract_last_uhf_deviation(str(d / "o2_triplet.out"))
    opi_s2 = out.get_s2()

    assert opi_s2 is not None
    expected, ideal = opi_s2
    assert abs(delfin_dev - (expected - ideal)) <= S2_TOL
    # the fixture itself: expectation 2.007206, ideal 2.0, dev 0.007206
    assert delfin_dev == pytest.approx(0.007206, abs=S2_TOL)
    assert expected == pytest.approx(2.007206, abs=S2_TOL)


@pytest.mark.skipif(not _OPI, reason="OPI (delfin[orca]) not installed")
@pytest.mark.parametrize("name", ["h2o_optfreq", "h2o_tddft"])
def test_s2_none_on_closed_shell(name):
    """Closed-shell runs have no UHF block: DELFIN None, OPI None -- the same."""
    d = _lc1_fixture_dir(f"tests/fixtures/orca611/{name}")
    delfin_dev = extract_last_uhf_deviation(str(d / f"{name}.out"))
    out = Output(name, working_dir=d, version_check=False)
    out.parse(read_gbw_json=False)

    assert delfin_dev is None
    # OPI returns no <S**2> tuple either -- closed shell, nothing to report.
    assert out.get_s2() is None


@pytest.mark.skipif(not _OPI, reason="OPI (delfin[orca]) not installed")
def test_not_converged_behaves_sensibly():
    """On a failed run both parsers return no values; only the flags say why."""
    d = _lc1_fixture_dir(f"tests/fixtures/orca611/{NOT_CONVERGED}")
    out = Output(NOT_CONVERGED, working_dir=d, version_check=False)
    out.parse(read_gbw_json=False)

    # the failure is signalled ONLY by these flags -- never by a silent value
    assert out.scf_converged() is False
    assert out.terminated_normally() is False
    assert out.geometry_optimization_converged() is False

    # no values emerge from a run that never converged
    assert out.get_s2() is None
    assert out.get_mult() is None
    assert out.get_charge() is None
    assert out.get_nelectrons() == (None, None)

    # DELFIN's text parser agrees: no UHF block in the failed .out
    delfin_dev = extract_last_uhf_deviation(str(d / f"{NOT_CONVERGED}.out"))
    assert delfin_dev is None

    # The .inp still carries the *intended* charge/multiplicity (input, not
    # output) -- DELFIN, unlike OPI, can still report what was requested.
    assert _delfin_charge_from_inp(d / f"{NOT_CONVERGED}.inp") == 0
    assert _delfin_mult_from_inp(d / f"{NOT_CONVERGED}.inp") == 3


# ---------------------------------------------------------------------------
# DELFIN's input-side readers (what it compares against OPI's output reads)
# ---------------------------------------------------------------------------


def _delfin_charge_from_inp(inp: Path):
    """DELFIN reads charge from the ``* xyz <charge> <mult>`` input header."""
    m = _coord_header(inp)
    return int(m.group(1)) if m else None


def _delfin_mult_from_inp(inp: Path):
    m = _coord_header(inp)
    return int(m.group(2)) if m else None


def _coord_header(inp: Path):
    import re

    text = inp.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"\*\s*(?:xyz|xyzfile)\s+(-?\d+)\s+(\d+)", text)
    return m


def _delfin_nelectrons_from_xyz(inp: Path, charge: int):
    """Total electrons = sum of atomic numbers minus charge (like OCCUPIER)."""
    import re

    _AT = {
        "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "NE": 10,
    }
    text = inp.read_text(encoding="utf-8", errors="replace")
    coord = re.search(r"\*\s*xyz\s+-?\d+\s+\d+\s*\n(.*?)\*\s*", text, re.DOTALL)
    if not coord:
        return None
    atoms = [ln.split()[0].upper() for ln in coord.group(1).splitlines() if ln.split()]
    total = sum(_AT.get(a, 0) for a in atoms)
    return total - charge
