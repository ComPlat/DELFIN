"""Every quantity falls back to the text parser when OPI is absent.

Runs under the CI interpreter (no OPI installed): :mod:`delfin.orca_read`
must import lazily and still answer, from DELFIN's own text readers, on
the same fixture runs the OPI tests use.
"""
from pathlib import Path

import pytest

from delfin.orca_read import Quantity, read

# The real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); every test here SKIPS without it.
from test_lc1_fixtures_helper import private_fixtures_or_skip  # noqa: E402

TOL = 1e-6


def fixture_copy(name: str, tmp_path: Path) -> Path:
    """Copy the fixture run's files from the private store into tmp."""
    src = private_fixtures_or_skip(f"tests/fixtures/orca611/{name}")
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        (dst / f.name).write_bytes(f.read_bytes())
    return dst


def test_opi_is_absent_here():
    """The test means anything only when OPI is really missing."""
    try:
        import opi  # noqa: F401
    except ImportError:
        return
    pytest.skip("OPI installed on this interpreter; this file pins the no-OPI path")


@pytest.mark.parametrize("name,final", [
    ("h2o_optfreq", -76.321269021507),
    ("o2_triplet", -150.051820514846),
    ("h2o_tddft", -76.041359207603),
])
def test_final_energy_from_text(name, final, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value == pytest.approx(final, abs=TOL)
    assert r.source == "text"


@pytest.mark.parametrize("name,g,zpe,h", [
    ("h2o_optfreq", -76.3177312, 0.02120238, -76.29628713),
    ("o2_triplet", None, None, None),
    ("h2o_tddft", None, None, None),
])
def test_thermochemistry_from_text(name, g, zpe, h, tmp_path):
    d = fixture_copy(name, tmp_path)
    for q, expect in ((Quantity.GIBBS, g), (Quantity.ZPE, zpe),
                      (Quantity.ENTHALPY, h)):
        r = read(d, name, q)
        if expect is None:
            assert r.value is None and r.source == "none", q
        else:
            assert r.value == pytest.approx(expect, abs=TOL)
            assert r.source == "text"


def test_inner_energy_from_text(tmp_path):
    d = fixture_copy("h2o_optfreq", tmp_path)
    r = read(d, "h2o_optfreq", Quantity.INNER_ENERGY)
    assert r.value == pytest.approx(-76.29723134, abs=TOL)
    assert r.source == "text"


def test_entropy_term_from_text(tmp_path):
    """DELFIN's text reader returns T*S (in Eh), not the molar entropy."""
    d = fixture_copy("h2o_optfreq", tmp_path)
    r = read(d, "h2o_optfreq", Quantity.ENTROPY)
    assert r.value == pytest.approx(0.02144407, abs=TOL)
    assert r.source == "text"
    assert r.unit == "Eh (T*S)"
    assert "T*S" in r.note or "naming" in r.note.lower()


def test_s2_deviation_from_text(tmp_path):
    d = fixture_copy("o2_triplet", tmp_path)
    r = read(d, "o2_triplet", Quantity.S2_DEVIATION)
    assert r.value == pytest.approx(0.007206, abs=TOL)
    assert r.source == "text"


def test_s2_deviation_none_on_closed_shell(tmp_path):
    for name in ("h2o_optfreq", "h2o_tddft"):
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.S2_DEVIATION)
        assert r.value is None
        assert r.source == "none"


def test_multiplicity_charge_nelectrons_from_input(tmp_path):
    cases = {
        "h2o_optfreq": (1, 0, 10),
        "o2_triplet": (3, 0, 16),
        "h2o_tddft": (1, 0, 10),
        "scf_not_converged": (3, 0, 16),
    }
    for name, (mult, charge, nelec) in cases.items():
        d = fixture_copy(name, tmp_path)
        assert read(d, name, Quantity.MULTIPLICITY).value == mult
        assert read(d, name, Quantity.CHARGE).value == charge
        assert read(d, name, Quantity.NUMBER_OF_ELECTRONS).value == nelec


def test_terminated_normally_and_failure_reason(tmp_path):
    for name, term in [("h2o_optfreq", True), ("o2_triplet", True),
                       ("h2o_tddft", True), ("scf_not_converged", False)]:
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.TERMINATED_NORMALLY)
        assert r.value is term
        assert r.source == "text"
    d = fixture_copy("scf_not_converged", tmp_path)
    r = read(d, "scf_not_converged", Quantity.FAILURE_REASON)
    assert r.source == "text"
    assert "error termination" in str(r.value).lower()


def test_geometry_convergence_none_for_non_opt_runs(tmp_path):
    for name, expect in [("h2o_optfreq", True), ("o2_triplet", None),
                         ("h2o_tddft", None), ("scf_not_converged", None)]:
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.GEOMETRY_CONVERGED)
        assert r.value is expect


def test_unconverged_run_never_raises(tmp_path):
    d = fixture_copy("scf_not_converged", tmp_path)
    for q in Quantity:
        r = read(d, "scf_not_converged", q)
        if q is Quantity.TERMINATED_NORMALLY:
            assert r.value is False
        elif q is Quantity.FAILURE_REASON:
            assert r.value
        elif q in (Quantity.MULTIPLICITY, Quantity.CHARGE,
                   Quantity.NUMBER_OF_ELECTRONS):
            # read from the .inp: the run's *intended* charge/multiplicity
            # stays reportable even when the run failed (phase-1 q04).
            assert r.source == "text"
            assert r.value is not None
        else:
            assert r.value is None, q


def test_missing_output_directory_never_raises(tmp_path):
    r = read(tmp_path / "does_not_exist", "nothing", Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value is None
    assert r.source == "none"
