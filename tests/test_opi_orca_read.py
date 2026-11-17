"""OPI-first reading through :mod:`delfin.orca_read` (requires OPI).

Every test works on a copy of the fixture run in ``tmp_path`` -- never on
``tests/fixtures/`` itself, because OPI's ``parse()`` may write JSON files
next to the output it reads.
"""
from pathlib import Path

import pytest

pytest.importorskip("opi", reason="OPI (orca-pi) not installed")

from delfin.orca_read import Quantity, read

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir  # noqa: E402

#: The .out prints energies/s2 to 8 resp. 6 decimals; OPI's JSON to full
#: double precision. The phase-1 equivalence tests used 1e-6 (S2) and the
#: text's printed precision (energies); the shadow check uses the same.
TOL = 1e-6


def fixture_copy(name: str, tmp_path: Path) -> Path:
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for src in lc1_fixture_dir(f"tests/fixtures/orca611/{name}").iterdir():
        (dst / src.name).write_bytes(src.read_bytes())
    return dst


def value(work_dir: Path, name: str, quantity: Quantity):
    r = read(work_dir, name, quantity)
    assert isinstance(r.value, (int, float)) or r.value is None
    return r


@pytest.mark.parametrize(
    "name,expect",
    [
        # (final E, G, ZPE, H, U) in Eh, from the fixture .out/property.json
        ("h2o_optfreq", (-76.32126902150691, -76.3177312027805,
                         0.021202376295716842, -76.29628713260958,
                         -76.29723134165184)),
        ("o2_triplet", (-150.05182051484584, None, None, None, None)),
        ("h2o_tddft", (-76.04135920760325, None, None, None, None)),
    ],
)
def test_energies_opi_first_and_equal_to_text(name, expect, tmp_path):
    """Energies route through OPI; value equals the text value within TOL."""
    d = fixture_copy(name, tmp_path)
    for quantity, expected in zip(
        (Quantity.FINAL_SINGLE_POINT_ENERGY, Quantity.GIBBS,
         Quantity.ZPE, Quantity.ENTHALPY, Quantity.INNER_ENERGY), expect):
        r = read(d, name, quantity)
        if expected is None:
            assert r.value is None
            assert r.source == "none"
            continue
        assert r.value == pytest.approx(expected, abs=TOL), quantity
        assert r.source == "opi", quantity
        assert r.unit == "Eh"


def test_s2_deviation_on_open_shell(tmp_path):
    """<S**2> deviation: OPI's (expected - ideal) == text's Deviation value."""
    d = fixture_copy("o2_triplet", tmp_path)
    r = read(d, "o2_triplet", Quantity.S2_DEVIATION)
    assert r.value == pytest.approx(0.007206, abs=TOL)
    assert r.source == "opi"


def test_s2_none_on_closed_shell(tmp_path):
    for name in ("h2o_optfreq", "h2o_tddft"):
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.S2_DEVIATION)
        assert r.value == 0.0
        assert r.source == "opi"


def test_multiplicity_charge_nelectrons(tmp_path):
    cases = {
        "h2o_optfreq": (1, 0, 10),
        "o2_triplet": (3, 0, 16),
        "h2o_tddft": (1, 0, 10),
    }
    for name, (mult, charge, nelec) in cases.items():
        d = fixture_copy(name, tmp_path)
        assert read(d, name, Quantity.MULTIPLICITY).value == mult
        assert read(d, name, Quantity.CHARGE).value == charge
        assert read(d, name, Quantity.NUMBER_OF_ELECTRONS).value == nelec


def test_terminated_normally(tmp_path):
    for name, expect in [("h2o_optfreq", True), ("o2_triplet", True),
                         ("h2o_tddft", True), ("scf_not_converged", False)]:
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.TERMINATED_NORMALLY)
        assert r.value is expect
        assert r.unit == ""


def test_unconverged_run_gives_none(tmp_path):
    """A failed run yields no values; flags stay honest, no exception."""
    d = fixture_copy("scf_not_converged", tmp_path)
    for q in Quantity:
        r = read(d, "scf_not_converged", q)
        if q is Quantity.TERMINATED_NORMALLY:
            assert r.value is False
        elif q is Quantity.FAILURE_REASON:
            assert r.source == "text" and r.value
        elif q in (Quantity.MULTIPLICITY, Quantity.CHARGE,
                   Quantity.NUMBER_OF_ELECTRONS):
            # input-side reads survive a failed run (intended values)
            assert r.source == "text" and r.value is not None
        else:
            assert r.value is None, q
            assert r.source in ("opi", "none"), q


def test_failure_reason_is_text(tmp_path):
    """Failure reason: OPI says THAT a run failed, not WHY -- text reader."""
    d = fixture_copy("scf_not_converged", tmp_path)
    r = read(d, "scf_not_converged", Quantity.FAILURE_REASON)
    assert r.source == "text"
    assert "error termination" in str(r.value).lower()


def test_geometry_convergence_semantics(tmp_path):
    """DELFIN keeps None (no optimization) where OPI answers False."""
    for name, kind in [("h2o_optfreq", "opt"), ("o2_triplet", "sp"),
                       ("scf_not_converged", "failed")]:
        d = fixture_copy(name, tmp_path)
        r = read(d, name, Quantity.GEOMETRY_CONVERGED)
        if kind == "opt":
            assert r.value is True
        else:
            assert r.value is None


def test_missing_property_files_fall_back_to_text(tmp_path):
    """OPI needs its JSON files; without them the text parser reads."""
    d = fixture_copy("o2_triplet", tmp_path)
    (d / "o2_triplet.property.json").unlink()
    (d / "o2_triplet.property.txt").unlink()
    r = read(d, "o2_triplet", Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value == pytest.approx(-150.05182051484584, abs=TOL)
    assert r.source == "text"
    assert "property" in r.note.lower()


def test_orca_60_header_falls_back_to_text(tmp_path):
    """ORCA 6.0.x output: OPI is out of its contract -- text reader."""
    d = fixture_copy("o2_triplet", tmp_path)
    out = d / "o2_triplet.out"
    out.write_text(out.read_text(errors="replace")
                   .replace("Program Version 6.1.1", "Program Version 6.0.1"))
    # property JSON still says 6.1.1 -- a 6.0.1 run converted silently keeps
    # its energies; here the .out header is the gate, so OPI must not read.
    r = read(d, "o2_triplet", Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.source == "text"
    assert r.value == pytest.approx(-150.05182051484584, abs=TOL)


def test_forced_mismatch_returns_text_value_and_warns(tmp_path, monkeypatch, caplog):
    """Shadow check: an OPI value disagreeing with the text is rejected."""
    import logging

    d = fixture_copy("o2_triplet", tmp_path)
    monkeypatch.setattr("delfin.orca_read._opi_values",
                        lambda *a, **k: {"final_energy": -123.0})
    with caplog.at_level(logging.WARNING, logger="delfin.orca_read"):
        r = read(d, "o2_triplet", Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value == pytest.approx(-150.05182051484584, abs=TOL)
    assert r.source == "text"
    assert "o2_triplet" in caplog.text
    assert "-123.0" in caplog.text
    # no absolute paths in the message
    for line in caplog.text.splitlines():
        if "orca_read" in line or "mismatch" in line.lower():
            assert "/home/" not in line


def test_setting_orca_reader_text_never_uses_opi(tmp_path, monkeypatch):
    """``orca.reader = "text"``: every quantity comes from the text parser."""
    from delfin import orca_read
    d = fixture_copy("h2o_optfreq", tmp_path)
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    for q in (Quantity.FINAL_SINGLE_POINT_ENERGY, Quantity.GIBBS, Quantity.ZPE,
              Quantity.MULTIPLICITY, Quantity.CHARGE):
        r = read(d, "h2o_optfreq", q)
        assert r.source == "text", q
        assert r.value is not None, q


def test_setting_orca_reader_text_on_failed_run(tmp_path, monkeypatch):
    from delfin import orca_read
    d = fixture_copy("scf_not_converged", tmp_path)
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    r = read(d, "scf_not_converged", Quantity.TERMINATED_NORMALLY)
    assert r.value is False
    assert r.source == "text"
