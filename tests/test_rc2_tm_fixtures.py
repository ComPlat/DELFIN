"""RC2 TM fixtures: the reference files exist and pin the intended physics.

Red control: on the previous commit these fail because no
``tests/fixtures/orca611/rc2_*/`` fixture exists yet. Each fixture is
accepted only after its run satisfies the physics gates (converged SCF +
geometry, 0 imaginary modes, intended spin).

Uses only Quantity members that exist in delfin/orca_read.py:67-95:
MULTIPLICITY, CHARGE, S2_DEVIATION, GEOMETRY_CONVERGED,
TERMINATED_NORMALLY, FINAL_SINGLE_POINT_ENERGY. Imaginary modes come from
delfin.imag.collect_imaginary_modes (imag.py:166), the existing reader for
ORCA's "***imaginary mode***" lines. Reviewer F1 (nacht-s4) caught the
phantom members of the first draft.
"""
from pathlib import Path

import pytest

from delfin.imag import collect_imaginary_modes
from delfin.orca_read import Quantity, read

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir


def _fixtures(name: str):
    """The specific orca611 subfolder: legacy path while its files are
    committed, else the private store (skip without the env var)."""
    return lc1_fixture_dir(f"tests/fixtures/orca611/{name}")

# name -> (charge, multiplicity)
SYSTEMS = {
    "rc2_fe_co5": (0, 1),
    "rc2_cu_nh3_4": (2, 2),
    "rc2_ru_o4": (0, 1),
}


def fixture_copy(name: str, tmp_path: Path) -> Path:
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for src in _fixtures(name).iterdir():
        (dst / src.name).write_bytes(src.read_bytes())
    return dst


def test_there_are_three_rc2_fixtures():
    for name in SYSTEMS:
        d = _fixtures(name)
        assert d.is_dir(), name
        assert (d / f"{name}.inp").is_file(), name
        assert (d / f"{name}.out").is_file(), name
        # ORCA 6.1.1 manual section 9.4.2: %output JSONPropFile True writes
        # basename.property.txt AND basename.property.json (a bare
        # basename.json comes only from the separate orca_2json OPI route,
        # which the rc2 practice runs did not use).
        assert (d / f"{name}.property.json").is_file(), name
        assert (d / "PROVENANCE.md").is_file(), name


@pytest.mark.parametrize("name,charge,mult", sorted((k, *v) for k, v in SYSTEMS.items()))
def test_charge_and_multiplicity_pinned(name, charge, mult, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.MULTIPLICITY)
    assert r.value == mult, (name, r.value)
    c = read(d, name, Quantity.CHARGE)
    assert c.value == charge, (name, c.value)


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_zero_imaginary_modes(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    modes = collect_imaginary_modes(str(d / f"{name}.out"))
    assert modes == [], (name, modes)


def test_cu_doublet_s2_near_ideal(tmp_path):
    """Cu(II) d9: pure doublet has <S^2> = 0.75; deviations are reported."""
    d = fixture_copy("rc2_cu_nh3_4", tmp_path)
    dev = read(d, "rc2_cu_nh3_4", Quantity.S2_DEVIATION)
    assert dev.value is not None, "S2_DEVIATION must be readable for a doublet"
    # <S^2> = 0.75 + deviation (orca_read.py:82-83 defines the deviation as
    # expected - ideal); a pure doublet has 0.75.
    s2 = 0.75 + dev.value
    assert s2 == pytest.approx(0.75, abs=0.05), f"<S^2> = {s2}"
    assert dev.value < 0.05, f"spin contamination: deviation {dev.value}"


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_terminated_normally(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.TERMINATED_NORMALLY)
    assert r.value is True, (name, r.value)


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_geometry_converged(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.GEOMETRY_CONVERGED)
    assert r.value is True, (name, r.value)


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_final_energy_present(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value is not None, name


def test_no_host_identity_in_rc2_fixtures():
    for name in SYSTEMS:
        for f in sorted(_fixtures(name).rglob("*")):
            if not f.is_file() or f.suffix == ".gbw":
                continue
            text = f.read_text(errors="replace")
            assert "/home/" not in text, f.name
            for line in text.splitlines():
                if "Host name:" in line:
                    assert line.split("Host name:", 1)[1].strip() == "fixture-host", f.name
