"""RC2 TM fixtures: the reference files exist and pin the intended physics.

Red control: on the previous commit these fail because no
``tests/fixtures/orca611/rc2_*/`` fixture exists yet. Each fixture is
accepted only after its run satisfies the physics gates (converged SCF +
geometry, 0 imaginary modes, intended spin).
"""
import re
from pathlib import Path

import pytest

from delfin.orca_read import Quantity, read

FIXTURES = Path(__file__).parent / "fixtures" / "orca611"

# name -> (charge, multiplicity)
SYSTEMS = {
    "rc2_fe_co5": (0, 1),
    "rc2_cu_nh3_4": (2, 2),
    "rc2_ru_o4": (0, 1),
}


def fixture_copy(name: str, tmp_path: Path) -> Path:
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for src in (FIXTURES / name).iterdir():
        (dst / src.name).write_bytes(src.read_bytes())
    return dst


def test_there_are_three_rc2_fixtures():
    for name in SYSTEMS:
        d = FIXTURES / name
        assert d.is_dir(), name
        assert (d / f"{name}.inp").is_file(), name
        assert (d / f"{name}.out").is_file(), name
        assert (d / f"{name}.json").is_file(), name
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
    r = read(d, name, Quantity.IMAGINARY_FREQ_COUNT)
    assert r.value == 0, (name, r.value)


def test_cu_doublet_s2_near_ideal(tmp_path):
    """Cu(II) d9: pure doublet has <S^2> = 0.75; deviations are reported."""
    d = fixture_copy("rc2_cu_nh3_4", tmp_path)
    r = read(d, "rc2_cu_nh3_4", Quantity.S2_EXPECTED)
    dev = read(d, "rc2_cu_nh3_4", Quantity.S2_DEVIATION)
    assert r.value is not None
    assert r.value == pytest.approx(0.75, abs=0.01)
    assert dev.value is not None
    # heavy spin contamination -> flag, do not hide
    assert dev.value < 0.05, f"spin contamination: deviation {dev.value}"


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_scf_and_geometry_converged(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    scf = read(d, name, Quantity.SCF_CONVERGED)
    opt = read(d, name, Quantity.OPT_CONVERGED)
    assert scf.value is True, (name, scf.value)
    assert opt.value is True, (name, opt.value)


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_final_energy_present(name, tmp_path):
    d = fixture_copy(name, tmp_path)
    r = read(d, name, Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert r.value is not None, name


def test_no_host_identity_in_rc2_fixtures():
    for name in SYSTEMS:
        for f in sorted((FIXTURES / name).rglob("*")):
            if not f.is_file() or f.suffix == ".gbw":
                continue
            text = f.read_text(errors="replace")
            assert "/home/" not in text, f.name
            for line in text.splitlines():
                if "Host name:" in line:
                    assert line.split("Host name:", 1)[1].strip() == "fixture-host", f.name
