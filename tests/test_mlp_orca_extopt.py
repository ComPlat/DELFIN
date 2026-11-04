"""Tests for the ORCA ExtOpt/MLP wrapper (delfin/mlp_tools/orca_extopt.py).

Every file here is synthetic (first line marker "# SYNTHETIC test fixture
written for DELFIN -- not ORCA output" where a file fixture is written);
the calculator is a fake ASE calculator, so this suite runs with no torch,
no GPU and no ORCA.  Physical conventions verified against the ORCA 6.1.1
manual, chapter 4.1.10 "Optimizing with External Methods" (manual text
structurereactivity/optimizations.txt lines 1537-1701).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from conftest import child_env

from delfin.mlp_tools.orca_extopt import (
    EV_TO_EH,
    EV_PER_ANG_TO_EH_PER_BOHR,
    parse_extinp,
    run_wrapper,
    write_engrad,
)

SYN = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"

#: CODATA 2018 anchor as an independent literal (not derived from the
#: module's factor, which used to hide a reciprocal-factor bug).
ANG_PER_BOHR = 0.529177210903


def test_unit_factors_are_consistent():
    # CODATA: 1 Eh = 27.211386245988 eV; 1 Bohr = 0.529177210903 Angstrom.
    assert EV_TO_EH == pytest.approx(1.0 / 27.211386245988, rel=1e-10)
    # Eh/Bohr per eV/Angstrom = Bohr/Eh = 0.0194469..., NOT ~51.42.
    assert EV_PER_ANG_TO_EH_PER_BOHR == pytest.approx(
        ANG_PER_BOHR / 27.211386245988, rel=1e-10
    )


class FakeCalculator:
    """Fake ASE calculator: E = sum(x^2) (eV), F = -2x (eV/Angstrom)."""

    def get_potential_energy(self, atoms=None):
        positions = atoms.positions
        return float(np.sum(positions**2))

    def get_forces(self, atoms=None):
        positions = atoms.positions
        return -2.0 * positions


def _write_xyz(tmp_path: Path, name: str) -> Path:
    # ORCA-written xyz files carry no comment marker; the "# SYNTHETIC"
    # marker rule applies to committed fixture files, not to files a test
    # generates at runtime in tmp.
    xyz = tmp_path / name
    xyz.write_text(
        "3\nsynthetic water, positions in Angstrom\n"
        + "O 0.0 0.0 0.0\n"
        + "H 0.757 0.586 0.0\n"
        + "H -0.757 0.586 0.0\n"
    )
    return xyz


def test_parse_extinp_reads_orca_written_file(tmp_path):
    """Format per manual 4.1.10: xyz name, charge, mult, ncores, do_gradient,
    optional pointcharge file; '#' comments ignored."""
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(
        SYN
        + "test_EXT.xyz # xyz filename\n"
        + "0 # charge\n"
        + "1 # multiplicity\n"
        + "2 # NCores\n"
        + "1 # do gradient\n"
        + "pointcharges.pc # optional point charges\n"
    )
    info = parse_extinp(extinp)
    assert info.xyz_filename == "test_EXT.xyz"
    assert info.charge == 0
    assert info.multiplicity == 1
    assert info.ncores == 2
    assert info.do_gradient is True
    assert info.pointcharge_filename == "pointcharges.pc"


def test_parse_extinp_without_gradient_and_pointcharges(tmp_path):
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(SYN + "test_EXT.xyz\n0\n1\n1\n0\n")
    info = parse_extinp(extinp)
    assert info.do_gradient is False
    assert info.pointcharge_filename is None


def test_parse_extinp_rejects_short_file(tmp_path):
    extinp = tmp_path / "bad.extinp"
    extinp.write_text(SYN + "only_one_line\n")
    with pytest.raises(ValueError):
        parse_extinp(extinp)


def test_write_engrad_format_matches_manual(tmp_path):
    """engrad file: natoms, energy in Eh, gradient in Eh/bohr flattened
    atom-major (Atom1X, Atom1Y, Atom1Z, Atom2X, ...), '#' comment blocks."""
    engrad = tmp_path / "test_EXT.engrad"
    energy_eh = -5.504066223730
    grad = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=float)
    write_engrad(engrad, n_atoms=2, energy_eh=energy_eh, gradient_eh_bohr=grad)
    text = engrad.read_text()
    # Comments present, values in Eh
    assert "# Number of atoms" in text
    assert "# The current total energy in Eh" in text
    assert "# The current gradient in Eh/bohr" in text
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert lines[0] == "2"
    assert float(lines[1]) == pytest.approx(energy_eh)
    values = [float(v) for v in lines[2:]]
    assert values == pytest.approx(list(grad))


def test_engrad_unit_conversion(tmp_path):
    """An energy given in eV (ASE convention) must land as Eh, and forces
    in eV/Angstrom as Eh/bohr, with the documented conversion factors."""
    energy_ev = -1.5
    forces_ev_ang = np.array([[0.1, -0.2, 0.3], [0.0, 0.0, 0.0]])
    write_engrad(
        tmp_path / "g.engrad",
        n_atoms=2,
        energy_eh=energy_ev * EV_TO_EH,
        gradient_eh_bohr=forces_ev_ang.reshape(-1) * EV_PER_ANG_TO_EH_PER_BOHR,
    )
    lines = [
        ln
        for ln in (tmp_path / "g.engrad").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    assert float(lines[1]) == pytest.approx(-1.5 * EV_TO_EH)
    assert float(lines[2]) == pytest.approx(0.1 * EV_PER_ANG_TO_EH_PER_BOHR)


class FakeMLPFactory:
    """Stands in for delfin.mlp_tools.calculators.create_calculator."""

    def __init__(self):
        self.calls = []

    def __call__(self, backend, *, device="cpu", charge=0, mult=1, **kw):
        self.calls.append((backend, device, charge, mult))
        return FakeCalculator()


def test_run_wrapper_full_round_trip(tmp_path, monkeypatch):
    """ORCA calls: script basename_EXT.extinp.tmp.  The wrapper reads the
    xyz, evaluates the fake calculator, writes basename_EXT.engrad."""
    xyz = _write_xyz(tmp_path, "test_EXT.xyz")
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(SYN + "test_EXT.xyz\n0\n1\n1\n1\n")

    fake = FakeMLPFactory()
    monkeypatch.setattr(
        "delfin.mlp_tools.orca_extopt.create_calculator", fake
    )
    run_wrapper(extinp, backend="ani2x", calculator_factory=fake)
    engrad = tmp_path / "test_EXT.engrad"
    assert engrad.exists()
    lines = [ln for ln in engrad.read_text().splitlines() if ln and not ln.startswith("#")]
    assert lines[0] == "3"
    # E = 0^2 + (0.757^2+0.586^2) * 2 (eV) -> Eh
    e_ev = 0.757**2 + 0.586**2 + (0.757**2 + 0.586**2)
    assert float(lines[1]) == pytest.approx(e_ev * EV_TO_EH, rel=1e-9)
    # gradient: dE/dx = -force = +2x per atom, atom-major, in Eh/bohr
    grad = [float(v) for v in lines[2:]]
    assert len(grad) == 9
    assert grad[0] == pytest.approx(0.0)
    assert grad[3] == pytest.approx(2.0 * 0.757 * EV_PER_ANG_TO_EH_PER_BOHR)
    # backend and thread limits reach create_calculator
    assert fake.calls == [("ani2x", "cpu", 0, 1)]


def test_run_wrapper_energy_only(tmp_path, monkeypatch):
    """do_gradient = 0: only energy is requested (manual 4.1.10: the
    gradient is only calculated when needed)."""
    _write_xyz(tmp_path, "test_EXT.xyz")
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(SYN + "test_EXT.xyz\n0\n1\n1\n0\n")
    fake = FakeMLPFactory()
    monkeypatch.setattr("delfin.mlp_tools.orca_extopt.create_calculator", fake)
    run_wrapper(extinp, backend="mace_off", calculator_factory=fake)
    lines = [
        ln
        for ln in (tmp_path / "test_EXT.engrad").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    assert lines[0] == "3"
    assert len(lines) == 2  # natoms, energy -- no gradient block


def test_run_wrapper_missing_xyz_fails_clearly(tmp_path, monkeypatch):
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(SYN + "test_EXT.xyz\n0\n1\n1\n1\n")
    monkeypatch.setattr(
        "delfin.mlp_tools.orca_extopt.create_calculator", FakeMLPFactory()
    )
    with pytest.raises(FileNotFoundError):
        run_wrapper(extinp, backend="ani2x")


def test_cli_entrypoint(tmp_path):
    """The wrapper must be callable as a command line script, exactly the
    way ORCA calls it: script <extinp> [args].  The fake calculator is
    injected inside the child process (no torch in this env)."""
    _write_xyz(tmp_path, "test_EXT.xyz")
    extinp = tmp_path / "test_EXT.extinp.tmp"
    extinp.write_text(SYN + "test_EXT.xyz\n0\n1\n1\n0\n")
    helper = tmp_path / "_fakecalc_cli.py"
    helper.write_text(
        "import numpy as np\n"
        "class FakeCalculator:\n"
        "    def get_potential_energy(self, atoms=None):\n"
        "        return float(np.sum(atoms.positions**2))\n"
        "    def get_forces(self, atoms=None):\n"
        "        return -2.0 * atoms.positions\n"
    )
    child = (
        "import sys\n"
        f"sys.path.insert(0, {str(tmp_path)!r})\n"
        "from _fakecalc_cli import FakeCalculator\n"
        f"import delfin.mlp_tools.orca_extopt as m\n"
        "m.create_calculator = lambda backend, **kw: FakeCalculator()\n"
        f"sys.exit(m.main([{str(extinp)!r}, '--backend', 'fake', '--opt', 'ignored']))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", child],
        capture_output=True,
        text=True,
        env=child_env(tmp_path),
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "test_EXT.engrad").exists()
