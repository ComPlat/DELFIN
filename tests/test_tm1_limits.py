"""Builder tests for package TM1 phase 1 (nacht-s20).

Physical limits and the reviewer findings, driven through the real
public path (registered adapter execute()) with recording fakes on PATH.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get


@pytest.fixture
def tm_fakes(tmp_path, monkeypatch):
    """Recording fakes for jobex AND aoforce on PATH.

    aoforce exits 0 and writes a minimal hessian file, like the real
    program; jobex exits 1 so the adapter reports FAILED without a real
    Turbomole.
    """
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    for name in ("ridft", "dscf", "jobex", "aoforce"):
        fake = bin_dir / name
        fake.write_text("#!/bin/sh\nexit 0\n")
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    (bin_dir / "aoforce").write_text(
        f"#!/bin/sh\necho '$hessian' > {work_dir}/hessian\nexit 0\n")
    for f in (bin_dir / "aoforce",):
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return work_dir


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


def _run(name, work_dir, **kwargs):
    adapter = get(name)
    assert adapter is not None, f"adapter {name} not registered"
    return adapter.execute(work_dir, geometry=_xyz(work_dir), cores=1, **kwargs)


class TestChargeMultLimits:
    """Physical limits are enforced, not described."""

    def test_non_integer_charge_rejected(self, tm_fakes):
        result = _run("turbomole_sp", tm_fakes, charge="0")
        assert result.status.name == "FAILED"
        assert "charge must be an integer" in (result.error or "")

    def test_mult_zero_rejected(self, tm_fakes):
        result = _run("turbomole_sp", tm_fakes, charge=0, mult=0)
        assert result.status.name == "FAILED"
        assert "mult must be >= 1" in (result.error or "")

    def test_negative_mult_rejected(self, tm_fakes):
        result = _run("turbomole_opt", tm_fakes, charge=0, mult=-2)
        assert result.status.name == "FAILED"

    def test_valid_neutral_accepted(self, tm_fakes):
        result = _run("turbomole_sp", tm_fakes, charge=0, mult=1)
        assert result.status.name == "SUCCESS"

    def test_max_cycles_string_rejected(self, tm_fakes):
        result = _run("turbomole_opt", tm_fakes, charge=0, max_cycles="200")
        assert result.status.name == "FAILED"
        assert "max_cycles" in (result.error or "")


class TestOpenShell:
    def test_triplet_writes_unpaired(self, tm_fakes):
        _run("turbomole_sp", tm_fakes, charge=0, mult=3)
        control = (tm_fakes / "control").read_text()
        assert "$unpaired 2" in control

    def test_singlet_writes_no_unpaired(self, tm_fakes):
        _run("turbomole_sp", tm_fakes, charge=0, mult=1)
        control = (tm_fakes / "control").read_text()
        assert "$unpaired" not in control


class TestFreqHessianSurfaced:
    def test_hessian_file_reported_in_data(self, tm_fakes):
        result = _run("turbomole_freq", tm_fakes, charge=0)
        assert result.status.name == "SUCCESS", result.error
        assert result.data.get("hessian") == str(tm_fakes / "hessian")

    def test_missing_hessian_fails_open_but_honest(self, tmp_path, monkeypatch):
        # aoforce exits 0 and writes NOTHING: adapter must fail, not lie.
        work_dir = tmp_path / "w2"
        work_dir.mkdir()
        bin_dir = tmp_path / "fb2"
        bin_dir.mkdir()
        for name in ("ridft", "aoforce"):
            f = bin_dir / name
            f.write_text("#!/bin/sh\nexit 0\n")
            f.chmod(f.stat().st_mode | stat.S_IEXEC)
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        result = _run("turbomole_freq", work_dir, charge=0)
        assert result.status.name == "FAILED"
        assert "no hessian" in (result.error or "")


class TestMissingBinary:
    def test_absent_binary_is_failed_not_raised(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "w3"
        work_dir.mkdir()
        bin_dir = tmp_path / "fb3"
        bin_dir.mkdir()  # empty PATH dir: no Turbomole at all
        monkeypatch.setenv("PATH", str(bin_dir))
        result = _run("turbomole_sp", work_dir, charge=0)
        assert result.status.name == "FAILED"
        assert "not found in PATH" in (result.error or "")
