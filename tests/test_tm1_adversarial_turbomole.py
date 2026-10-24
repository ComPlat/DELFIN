"""Adversarial tests for package TM1 (reviewer nacht-s25).

Drives the real public path: the registered Turbomole adapters' execute()
and the scientific_lint entry point. Turbomole itself is never needed —
where a binary is required, a recording fake is put on PATH.

Findings here are defects until the builder proves otherwise.
"""

from __future__ import annotations

import json
import os
import stat
import textwrap
from pathlib import Path

import pytest

from delfin.tools._registry import get
from delfin.tools._science import scientific_lint


@pytest.fixture
def tm_env(tmp_path, monkeypatch):
    """A work_dir plus a recording fake `jobex` first on PATH.

    Phase-1 note: the adapter now calls ridft (SP) and aoforce (freq) too,
    so the env provides recording fakes for ridft and aoforce. ridft writes
    an `energy` file (success); aoforce writes a `hessian` file (success).
    """
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    log = tmp_path / "fake_jobex.log"
    fake = bin_dir / "jobex"
    fake.write_text(
        f"#!/bin/sh\nprintf '%s\\n' \"$@\" > {log}\nexit 1\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    ridft = bin_dir / "ridft"
    ridft.write_text(f"#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > {work_dir}/energy\n")
    ridft.chmod(ridft.stat().st_mode | stat.S_IEXEC)
    aoforce = bin_dir / "aoforce"
    aoforce.write_text(f"#!/bin/sh\nprintf '$hessian\\n1\\n$end\\n' > {work_dir}/hessian\n")
    aoforce.chmod(aoforce.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return work_dir, log


def _xyz(work_dir: Path, n: int = 3) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


def _run(name: str, work_dir: Path, **kwargs):
    adapter = get(name)
    assert adapter is not None, f"adapter {name} not registered"
    return adapter.execute(work_dir, geometry=_xyz(work_dir), cores=1, **kwargs)


class TestOpenShellMultiplicity:
    """mult != 1 must not be silently dropped."""

    def test_triplet_control_declares_unpaired_electrons(self, tm_env):
        """A triplet SP must declare its spin in the Turbomole control file.

        Current manual control writer (delfin/tools/adapters/turbomole.py,
        _setup_control) never writes $unpaired — an open-shell request
        silently computes a closed-shell singlet. Physical-realism defect.
        """
        work_dir, _ = tm_env
        _run("turbomole_sp", work_dir, charge=0, mult=3)
        control = (work_dir / "control").read_text()
        assert "$unpaired" in control, (
            "triplet requested but control has no $unpaired — "
            "open shell silently computed as singlet"
        )


class TestOptArgumentSafety:
    """max_cycles is interpolated into a shell command — no injection."""

    def test_semicolon_in_max_cycles_is_not_executed(self, tm_env):
        """`jobex -ri -c {max_cycles}` runs with shell=True (turbomole.py:262).

        A hostile max_cycles must never reach a shell as a command.
        The fake jobex records what it saw; the injected payload must not run.
        """
        work_dir, log = tm_env
        injected = work_dir / "pwned"
        payload = f"200; touch {injected}"
        result = _run("turbomole_opt", work_dir, charge=0, max_cycles=payload)
        assert not injected.exists(), (
            f"max_cycles was passed to a shell verbatim; the payload ran: "
            f"{injected} created. argv seen: {log.read_text() if log.exists() else '(none)'}"
        )
        assert result.status.name == "FAILED"  # fake exits 1; run itself is contained

    def test_max_cycles_still_numeric_when_clean(self, tm_env):
        """A clean integer keeps the real semantics (allow-list, not deny-list)."""
        work_dir, log = tm_env
        _run("turbomole_opt", work_dir, charge=0, max_cycles=37)
        argv = log.read_text().splitlines()
        assert argv == ["-ri", "-c", "37"], f"argv seen by jobex: {argv}"


class TestFreqHessian:
    """turbomole_freq declares `produces ("qc_output", "hessian")`.

    Fixed in c64fb267: the adapter surfaces the aoforce hessian file in
    data["hessian"], or fails honestly when aoforce exits 0 without one.
    Both directions tested; the negative case builds its own bin dir whose
    aoforce writes nothing.
    """

    def test_freq_produces_declared_hessian(self, tm_env):
        work_dir, _ = tm_env
        result = _run("turbomole_freq", work_dir, charge=0)
        assert result.status.name == "SUCCESS", result.error
        assert result.data.get("hessian"), (
            "adapter declares a hessian product but execute() produced none"
        )

    def test_aoforce_silent_success_is_reported_failed(self, tmp_path, monkeypatch):
        """aoforce exiting 0 with NO hessian must fail honestly: SUCCESS
        requires the declared product to exist, not just rc == 0."""
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir()
        ridft = bin_dir / "ridft"
        ridft.write_text(f"#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > {work_dir}/energy\n")
        ridft.chmod(ridft.stat().st_mode | stat.S_IEXEC)
        aoforce = bin_dir / "aoforce"
        aoforce.write_text("#!/bin/sh\nexit 0\n")
        aoforce.chmod(aoforce.stat().st_mode | stat.S_IEXEC)
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        result = _run("turbomole_freq", work_dir, charge=0)
        assert result.status.name == "FAILED"
        assert "no hessian" in result.error


class TestScientificLintParity:
    """WB1 rule 10 must reject a physically impossible TM spec."""

    def test_impossible_multiplicity_flagged(self):
        """Physics: Nα+Nβ and Nα−Nβ share parity, so nelec ≡ mult−1 (mod 2).

        mult=3 (odd) needs an EVEN electron count; 9 electrons with a
        triplet is impossible. Verified direction with a probe run first
        (test_zz_probe_parity.py): 3/10 and 4/9 pass, 3/9 and 4/10 flag.
        """
        findings = scientific_lint(
            {"steps": [{"label": "sp", "kwargs": {
                "method": "b3-lyp", "basis": "def2-SVP",
                "mult": 3, "electrons": 9, "program": "turbomole"}}]}
        )
        assert any(f.rule == "electron-multiplicity-parity" for f in findings), (
            [f.rule for f in findings]
        )

    def test_valid_triplet_not_flagged(self):
        findings = scientific_lint(
            {"steps": [{"label": "sp", "kwargs": {
                "method": "b3-lyp", "basis": "def2-SVP",
                "mult": 3, "electrons": 10, "program": "turbomole"}}]}
        )
        assert not any(f.rule == "electron-multiplicity-parity" for f in findings), (
            [f.rule for f in findings]
        )
