"""Site-install defects from the TM2 real run #608 (operator, 08.11.).

Two findings, each with a red control driven through the public path:

1. The ``define`` dialog can break (non-zero exit) — the old
   ``_setup_control`` never checked its returncode, so the adapter went on
   to run ridft on a missing/partial ``control`` and the failure surfaced
   as a misleading SCF/convergence error. The adapter must report the
   ``define`` failure itself, before any SCF program starts.

2. ``detect_installed_version`` returned None on the site install: the
   TURBODIR path carries no version digits there and no explicit
   module name is passed. The environment module system leaves the
   resolved module name behind in ``LOADEDMODULES``/``_LMFILES_`` —
   detection must read it instead of giving up (never assume a version,
   but also never report unknown when the evidence is present).
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get
from delfin.turbomole_version import detect_installed_version


def _make_bin(bin_dir: Path, name: str, body: str) -> Path:
    p = bin_dir / name
    p.write_text(body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


class TestDefineFailureIsSurfaced:
    """A broken define dialog is a define failure, not an SCF failure."""

    def test_define_nonzero_exit_fails_before_ridft(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir()
        # define breaks (real site: dialog aborts, non-zero exit, stderr).
        _make_bin(bin_dir, "define",
                  "#!/bin/sh\necho 'define: input error' >&2\nexit 1\n")
        # If the adapter wrongly continues, ridft would "converge" quietly:
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        adapter = get("turbomole_sp")
        assert adapter is not None
        result = adapter.execute(work_dir, geometry=_xyz(work_dir),
                                 cores=1, charge=0)
        assert result.status.name == "FAILED"
        # The error names define, not ridft/SCF (the old misreport).
        assert "define" in (result.error or "").lower(), result.error
        # ridft must never have run: no energy file was produced by it.
        assert not (work_dir / "energy").is_file(), (
            "the adapter ran ridft although define failed")

    def test_define_ok_still_runs_ridft(self, tmp_path, monkeypatch):
        """Nothing new is refused: a working define path still succeeds."""
        work_dir = tmp_path / "work2"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin2"
        bin_dir.mkdir()
        _make_bin(bin_dir, "define",
                  f"#!/bin/sh\nprintf 'x\\n' > {work_dir}/control\n")
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        adapter = get("turbomole_sp")
        result = adapter.execute(work_dir, geometry=_xyz(work_dir),
                                 cores=1, charge=0)
        assert result.status.name == "SUCCESS", result.error


class TestVersionFromLoadedModules:
    """The module system's own records are a version source."""

    def test_loadedmodules_env_is_read(self, monkeypatch):
        # Site reality: TURBODIR without digits, but the module is loaded.
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("LOADEDMODULES", "chem/turbomole/7.9:intel/2022")
        ver = detect_installed_version()
        assert ver is not None, (
            "detect_installed_version must read LOADEDMODULES on the site "
            "install instead of returning None")
        assert str(ver) == "7.9"
        assert ver.source == "module"

    def test_lmfiles_env_is_read(self, monkeypatch):
        monkeypatch.delenv("LOADEDMODULES", raising=False)
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("_LMFILES_", "/opt/modules/chem/turbomole/7.9")
        ver = detect_installed_version()
        assert ver is not None
        assert str(ver) == "7.9"

    def test_explicit_module_name_still_wins_over_env(self, monkeypatch):
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("LOADEDMODULES", "chem/turbomole/7.9")
        ver = detect_installed_version(module_name="chem/turbomole/8.1")
        assert str(ver) == "8.1"

    def test_no_evidence_still_returns_none(self, monkeypatch):
        for name in ("TURBODIR", "LOADEDMODULES", "_LMFILES_"):
            monkeypatch.delenv(name, raising=False)
        assert detect_installed_version() is None
