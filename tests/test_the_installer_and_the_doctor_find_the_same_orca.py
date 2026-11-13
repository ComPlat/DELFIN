"""The installer and `delfin doctor` answer "where is ORCA" with one function.

Measured on a fresh install (2026-11-13): the installer said "ORCA not found"
while `delfin doctor` on the same account reported ORCA at /opt/orca_6_1_1.
The installer looked at its settings, PATH and orca_* directories in HOME and
the prefix; the doctor asked DELFIN's tool resolver (ORCA_BINARY/ORCA_PATH,
PATH, the system directories, modules). Both now ask
delfin.qm_runtime.find_orca_dir.
"""

import os
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
INSTALLER = REPO / "delfin" / "installers" / "install_delfin.sh"


def _orca_dir(path: pathlib.Path) -> pathlib.Path:
    path.mkdir(parents=True)
    for name in ("orca", "orca_plot"):
        (path / name).write_text("#!/bin/sh\nexit 0\n")
        (path / name).chmod(0o755)
    return path


def _nothing_resolved(monkeypatch):
    from delfin import qm_runtime
    monkeypatch.setattr(qm_runtime, "resolve_tool", lambda name: None)


def test_an_unpacked_orca_under_the_given_roots_is_found(tmp_path, monkeypatch):
    from delfin import qm_runtime
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    _nothing_resolved(monkeypatch)
    older = _orca_dir(tmp_path / "software" / "orca_6_0_1_linux_x86-64_shared_openmpi416")
    newer = _orca_dir(tmp_path / "software" / "orca_6_1_1_linux_x86-64_shared_openmpi418_avx2")

    assert qm_runtime.find_orca_dir([tmp_path / "software"]) == str(newer)
    assert str(older) != str(newer)


def test_what_the_resolver_finds_is_found(tmp_path, monkeypatch):
    from delfin import qm_runtime
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    site = _orca_dir(tmp_path / "opt" / "orca_6_1_1")
    monkeypatch.setattr(qm_runtime, "resolve_tool", lambda name: qm_runtime.ResolvedTool(
        requested_name=name, canonical_name=name, path=str(site / "orca"), source="system"))

    assert qm_runtime.find_orca_dir() == str(site)


def test_a_program_merely_called_orca_is_not_orca(tmp_path, monkeypatch):
    """/usr/bin/orca is a screen reader on many desktops."""
    from delfin import qm_runtime
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    reader = tmp_path / "bin" / "orca"
    reader.parent.mkdir()
    reader.write_text("#!/bin/sh\nexit 0\n")
    reader.chmod(0o755)
    monkeypatch.setattr(qm_runtime, "resolve_tool", lambda name: qm_runtime.ResolvedTool(
        requested_name=name, canonical_name=name, path=str(reader), source="PATH"))

    assert qm_runtime.find_orca_dir() == ""


def test_the_doctor_reports_the_orca_the_installer_would_wire(tmp_path, monkeypatch):
    from delfin import doctor, qm_runtime
    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    for name in ("ORCA_BINARY", "ORCA_PATH", "DELFIN_ORCA_BINARY"):
        monkeypatch.delenv(name, raising=False)
    unpacked = _orca_dir(home / "software" / "orca_6_1_1_linux_x86-64_shared_openmpi418_avx2")
    real_resolve = qm_runtime.resolve_tool

    def resolve(name):
        # The host's own ORCA (if any) is not this test's subject; what the
        # environment names is resolved as before.
        if os.environ.get("ORCA_BINARY"):
            return real_resolve(name)
        return None

    monkeypatch.setattr(qm_runtime, "resolve_tool", resolve)
    qm_runtime.clear_resolver_cache()

    result = doctor.check_orca()

    assert result.status == doctor.OK, result
    assert str(unpacked) in result.detail
    assert "ORCA_BINARY" not in os.environ


def test_the_installer_wires_the_orca_the_resolver_finds(tmp_path):
    """ORCA_BINARY is one of the places the doctor's resolver reads."""
    site = _orca_dir(tmp_path / "site" / "orca_6_1_1_linux_x86-64_shared_openmpi418_avx2")
    home = tmp_path / "home"
    home.mkdir()
    done = subprocess.run(
        ["bash", str(INSTALLER), "--dry-run", "--profile", "core", "--no-openmpi"],
        capture_output=True, text=True, timeout=120,
        env={"HOME": str(home), "DELFIN_REPO": str(REPO), "DELFIN_PYTHON": sys.executable,
             "PATH": os.pathsep.join([str(pathlib.Path(sys.executable).parent), "/usr/bin", "/bin"]),
             "ORCA_BINARY": str(site / "orca")})

    assert done.returncode == 0, done.stderr
    assert f"ORCA     {site}" in done.stdout, done.stdout
