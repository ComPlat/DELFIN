"""One OPI everywhere DELFIN runs.

ORCA's Python interface (OPI, PyPI ``orca-pi``) is how DELFIN builds ORCA
inputs and reads ORCA outputs, so every installation must carry the same
version: an input another version would build differently, or an output it
reads differently, is not the same calculation.

* The ``orca`` extra pins one exact version.
* The Pi installs it with every other extra (OPI is pure Python; the Pi
  builds and checks inputs and reads results even where ORCA itself does
  not run).
* ``delfin doctor`` names the OPI it finds: the pinned one is ok, none is
  "missing" with the install line, another version is "broken".
"""

from __future__ import annotations

import re
from pathlib import Path


from delfin import doctor

ROOT = Path(__file__).resolve().parents[1]


def _pin() -> str:
    import tomllib
    extra = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"][
        "optional-dependencies"]["orca"]
    (spec,) = [d for d in extra if d.startswith("orca-pi")]
    m = re.fullmatch(r"orca-pi==(\d+\.\d+\.\d+)", spec)
    assert m, f"the orca extra must pin one exact OPI version, not {spec!r}"
    return m.group(1)


def test_the_extra_pins_one_exact_version():
    assert _pin() == doctor.OPI_VERSION


def test_the_pi_installs_the_orca_extra():
    text = (ROOT / "scripts" / "pi" / "pi_install.sh").read_text()
    m = re.search(r'SKIP_EXTRAS="\$\{DELFIN_PI_SKIP_EXTRAS:-([^}]*)\}"', text)
    assert m and "orca" not in m.group(1).split()


def _fake_opi(monkeypatch, version):
    from importlib import metadata

    real = metadata.version

    def fake(name):
        if name == "orca-pi":
            if version is None:
                raise metadata.PackageNotFoundError(name)
            return version
        return real(name)
    monkeypatch.setattr(metadata, "version", fake)


def test_the_doctor_names_the_pinned_opi_as_ok(monkeypatch):
    _fake_opi(monkeypatch, doctor.OPI_VERSION)
    result = doctor.check_opi()
    assert result.status == doctor.OK and doctor.OPI_VERSION in result.detail


def test_without_opi_the_doctor_says_how_to_install_it(monkeypatch):
    _fake_opi(monkeypatch, None)
    result = doctor.check_opi()
    assert result.status == doctor.MISSING and "[orca]" in result.fix_hint


def test_another_opi_version_is_broken(monkeypatch):
    _fake_opi(monkeypatch, "1.0.1")
    result = doctor.check_opi()
    assert result.status == doctor.BROKEN and doctor.OPI_VERSION in result.fix_hint


def test_the_doctor_runs_the_check():
    names = {r.name for r in doctor.run_all()}
    assert "opi" in names
