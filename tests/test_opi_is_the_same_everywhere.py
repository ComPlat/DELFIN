"""One OPI everywhere DELFIN runs.

ORCA's Python interface (OPI, PyPI ``orca-pi``) is how DELFIN builds ORCA
inputs and reads ORCA outputs, so every installation must carry the same
version: an input another version would build differently, or an output it
reads differently, is not the same calculation.

* The ``orca`` extra pins one exact version.
* The Pi installs it with every other extra (OPI is pure Python; the Pi
  builds and checks inputs and reads results even where ORCA itself does
  not run), through a numpy override of its own (scripts/pi/pi-overrides.txt).
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


def test_the_pi_installs_the_orca_extra_and_says_why_it_resolves():
    """The Pi installs the orca extra. It was skipped while orca-pi 2.0.0
    (numpy>=2.2.6, rdkit<2026) and mendeleev (numpy<2 declared on Python
    < 3.12) had no common solution on the Pi's Python 3.11. The Pi's own
    override file replaces the numpy requirements with a range inside the
    core's; the script must name that reason next to the skip list, so the
    override goes when mendeleev's metadata allows numpy 2 on 3.11."""
    text = (ROOT / "scripts" / "pi" / "pi_install.sh").read_text()
    m = re.search(r'SKIP_EXTRAS="\$\{DELFIN_PI_SKIP_EXTRAS:-([^}]*)\}"', text)
    assert m and "orca" not in m.group(1).split()
    assert "orca-pi" in text and "numpy>=2.2.6" in text and "mendeleev" in text
    assert 'PI_OVERRIDES="$ROOT/src/scripts/pi/pi-overrides.txt"' in text
    assert '--overrides "$PI_OVERRIDES"' in text


def _numpy_range(spec: str) -> tuple:
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version
    s = SpecifierSet(spec.split("numpy", 1)[1].split(";")[0].strip())
    lows = [Version(x.version) for x in s if x.operator == ">="]
    highs = [Version(x.version) for x in s if x.operator == "<"]
    return max(lows), min(highs)


def test_the_pi_numpy_override_stays_inside_the_core_range():
    """An override replaces every numpy requirement, the core's included,
    so its range must lie inside the core's: floor at or above, ceiling
    equal. Raising the core ceiling fails here until the override follows,
    and the override must still satisfy orca-pi (numpy>=2.2.6)."""
    import tomllib
    from packaging.version import Version
    core = [d for d in tomllib.loads((ROOT / "pyproject.toml").read_text())[
        "project"]["dependencies"] if re.match(r"numpy\b", d)]
    (core_spec,) = core
    lines = [ln.strip() for ln in (ROOT / "scripts" / "pi" / "pi-overrides.txt")
             .read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    assert len(lines) == 1 and lines[0].startswith("numpy"), lines
    lo, hi = _numpy_range(lines[0])
    core_lo, core_hi = _numpy_range(core_spec)
    assert lo >= core_lo and hi == core_hi
    assert lo >= Version("2.2.6")


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
    assert result.status == doctor.MISSING and f"orca-pi=={doctor.OPI_VERSION}" in result.fix_hint


def test_another_opi_version_is_broken(monkeypatch):
    _fake_opi(monkeypatch, "1.0.1")
    result = doctor.check_opi()
    assert result.status == doctor.BROKEN and doctor.OPI_VERSION in result.fix_hint


def test_the_doctor_runs_the_check():
    names = {r.name for r in doctor.run_all()}
    assert "opi" in names
