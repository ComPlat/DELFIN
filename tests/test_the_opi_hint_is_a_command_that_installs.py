"""delfin doctor's OPI hint is a command that installs OPI into DELFIN's Python.

The hint was `pip install "delfin-complat[orca]"`. On the installer's Python
3.11 that set has no solution (orca-pi 2.0.0 wants numpy>=2.2.6 and
rdkit<2026; mendeleev declares numpy<2 below Python 3.12), and from a git
checkout it would fetch DELFIN from PyPI over the editable install. Measured
2026-11-13 in a Python 3.11 venv with DELFIN's mendeleev/numpy/rdkit: a
plain `pip install orca-pi==2.0.0` would pull numpy 2.4.6 (outside DELFIN's
<2.3); with numpy and rdkit bounded it installs numpy 2.2.6, rdkit 2025.9.6
and orca-pi 2.0.0, pip names mendeleev's declared numpy<2, and mendeleev's
element data reads as before (Fe 55.845).
"""

import sys

from delfin import doctor


def _no_opi(monkeypatch):
    from importlib import metadata
    real = metadata.version

    def fake(name):
        if name == "orca-pi":
            raise metadata.PackageNotFoundError(name)
        return real(name)
    monkeypatch.setattr(metadata, "version", fake)


def test_the_hint_installs_opi_into_this_python_within_delfins_bounds(monkeypatch):
    _no_opi(monkeypatch)
    hint = doctor.check_opi().fix_hint

    assert f"{sys.executable} -m pip install" in hint
    assert f'"orca-pi=={doctor.OPI_VERSION}"' in hint
    assert '"numpy>=2.2.6,<2.3"' in hint
    assert '"rdkit>=2025.3.2,<2026"' in hint
    assert "delfin-complat[orca]" not in hint
