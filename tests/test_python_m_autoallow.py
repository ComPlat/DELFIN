"""`python -m <module>` (a user's own CLI) is auto-allowed like `python x.py`,
but `python -m pip install/uninstall/download` stays behind the confirm gate.

bug 20260718-193300: the agent built a `molkit` package but could not run its
own CLI (`python -m molkit "Fe2(SO4)3"`) because only a fixed module whitelist
was auto-allowed — even though running a script (equally arbitrary code) was.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path) -> KitToolPermissions:
    return KitToolPermissions(workspace=tmp_path)


@pytest.mark.parametrize("cmd", [
    'python -m molkit "Fe2(SO4)3"',
    "python3 -m molkit H2O",
    "python -m mypackage.cli --flag",
    "python3.11 -m app",
    "python -m pytest tests/ -v",        # existing whitelist still works
    "python -m http.server",
    "python -m pip show numpy",          # read-only pip stays allowed
])
def test_python_m_module_auto_allowed(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is True, cmd


@pytest.mark.parametrize("cmd", [
    "python -m pip install requests",
    "python -m pip uninstall numpy",
    "python -m pip download scipy",
    "python3 -m pip install -r requirements.txt",
])
def test_pip_install_still_gated(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is False, cmd
