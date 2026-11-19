"""Every place that builds a py3Dmol viewer gets the </script> rule.

py3Dmol writes viewer arguments into an inline <script> with json.dumps,
which leaves ``</script>`` intact; molecule_viewer._harden_py3dmol
rewrites them as unicode escapes. Since the lazy import (5e79c6a87) the
rule is installed by the accessor ``molecule_viewer.py3dmol()`` only,
while tab_chemdarwin and tab_turbomole_builder import py3Dmol
themselves: a viewer built there before any accessor call carried the
text unescaped. Input: a fresh py3Dmol (reloaded, so no rule is left
from an earlier test) and an XYZ comment line holding </script>.
Output: no raw script end in the generated HTML.
"""
from __future__ import annotations

import importlib
import pathlib
import re

import pytest

py3Dmol = pytest.importorskip("py3Dmol")

MARKER = "</script><b>marker</b>"
_ROOT = pathlib.Path(__file__).resolve().parents[1] / "delfin"


def _fresh():
    mod = importlib.reload(py3Dmol)
    assert not getattr(mod.tostr, "_delfin_script_safe", False)
    return mod


def _html(mod):
    view = mod.view(width=100, height=100)
    view.addModel(f"1\n{MARKER}\nH 0 0 0\n", "xyz")
    return view._make_html()


@pytest.mark.parametrize("module", [
    "delfin.dashboard.tab_chemdarwin",
    "delfin.dashboard.tab_turbomole_builder",
])
def test_a_tab_that_imports_py3dmol_itself_installs_the_rule(module):
    mod = _fresh()
    importlib.reload(importlib.import_module(module))
    assert "<b>marker</b>" not in _html(mod)


def test_the_accessor_hardens_a_module_that_was_reloaded():
    from delfin.dashboard import molecule_viewer
    molecule_viewer.py3dmol()
    mod = _fresh()
    assert molecule_viewer.py3dmol() is mod
    assert "<b>marker</b>" not in _html(mod)


def test_no_module_imports_py3dmol_without_the_rule():
    """A new direct import must harden what it imported."""
    offenders = []
    for path in _ROOT.rglob("*.py"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if not re.search(r"^\s*import py3Dmol\b|^\s*from py3Dmol\b", text, re.M):
            continue
        if "_harden_py3dmol(" not in text:
            offenders.append(str(path.relative_to(_ROOT)))
    assert offenders == [], offenders
