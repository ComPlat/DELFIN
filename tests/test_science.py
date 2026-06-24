"""Tests for the scientific-plausibility lint (Track 5)."""

from __future__ import annotations

from delfin.tools import platform
from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _step(label, **kw):
    return {"name": "s", "steps": [{"step": "orca_sp", "label": label, "kwargs": kw}]}


def test_composite_method_flags_redundant_basis_and_dispersion():
    out = _rules(_step("sp", charge=0, method="r2scan-3c", basis="def2-TZVP", dispersion="D4"))
    assert ("composite-with-basis", "sp", "warning") in out
    assert ("composite-with-dispersion", "sp", "warning") in out


def test_relativity_requires_relativistic_basis():
    bad = _rules(_step("opt", charge=0, relativity="ZORA", basis="def2-TZVP"))
    assert ("relativity-without-rel-basis", "opt", "warning") in bad
    # a properly recontracted basis passes
    assert _rules(_step("opt", charge=0, relativity="ZORA", basis="ZORA-def2-TZVP")) == set()


def test_solvent_model_without_solvent():
    assert ("solvent-model-without-solvent", "sp", "warning") in _rules(
        _step("sp", charge=0, solvent_model="SMD"))
    # with a solvent it is fine
    assert _rules(_step("sp", charge=0, solvent_model="SMD", solvent="water")) == set()


def test_invalid_multiplicity_is_an_error():
    assert ("invalid-multiplicity", "sp", "error") in _rules(_step("sp", charge=0, mult=0))
    assert _rules(_step("sp", charge=0, mult=1)) == set()


def test_clean_spec_and_placeholders_yield_nothing():
    assert _rules(_step("sp", charge=0, method="B3LYP", basis="def2-SVP")) == set()
    # unresolved template placeholders are skipped (checked after build)
    assert _rules(_step("sp", charge=0, method="{method}", basis="{basis}")) == set()


def test_lint_walks_branches_and_application_specs():
    spec = {
        "name": "b", "steps": [],
        "branches": {"ox": [{"step": "orca_sp", "label": "ox_sp",
                             "kwargs": {"charge": 1, "mult": 0}}]},
    }
    assert ("invalid-multiplicity", "ox_sp", "error") in _rules(spec)
    # application dict ({"spec": {...}}) is unwrapped
    app = {"name": "a", "spec": _step("sp", charge=0, mult=-1)}
    assert ("invalid-multiplicity", "sp", "error") in _rules(app)


def test_platform_scientific_lint_returns_dicts():
    findings = platform.scientific_lint(_step("sp", charge=0, mult=0))
    assert findings and findings[0]["rule"] == "invalid-multiplicity"
    assert set(findings[0]) == {"step", "level", "rule", "message"}
