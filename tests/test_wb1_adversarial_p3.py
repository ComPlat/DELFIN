"""Adversarial P3 tests (reviewer nacht-s15) for WB1 rules 5/6/7.

All probes go through the real entry point scientific_lint. Each case
states the physically right answer; a failure means the rule mis-fires
or misses.
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _spec(steps):
    return {"name": "s", "steps": steps}


def _st(step, label=None, **kw):
    d = {"step": step, "kwargs": kw}
    if label:
        d["label"] = label
    return d


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


# ------------------------------------------------------------ rule 5
def test_rule5_fires_for_noncovalent_without_dispersion():
    spec = _spec([_st("orca_sp", label="a", system="noncovalent",
                      method="B3LYP", basis="def2-TZVP")])
    assert ("dispersion-recommended", "a", "warning") in _rules(spec)


def test_rule5_silent_with_explicit_dispersion():
    spec = _spec([_st("orca_sp", label="a", system="noncovalent",
                      method="B3LYP", basis="def2-TZVP", dispersion="D4")])
    assert "dispersion-recommended" not in {r for r, _, _ in _rules(spec)}


def test_rule5_silent_for_composite_method():
    spec = _spec([_st("orca_sp", label="a", system="flexible",
                      method="r2scan-3c", basis="def2-mTZVP")])
    assert "dispersion-recommended" not in {r for r, _, _ in _rules(spec)}


def test_rule5_uppercase_method_still_recognised_as_composite():
    # users write R2SCAN-3c; the composite set is matched lowercased.
    spec = _spec([_st("orca_sp", label="a", system="flexible",
                      method="R2SCAN-3C")])
    assert "dispersion-recommended" not in {r for r, _, _ in _rules(spec)}


def test_rule5_silent_on_non_orca_step():
    # a crest step with a noncovalent hint must not demand dispersion
    spec = _spec([_st("crest_conformers", label="c1", system="noncovalent")])
    assert "dispersion-recommended" not in {r for r, _, _ in _rules(spec)}


def test_rule5_silent_when_system_hint_absent():
    spec = _spec([_st("orca_sp", label="a", method="B3LYP")])
    assert "dispersion-recommended" not in {r for r, _, _ in _rules(spec)}


# ------------------------------------------------------------ rule 6
def test_rule6_fires_for_flexible_without_conformer_search():
    spec = _spec([_st("orca_opt", label="o1", system="flexible",
                      method="B3LYP", basis="def2-SVP", dispersion="D4")])
    assert ("conformer-search-missing", "o1", "warning") in _rules(spec)


def test_rule6_silent_when_crest_ran_earlier():
    spec = _spec([_st("crest_conformers", label="c1"),
                  _st("orca_opt", label="o1", system="flexible",
                      method="B3LYP", basis="def2-SVP", dispersion="D4")])
    rules = _rules(spec)
    assert "conformer-search-missing" not in {r for r, _, _ in rules}


def test_rule6_finding_F6_crest_after_the_orca_step_does_not_close_it():
    # A conformer search satisfies rule 6 only if it came BEFORE the DFT
    # step. If crest runs after, the DFT geometry was still picked from
    # an unsearched start — the finding must still fire on the orca step.
    spec = _spec([_st("orca_opt", label="o1", system="flexible",
                      method="B3LYP", basis="def2-SVP", dispersion="D4"),
                  _st("crest_conformers", label="c1")])
    assert ("conformer-search-missing", "o1", "warning") in _rules(spec)


# ------------------------------------------------------------ rule 7
def test_rule7_fires_for_3d_complex_without_scan_or_bs():
    spec = _spec([_st("orca_sp", label="a", system="3d_complex",
                      metals=["Fe"], method="B3LYP", basis="def2-SVP")])
    assert ("spin-state-scan-missing", "a", "warning") in _rules(spec)


def test_rule7_silent_with_broken_symmetry_kwarg():
    spec = _spec([_st("orca_sp", label="a", system="3d_complex",
                      metals=["Fe"], method="B3LYP", basis="def2-SVP",
                      broken_symmetry=True)])
    assert "spin-state-scan-missing" not in {r for r, _, _ in _rules(spec)}


def test_rule7_silent_after_occupier_scan_earlier():
    spec = _spec([_st("occupier", label="scan1"),
                  _st("orca_sp", label="a", system="3d_complex",
                      metals=["Fe"], method="B3LYP", basis="def2-SVP")])
    rules = _rules(spec)
    assert "spin-state-scan-missing" not in {r for r, _, _ in rules}


def test_rule7_silent_without_metals_kwarg():
    # the hint alone must not invent chemistry (documented design)
    spec = _spec([_st("orca_sp", label="a", system="3d_complex",
                      method="B3LYP", basis="def2-SVP")])
    assert "spin-state-scan-missing" not in {r for r, _, _ in _rules(spec)}


# -------------------------------------------------------- robustness
def test_robustness_system_hint_shapes():
    # hint on a branch step, hint=None, hint unicode — nothing crashes
    assert _rules({"steps": [], "branches": {
        "a": [_st("orca_sp", label="b1", system="flexible",
                  method="B3LYP", dispersion="D4")]}}) == set()
    assert _rules(_spec([_st("orca_sp", label="a", system=None)])) == set()
    assert _rules(_spec([_st("orca_sp", label="a",
                             system="спутник")])) == set()
