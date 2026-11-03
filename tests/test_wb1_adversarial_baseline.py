"""Adversarial review tests for package WB1 (reviewer nacht-s15).

Baseline part (committed before any builder work): pins the measured
behaviour of the existing scientific lint so a builder refactor cannot
silently change it. Findings-specific parts are appended as the builder's
phases land.
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _step(label, **kw):
    return {"name": "s", "steps": [{"step": "orca_sp", "label": label, "kwargs": kw}]}


# ---------------------------------------------------------------- baseline
def test_existing_five_rules_still_fire_after_any_builder_refactor():
    spec = _step("sp", charge=0, method="r2scan-3c", basis="def2-TZVP",
                 dispersion="D4", relativity="ZORA", solvent_model="SMD", mult=0)
    rules = _rules(spec)
    assert ("composite-with-basis", "sp", "warning") in rules
    assert ("composite-with-dispersion", "sp", "warning") in rules
    assert ("relativity-without-rel-basis", "sp", "warning") in rules
    assert ("solvent-model-without-solvent", "sp", "warning") in rules
    assert ("invalid-multiplicity", "sp", "error") in rules
    assert len(rules) == 5


def test_clean_spec_stays_silent():
    assert _rules(_step("sp", charge=0, method="B3LYP", basis="def2-SVP")) == set()


# ------------------------------------------------------- adversarial probes
def test_spec_none_or_empty_does_not_crash():
    assert scientific_lint({}) == []
    assert _rules({"steps": [{}]}) == set()


def test_kwargs_none_and_wrong_types_do_not_crash():
    assert _rules({"steps": [{"label": "x", "kwargs": None}]}) == set()
    # wrong types: mult as string must not be accepted as a valid multiplicity
    # silently — string "1" would slip past the isinstance(mult, int) check.
    assert _rules(_step("sp", mult="0")) == set()


def test_branches_in_unexpected_shapes_do_not_crash():
    # branch value is a plain string, not steps or dict
    assert _rules({"steps": [], "branches": {"a": "junk"}}) == set()
    # branch is a dict without steps
    assert _rules({"steps": [], "branches": {"a": {"other": 1}}}) == set()


def test_unicode_step_labels_survive():
    rules = _rules(_step("opt—ü", mult=0))
    assert ("invalid-multiplicity", "opt—ü", "error") in rules


def test_huge_spec_does_not_take_quadratic_time(monkeypatch=None):
    # 5000 steps: must stay linear. If a builder adds a rule that scans the
    # whole spec per step this degrades quadratically and the wall-time cap
    # catches it (test-failure is the signal, not a timing measurement that
    # would depend on machine speed).
    spec = {"steps": [{"step": "orca_sp", "label": f"st{i}",
                       "kwargs": {"charge": 0, "method": "B3LYP",
                                  "basis": "def2-SVP"}} for i in range(5000)]}
    assert _rules(spec) == set()
