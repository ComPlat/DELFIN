"""WB1 package tests — physical realism rules + baseline protection.

Part 1 (baseline, from reviewer s15 finding F1): the lint must not crash on
unexpected spec shapes; a malformed branch value is skipped, not fatal.
Part 2: the existing five rules keep firing through the refactor.
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _step(label, **kw):
    return {"name": "s", "steps": [{"step": "orca_sp", "label": label, "kwargs": kw}]}


# ── F1: unexpected branch shapes must not crash the lint ─────────────────

def test_branches_in_unexpected_shapes_do_not_crash():
    # a branch value that is a bare string, not a steps list nor a dict
    spec = _step("sp", charge=0, mult=0)
    spec["branches"] = {"a": "junk"}
    assert ("invalid-multiplicity", "sp", "error") in _rules(spec)

    # a branch value that is a dict without "steps"
    spec2 = _step("sp", charge=0, mult=0)
    spec2["branches"] = {"a": {"other": 1}}
    assert ("invalid-multiplicity", "sp", "error") in _rules(spec2)

    # steps that are not dicts at all
    spec3 = {"name": "s", "steps": ["junk", None], "branches": {}}
    assert _rules(spec3) == set()


# ── baseline: the existing five rules keep firing through the refactor ───

def test_existing_five_rules_still_fire_after_any_builder_refactor():
    out = _rules(_step("sp", charge=0, method="r2scan-3c", basis="def2-TZVP",
                       dispersion="D4", solvent_model="SMD"))
    assert ("composite-with-basis", "sp", "warning") in out
    assert ("composite-with-dispersion", "sp", "warning") in out
    assert ("solvent-model-without-solvent", "sp", "warning") in out
    assert ("relativity-without-rel-basis", "opt", "warning") in _rules(
        _step("opt", charge=0, relativity="ZORA", basis="def2-TZVP"))
    assert ("invalid-multiplicity", "sp", "error") in _rules(
        _step("sp", charge=0, mult=0))
