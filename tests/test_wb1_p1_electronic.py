"""WB1 P1 tests — electronic rules: diffuse basis, relativistic elements,
electron-count/multiplicity parity.

Sources (verified this session):
- Diffuse functions for anions: J. Phys. Chem. A 2013, 117, 10.1021/jp312755z;
  minimally augmented sets 10.1021/ct900566x; dipole/Rydberg states
  10.1007/s00214-015-1647-1.
- Parity: DELFIN's own rule — delfin/config.py:709-764 (ORCA stops a run with
  "multiplicity (4) is even and number of electrons (242) is even ->
  impossible"); parity-derived multiplicity in delfin/utils.py:414-428.
- Relativistic elements: delfin/utils.py:203 `_should_use_rel` (4d/5d/4f/5f/mixed).
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint
from delfin.utils import _should_use_rel


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _step(label, **kw):
    return {"name": "s", "steps": [{"step": "orca_sp", "label": label, "kwargs": kw}]}


# ── Rule 3: anions / Rydberg / diffuse excited states need a diffuse basis ──

def test_anion_with_non_diffuse_basis_is_flagged():
    out = _rules(_step("sp", charge=-1, mult=1, basis="def2-TZVP"))
    assert ("diffuse-basis-required", "sp", "warning") in out


def test_anion_with_diffuse_basis_is_silent():
    # ma- (minimally augmented) and aug- count as diffuse
    assert _rules(_step("sp", charge=-1, mult=1, basis="ma-def2-TZVP")) == set()
    assert _rules(_step("sp", charge=-1, mult=1, basis="aug-cc-pVTZ")) == set()


def test_rydberg_or_excited_states_require_diffuse_basis():
    for marker in ("rydberg", "diffuse_excited"):
        out = _rules(_step("es", charge=0, mult=1, basis="def2-SVP",
                           excited_state=marker))
        assert ("diffuse-basis-required", "es", "warning") in out, marker
    assert _rules(_step("es", charge=0, mult=1, basis="aug-cc-pVDZ",
                        excited_state="rydberg")) == set()


def test_composite_method_on_anion_is_not_flagged_for_basis():
    # hf-3c defines its own basis; the rule cannot fix it and stays silent
    assert _rules(_step("sp", charge=-1, mult=1, method="hf-3c")) == set()


# ── Rule 4: 4d/5d/4f/5f elements need a relativistic Hamiltonian + basis ──

def test_heavy_metal_without_relativity_is_flagged():
    # 4d element (Pd): _should_use_rel says True
    assert _should_use_rel(["Pd"])
    out = _rules(_step("sp", charge=0, mult=1, basis="def2-TZVP", metals=["Pd"]))
    assert ("relativistic-required", "sp", "warning") in out


def test_heavy_metal_with_relativity_is_silent():
    assert _rules(_step("sp", charge=0, mult=1, basis="ZORA-def2-TZVP",
                        metals=["Pd"], relativity="ZORA")) == set()


def test_light_metal_is_silent():
    # 3d only: no relativistic treatment required
    assert not _should_use_rel(["Fe"])
    assert _rules(_step("sp", charge=0, mult=1, basis="def2-SVP",
                        metals=["Fe"])) == set()


def test_rule4_decision_pin_agrees_with_should_use_rel():
    # the pin: the rule's decision must track utils._should_use_rel on the
    # same element set, so a change of one cannot silently diverge
    for metals, expect_findings in (
        (["Ag"], True),        # 4d
        (["Pt"], True),        # 5d
        (["Eu"], True),        # 4f
        (["Fe", "Pt"], True),  # mixed
        (["Fe"], False),       # 3d
        (["Zn"], False),       # 3d
    ):
        assert _should_use_rel(metals) is expect_findings, metals
        out = _rules(_step("sp", charge=0, mult=1, basis="def2-SVP",
                           metals=metals))
        found = any(r == "relativistic-required" for r, _s, _l in out)
        assert found is expect_findings, (metals, out)


# ── Rule 10: electron count compatible with multiplicity at every step ──

def test_even_electron_count_with_even_multiplicity_is_flagged():
    # 4 electrons, mult 4: ORCA stops ("multiplicity 4 is even and number of
    # electrons 4 is even -> impossible")
    assert ("electron-multiplicity-parity", "sp", "error") in _rules(
        _step("sp", charge=0, mult=4, electrons=4))


def test_odd_electron_count_with_even_multiplicity_is_silent():
    assert _rules(_step("sp", charge=0, mult=2, electrons=3)) == set()


def test_parity_mismatch_is_an_error_at_every_step():
    # even multiplicity needs an odd electron count (config.py:741: even
    # electrons -> odd multiplicity); mult 4 with 4 e- and mult 2 with 2 e-
    # are both "multiplicity ... -> impossible" for ORCA
    spec = {"name": "s", "steps": [
        {"step": "orca_opt", "label": "opt", "kwargs": {"mult": 4, "electrons": 4}},
        {"step": "orca_freq", "label": "freq", "kwargs": {"mult": 2, "electrons": 2}},
    ]}
    out = _rules(spec)
    assert ("electron-multiplicity-parity", "opt", "error") in out
    assert ("electron-multiplicity-parity", "freq", "error") in out


def test_no_electrons_kwarg_means_rule_is_silent():
    # the rule fires only on data it actually has (static lint, no parsing)
    assert _rules(_step("sp", charge=0, mult=4)) == set()
