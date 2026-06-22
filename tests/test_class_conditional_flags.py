"""Tests for ``_class_conditional_flag`` per-class wrappers (Phase 3B + Wave-5).

Phase 3B (commit 5dd520a) introduced the generic helper that combines:

  1. ``DELFIN_<NAME>_CLASSES="cls1,cls2"`` — flag True iff
     ``_classify_complex_class(mol)`` is in that list.
  2. ``DELFIN_<NAME>=1/0`` — scalar fallback (default 0 unless caller passes
     ``default=1``).

This Wave extends the same helper to two more existing scalar flags:

  - ``DELFIN_UFF_SOFT_DONORS`` (consumed in ``_optimize_xyz_openbabel``).
  - ``DELFIN_B5_RIGID_H``      (consumed in ``_apply_baustein5_if_enabled``).

Pool-verdict use case: if soft-donor mode helps only the sigma class while
hurting hapto, an operator can ship

    export DELFIN_UFF_SOFT_DONORS_CLASSES="sigma"

without touching the rest of the pipeline.  Empty ``_CLASSES`` env (default)
is bit-identical to the pre-wave scalar behaviour.

The tests below cover the contract for *both* call-sites:

  * env unset            → uses scalar (back-compat, default OFF)
  * scalar=1             → enabled for every class
  * CLASSES="sigma"      → enabled only for sigma
  * CLASSES="a,b"        → enabled for both listed classes
  * CLASSES set + mol=None / unknown class → safe fallback (disabled)
  * malformed CLASSES (whitespace / trailing commas / unknown labels)
    → safe fallback (disabled, no crash)

All tests skip cleanly when RDKit is not importable so the file is harmless
in minimal CI environments.
"""
from __future__ import annotations

import os
from typing import Iterable, Optional

import pytest

pytest.importorskip("rdkit", reason="RDKit required for class-conditional tests")
from rdkit import Chem  # noqa: E402 — after importorskip guard

from delfin.smiles_converter import (  # noqa: E402
    _class_conditional_flag,
    _classify_complex_class,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

#: SMILES → expected ``_classify_complex_class`` label.  Chosen so the
#: test-suite exercises sigma, multi_sigma, hapto and no_metal without
#: depending on hard-to-parse cyclopentadienyl SMILES (which the classifier
#: would still label as hapto via ``_find_hapto_groups``).
_SMILES_FIXTURES = {
    "sigma":       "[Pt](Cl)(Cl)(N)N",
    "multi_sigma": "[Mn]([C]#O)([C]#O)[Mn]([C]#O)([C]#O)",
    "no_metal":    "CC(=O)O",
}


def _mol(smi: str):
    """Parse SMILES with sanitize=False to keep dative/bracket forms intact."""
    m = Chem.MolFromSmiles(smi, sanitize=False)
    assert m is not None, f"SMILES did not parse: {smi}"
    return m


def _scrub_env(*names: Iterable[str]) -> None:
    """Remove the given env-vars from ``os.environ`` (no-op if absent)."""
    for n in names:
        os.environ.pop(n, None)


@pytest.fixture(autouse=True)
def _clean_env():
    """Strip all DELFIN_TEST_* vars before/after every test to keep isolation."""
    targets = [
        "DELFIN_TEST_FLAG",
        "DELFIN_TEST_FLAG_CLASSES",
        "DELFIN_UFF_SOFT_DONORS",
        "DELFIN_UFF_SOFT_DONORS_CLASSES",
        "DELFIN_B5_RIGID_H",
        "DELFIN_B5_RIGID_H_CLASSES",
    ]
    _scrub_env(*targets)
    yield
    _scrub_env(*targets)


# ----------------------------------------------------------------------
# 1. env unset → uses scalar default (back-compat)
# ----------------------------------------------------------------------

class TestEnvUnsetScalarBackCompat:
    """When neither ``DELFIN_<name>`` nor ``DELFIN_<name>_CLASSES`` is set the
    helper must return ``bool(default)`` — and the production default is 0
    (= disabled), so the call must be bit-identical to the pre-wave behaviour.
    """

    def test_default_zero_means_disabled(self):
        mol = _mol(_SMILES_FIXTURES["sigma"])
        assert _class_conditional_flag("DELFIN_TEST_FLAG", mol) is False

    def test_default_one_means_enabled(self):
        # ``default=1`` simulates a flag that ships ON by default
        mol = _mol(_SMILES_FIXTURES["sigma"])
        assert _class_conditional_flag(
            "DELFIN_TEST_FLAG", mol, default=1
        ) is True

    def test_uff_soft_donors_default_off(self):
        """Production wiring: DELFIN_UFF_SOFT_DONORS unset → False on every class."""
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag("DELFIN_UFF_SOFT_DONORS", mol) is False, (
                f"unexpected True on class={label}"
            )

    def test_b5_rigid_h_default_off(self):
        """Production wiring: DELFIN_B5_RIGID_H unset → False on every class."""
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag("DELFIN_B5_RIGID_H", mol) is False, (
                f"unexpected True on class={label}"
            )


# ----------------------------------------------------------------------
# 2. scalar=1 → enabled for every class
# ----------------------------------------------------------------------

class TestScalarOneEnablesEverywhere:
    def test_scalar_one_uff_soft_donors_enables_all_classes(self):
        os.environ["DELFIN_UFF_SOFT_DONORS"] = "1"
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_UFF_SOFT_DONORS", mol
            ) is True, f"scalar=1 must enable class={label}"

    def test_scalar_one_b5_rigid_h_enables_all_classes(self):
        os.environ["DELFIN_B5_RIGID_H"] = "1"
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_B5_RIGID_H", mol
            ) is True, f"scalar=1 must enable class={label}"


# ----------------------------------------------------------------------
# 3. CLASSES="sigma" → enabled only for sigma
# ----------------------------------------------------------------------

class TestSingleClassGate:
    def test_uff_soft_donors_sigma_only(self):
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = "sigma"
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is False
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["no_metal"])
        ) is False

    def test_b5_rigid_h_sigma_only(self):
        os.environ["DELFIN_B5_RIGID_H_CLASSES"] = "sigma"
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is False

    def test_classes_overrides_scalar_zero(self):
        """When ``_CLASSES`` is non-empty the scalar is ignored entirely."""
        os.environ["DELFIN_UFF_SOFT_DONORS"] = "0"
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = "sigma"
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["sigma"])
        ) is True

    def test_classes_overrides_scalar_one(self):
        """``_CLASSES`` precedence still holds when scalar would say True."""
        os.environ["DELFIN_B5_RIGID_H"] = "1"
        os.environ["DELFIN_B5_RIGID_H_CLASSES"] = "sigma"
        # multi_sigma is *not* in the list → must be False even though
        # the scalar alone would have returned True.
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is False


# ----------------------------------------------------------------------
# 4. CLASSES="sigma,multi_sigma" → enabled for both classes
# ----------------------------------------------------------------------

class TestMultipleClassGate:
    def test_uff_soft_donors_sigma_and_multi_sigma(self):
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = "sigma,multi_sigma"
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["no_metal"])
        ) is False

    def test_b5_rigid_h_sigma_and_multi_sigma(self):
        os.environ["DELFIN_B5_RIGID_H_CLASSES"] = "sigma,multi_sigma"
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["no_metal"])
        ) is False


# ----------------------------------------------------------------------
# 5. Unknown class label / mol=None → safe fallback (disabled)
# ----------------------------------------------------------------------

class TestUnknownClassSafeFallback:
    def test_unknown_class_in_env_disables(self):
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = "unicorn"
        # No class label produced by _classify_complex_class matches the bogus
        # token, so every real fixture must come out False — no crash.
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_UFF_SOFT_DONORS", mol
            ) is False, f"unknown class token should disable class={label}"

    def test_mol_none_with_classes_env_disables(self):
        """``_CLASSES`` is non-empty but ``mol`` is None → classifier yields
        'no_metal'; 'no_metal' is not in {'sigma'} so the flag is False.  No
        AttributeError / NoneType crash."""
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = "sigma"
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", None
        ) is False

    def test_mol_none_with_scalar_one_still_enables(self):
        """``mol=None`` falls through to the scalar path → scalar wins."""
        os.environ["DELFIN_B5_RIGID_H"] = "1"
        assert _class_conditional_flag("DELFIN_B5_RIGID_H", None) is True

    def test_class_label_no_metal_when_listed(self):
        """Explicitly listing 'no_metal' enables that class only."""
        os.environ["DELFIN_B5_RIGID_H_CLASSES"] = "no_metal"
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["no_metal"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["sigma"])
        ) is False


# ----------------------------------------------------------------------
# 6. Malformed env → safe fallback (disabled)
# ----------------------------------------------------------------------

class TestMalformedEnvSafeFallback:
    @pytest.mark.parametrize("raw", [
        "   ",          # whitespace only — parses to empty set → scalar path
        ",,,",          # trailing commas only
        ", , , ,",      # whitespace + commas
    ])
    def test_whitespace_or_empty_classes_falls_back_to_scalar(self, raw):
        """``_CLASSES`` that parses to an empty set is treated as unset, so the
        scalar (default 0) path applies → False without raising."""
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = raw
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_UFF_SOFT_DONORS", mol
            ) is False, f"malformed env should fall back on class={label}"

    def test_malformed_classes_with_scalar_one_uses_scalar(self):
        """Empty-parsed CLASSES + scalar=1 → scalar wins (back-compat path)."""
        os.environ["DELFIN_UFF_SOFT_DONORS"] = "1"
        os.environ["DELFIN_UFF_SOFT_DONORS_CLASSES"] = ", ,"
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", _mol(_SMILES_FIXTURES["sigma"])
        ) is True

    def test_classes_with_extra_whitespace_still_parses(self):
        """``_CLASSES`` tokens are stripped — ' sigma , multi_sigma ' must
        still gate on {sigma, multi_sigma}."""
        os.environ["DELFIN_B5_RIGID_H_CLASSES"] = "  sigma ,  multi_sigma  "
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["multi_sigma"])
        ) is True
        assert _class_conditional_flag(
            "DELFIN_B5_RIGID_H", _mol(_SMILES_FIXTURES["no_metal"])
        ) is False

    def test_scalar_non_integer_falls_back_to_default(self):
        """Non-int scalar (``_delfin_env_int`` swallow path) must not crash."""
        os.environ["DELFIN_UFF_SOFT_DONORS"] = "definitely-not-an-int"
        mol = _mol(_SMILES_FIXTURES["sigma"])
        assert _class_conditional_flag(
            "DELFIN_UFF_SOFT_DONORS", mol
        ) is False


# ----------------------------------------------------------------------
# 7. Integration sanity: classification labels we depend on
# ----------------------------------------------------------------------

class TestFixtureLabelsAreStable:
    """If ``_classify_complex_class`` ever changes its labels for the fixture
    SMILES the per-class env-gate would silently break.  Pin the labels."""

    @pytest.mark.parametrize("label,smi", list(_SMILES_FIXTURES.items()))
    def test_fixture_label(self, label, smi):
        assert _classify_complex_class(_mol(smi)) == label
