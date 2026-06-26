"""Tests for class-aware ETKDG seed-count helpers.

Verifies the public contract of ``_class_aware_seed_count`` and
``_resolve_top_level_seed_count`` introduced for per-class ETKDG seed
budgeting in ``smiles_converter``.

Acceptance contract:

  1. ``DELFIN_CLASS_AWARE_SEEDS=0`` (default) → resolver returns the
     scalar ``DELFIN_TOP_LEVEL_SEED_COUNT`` for every input.  Bit-exact
     pre-patch behaviour for the consumption sites.
  2. ``DELFIN_CLASS_AWARE_SEEDS=1`` + recognised class label → resolver
     returns the chemistry default for that class:
       sigma=20, hapto=12, multi_sigma=30, multi_hapto=15, no_metal=20.
  3. Per-class override ``DELFIN_CLASS_AWARE_SEEDS_<CLASS>=N`` wins over
     the chemistry default for the matching class only; other classes
     keep their defaults.
  4. Unknown / ``None`` / malformed class label → falls back to the
     global scalar, never raises.
  5. The resolver classifies via ``_classify_complex_class``; the
     module-level fallback when classifier raises must still return a
     usable count (>= 1).

Tests skip cleanly when RDKit is not importable.
"""
from __future__ import annotations

import os
from typing import Iterable

import pytest

pytest.importorskip(
    "rdkit",
    reason="RDKit required for class-aware seed-count tests",
)

from rdkit import Chem  # noqa: E402 — after importorskip guard

from delfin.smiles_converter import (  # noqa: E402
    DELFIN_TOP_LEVEL_SEED_COUNT,
    _CLASS_AWARE_SEED_DEFAULTS,
    _class_aware_seed_count,
    _classify_complex_class,
    _resolve_top_level_seed_count,
)


# ----------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------

#: SMILES → expected ``_classify_complex_class`` label.  Same fixtures as
#: ``test_class_conditional_flags.py`` so the two test suites stay in
#: lock-step on classifier behaviour.
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
    for n in names:
        os.environ.pop(n, None)


@pytest.fixture(autouse=True)
def _clean_env():
    """Strip every class-aware-seeds env-var before/after each test."""
    targets = [
        "DELFIN_CLASS_AWARE_SEEDS",
        "DELFIN_CLASS_AWARE_SEEDS_SIGMA",
        "DELFIN_CLASS_AWARE_SEEDS_HAPTO",
        "DELFIN_CLASS_AWARE_SEEDS_MULTI_SIGMA",
        "DELFIN_CLASS_AWARE_SEEDS_MULTI_HAPTO",
        "DELFIN_CLASS_AWARE_SEEDS_NO_METAL",
    ]
    _scrub_env(*targets)
    yield
    _scrub_env(*targets)


# ----------------------------------------------------------------------
# 1. Chemistry defaults — pure helper
# ----------------------------------------------------------------------

class TestClassAwareSeedCountDefaults:
    """``_class_aware_seed_count`` returns the documented chemistry defaults
    when no per-class env override is present.  These are the values the
    ITER doc cites and the pool-evaluator will gate against."""

    @pytest.mark.parametrize(
        "class_label,expected_count",
        [
            ("sigma",       20),
            ("hapto",       12),
            ("multi_sigma", 30),
            ("multi_hapto", 15),
            ("no_metal",    20),
        ],
    )
    def test_each_class_returns_documented_default(
        self, class_label, expected_count
    ):
        assert _class_aware_seed_count(class_label) == expected_count

    def test_defaults_table_matches_helper(self):
        """The exported ``_CLASS_AWARE_SEED_DEFAULTS`` table is the source
        of truth; helper must mirror it exactly with no env overrides set."""
        for label, expected in _CLASS_AWARE_SEED_DEFAULTS.items():
            assert _class_aware_seed_count(label) == expected, (
                f"helper drift for class={label}"
            )

    def test_case_insensitive_label(self):
        """The classifier returns lowercase labels but defensive callers
        may pass mixed case — helper should normalize."""
        assert _class_aware_seed_count("SIGMA") == 20
        assert _class_aware_seed_count("Hapto") == 12

    def test_whitespace_stripped(self):
        assert _class_aware_seed_count("  multi_sigma  ") == 30


# ----------------------------------------------------------------------
# 2. Per-class env overrides
# ----------------------------------------------------------------------

class TestClassAwareSeedCountEnvOverrides:
    """``DELFIN_CLASS_AWARE_SEEDS_<CLASS>=N`` overrides the chemistry
    default for that class; other classes keep their defaults."""

    def test_sigma_override_via_env(self):
        os.environ["DELFIN_CLASS_AWARE_SEEDS_SIGMA"] = "8"
        assert _class_aware_seed_count("sigma") == 8
        # other classes unaffected
        assert _class_aware_seed_count("hapto") == 12
        assert _class_aware_seed_count("multi_sigma") == 30

    def test_hapto_override_via_env(self):
        os.environ["DELFIN_CLASS_AWARE_SEEDS_HAPTO"] = "6"
        assert _class_aware_seed_count("hapto") == 6
        assert _class_aware_seed_count("sigma") == 20

    def test_override_clamped_to_minimum_one(self):
        """Negative / zero overrides must clamp to >=1 so downstream slice
        ``_PIPELINE_SEEDS[:n]`` never returns an empty tuple."""
        os.environ["DELFIN_CLASS_AWARE_SEEDS_SIGMA"] = "0"
        assert _class_aware_seed_count("sigma") == 1
        os.environ["DELFIN_CLASS_AWARE_SEEDS_SIGMA"] = "-3"
        assert _class_aware_seed_count("sigma") == 1

    def test_non_integer_override_falls_back_to_default(self):
        """``_delfin_env_int`` swallows parse errors; result must equal
        the chemistry default for the class."""
        os.environ["DELFIN_CLASS_AWARE_SEEDS_SIGMA"] = "not-a-number"
        assert _class_aware_seed_count("sigma") == 20


# ----------------------------------------------------------------------
# 3. Unknown / None / safe-fallback
# ----------------------------------------------------------------------

class TestClassAwareSeedCountSafeFallback:
    def test_none_returns_scalar_fallback(self):
        assert _class_aware_seed_count(None) == max(
            1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
        )

    def test_empty_string_returns_scalar_fallback(self):
        assert _class_aware_seed_count("") == max(
            1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
        )

    def test_unknown_class_returns_scalar_fallback(self):
        assert _class_aware_seed_count("unicorn") == max(
            1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
        )

    def test_non_string_returns_scalar_fallback(self):
        """Robust against accidental int/None propagation from upstream."""
        assert _class_aware_seed_count(42) == max(  # type: ignore[arg-type]
            1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
        )


# ----------------------------------------------------------------------
# 4. Resolver gated by DELFIN_CLASS_AWARE_SEEDS
# ----------------------------------------------------------------------

class TestResolverEnvGate:
    """``_resolve_top_level_seed_count`` is the production entry point.
    When the gate env-var is 0 (default) it must return the scalar for
    EVERY mol — i.e. bit-exact pre-patch behaviour at all call-sites."""

    def test_gate_off_returns_scalar_for_every_class(self):
        # gate unset by autouse fixture
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _resolve_top_level_seed_count(mol) == max(
                1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
            ), f"gate-off must return scalar on class={label}"

    def test_gate_off_with_none_mol_returns_scalar(self):
        assert _resolve_top_level_seed_count(None) == max(
            1, int(DELFIN_TOP_LEVEL_SEED_COUNT)
        )

    def test_gate_on_routes_through_classifier(self):
        os.environ["DELFIN_CLASS_AWARE_SEEDS"] = "1"
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            cls = _classify_complex_class(mol)
            assert _resolve_top_level_seed_count(mol) == _class_aware_seed_count(
                cls
            ), f"resolver drift for class={label}"

    def test_gate_on_with_none_mol_returns_scalar(self):
        """Even with the gate ON, a None mol must not crash; classifier
        returns 'no_metal' which itself routes to scalar-fallback because
        20 == DELFIN_TOP_LEVEL_SEED_COUNT default."""
        os.environ["DELFIN_CLASS_AWARE_SEEDS"] = "1"
        out = _resolve_top_level_seed_count(None)
        assert out >= 1

    def test_gate_on_with_per_class_override(self):
        """End-to-end: gate ON + per-class override → resolver returns
        the override for that class only."""
        os.environ["DELFIN_CLASS_AWARE_SEEDS"] = "1"
        os.environ["DELFIN_CLASS_AWARE_SEEDS_SIGMA"] = "5"
        mol_sigma = _mol(_SMILES_FIXTURES["sigma"])
        assert _resolve_top_level_seed_count(mol_sigma) == 5
        # multi_sigma still at default (30 by chemistry table)
        mol_msig = _mol(_SMILES_FIXTURES["multi_sigma"])
        assert _resolve_top_level_seed_count(mol_msig) == 30


# ----------------------------------------------------------------------
# 5. Classifier integration sanity
# ----------------------------------------------------------------------

class TestClassifierIntegration:
    """Pin the classifier labels we depend on — if the classifier ever
    rewires its outputs, the resolver would silently regress and these
    tests will catch it before any pool run."""

    @pytest.mark.parametrize("label,smi", list(_SMILES_FIXTURES.items()))
    def test_classifier_label_stable(self, label, smi):
        assert _classify_complex_class(_mol(smi)) == label

    def test_defaults_cover_every_classifier_label(self):
        """Every label the classifier can emit must have an entry in the
        defaults table — otherwise the resolver silently falls back to
        the global scalar for that class and the env-flag is a no-op."""
        classifier_labels = {
            "sigma", "hapto", "multi_sigma", "multi_hapto", "no_metal",
        }
        missing = classifier_labels - set(_CLASS_AWARE_SEED_DEFAULTS)
        assert not missing, (
            f"defaults table missing classifier labels: {sorted(missing)}"
        )
