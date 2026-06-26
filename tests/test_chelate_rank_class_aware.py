"""Tests for the class-aware chelate-rank patch (ITER-chelate_rank_class_aware).

The patch adds a default-OFF env-flag
``DELFIN_CHELATE_RANK_CLASS_AWARE`` that augments the chelate-conformer
ranking inside ``_chelate_conformer_candidates`` with a per-class
secondary score:

  composite = delta + alpha * (element_weight + spread_factor * pucker)

where

  * ``element_weight``  — donor-element-table value averaged across the
                          fragment's donor atoms (constant per fragment),
  * ``pucker``          — per-conformer RMS donor out-of-plane deviation,
  * ``spread_factor``   — per-class multiplier (negative = reward pucker).

The composite re-orders the survivor list when the env-flag is on, and
is bit-exact with HEAD otherwise.  This file covers the unit-level
contract of the new public surface:

  1. Env flag default OFF — helper untouched, weights table populated.
  2. Element weights table coverage — every class has a sane dict.
  3. ``_chelate_class_donor_penalty`` — returns 0 for unknown class,
     empty donor list, ``mol=None``, ``no_metal`` class; returns the
     mean of weights for known donors; ignores unknown donor elements.
  4. ``_class_conditional_flag`` integration — per-class env gate works.
  5. End-to-end re-ranking — the composite key changes the order of
     synthetic ``(delta, coords, pucker)`` accepted tuples consistent
     with the flag and class.

The integration test uses a tiny stub of the sort loop because invoking
the full ``_chelate_conformer_candidates`` requires RDKit ETKDG which is
heavy and stochastic; the math under test is the composite-key sort.
"""
from __future__ import annotations

import os
from typing import List, Tuple

import pytest

pytest.importorskip("rdkit", reason="RDKit required for chelate-rank tests")
from rdkit import Chem  # noqa: E402

from delfin.smiles_converter import (  # noqa: E402
    DELFIN_CHELATE_RANK_CLASS_ALPHA,
    DELFIN_CHELATE_RANK_CLASS_AWARE,
    _CHELATE_CLASS_DONOR_WEIGHTS,
    _chelate_class_donor_penalty,
    _class_conditional_flag,
    _classify_complex_class,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

_SMILES_FIXTURES = {
    "sigma":       "[Pt](Cl)(Cl)(N)N",
    "multi_sigma": "[Mn]([C]#O)([C]#O)[Mn]([C]#O)([C]#O)",
    "no_metal":    "CC(=O)O",
}


def _mol(smi: str):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    assert m is not None, f"SMILES did not parse: {smi}"
    return m


@pytest.fixture(autouse=True)
def _clean_env():
    targets = [
        "DELFIN_CHELATE_RANK_CLASS_AWARE",
        "DELFIN_CHELATE_RANK_CLASS_AWARE_CLASSES",
        "DELFIN_CHELATE_RANK_CLASS_ALPHA",
    ]
    for n in targets:
        os.environ.pop(n, None)
    yield
    for n in targets:
        os.environ.pop(n, None)


# ----------------------------------------------------------------------
# 1. Module-level constants — sane defaults
# ----------------------------------------------------------------------

class TestDefaultsAndTableShape:
    """The default state must be the bit-exact HEAD path."""

    def test_default_flag_is_off(self):
        assert DELFIN_CHELATE_RANK_CLASS_AWARE == 0

    def test_default_alpha_is_small_positive(self):
        # Alpha must be small enough that the secondary score only
        # reorders near-ties.  >0 so the score actually mixes; <0.1 so
        # large delta differences always dominate.
        assert 0.0 < DELFIN_CHELATE_RANK_CLASS_ALPHA < 0.10

    @pytest.mark.parametrize("cls", [
        "sigma", "hapto", "multi_sigma", "multi_hapto", "no_metal",
    ])
    def test_every_class_present_in_weights_table(self, cls):
        assert cls in _CHELATE_CLASS_DONOR_WEIGHTS

    def test_no_metal_weights_empty(self):
        """no_metal must never re-rank — empty dict guarantees that."""
        assert _CHELATE_CLASS_DONOR_WEIGHTS["no_metal"] == {}

    def test_weight_magnitudes_are_bounded(self):
        """Per-donor weights must be bounded so the per-fragment mean
        stays well within ``alpha`` after multiplication."""
        for cls, table in _CHELATE_CLASS_DONOR_WEIGHTS.items():
            for elem, w in table.items():
                assert -1.0 < w < 1.0, (
                    f"class={cls} elem={elem} weight={w} out of bounds"
                )

    def test_sigma_table_prefers_n_over_p(self):
        """Sigma class: N donor must outrank P donor (lower weight)."""
        sig = _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]
        assert sig["N"] < sig["P"]

    def test_hapto_table_prefers_c_donor(self):
        """Hapto class: carbon (η-anchor) must score better than P/Cl."""
        hap = _CHELATE_CLASS_DONOR_WEIGHTS["hapto"]
        assert hap["C"] < hap["P"]
        assert hap["C"] < hap["Cl"]


# ----------------------------------------------------------------------
# 2. _chelate_class_donor_penalty — pure-function contract
# ----------------------------------------------------------------------

class TestDonorPenaltyHelper:
    def test_empty_donor_list_zero(self):
        mol = _mol(_SMILES_FIXTURES["sigma"])
        assert _chelate_class_donor_penalty(mol, [], "sigma") == 0.0

    def test_mol_none_zero(self):
        assert _chelate_class_donor_penalty(None, [0, 1], "sigma") == 0.0

    def test_unknown_class_zero(self):
        mol = _mol(_SMILES_FIXTURES["sigma"])
        # Donor indices that exist in mol — but class label is unknown.
        donors = [
            i for i, a in enumerate(mol.GetAtoms())
            if a.GetSymbol() in ("N", "Cl")
        ]
        assert _chelate_class_donor_penalty(mol, donors, "unicorn") == 0.0

    def test_no_metal_class_zero(self):
        mol = _mol(_SMILES_FIXTURES["sigma"])
        donors = [
            i for i, a in enumerate(mol.GetAtoms())
            if a.GetSymbol() == "N"
        ]
        # ``no_metal`` table is empty → penalty stays 0 regardless of
        # which donor atoms are passed in.
        assert _chelate_class_donor_penalty(mol, donors, "no_metal") == 0.0

    def test_mean_of_known_weights(self):
        """Two donors of the same element → penalty == that element's
        weight (mean of identical values)."""
        mol = _mol(_SMILES_FIXTURES["sigma"])
        n_donors = [
            i for i, a in enumerate(mol.GetAtoms())
            if a.GetSymbol() == "N"
        ]
        assert len(n_donors) >= 2
        pen = _chelate_class_donor_penalty(mol, n_donors[:2], "sigma")
        expected = _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]["N"]
        assert pen == pytest.approx(expected, abs=1e-9)

    def test_mixed_donor_set_mean(self):
        """Mixed N/Cl donor set → penalty == arithmetic mean of the two
        element weights."""
        mol = _mol(_SMILES_FIXTURES["sigma"])
        idx_by_sym = {sym: [] for sym in ("N", "Cl")}
        for i, a in enumerate(mol.GetAtoms()):
            if a.GetSymbol() in idx_by_sym:
                idx_by_sym[a.GetSymbol()].append(i)
        donors = [idx_by_sym["N"][0], idx_by_sym["Cl"][0]]
        pen = _chelate_class_donor_penalty(mol, donors, "sigma")
        wN = _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]["N"]
        wCl = _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]["Cl"]
        assert pen == pytest.approx((wN + wCl) / 2.0, abs=1e-9)

    def test_unknown_element_counts_neutral(self):
        """An exotic donor element (Se) contributes 0 to the sum but is
        counted, so the mean is reduced in magnitude."""
        # Build a tiny mol with N (known) + Se (unknown to sigma table).
        m = Chem.MolFromSmiles("N[Se]", sanitize=False)
        assert m is not None
        donors = list(range(m.GetNumAtoms()))
        pen = _chelate_class_donor_penalty(m, donors, "sigma")
        wN = _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]["N"]
        # Mean of (wN, 0.0) == wN / 2
        assert pen == pytest.approx(wN / 2.0, abs=1e-9)

    def test_invalid_donor_index_skipped(self):
        """Out-of-range donor indices are silently skipped (no crash)."""
        mol = _mol(_SMILES_FIXTURES["sigma"])
        # Mix a valid N donor with a garbage index.
        n_idx = next(
            i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "N"
        )
        pen = _chelate_class_donor_penalty(
            mol, [n_idx, 9999], "sigma"
        )
        # Garbage index ignored → penalty == weight(N) (one donor counted).
        assert pen == pytest.approx(
            _CHELATE_CLASS_DONOR_WEIGHTS["sigma"]["N"], abs=1e-9
        )


# ----------------------------------------------------------------------
# 3. _class_conditional_flag integration with the new env-flag
# ----------------------------------------------------------------------

class TestClassConditionalFlagIntegration:
    def test_flag_default_off_every_class(self):
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_CHELATE_RANK_CLASS_AWARE", mol
            ) is False, f"unexpected True on class={label}"

    def test_scalar_one_enables_every_class(self):
        os.environ["DELFIN_CHELATE_RANK_CLASS_AWARE"] = "1"
        for label, smi in _SMILES_FIXTURES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_CHELATE_RANK_CLASS_AWARE", mol
            ) is True, f"scalar=1 must enable class={label}"

    def test_classes_env_gates_only_listed_class(self):
        os.environ["DELFIN_CHELATE_RANK_CLASS_AWARE_CLASSES"] = "sigma"
        assert _class_conditional_flag(
            "DELFIN_CHELATE_RANK_CLASS_AWARE",
            _mol(_SMILES_FIXTURES["sigma"]),
        ) is True
        assert _class_conditional_flag(
            "DELFIN_CHELATE_RANK_CLASS_AWARE",
            _mol(_SMILES_FIXTURES["multi_sigma"]),
        ) is False


# ----------------------------------------------------------------------
# 4. End-to-end composite-key re-ranking (math under test)
# ----------------------------------------------------------------------

def _composite_sort(
    accepted: List[Tuple[float, str, float]],
    *,
    enabled: bool,
    penalty_const: float,
    spread_factor: float,
    alpha: float = DELFIN_CHELATE_RANK_CLASS_ALPHA,
) -> List[str]:
    """Mirror of the sort logic in ``_chelate_conformer_candidates``.

    ``accepted`` is a list of ``(delta, tag, pucker)``; the ``tag``
    stands in for the coords-map dict so the test can assert order.
    Returns the ordered list of tags.
    """
    if enabled:
        sorted_ = sorted(
            accepted,
            key=lambda item: (
                item[0]
                + alpha * (penalty_const + spread_factor * item[2])
            ),
        )
    else:
        sorted_ = sorted(accepted, key=lambda item: item[0])
    return [tag for _d, tag, _p in sorted_]


class TestCompositeReranking:
    def test_disabled_path_uses_delta_only(self):
        """When ``enabled=False`` the order must be by ``delta`` only —
        bit-exact HEAD."""
        cands = [
            (0.05, "A", 0.30),  # better fit, less pucker
            (0.08, "B", 0.50),  # worse fit, more pucker
            (0.03, "C", 0.10),  # best fit, least pucker
        ]
        order = _composite_sort(
            cands, enabled=False,
            penalty_const=0.0, spread_factor=0.0,
        )
        assert order == ["C", "A", "B"]

    def test_enabled_zero_penalty_zero_factor_matches_disabled(self):
        cands = [
            (0.05, "A", 0.30),
            (0.08, "B", 0.50),
            (0.03, "C", 0.10),
        ]
        order_off = _composite_sort(
            cands, enabled=False,
            penalty_const=0.0, spread_factor=0.0,
        )
        order_on = _composite_sort(
            cands, enabled=True,
            penalty_const=0.0, spread_factor=0.0,
        )
        assert order_off == order_on

    def test_sigma_rewards_pucker_when_deltas_close(self):
        """Sigma class: spread_factor=-0.50 should promote the more
        puckered conformer when delta values are within alpha."""
        # delta gap = 0.005 (< alpha=0.04) so secondary term can flip
        # the order if pucker contributions differ enough.
        cands = [
            (0.030, "low_pucker",  0.05),  # best delta but flat
            (0.035, "high_pucker", 0.50),  # slightly worse delta, puckered
        ]
        order = _composite_sort(
            cands, enabled=True,
            penalty_const=0.0, spread_factor=-0.50,
        )
        # composite_low  = 0.030 + 0.04*(0 + (-0.50)*0.05) = 0.029
        # composite_high = 0.035 + 0.04*(0 + (-0.50)*0.50) = 0.025
        # high_pucker wins.
        assert order[0] == "high_pucker"

    def test_large_delta_difference_dominates(self):
        """Even with strong pucker reward, a clearly better delta wins —
        the secondary term must NOT flip far-apart fits."""
        cands = [
            (0.01, "tight_flat",   0.05),  # very tight, flat
            (0.20, "loose_pucker", 0.80),  # very loose, very puckered
        ]
        order = _composite_sort(
            cands, enabled=True,
            penalty_const=0.0, spread_factor=-0.50,
        )
        # composite_tight = 0.01 + 0.04*(-0.025) = 0.009
        # composite_loose = 0.20 + 0.04*(-0.40)  = 0.184
        assert order[0] == "tight_flat"

    def test_hapto_neutral_factor_minimal_reorder(self):
        """Hapto class: spread_factor=-0.10 — much weaker reward than
        sigma; same ``cands`` as the sigma-rewards test must stay in
        delta order."""
        cands = [
            (0.030, "low_pucker",  0.05),
            (0.035, "high_pucker", 0.50),
        ]
        order = _composite_sort(
            cands, enabled=True,
            penalty_const=0.0, spread_factor=-0.10,
        )
        # composite_low  = 0.030 + 0.04*(-0.005) = 0.0298
        # composite_high = 0.035 + 0.04*(-0.050) = 0.0330
        # low_pucker (= the original delta-best) wins.
        assert order[0] == "low_pucker"

    def test_constant_penalty_alone_does_not_reorder(self):
        """The element_weight term is per-fragment-constant, so a non-zero
        penalty_const with spread_factor=0 must NOT reorder the list."""
        cands = [
            (0.05, "A", 0.30),
            (0.03, "B", 0.10),
            (0.08, "C", 0.50),
        ]
        order_off = _composite_sort(
            cands, enabled=False,
            penalty_const=0.0, spread_factor=0.0,
        )
        order_on = _composite_sort(
            cands, enabled=True,
            penalty_const=-0.4, spread_factor=0.0,
        )
        # Adding a constant to every composite key cannot change order.
        assert order_off == order_on


# ----------------------------------------------------------------------
# 5. Integration with _chelate_conformer_candidates signature (sanity)
# ----------------------------------------------------------------------

class TestChelateConformerCandidatesUnchangedDefault:
    """The function signature must be preserved (no new required args).
    Default-OFF path must remain bit-exact: this test passes a fragment
    that triggers the early-return (no donors) and asserts an empty
    list — proves the new code does not crash before the early-out.
    """

    def test_empty_donor_list_returns_empty(self):
        from delfin.smiles_converter import _chelate_conformer_candidates

        mol = _mol(_SMILES_FIXTURES["sigma"])
        # No donors at all → early return [].
        result = _chelate_conformer_candidates(
            mol,
            frag_atom_indices=set(range(mol.GetNumAtoms())),
            donor_atom_indices=[],
            target_bite=2.0,
        )
        assert result == []

    def test_single_donor_returns_empty(self):
        from delfin.smiles_converter import _chelate_conformer_candidates

        mol = _mol(_SMILES_FIXTURES["sigma"])
        n_idx = next(
            i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "N"
        )
        result = _chelate_conformer_candidates(
            mol,
            frag_atom_indices=set(range(mol.GetNumAtoms())),
            donor_atom_indices=[n_idx],
            target_bite=2.0,
        )
        # Only one donor → < 2 → early return [].
        assert result == []

    def test_enabled_flag_does_not_crash_on_empty_donors(self):
        """Even with the flag ON the early-return for n_donors<2 must
        execute without computing pucker / penalty (no IndexError)."""
        os.environ["DELFIN_CHELATE_RANK_CLASS_AWARE"] = "1"
        from delfin.smiles_converter import _chelate_conformer_candidates

        mol = _mol(_SMILES_FIXTURES["sigma"])
        result = _chelate_conformer_candidates(
            mol,
            frag_atom_indices=set(range(mol.GetNumAtoms())),
            donor_atom_indices=[],
            target_bite=2.0,
        )
        assert result == []
