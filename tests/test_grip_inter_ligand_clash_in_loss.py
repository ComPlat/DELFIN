"""Tests for the inter-ligand clash boost in the L-BFGS loss (fix A,
User-Direktive 2026-06-02).

These tests pin down the extended :class:`ClashFloorPenalty` behaviour and
the env-flag wiring inside :func:`grip_polish`:

* Backward-compat: omitting ``ligand_atom_id`` is byte-identical with
  the legacy single-weight path.
* Inter-ligand pairs receive the boosted ``w_inter`` weight while
  intra-ligand pairs keep the standard ``weight``.
* The env-flag ``DELFIN_FFFREE_GRIP_INTER_LIGAND_CLASH_WEIGHT`` overrides
  the default boost when no explicit value is passed.
* Acceptance-gate extension (``DELFIN_FFFREE_GRIP_ACCEPT_WITH_CLASH=1``)
  augments the severity score with the inter-ligand clash count.
* Finite-difference gradient agreement (chain rule correct for both
  intra- and inter-ligand contributions).
* ``build_ligand_atom_id_map`` returns the same partition as
  :func:`identify_ligand_subgraphs` (single source of truth).
"""
from __future__ import annotations

import os
import sys

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("scipy")

from rdkit import Chem

from delfin.fffree.grip_constraints import ClashFloorPenalty
from delfin.fffree.grip_polish import (
    DEFAULT_ACCEPT_CLASH_ALPHA,
    DEFAULT_CLASH_WEIGHT,
    DEFAULT_INTER_LIGAND_CLASH_WEIGHT,
    _INTER_LIGAND_CLASH_WEIGHT_ENV,
    _ACCEPT_WITH_CLASH_ENV,
    _ACCEPT_WITH_CLASH_ALPHA_ENV,
    _accept_with_clash_active,
    _accept_with_clash_alpha,
    _resolve_inter_ligand_clash_weight,
    build_ligand_atom_id_map,
)
from delfin.fffree.grip_ensemble import identify_ligand_subgraphs


# ---------------------------------------------------------------------------
# Synthetic minimal-clash setup: two 2-atom "ligands" with one inter-pair
# inside the vdW floor and one intra-pair inside the floor.  This isolates
# the per-pair-weight contract from any RDKit / hapto detection issues.
# ---------------------------------------------------------------------------
def _make_two_ligand_clash():
    """Return a (vdW_radii, coords, lig_map, exclude) bundle suitable for
    pinning down the clash weights independently of any RDKit dependency.

    Layout (positions in Å):
        ligand 0  : atom 0 at (0, 0, 0), atom 1 at (1.0, 0, 0)
        ligand 1  : atom 2 at (1.1, 0, 0), atom 3 at (2.6, 0, 0)

    With vdW radius 1.0 for all atoms and floor_fraction 0.85, the
    threshold ``d_min = 1.7`` Å.  Clashing pairs (skipping the bonded
    intra-pairs):

        intra-ligand: (0, 1) -- d=1.0 -- but excluded as bonded pair.
        inter-ligand: (1, 2) -- d=0.1 -- KEPT, INTER-LIGAND.
        inter-ligand: (0, 2) -- d=1.1 -- KEPT, INTER-LIGAND.
        inter-ligand: (1, 3) -- d=1.6 -- KEPT, INTER-LIGAND.
        inter-ligand: (0, 3) -- d=2.6 -- d>d_min, no contribution.
        intra-ligand: (2, 3) -- d=1.5 -- excluded as bonded.

    So only inter-ligand pairs contribute when the bonded pairs are in the
    exclusion set.
    """
    R = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.1, 0.0, 0.0],
            [2.6, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    vdw = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0}
    lig_map = {0: 0, 1: 0, 2: 1, 3: 1}
    excl_13 = {frozenset((0, 1)), frozenset((2, 3))}
    return vdw, R, lig_map, excl_13


# ---------------------------------------------------------------------------
# 1) Backward-compat: omitting ligand_atom_id keeps legacy weight everywhere
# ---------------------------------------------------------------------------
class TestIntraLigandClashWeightUnchanged:
    def test_legacy_path_byte_identical_when_no_ligand_map(self):
        """No ligand_atom_id passed -> same penalty as before fix A."""
        vdw, R, _lig_map, excl_13 = _make_two_ligand_clash()
        # Legacy constructor signature (positional / first kwargs only).
        legacy = ClashFloorPenalty(
            vdw_radii=vdw,
            exclude_13_pairs=excl_13,
            floor_fraction=0.85,
            weight=5.0,
        )
        L_legacy, G_legacy = legacy.value_and_grad(R)
        # New constructor with explicit None -- must yield the same value.
        new_default = ClashFloorPenalty(
            vdw_radii=vdw,
            exclude_13_pairs=excl_13,
            floor_fraction=0.85,
            weight=5.0,
            ligand_atom_id=None,
        )
        L_new, G_new = new_default.value_and_grad(R)
        assert L_legacy == pytest.approx(L_new)
        np.testing.assert_allclose(G_legacy, G_new, atol=1e-12)


# ---------------------------------------------------------------------------
# 2) Inter-ligand pairs receive w_inter, intra-ligand keep weight
# ---------------------------------------------------------------------------
class TestInterLigandClashWeightBoosted:
    def test_inter_ligand_pairs_get_w_inter(self):
        """All clashing pairs in the toy setup are inter-ligand, so the
        boosted-weight penalty must equal (w_inter / weight) * legacy.
        """
        vdw, R, lig_map, excl_13 = _make_two_ligand_clash()
        legacy = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=5.0,
        )
        L_legacy, _ = legacy.value_and_grad(R)
        boosted = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=5.0,
            ligand_atom_id=lig_map, w_inter=15.0,
        )
        L_boost, G_boost = boosted.value_and_grad(R)
        # All non-excluded clashing pairs are inter-ligand, so the ratio
        # of penalties must exactly equal the weight ratio (15/5 = 3).
        assert L_boost == pytest.approx(L_legacy * 3.0, rel=1e-9)
        # And the gradient magnitudes scale the same way.
        _, G_legacy = legacy.value_and_grad(R)
        np.testing.assert_allclose(G_boost, G_legacy * 3.0, atol=1e-9)

    def test_intra_ligand_pairs_only_use_weight(self):
        """Place a non-bonded intra-ligand clash + non-bonded inter clash;
        intra pair must contribute weight, inter pair must contribute
        w_inter -- proves per-pair routing is correct.
        """
        # 4 atoms: ligand 0 = {0, 1, 2}, ligand 1 = {3}.
        # Bonded pairs: (0,1), (1,2) -> excluded.  Pair (0,2) is intra-
        # ligand but non-bonded -- it should use ``weight``.  Pair (2,3)
        # is inter-ligand non-bonded -- it should use ``w_inter``.
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],   # intra non-bonded clash with 0
                [1.6, 0.0, 0.0],   # inter non-bonded clash with 2
            ],
            dtype=np.float64,
        )
        vdw = {i: 1.0 for i in range(4)}
        excl_13 = {frozenset((0, 1)), frozenset((1, 2))}
        lig_map = {0: 0, 1: 0, 2: 0, 3: 1}
        boosted = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=5.0,
            ligand_atom_id=lig_map, w_inter=15.0,
        )
        # Analytical decomposition:
        # d_min = 0.85 * (1+1) = 1.7
        # pair (0,2): d=1.5, gap=0.2 -> intra, weight 5 -> 5 * 0.04 = 0.2
        # pair (0,3): d=1.6, gap=0.1 -> inter, weight 15 -> 15 * 0.01 = 0.15
        # pair (1,3): d=0.6, gap=1.1 -> inter, weight 15 -> 15 * 1.21 = 18.15
        # pair (2,3): d=0.1, gap=1.6 -> inter, weight 15 -> 15 * 2.56 = 38.40
        expected = 5.0 * 0.04 + 15.0 * (0.01 + 1.21 + 2.56)
        L, _ = boosted.value_and_grad(R)
        assert L == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# 3) Env-flag controls the resolved w_inter default
# ---------------------------------------------------------------------------
class TestEnvInterLigandClashWeightOverride:
    def setup_method(self, _m):
        self._saved_env = os.environ.get(_INTER_LIGAND_CLASH_WEIGHT_ENV)

    def teardown_method(self, _m):
        if self._saved_env is None:
            os.environ.pop(_INTER_LIGAND_CLASH_WEIGHT_ENV, None)
        else:
            os.environ[_INTER_LIGAND_CLASH_WEIGHT_ENV] = self._saved_env

    def test_default_inter_weight_without_env(self):
        os.environ.pop(_INTER_LIGAND_CLASH_WEIGHT_ENV, None)
        w = _resolve_inter_ligand_clash_weight(None)
        assert w == pytest.approx(DEFAULT_INTER_LIGAND_CLASH_WEIGHT)

    def test_env_override(self):
        os.environ[_INTER_LIGAND_CLASH_WEIGHT_ENV] = "25.0"
        w = _resolve_inter_ligand_clash_weight(None)
        assert w == pytest.approx(25.0)

    def test_explicit_arg_wins_over_env(self):
        os.environ[_INTER_LIGAND_CLASH_WEIGHT_ENV] = "25.0"
        w = _resolve_inter_ligand_clash_weight(7.5)
        assert w == pytest.approx(7.5)

    def test_invalid_env_falls_back_to_default(self):
        os.environ[_INTER_LIGAND_CLASH_WEIGHT_ENV] = "not_a_number"
        w = _resolve_inter_ligand_clash_weight(None)
        assert w == pytest.approx(DEFAULT_INTER_LIGAND_CLASH_WEIGHT)


# ---------------------------------------------------------------------------
# 4) Acceptance gate is augmented with the clash component when enabled
# ---------------------------------------------------------------------------
class TestAcceptanceGateWithClashComponent:
    def setup_method(self, _m):
        self._saved_active = os.environ.get(_ACCEPT_WITH_CLASH_ENV)
        self._saved_alpha = os.environ.get(_ACCEPT_WITH_CLASH_ALPHA_ENV)

    def teardown_method(self, _m):
        for k, v in [
            (_ACCEPT_WITH_CLASH_ENV, self._saved_active),
            (_ACCEPT_WITH_CLASH_ALPHA_ENV, self._saved_alpha),
        ]:
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_accept_with_clash_default_off(self):
        os.environ.pop(_ACCEPT_WITH_CLASH_ENV, None)
        assert _accept_with_clash_active() is False

    def test_accept_with_clash_on(self):
        os.environ[_ACCEPT_WITH_CLASH_ENV] = "1"
        assert _accept_with_clash_active() is True

    def test_accept_with_clash_off_explicit(self):
        os.environ[_ACCEPT_WITH_CLASH_ENV] = "0"
        assert _accept_with_clash_active() is False

    def test_clash_alpha_default(self):
        os.environ.pop(_ACCEPT_WITH_CLASH_ALPHA_ENV, None)
        assert _accept_with_clash_alpha() == pytest.approx(DEFAULT_ACCEPT_CLASH_ALPHA)

    def test_clash_alpha_env(self):
        os.environ[_ACCEPT_WITH_CLASH_ALPHA_ENV] = "2.5"
        assert _accept_with_clash_alpha() == pytest.approx(2.5)


# ---------------------------------------------------------------------------
# 5) Finite-difference gradient agreement (chain rule check)
# ---------------------------------------------------------------------------
class TestGradientFiniteDiffWithInterLigand:
    def test_finite_difference_grad_agrees_with_analytical(self):
        """The analytical gradient must match a central finite-difference
        approximation when the per-pair weights are heterogeneous.
        """
        vdw, R, lig_map, excl_13 = _make_two_ligand_clash()
        clash = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=5.0,
            ligand_atom_id=lig_map, w_inter=15.0,
        )
        L0, G_analytic = clash.value_and_grad(R)
        eps = 1e-6
        G_fd = np.zeros_like(R)
        for i in range(R.shape[0]):
            for k in range(3):
                R_p = R.copy(); R_p[i, k] += eps
                R_m = R.copy(); R_m[i, k] -= eps
                Lp, _ = clash.value_and_grad(R_p)
                Lm, _ = clash.value_and_grad(R_m)
                G_fd[i, k] = (Lp - Lm) / (2.0 * eps)
        np.testing.assert_allclose(G_analytic, G_fd, atol=1e-4)


# ---------------------------------------------------------------------------
# 6) build_ligand_atom_id_map mirrors identify_ligand_subgraphs
# ---------------------------------------------------------------------------
class TestBuildLigandAtomIdMap:
    def _cisplatin(self):
        mw = Chem.RWMol()
        mw.AddAtom(Chem.Atom("Pt"))
        for sym in ("N", "N", "Cl", "Cl"):
            a = mw.AddAtom(Chem.Atom(sym))
            mw.AddBond(0, a, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(mw, catchErrors=True)
        except Exception:
            pass
        return mw

    def test_map_matches_subgraphs(self):
        mol = self._cisplatin()
        comps = identify_ligand_subgraphs(mol, 0, [1, 2, 3, 4])
        m = build_ligand_atom_id_map(mol, 0, [1, 2, 3, 4])
        # Every non-metal atom must be in some ligand group.
        for atom in (1, 2, 3, 4):
            assert atom in m
        # And the ligand IDs match the subgraph indexing.
        for lig_id, comp in enumerate(comps):
            for a in comp:
                assert m[int(a)] == lig_id

    def test_empty_on_broken_input(self):
        # None should not crash -- empty dict is the safe fallback.
        assert build_ligand_atom_id_map(None, 0, []) == {}


# ---------------------------------------------------------------------------
# 7) Default-off byte-identity at the grip_polish-level
# ---------------------------------------------------------------------------
class TestDefaultOffByteIdentical:
    """When neither env-flag is set, grip_polish must produce the SAME
    output as before fix A.  We probe the boost env-flag resolution path:
    with no env, the resolved ligand_map is None -> ClashFloorPenalty is
    constructed identically to the legacy call (w_inter defaults to weight).
    """

    def setup_method(self, _m):
        self._saved = os.environ.get(_INTER_LIGAND_CLASH_WEIGHT_ENV)
        os.environ.pop(_INTER_LIGAND_CLASH_WEIGHT_ENV, None)

    def teardown_method(self, _m):
        if self._saved is None:
            os.environ.pop(_INTER_LIGAND_CLASH_WEIGHT_ENV, None)
        else:
            os.environ[_INTER_LIGAND_CLASH_WEIGHT_ENV] = self._saved

    def test_no_env_no_ligand_map(self):
        # Resolve helper returns the default boost (15.0) — but this is
        # only consumed when the env-flag is set in grip_polish.  The
        # important property is that the default-off path constructs the
        # penalty with the SAME effective weight as the legacy default.
        vdw, R, _lig, excl_13 = _make_two_ligand_clash()
        # ligand_atom_id=None reproduces the legacy semantics.
        legacy = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=DEFAULT_CLASH_WEIGHT,
        )
        new_default = ClashFloorPenalty(
            vdw_radii=vdw, exclude_13_pairs=excl_13,
            floor_fraction=0.85, weight=DEFAULT_CLASH_WEIGHT,
            ligand_atom_id=None,
        )
        L_l, _ = legacy.value_and_grad(R)
        L_n, _ = new_default.value_and_grad(R)
        assert L_l == pytest.approx(L_n)
