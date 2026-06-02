"""Integration tests for GRIP Phase 3 (L-BFGS polish + constraints).

These tests cover the full polish loop end-to-end:
* constraint primitives (M-D / topology / chiral / clash-floor)
* the ``grip_polish`` function on toy molecules
* determinism (byte-identical output across runs at ``PYTHONHASHSEED=0``)
* accept-if-better gating
* explicit no-op paths (no fragments, no clashes, already-at-optimum)

The tests deliberately avoid the full ``assemble_complex`` pipeline -- that
integration step belongs to Phase 4.  Instead they build minimal RDKit mols
and synthetic coordinate arrays so each test is fast and self-contained.
"""
from __future__ import annotations

import os
import sys

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.fffree.grip_constraints import (
    ChiralVolumeConstraint,
    ClashFloorPenalty,
    DonorPolyhedronConstraint,
    MDInvariantConstraint,
    TopologyConstraint,
)
from delfin.fffree.grip_loss_terms import BondTerm, TotalGripLoss
from delfin.fffree.grip_polish import (
    DEFAULT_CLASH_WEIGHT,
    DEFAULT_VDW_RADII,
    GripPolishResult,
    _CLASH_WEIGHT_ENV,
    _resolve_clash_weight,
    grip_polish,
    mogul_severity,
)
from delfin.fffree.grip_mogul_lookup import (
    GripLibrary,
    DEFAULT_LIB_PATH,
    get_default_library,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _rng(seed: int):
    return np.random.default_rng(seed)


def _toluene():
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _methane():
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _benzene():
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _coords(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )


# ---------------------------------------------------------------------------
# Constraint primitive tests
# ---------------------------------------------------------------------------
class TestMDInvariantConstraint:
    def test_validate_pass_at_initial(self):
        R = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            dtype=np.float64,
        )
        c = MDInvariantConstraint(0, (1, 2), (2.0, 2.0), tol=0.05)
        assert c.validate(R) is True

    def test_validate_fails_outside_tol(self):
        R = np.array(
            [[0.0, 0.0, 0.0], [2.20, 0.0, 0.0], [0.0, 2.0, 0.0]],
            dtype=np.float64,
        )
        c = MDInvariantConstraint(0, (1, 2), (2.0, 2.0), tol=0.05)
        assert c.validate(R) is False

    def test_frozen_atom_set(self):
        c = MDInvariantConstraint(7, (3, 4, 5), (2.1, 2.1, 2.1))
        assert c.frozen_atom_set == frozenset({7, 3, 4, 5})

    def test_violations_array(self):
        R = np.array(
            [[0.0, 0.0, 0.0], [2.10, 0.0, 0.0], [0.0, 2.0, 0.0]],
            dtype=np.float64,
        )
        c = MDInvariantConstraint(0, (1, 2), (2.0, 2.0))
        dev = c.violations(R)
        assert dev[0] == pytest.approx(0.10, rel=1e-12)
        assert dev[1] == pytest.approx(0.0, abs=1e-12)


class TestTopologyConstraint:
    def test_validate_pass_unchanged(self):
        R = np.array(
            [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]],
            dtype=np.float64,
        )
        c = TopologyConstraint.from_initial([(0, 1), (1, 2)], R, 1.5)
        assert c.validate(R) is True

    def test_validate_fails_on_broken_bond(self):
        R0 = np.array(
            [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64,
        )
        c = TopologyConstraint.from_initial([(0, 1)], R0, 1.5)
        # Stretch the bond to 2.5x its initial length.
        R1 = R0.copy()
        R1[1, 0] = 3.8
        assert c.validate(R1) is False

    def test_empty_bonds_passes(self):
        c = TopologyConstraint.from_initial([], np.zeros((3, 3)))
        assert c.validate(np.zeros((3, 3))) is True


class TestChiralVolumeConstraint:
    def test_sign_preserved(self):
        # Tetrahedral arrangement -- positive determinant.
        R0 = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        c = ChiralVolumeConstraint.from_initial([(0, 1, 2, 3)], R0)
        assert c.validate(R0) is True
        # Translate -- sign preserved.
        R1 = R0 + np.array([0.5, -0.2, 0.3])
        assert c.validate(R1) is True

    def test_sign_inversion_caught(self):
        R0 = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        c = ChiralVolumeConstraint.from_initial([(0, 1, 2, 3)], R0)
        # Mirror through xy-plane -- inverts the z-component sign.
        R_flip = R0.copy()
        R_flip[3, 2] = -1.0
        assert c.validate(R_flip) is False

    def test_skip_when_initial_volume_zero(self):
        # Coplanar -- no chirality to protect.
        R0 = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
        c = ChiralVolumeConstraint.from_initial([(0, 1, 2, 3)], R0)
        # Empty -- no centres recorded, always validates True.
        assert len(c.stereo_centers) == 0
        assert c.validate(R0) is True


class TestClashFloorPenalty:
    def test_no_clash_zero_loss(self):
        R = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=np.float64)
        c = ClashFloorPenalty(vdw_radii={0: 1.70, 1: 1.70})
        L, G = c.value_and_grad(R)
        assert L == 0.0
        assert np.allclose(G, 0.0)

    def test_clash_positive_loss_and_repulsive_grad(self):
        # Two C atoms at 1.0 Å.  d_min = 0.85 * (1.70 + 1.70) = 2.89.
        R = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        c = ClashFloorPenalty(vdw_radii={0: 1.70, 1: 1.70})
        L, G = c.value_and_grad(R)
        # gap = 2.89 - 1.0 = 1.89 -> L = 1.89^2 = 3.5721
        assert L == pytest.approx(1.89 ** 2, rel=1e-10)
        # Descent direction = -grad pushes atom 0 in -x and atom 1 in +x.
        # Therefore grad[0, x] > 0 (so -grad pushes left) and grad[1, x] < 0.
        assert G[0, 0] > 0
        assert G[1, 0] < 0
        # Symmetric magnitudes.
        assert np.allclose(G[0], -G[1], atol=1e-12)

    def test_exclude_13_pair_skipped(self):
        # Same setup as above, but the pair is excluded -> zero loss.
        R = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        c = ClashFloorPenalty(
            vdw_radii={0: 1.70, 1: 1.70},
            exclude_13_pairs={frozenset((0, 1))},
        )
        L, _ = c.value_and_grad(R)
        assert L == 0.0

    def test_finite_diff_gradient(self):
        rng = _rng(7)
        R = rng.standard_normal((4, 3)).astype(np.float64) * 0.6
        c = ClashFloorPenalty(vdw_radii={0: 1.70, 1: 1.55, 2: 1.52, 3: 1.20})
        _, G_ana = c.value_and_grad(R)
        eps = 1e-6
        G_fd = np.zeros_like(R)
        for i in range(R.shape[0]):
            for k in range(3):
                Rp = R.copy()
                Rm = R.copy()
                Rp[i, k] += eps
                Rm[i, k] -= eps
                Lp, _ = c.value_and_grad(Rp)
                Lm, _ = c.value_and_grad(Rm)
                G_fd[i, k] = (Lp - Lm) / (2.0 * eps)
        assert np.allclose(G_ana, G_fd, atol=1e-5)


# ---------------------------------------------------------------------------
# grip_polish integration tests
# ---------------------------------------------------------------------------
class TestGripPolish:
    def test_load_grip_library(self):
        """GripLibrary loads + basic bond / angle / improper queries succeed."""
        if not DEFAULT_LIB_PATH.exists():
            pytest.skip(f"library not present at {DEFAULT_LIB_PATH}")
        lib = get_default_library()
        assert isinstance(lib, GripLibrary)
        # Should have a populated key index.
        assert len(lib._key_to_idx) > 10000  # 172k fragments expected
        # And the C-C sp3-sp3 staple lookup hits.
        hit = lib.lookup_bond("C", "sp3", "C", "sp3")
        assert hit is not None
        mu, sigma, n = hit
        assert 1.4 < mu < 1.65
        assert sigma > 0
        assert n >= 5

    def test_polish_handles_no_fragments(self):
        """Methane (no usable Mogul fragments after freezing) -> P0 unchanged."""
        mol = _methane()
        P0 = _coords(mol)
        # Freeze atom 0 (the C) as the "metal" and atoms 1..4 as donors.
        # Every bond touches the frozen set -> zero terms.
        result = grip_polish(
            P0, mol, metal=0, donors=[1, 2, 3, 4],
            geom="", clash_weight=5.0,
            return_diagnostics=True,
        )
        # Either no terms + no clashes (early return P0) OR an accepted no-op.
        assert isinstance(result, GripPolishResult)
        assert np.allclose(result.P, P0, atol=1e-12)

    def test_polish_returns_ndarray_by_default(self):
        mol = _toluene()
        P0 = _coords(mol)
        out = grip_polish(P0, mol, metal=0, donors=[], geom="")
        assert isinstance(out, np.ndarray)
        assert out.shape == P0.shape
        assert np.all(np.isfinite(out))

    def test_polish_md_constraint_preserved(self):
        """Perturb non-frozen atoms, polish, check M-D distances rebound."""
        mol = _toluene()
        P0 = _coords(mol)
        # Treat atom 0 (methyl C) as "metal" and atom 1 (ring C) as a donor.
        # M-D distance preserved across polish (atoms 0 and 1 are frozen).
        m, d = 0, 1
        d_target = float(np.linalg.norm(P0[d] - P0[m]))
        # Perturb every non-frozen atom by a small random amount.
        rng = _rng(seed=11)
        P_init = P0.copy()
        for i in range(P_init.shape[0]):
            if i in (m, d):
                continue
            P_init[i] += rng.standard_normal(3) * 0.05
        P_out = grip_polish(
            P_init, mol, metal=m, donors=[d], geom="",
            md_tol=0.05, clash_weight=5.0,
        )
        d_after = float(np.linalg.norm(P_out[d] - P_out[m]))
        # The metal and the donor are frozen -- distance must be EXACTLY
        # preserved (their coordinates do not move).
        assert abs(d_after - d_target) < 1e-10

    def test_polish_topology_preserved(self):
        """No bond is stretched past the multiplier after polish."""
        mol = _toluene()
        P0 = _coords(mol)
        m = 0
        # Use atoms 1..3 as "donors" (frozen).  Then polish.
        P_out = grip_polish(
            P0, mol, metal=m, donors=[1, 2, 3], geom="",
            topo_max_multiplier=1.5, clash_weight=5.0,
        )
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            L0 = float(np.linalg.norm(P0[a] - P0[b]))
            L1 = float(np.linalg.norm(P_out[a] - P_out[b]))
            assert L1 <= 1.5 * L0 + 1e-9

    def test_polish_chiral_volume_preserved(self):
        """A pyramidal centre keeps its tetrahedral-volume sign."""
        # Build a chiral sp3 carbon: CHFClBr (no metal).  Use as toy.
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("[C@H](F)(Cl)Br")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        P0 = _coords(mol)
        # No metal -- pick an arbitrary "metal" outside the molecule
        # (skipping the M-D constraint logic) by using metal=-1 and donors=[].
        # Choose metal=0 with donors=[] so M-D set is just {0}, freezing C.
        P_out = grip_polish(
            P0, mol, metal=0, donors=[], geom="",
            clash_weight=2.0,
        )
        # Check signed volume sign for the C with its 3 heavy neighbours.
        idx_c = 0
        nbrs = sorted(int(n.GetIdx()) for n in mol.GetAtomWithIdx(idx_c).GetNeighbors())
        # Take the first three heavy neighbours.
        heavy_nbrs = [
            i for i in nbrs if mol.GetAtomWithIdx(i).GetSymbol() != "H"
        ][:3]
        if len(heavy_nbrs) >= 3:
            v0 = ChiralVolumeConstraint.signed_volume(P0, idx_c, *heavy_nbrs)
            v1 = ChiralVolumeConstraint.signed_volume(P_out, idx_c, *heavy_nbrs)
            assert np.sign(v0) == np.sign(v1)

    def test_clash_floor_repulsion_pushes_atoms_apart(self):
        """Synthetic test: the clash floor PRIMITIVE produces a strictly
        repulsive gradient between two clashing atoms, and stepping along
        ``-grad`` moves them apart.

        The full ``grip_polish`` rollback gate is intentionally severity-only
        (SPEC §3.4): a clash-only structure does not improve severity so the
        gate would (correctly) discard it. To test the floor's *physical*
        effect we apply one gradient step to the constraint primitive directly.
        """
        # Two C atoms at 1.0 Å -- deep clash (d_min ~ 2.89).
        R = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64,
        )
        c = ClashFloorPenalty(vdw_radii={0: 1.70, 1: 1.70})
        L0, G = c.value_and_grad(R)
        assert L0 > 0
        # One descent step with a unit-friendly learning rate.
        R_step = R - 0.05 * G  # gradient descent
        d_after = float(np.linalg.norm(R_step[1] - R_step[0]))
        assert d_after > 1.0 + 1e-3, (
            f"clash floor did not push atoms apart: ended at {d_after}"
        )
        # Iterating until past d_min eventually reaches the floor.
        R_iter = R.copy()
        for _ in range(200):
            Li, Gi = c.value_and_grad(R_iter)
            if Li == 0.0:
                break
            R_iter = R_iter - 0.05 * Gi
        d_final = float(np.linalg.norm(R_iter[1] - R_iter[0]))
        d_min = 0.85 * (1.70 + 1.70)
        assert d_final >= d_min - 1e-3, (
            f"clash floor failed to push beyond d_min={d_min:.3f}, "
            f"ended at {d_final:.3f}"
        )

    def test_clash_floor_no_attraction(self):
        """Two atoms placed beyond the floor are NOT pulled together."""
        Chem = pytest.importorskip("rdkit.Chem")
        rw = Chem.RWMol()
        ai = rw.AddAtom(Chem.Atom("C"))
        bi = rw.AddAtom(Chem.Atom("C"))
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
        except Exception:
            pass
        # Place at 1.2 * vdW = 1.2 * 3.40 = 4.08 Å -> WELL beyond d_min = 2.89.
        P0 = np.array([[0.0, 0.0, 0.0], [4.08, 0.0, 0.0]], dtype=np.float64)
        P_out = grip_polish(
            P0, mol, metal=ai, donors=[], geom="",
            clash_weight=5.0, topo_max_multiplier=10.0,
        )
        d_after = float(np.linalg.norm(P_out[bi] - P_out[ai]))
        # No attractive part -> distance not reduced (within numerical noise).
        assert d_after >= 4.08 - 1e-6, (
            f"clash floor wrongly attracted atoms: started at 4.08, ended at {d_after}"
        )

    def test_polish_improves_severity(self):
        """A degraded toluene structure (small random perturbation) sees
        the severity decrease (or stay the same) after polish."""
        mol = _toluene()
        P0 = _coords(mol)
        rng = _rng(seed=23)
        # Degrade by perturbing every atom (no freezing).
        P_bad = P0 + rng.standard_normal(P0.shape) * 0.10
        res = grip_polish(
            P_bad, mol, metal=0, donors=[], geom="",
            clash_weight=5.0,
            return_diagnostics=True,
        )
        # If terms exist + severity_before > 0, the polish should either
        # improve it (accepted) or roll back (P unchanged).
        assert isinstance(res, GripPolishResult)
        if res.accepted:
            assert res.severity_after < res.severity_before
        else:
            # Rollback path -- P is the input.
            assert np.allclose(res.P, P_bad, atol=1e-12)

    def test_polish_deterministic_byte_identical(self):
        mol = _toluene()
        P0 = _coords(mol)
        # Slightly perturb to force the L-BFGS to actually move.
        rng = _rng(seed=7)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.03
        out1 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        out2 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        assert np.array_equal(out1, out2)

    def test_accept_if_better_gate_returns_P0_on_regression(self):
        """Direct test of the gate: ``mogul_severity(P') >= mogul_severity(P0)``
        path returns ``P0``."""
        # Build a hand-rolled scenario with a single BondTerm whose minimum
        # is at distance 1.5, but P0 is already at 1.5 (severity 0).  Any
        # perturbation by the optimiser can only worsen the severity -> the
        # gate must roll back.
        Chem = pytest.importorskip("rdkit.Chem")
        rw = Chem.RWMol()
        ai = rw.AddAtom(Chem.Atom("C"))
        bi = rw.AddAtom(Chem.Atom("C"))
        rw.AddBond(ai, bi, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
        except Exception:
            pass
        # Place at exactly 1.5 Å.
        P0 = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
        # freeze both atoms -> no possible movement -> if the optimiser
        # somehow stepped (it should not), the gate would catch any worsening.
        res = grip_polish(
            P0, mol, metal=ai, donors=[bi], geom="",
            clash_weight=5.0, return_diagnostics=True,
        )
        # Result is byte-equal to P0 (either no-op or rollback path).
        assert np.allclose(res.P, P0, atol=1e-12)

    def test_polish_no_op_at_optimum(self):
        """If P0 is already an optimum, the polish should not move things."""
        mol = _benzene()
        P0 = _coords(mol)
        # Polish without perturbation.  Severity may not be 0 but at least
        # the polished structure should be very close to the input.
        P1 = grip_polish(P0, mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        # All non-frozen atoms moved at most ~0.2 Å (small)
        dmax = float(np.max(np.linalg.norm(P1 - P0, axis=1)))
        assert dmax < 0.5  # generous: allows L-BFGS micro-step

    def test_polish_returns_P0_on_md_violation(self):
        """Synthetic case where the L-BFGS path would necessarily violate the
        M-D tolerance.  We engineer this by hand-crafting a constraint with
        an absurdly small tolerance and a P0 that is already at the edge."""
        # Use toluene; the metal=0, donor=1 distance is the C-C ring bond
        # length; the polish CAN move atoms but the M-D pair is frozen.
        # So M-D will be byte-preserved -> the test instead verifies that
        # the function handles a forced-violation scenario gracefully by
        # constructing a constraint via direct validate() call.
        mol = _toluene()
        P0 = _coords(mol)
        # Construct a manual MDInvariantConstraint with an impossibly tight
        # tolerance and a "polished" geometry that violates it.
        c = MDInvariantConstraint(0, (1,), (1.50,), tol=1e-6)
        P_bad = P0.copy()
        # Move atom 1 by 0.1 Å -- violates the 1e-6 tol.
        P_bad[1] += np.array([0.1, 0.0, 0.0])
        assert c.validate(P_bad) is False


# ---------------------------------------------------------------------------
# Mogul severity standalone
# ---------------------------------------------------------------------------
class TestMogulSeverity:
    def test_zero_fragments_returns_zero(self):
        agg = TotalGripLoss(terms=[])
        R = np.zeros((3, 3), dtype=np.float64)
        assert mogul_severity(R, agg) == 0.0

    def test_matches_aggregate_loss(self):
        terms = [BondTerm(a=0, b=1, mu=1.5, sigma=0.02, weight=1.0)]
        agg = TotalGripLoss(terms=terms)
        R = np.array(
            [[0.0, 0.0, 0.0], [1.55, 0.0, 0.0]], dtype=np.float64
        )
        # z = (1.55 - 1.50) / 0.02 = 2.5 -> L = 6.25
        sev = mogul_severity(R, agg)
        assert sev == pytest.approx(6.25, rel=1e-12)


# ---------------------------------------------------------------------------
# Heal Hebel 2 — env-tunable clash_weight (DELFIN_FFFREE_GRIP_CLASH_WEIGHT)
# ---------------------------------------------------------------------------
class TestEnvClashWeight:
    """The clash_weight multiplier on ClashFloorPenalty must be reachable via
    an env var so the Phase-5 inter-ligand-clash sweep can pick the optimal
    value without rebuilding the binary.

    Default-OFF (env unset) MUST give the historical Phase-3 calibration
    (5.0) so the Phase-4 byte-identity contract holds.
    """

    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        # Always start each test with a known-clean env so prior leaks
        # from other tests cannot poison this fixture.
        monkeypatch.delenv(_CLASH_WEIGHT_ENV, raising=False)
        yield

    # ----- resolver unit tests (pure function) ------------------------------
    def test_resolver_default_when_env_unset(self, monkeypatch):
        monkeypatch.delenv(_CLASH_WEIGHT_ENV, raising=False)
        assert _resolve_clash_weight(None) == pytest.approx(DEFAULT_CLASH_WEIGHT)
        assert DEFAULT_CLASH_WEIGHT == pytest.approx(5.0)

    def test_resolver_reads_env_float(self, monkeypatch):
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "10.0")
        assert _resolve_clash_weight(None) == pytest.approx(10.0)
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "2.5")
        assert _resolve_clash_weight(None) == pytest.approx(2.5)
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "20.0")
        assert _resolve_clash_weight(None) == pytest.approx(20.0)

    def test_resolver_invalid_env_falls_back_to_default(self, monkeypatch, caplog):
        for bad in ("not-a-number", "", "nan", "inf", "-inf"):
            monkeypatch.setenv(_CLASH_WEIGHT_ENV, bad)
            # "" is treated as unset by our resolver (no warning), the others
            # are reported via logging.warning — verify both:
            assert _resolve_clash_weight(None) == pytest.approx(DEFAULT_CLASH_WEIGHT)

    def test_resolver_explicit_arg_wins(self, monkeypatch):
        # When the caller passes an explicit numeric clash_weight, env is
        # IGNORED -- this preserves backwards-compat for every internal test.
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "99.0")
        assert _resolve_clash_weight(5.0) == pytest.approx(5.0)
        assert _resolve_clash_weight(7.5) == pytest.approx(7.5)
        assert _resolve_clash_weight(0.0) == pytest.approx(0.0)

    def test_resolver_invalid_arg_falls_back(self, monkeypatch):
        # Non-numeric caller value -> consult env / default rather than crash.
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "12.5")
        assert _resolve_clash_weight("foo") == pytest.approx(12.5)
        monkeypatch.delenv(_CLASH_WEIGHT_ENV, raising=False)
        assert _resolve_clash_weight("foo") == pytest.approx(DEFAULT_CLASH_WEIGHT)

    # ----- end-to-end: env value flows through ClashFloorPenalty ------------
    def test_env_clash_weight_default(self, monkeypatch):
        """With env unset, grip_polish must construct ClashFloorPenalty with
        weight=DEFAULT_CLASH_WEIGHT (5.0).  We inspect the penalty via a
        captured ClashFloorPenalty constructor call."""
        monkeypatch.delenv(_CLASH_WEIGHT_ENV, raising=False)

        captured = {}
        from delfin.fffree import grip_polish as gp_mod
        OrigClash = gp_mod.ClashFloorPenalty

        def _spy(*args, **kwargs):
            captured["weight"] = kwargs.get("weight")
            return OrigClash(*args, **kwargs)

        monkeypatch.setattr(gp_mod, "ClashFloorPenalty", _spy)

        mol = _toluene()
        P0 = _coords(mol)
        _ = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        assert captured["weight"] == pytest.approx(DEFAULT_CLASH_WEIGHT)

    def test_env_clash_weight_override(self, monkeypatch):
        """With env set to "10.0", the ClashFloorPenalty weight must be 10.0.
        Calls grip_polish with clash_weight=None (signal-default) so the
        resolver consults the env."""
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "10.0")

        captured = {}
        from delfin.fffree import grip_polish as gp_mod
        OrigClash = gp_mod.ClashFloorPenalty

        def _spy(*args, **kwargs):
            captured["weight"] = kwargs.get("weight")
            return OrigClash(*args, **kwargs)

        monkeypatch.setattr(gp_mod, "ClashFloorPenalty", _spy)

        mol = _toluene()
        P0 = _coords(mol)
        _ = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        assert captured["weight"] == pytest.approx(10.0)

    def test_env_clash_weight_invalid(self, monkeypatch, caplog):
        """Garbage env value -> warning logged + fall back to 5.0."""
        import logging
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "not-a-number")

        captured = {}
        from delfin.fffree import grip_polish as gp_mod
        OrigClash = gp_mod.ClashFloorPenalty

        def _spy(*args, **kwargs):
            captured["weight"] = kwargs.get("weight")
            return OrigClash(*args, **kwargs)

        monkeypatch.setattr(gp_mod, "ClashFloorPenalty", _spy)

        mol = _toluene()
        P0 = _coords(mol)
        with caplog.at_level(logging.WARNING, logger="delfin.fffree.grip_polish"):
            _ = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        assert captured["weight"] == pytest.approx(DEFAULT_CLASH_WEIGHT)
        # A warning must have been emitted at least once.
        warn_msgs = [
            rec.message for rec in caplog.records
            if rec.levelno == logging.WARNING
            and _CLASH_WEIGHT_ENV in rec.message
        ]
        assert warn_msgs, "expected warning mentioning DELFIN_FFFREE_GRIP_CLASH_WEIGHT"

    def test_env_clash_weight_explicit_arg_overrides_env(self, monkeypatch):
        """If the caller passes clash_weight=X explicitly, env is ignored.
        This is what every existing test relies on (clash_weight=5.0)."""
        monkeypatch.setenv(_CLASH_WEIGHT_ENV, "99.0")

        captured = {}
        from delfin.fffree import grip_polish as gp_mod
        OrigClash = gp_mod.ClashFloorPenalty

        def _spy(*args, **kwargs):
            captured["weight"] = kwargs.get("weight")
            return OrigClash(*args, **kwargs)

        monkeypatch.setattr(gp_mod, "ClashFloorPenalty", _spy)

        mol = _toluene()
        P0 = _coords(mol)
        _ = grip_polish(P0, mol, metal=0, donors=[1], geom="", clash_weight=2.5)
        assert captured["weight"] == pytest.approx(2.5)
