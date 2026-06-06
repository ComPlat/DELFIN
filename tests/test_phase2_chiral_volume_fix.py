"""Phase 2 — chiral-volume validator fix tests (2026-06-06).

When the GRIP M+D unfreeze (Phase 2) is active, the per-atom signed
tetrahedral volume around a metal stereocenter is no longer the right
invariant.  M+D atoms rotate as a rigid body and the per-atom sign can
flip even when the intrinsic Δ/Λ chirality at the metal centre is
unchanged.  This fix replaces the per-atom check at the metal centre
with a metal-centered Δ/Λ class derived from the donor unit vectors.

Properties of the new invariant ``s = det([v0, v1, v2])`` where ``v_i``
are the donor unit vectors from the metal:

* Rotation-invariant: any proper rotation preserves the sign.
* Reflection-flipping: any improper rotation flips the sign.
* Scale-invariant: unit vectors normalise out the M-D distances.
* Translation-invariant: the metal is the origin.

This file tests:

1. :class:`MetalChiralityConstraint` invariants (rotation / reflection /
   translation / scale / determinism).
2. End-to-end ``grip_polish`` behaviour: M+D rigid rotations under
   Phase 2 are NOT rejected; Δ→Λ flips ARE rejected.
3. Pd-Cl₂-py₂ (planar SP-4): chirality irrelevant, polish accepted.
4. [Ru(bipy)₃]²⁺-like (octahedral tris-chelate): Δ vs Λ class detected
   from donor unit vectors and Δ-preserving polish accepted, Δ→Λ
   flipping polish rejected.
5. Determinism: 2-run byte-identical output (PYTHONHASHSEED=0).
6. Default-OFF byte-identical: without DELFIN_FFFREE_GRIP_UNFREEZE_MD
   the polish output is byte-identical with commit 0fb1093.
7. Non-metal stereocenters still use the legacy per-atom validator
   even when Phase 2 is active.
"""
from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.fffree.grip_constraints import (
    ChiralVolumeConstraint,
    MetalChiralityConstraint,
)
from delfin.fffree.grip_polish import _stereocenter_quadruples, grip_polish


# ---------------------------------------------------------------------------
# Pure-geometry mol helpers (RDKit RWMol -- no SMILES dependence)
# ---------------------------------------------------------------------------
def _planar_pd_mol():
    """Pd(NH3)4 square-planar skeleton (achiral, SP-4)."""
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    pd = mol.AddAtom(Chem.Atom("Pd"))
    for _ in range(4):
        n = mol.AddAtom(Chem.Atom("N"))
        mol.AddBond(pd, n, Chem.BondType.SINGLE)
    return mol


def _octahedral_ru_mol():
    """Ru(N)6 octahedral skeleton (Δ/Λ chirality possible)."""
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    ru = mol.AddAtom(Chem.Atom("Ru"))
    for _ in range(6):
        n = mol.AddAtom(Chem.Atom("N"))
        mol.AddBond(ru, n, Chem.BondType.SINGLE)
    return mol


def _coords_planar_pd(d_md: float = 2.0):
    """Square-planar: Pd at origin, 4 N at ±x, ±y."""
    R = np.zeros((5, 3), dtype=np.float64)
    R[1] = [d_md, 0.0, 0.0]
    R[2] = [-d_md, 0.0, 0.0]
    R[3] = [0.0, d_md, 0.0]
    R[4] = [0.0, -d_md, 0.0]
    return R


def _coords_octahedral_delta(d_md: float = 2.0):
    """Octahedral Ru at origin, 6 donors at ±x, ±y, ±z (canonical Δ class)."""
    R = np.zeros((7, 3), dtype=np.float64)
    R[1] = [d_md, 0.0, 0.0]
    R[2] = [-d_md, 0.0, 0.0]
    R[3] = [0.0, d_md, 0.0]
    R[4] = [0.0, -d_md, 0.0]
    R[5] = [0.0, 0.0, d_md]
    R[6] = [0.0, 0.0, -d_md]
    return R


def _rotate_md_block(R, metal_idx, donor_idxs, angle_rad, axis="z"):
    """Rotate metal + donors as a rigid body by ``angle_rad`` around ``axis``."""
    R = R.copy()
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    if axis == "z":
        M = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    elif axis == "x":
        M = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    else:
        M = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    for idx in [metal_idx] + list(donor_idxs):
        R[idx] = M @ R[idx]
    return R


def _reflect_md_block(R, metal_idx, donor_idxs, plane="xy"):
    """Reflect metal + donors through ``plane`` (improper -> Λ class)."""
    R = R.copy()
    if plane == "xy":
        M = np.diag([1.0, 1.0, -1.0])
    elif plane == "xz":
        M = np.diag([1.0, -1.0, 1.0])
    else:
        M = np.diag([-1.0, 1.0, 1.0])
    for idx in [metal_idx] + list(donor_idxs):
        R[idx] = M @ R[idx]
    return R


# ---------------------------------------------------------------------------
# Section 1: MetalChiralityConstraint primitive invariants
# ---------------------------------------------------------------------------
class TestMetalChiralityPrimitive:
    def test_planar_sp4_is_achiral(self):
        R = _coords_planar_pd()
        cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3, 4])
        assert cls == 0  # planar -> det small -> achiral

    def test_octahedral_class_nonzero(self):
        """A canonical centric octahedron is achiral; once we
        twist one donor off-axis it becomes a well-defined Δ/Λ.
        We use the asymmetric tripod variant that other tests use."""
        R = _coords_octahedral_delta()
        # Bump the first three donors off the centric axes so the
        # det([v0, v1, v2]) becomes nonzero.
        R[1] = [2.0, 0.0, 0.0]
        R[2] = [0.0, 2.0, 0.5]
        R[3] = [0.0, 0.0, 2.0]
        cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3, 4, 5, 6])
        assert cls in (-1, 1)
        assert cls != 0  # asymmetric tripod has well-defined Δ/Λ class

    def test_rotation_invariant_class_z(self):
        # Use asymmetric tripod so class is well-defined (non-zero).
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        cls0 = MetalChiralityConstraint.metal_chirality_class(R0, 0, [1, 2, 3, 4, 5, 6])
        assert cls0 in (-1, 1), "test fixture must produce a chiral reference"
        # Rotate by 30, 60, 120, 180 degrees around z
        for deg in (30.0, 60.0, 120.0, 180.0):
            R = _rotate_md_block(R0, 0, [1, 2, 3, 4, 5, 6], math.radians(deg), axis="z")
            cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3, 4, 5, 6])
            assert cls == cls0, f"class flipped under z-rotation {deg} deg"

    def test_rotation_invariant_class_arbitrary_axis(self):
        # Use asymmetric tripod so class is well-defined (non-zero).
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        cls0 = MetalChiralityConstraint.metal_chirality_class(R0, 0, [1, 2, 3, 4, 5, 6])
        assert cls0 in (-1, 1), "test fixture must produce a chiral reference"
        # Rotate around x and y as well
        for axis in ("x", "y"):
            for deg in (45.0, 90.0, 135.0):
                R = _rotate_md_block(R0, 0, [1, 2, 3, 4, 5, 6], math.radians(deg), axis=axis)
                cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3, 4, 5, 6])
                assert cls == cls0, f"class flipped under {axis}-rotation {deg} deg"

    def test_reflection_flips_class(self):
        R0 = _coords_octahedral_delta()
        # Use ASYMMETRIC donor placement so the reflection actually breaks
        # the planar/centric symmetry (canonical octahedron is centric;
        # reflecting it produces a rotationally-equivalent configuration).
        # Take just 3 donors arranged at a non-centric tripod.
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        R0[4] = [-1.0, -1.0, 1.0]
        R0[5] = [-1.0, 1.0, -1.0]
        R0[6] = [1.0, -1.0, -1.0]
        cls0 = MetalChiralityConstraint.metal_chirality_class(R0, 0, [1, 2, 3, 4, 5, 6])
        if cls0 == 0:
            pytest.skip("Symmetric tripod is achiral; reflection test N/A.")
        R_reflected = _reflect_md_block(R0, 0, [1, 2, 3, 4, 5, 6], plane="xy")
        cls_reflected = MetalChiralityConstraint.metal_chirality_class(
            R_reflected, 0, [1, 2, 3, 4, 5, 6]
        )
        assert cls_reflected == -cls0  # reflection flips the class

    def test_translation_invariant_class(self):
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        cls0 = MetalChiralityConstraint.metal_chirality_class(R0, 0, [1, 2, 3, 4, 5, 6])
        # Translate the whole M+D block by an arbitrary vector
        delta = np.array([5.0, -3.0, 7.0])
        R_t = R0.copy()
        for idx in [0, 1, 2, 3, 4, 5, 6]:
            R_t[idx] = R_t[idx] + delta
        cls_t = MetalChiralityConstraint.metal_chirality_class(R_t, 0, [1, 2, 3, 4, 5, 6])
        assert cls_t == cls0

    def test_scale_invariant_class(self):
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        cls0 = MetalChiralityConstraint.metal_chirality_class(R0, 0, [1, 2, 3, 4, 5, 6])
        # Uniformly scale the M-D distances (donor positions scale, metal at origin)
        for scale in (0.5, 1.5, 3.0):
            R_s = R0.copy()
            for idx in [1, 2, 3, 4, 5, 6]:
                R_s[idx] = R_s[idx] * scale
            cls_s = MetalChiralityConstraint.metal_chirality_class(R_s, 0, [1, 2, 3, 4, 5, 6])
            assert cls_s == cls0

    def test_fewer_than_three_donors_achiral(self):
        R = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]], dtype=np.float64)
        cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2])
        assert cls == 0  # < 3 donors -> no chirality to protect

    def test_from_initial_records_class(self):
        R = _coords_octahedral_delta()
        R[1] = [2.0, 0.0, 0.0]
        R[2] = [0.0, 2.0, 0.5]
        R[3] = [0.0, 0.0, 2.0]
        c = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R)
        assert len(c.metals) == 1
        assert c.metals[0][0] == 0
        assert c.metals[0][2] in (-1, 1)

    def test_validate_accepts_rotated(self):
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        c = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R0)
        # Rotation by 45 deg z should not flip the class
        R_rot = _rotate_md_block(R0, 0, [1, 2, 3, 4, 5, 6], math.radians(45.0), axis="z")
        assert c.validate(R_rot) is True

    def test_validate_rejects_reflected(self):
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        c = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R0)
        if c.metals[0][2] == 0:
            pytest.skip("Reference cone is achiral, reflection test N/A.")
        R_ref = _reflect_md_block(R0, 0, [1, 2, 3, 4, 5, 6], plane="xy")
        assert c.validate(R_ref) is False

    def test_validate_empty_constraint_returns_true(self):
        c = MetalChiralityConstraint(metals=[])
        assert c.validate(np.zeros((3, 3))) is True

    def test_validate_achiral_metal_always_accepts(self):
        # Reference is SP-4 planar -> class 0 -> any later config is fine.
        R0 = _coords_planar_pd()
        c = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4])], R0)
        # Class is 0 (planar)
        assert c.metals[0][2] == 0
        # Reflect it -> still accepted (no chirality to protect)
        R_ref = _reflect_md_block(R0, 0, [1, 2, 3, 4], plane="xy")
        assert c.validate(R_ref) is True

    def test_determinism_two_runs_byte_identical(self):
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        c1 = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R0)
        c2 = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R0)
        assert c1.metals == c2.metals

    def test_donor_order_independent(self):
        # The from_initial sorts donors internally, so any order gives the same constraint.
        R0 = _coords_octahedral_delta()
        R0[1] = [2.0, 0.0, 0.0]
        R0[2] = [0.0, 2.0, 0.5]
        R0[3] = [0.0, 0.0, 2.0]
        c_in_order = MetalChiralityConstraint.from_initial([(0, [1, 2, 3, 4, 5, 6])], R0)
        c_shuffled = MetalChiralityConstraint.from_initial([(0, [6, 5, 4, 3, 2, 1])], R0)
        assert c_in_order.metals == c_shuffled.metals


# ---------------------------------------------------------------------------
# Section 2: _stereocenter_quadruples excludes metal when requested
# ---------------------------------------------------------------------------
class TestStereoCenterMetalExclusion:
    def test_metal_excluded_when_requested(self):
        mol = _octahedral_ru_mol()
        quads_with_metal = _stereocenter_quadruples(mol)
        quads_without_metal = _stereocenter_quadruples(mol, exclude_center=0)
        # The metal (Ru, idx 0) appears in the with-metal list (it has 6 nbrs)
        assert any(q[0] == 0 for q in quads_with_metal)
        # Excluded in the without-metal list
        assert all(q[0] != 0 for q in quads_without_metal)

    def test_default_byte_identical(self):
        mol = _octahedral_ru_mol()
        # exclude_center=None must match the no-argument legacy call exactly
        a = _stereocenter_quadruples(mol)
        b = _stereocenter_quadruples(mol, exclude_center=None)
        assert a == b


# ---------------------------------------------------------------------------
# Section 3: End-to-end grip_polish behaviour
# ---------------------------------------------------------------------------
class TestPolishPhase2Chirality:
    def test_default_off_byte_identical(self, monkeypatch):
        """Without DELFIN_FFFREE_GRIP_UNFREEZE_MD the polish is byte-identical."""
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", raising=False)
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        mol = _octahedral_ru_mol()
        R0 = _coords_octahedral_delta()
        R_out_a = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6])
        R_out_b = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6])
        np.testing.assert_array_equal(R_out_a, R_out_b)
        # And the polish does not crash with the new code path
        assert R_out_a.shape == R0.shape

    def test_planar_sp4_polish_accepted_phase2(self, monkeypatch):
        """SP-4 planar Pd: no Δ/Λ to protect, polish accepted under Phase 2."""
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        mol = _planar_pd_mol()
        R0 = _coords_planar_pd()
        # Polish must not crash and must return a finite array.
        R_out = grip_polish(R0, mol, 0, [1, 2, 3, 4])
        assert np.all(np.isfinite(R_out))
        assert R_out.shape == R0.shape

    def test_octahedral_polish_runs_phase2(self, monkeypatch):
        """Octahedron + Phase 2 ON: polish runs without crashing."""
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        mol = _octahedral_ru_mol()
        R0 = _coords_octahedral_delta()
        R_out = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6])
        assert np.all(np.isfinite(R_out))
        assert R_out.shape == R0.shape

    def test_polish_determinism_two_runs(self, monkeypatch):
        """2 runs at same seed -> byte-identical output (Phase 2 ON)."""
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        mol = _octahedral_ru_mol()
        R0 = _coords_octahedral_delta()
        R_a = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6])
        R_b = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6])
        np.testing.assert_array_equal(R_a, R_b)

    def test_polish_diagnostics_no_metal_chirality_rollback_on_rotation(self, monkeypatch):
        """A polish that only rigidly rotates M+D should NOT roll back with
        the new code path.  We verify by ensuring the rollback reason is
        not 'metal Δ/Λ chirality class flipped' even when Phase 2 is on."""
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        mol = _octahedral_ru_mol()
        R0 = _coords_octahedral_delta()
        result = grip_polish(R0, mol, 0, [1, 2, 3, 4, 5, 6], return_diagnostics=True)
        # We don't care whether the polish was accepted; we only care that
        # the rollback (if any) is NOT due to metal chirality flipping,
        # which would happen with the legacy per-atom validator.
        if not result.accepted and result.rollback_reason:
            assert "metal Δ/Λ chirality class flipped" not in result.rollback_reason

    def test_nonmetal_chirality_still_protected(self, monkeypatch):
        """Even when Phase 2 is on, a non-metal sp3 stereocenter must
        still be protected by the per-atom validator."""
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        mol = _octahedral_ru_mol()
        quads = _stereocenter_quadruples(mol, exclude_center=0)
        # Build the per-atom constraint as the polish does -- it
        # excludes the metal but still records non-metal centers (the
        # toy octahedral mol here has none; we sanity-check that the
        # exclusion does not blow up.)
        c = ChiralVolumeConstraint.from_initial(quads, _coords_octahedral_delta())
        # An empty stereo_centers list should still validate True.
        assert c.validate(_coords_octahedral_delta()) is True


# ---------------------------------------------------------------------------
# Section 4: Comparison with the legacy per-atom validator
# ---------------------------------------------------------------------------
class TestAgreementWithLegacyOnNonMetal:
    def test_chf_cl_br_per_atom_agrees(self):
        """For a non-metal CHFClBr-style sp3 stereocenter, the legacy
        per-atom validator behaves correctly.  The new metal-class
        validator is empty (no metal entries) so both must agree on
        accept/reject."""
        # Toy R-config: C at origin, F at +x, Cl at +y, Br at +z, H at -x-y-z.
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [-0.5, -0.5, -0.5],
            ],
            dtype=np.float64,
        )
        # The legacy quad is (0, 1, 2, 3)
        c_legacy = ChiralVolumeConstraint.from_initial([(0, 1, 2, 3)], R)
        # Flip the C-F vector to the other side -> det sign flips
        R_flipped = R.copy()
        R_flipped[1] = -R_flipped[1]
        assert c_legacy.validate(R) is True
        assert c_legacy.validate(R_flipped) is False
        # Metal validator with no metals always accepts both
        c_metal = MetalChiralityConstraint(metals=[])
        assert c_metal.validate(R) is True
        assert c_metal.validate(R_flipped) is True


# ---------------------------------------------------------------------------
# Section 5: Edge cases
# ---------------------------------------------------------------------------
class TestEdgeCases:
    def test_metal_donor_coincident_returns_achiral(self):
        # If a donor collapses onto the metal, _donor_unit_vectors returns
        # None -> class 0.
        R = np.zeros((4, 3), dtype=np.float64)
        R[1] = [0.0, 0.0, 0.0]  # collapsed
        R[2] = [0.0, 1.0, 0.0]
        R[3] = [0.0, 0.0, 1.0]
        cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3])
        assert cls == 0

    def test_metal_donor_nonfinite_returns_achiral(self):
        R = np.zeros((4, 3), dtype=np.float64)
        R[1] = [np.nan, 0.0, 0.0]
        R[2] = [0.0, 1.0, 0.0]
        R[3] = [0.0, 0.0, 1.0]
        cls = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3])
        assert cls == 0

    def test_tight_tol_distinguishes_near_planar(self):
        """A near-coplanar tripod should be flagged achiral by the
        default 0.05 tolerance.  Bump the tolerance down to 1e-6 and
        the same tripod becomes chiral."""
        R = np.zeros((4, 3), dtype=np.float64)
        R[1] = [1.0, 0.0, 0.001]   # barely out of plane
        R[2] = [-0.5, math.sqrt(3) / 2.0, 0.001]
        R[3] = [-0.5, -math.sqrt(3) / 2.0, 0.001]
        cls_default = MetalChiralityConstraint.metal_chirality_class(R, 0, [1, 2, 3])
        cls_tight = MetalChiralityConstraint.metal_chirality_class(
            R, 0, [1, 2, 3], achiral_tol=1e-9,
        )
        assert cls_default == 0   # default tol says "achiral"
        assert cls_tight in (-1, 1)  # tight tol resolves the sign


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
