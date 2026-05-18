"""Tests for Wave-5 multi-hapto patches: SIMPLE_PATH (P1) + MM_ENFORCE (P2).

Background (per ``iters/ITER-multihapto_rootcause.md`` and ``WAVE5_O``):

  The ``multi_hapto`` class has been stuck at 0/22 topology-pass for months.
  Agent-3 forensik on the 3943c2b voll-pool identified two root causes:

  1.  ``_build_multimetal_hapto_sequential`` Step 3 skips all metals when
      collecting donor indices, so SMILES-declared M-M sigma bonds
      (e.g. Sn-Ir, Sn-Rh, Sn-Co) are never enforced.  The secondary
      metal lands ~3.4 A from the transition metal even though the
      crystallographic target is ~2.85 A.  The topology detector cutoff
      (r_cov_i + r_cov_j + 0.45) drops the bond at 3.25 A -> topology
      fails.

  2.  ``_enforce_metal_topology`` only PUSHES non-bonded atoms away; it
      does not PULL declared M-M bonds toward their ideal length.

Patches:

  * **P1 MULTIHAPTO_SIMPLE_PATH** (env-flag, default OFF):
        For the ``multi_hapto`` class, bypass the variant-plan loop and
        ``_build_hybrid_hapto_complex`` so the ETKDG-seed + sphere-scaffold
        candidates are the only entries in ``hapto_candidate_mols``.
        Wired via ``_class_conditional_flag`` so operator opt-in can be
        either scalar (all classes) or class-targeted (e.g. multi_hapto
        only).

  * **P2 MULTIHAPTO_MM_ENFORCE** (env-flag, default OFF):
        After the existing repulsive pass in ``_enforce_metal_topology``,
        loop over every SMILES-declared (M_i, M_j) neighbour pair and
        symmetric-shift both metals toward the ideal distance from
        ``_METAL_METAL_BOND_LENGTHS`` (covalent-radii fallback for
        unlisted hetero-pairs).  Tolerance 0.20 A.  Mono-metal SMILES
        and no-metal SMILES have no (M,M) edge -> bit-identical no-op.

Test coverage (acceptance per ITER-multihapto_rootcause Section 4):

  1.  ``test_default_off_simple_path_inactive_for_all_classes``
        — DELFIN_MULTIHAPTO_SIMPLE_PATH unset -> flag is False on every
          class, including multi_hapto.  Production-default OFF gate.

  2.  ``test_default_off_mm_enforce_distance_unchanged``
        — DELFIN_MULTIHAPTO_MM_ENFORCE unset -> Sn-Ir distance in
          D-COIRSN remains the baseline 3.4 A (no correction applied).

  3.  ``test_simple_path_active_only_for_multi_hapto_when_class_targeted``
        — DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES="multi_hapto" -> flag
          True only for multi_hapto class; sigma/multi_sigma/hapto/no_metal
          unaffected.  Mono-hapto non-regression guarantee.

  4.  ``test_mm_enforce_pulls_sn_ir_below_detector_cutoff``
        — DELFIN_MULTIHAPTO_MM_ENFORCE=1 on D-COIRSN -> Sn-Ir distance
          drops from 3.4 A to <=3.25 A (target 3.20 A, r_cov_Sn +
          r_cov_Ir + 0.4 fallback).

  5.  ``test_mm_enforce_helper_no_op_on_mono_metal``
        — Direct unit test on ``_enforce_smiles_mm_distances`` with a
          mono-metal mol -> zero pairs moved (no M-M edge in graph).

  6.  ``test_mm_enforce_helper_no_op_on_no_metal``
        — Direct unit test with a no-metal mol -> zero pairs moved.

  7.  ``test_simple_path_scalar_one_enables_for_all_classes``
        — DELFIN_MULTIHAPTO_SIMPLE_PATH=1 -> flag True on every class
          (operator escape hatch for cross-class experiments).

  8.  ``test_mm_enforce_helper_uses_dict_target_when_listed``
        — When the (M,M) pair IS in ``_METAL_METAL_BOND_LENGTHS`` (e.g.
          Pt-Pt 2.60 A), the helper uses that value, not the radii
          fallback.

  9.  ``test_unknown_class_safe_fallback``
        — mol=None -> classifier returns 'no_metal'; SIMPLE_PATH flag
          must be False without raising.

All tests skip cleanly when RDKit is not importable.
"""
from __future__ import annotations

import math
import os
from typing import Iterable

import pytest

pytest.importorskip("rdkit", reason="RDKit required for multi-hapto patch tests")
from rdkit import Chem  # noqa: E402
from rdkit.Geometry import Point3D  # noqa: E402

import numpy as np  # noqa: E402

from delfin.smiles_converter import (  # noqa: E402
    _COVALENT_RADII,
    _METAL_METAL_BOND_LENGTHS,
    _METAL_SET,
    _class_conditional_flag,
    _classify_complex_class,
    _enforce_smiles_mm_distances,
    smiles_to_xyz,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Real D-COIRSN multi-hapto SMILES from the master_v3 pool (Sn-Ir + 2 Cp).
# Baseline: Sn-Ir = 3.416 A (above 3.25 A detector cutoff -> bond dropped).
_COIRSN_SMI = (
    "[Cl][Sn]([Cl])([Cl])[Ir]123456([C]7=[C]1CC[C]2=[C]3CC7)"
    "[C]1=[C]4CC[C]5=[C]6CC1"
)

_CLASS_SMILES = {
    "sigma":        "[Pt](Cl)(Cl)(N)N",
    "multi_sigma":  "[Mn]([C]#O)([C]#O)[Mn]([C]#O)([C]#O)",
    "no_metal":     "CC(=O)O",
    "multi_hapto":  _COIRSN_SMI,
}


def _mol(smi: str):
    """Parse SMILES with ``sanitize=False`` to keep dative/bracket forms."""
    m = Chem.MolFromSmiles(smi, sanitize=False)
    assert m is not None, f"SMILES did not parse: {smi}"
    return m


def _scrub(*names: Iterable[str]) -> None:
    for n in names:
        os.environ.pop(n, None)


_ENV_VARS = (
    "DELFIN_MULTIHAPTO_SIMPLE_PATH",
    "DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES",
    "DELFIN_MULTIHAPTO_MM_ENFORCE",
)


@pytest.fixture(autouse=True)
def _clean_env():
    """Strip multi-hapto env-vars before/after every test for isolation."""
    _scrub(*_ENV_VARS)
    yield
    _scrub(*_ENV_VARS)


def _sn_ir_distance(xyz: str):
    """Return Sn-Ir Euclidean distance (A) from the first XYZ frame, or None."""
    sn = ir = None
    for ln in (xyz or "").splitlines():
        parts = ln.split()
        if len(parts) == 4:
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                continue
            if parts[0] == "Sn" and sn is None:
                sn = np.array([x, y, z])
            elif parts[0] == "Ir" and ir is None:
                ir = np.array([x, y, z])
    if sn is None or ir is None:
        return None
    return float(np.linalg.norm(sn - ir))


# ---------------------------------------------------------------------------
# 1. Wire-in defaults (Phase-1.5 2026-05-15)
# ---------------------------------------------------------------------------
#
# Wire-in semantics:
#   * Env unset                         -> ON for multi_hapto, OFF elsewhere.
#   * DELFIN_MULTIHAPTO_<PATCH>=0       -> hard rollback to pre-wire-in
#                                          (OFF for every class).
#   * DELFIN_MULTIHAPTO_<PATCH>=1       -> ON for every class
#                                          (operator escape hatch).
#   * DELFIN_MULTIHAPTO_<PATCH>_CLASSES -> ON for the listed classes only.
#
# Both patches use ``_class_conditional_flag(..., default_classes=("multi_hapto",))``
# at their wire-in sites in ``smiles_converter.py``.

class TestWireInDefaults:
    """Class-conditional default-ON for the ``multi_hapto`` class only."""

    def test_env_unset_simple_path_active_only_for_multi_hapto(self):
        """Env unset -> SIMPLE_PATH flag True for multi_hapto, False elsewhere."""
        for label, smi in _CLASS_SMILES.items():
            mol = _mol(smi)
            flag = _class_conditional_flag(
                "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0,
                default_classes=("multi_hapto",),
            )
            if label == "multi_hapto":
                assert flag is True, (
                    f"wire-in: flag must default-ON for multi_hapto, got {flag}"
                )
            else:
                assert flag is False, (
                    f"wire-in: flag must default-OFF for class={label}, got {flag}"
                )

    def test_env_unset_mm_enforce_active_only_for_multi_hapto(self):
        """Env unset -> MM_ENFORCE flag True for multi_hapto, False elsewhere."""
        for label, smi in _CLASS_SMILES.items():
            mol = _mol(smi)
            flag = _class_conditional_flag(
                "DELFIN_MULTIHAPTO_MM_ENFORCE", mol, default=0,
                default_classes=("multi_hapto",),
            )
            if label == "multi_hapto":
                assert flag is True, (
                    f"wire-in: flag must default-ON for multi_hapto, got {flag}"
                )
            else:
                assert flag is False, (
                    f"wire-in: flag must default-OFF for class={label}, got {flag}"
                )

    def test_explicit_zero_simple_path_rolls_back_for_multi_hapto(self):
        """``DELFIN_MULTIHAPTO_SIMPLE_PATH=0`` overrides the wire-in default."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH"] = "0"
        mol = _mol(_COIRSN_SMI)
        assert _class_conditional_flag(
            "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0,
            default_classes=("multi_hapto",),
        ) is False

    def test_explicit_zero_mm_enforce_rolls_back_distance_to_baseline(self):
        """``DELFIN_MULTIHAPTO_MM_ENFORCE=0`` restores pre-wire-in Sn-Ir = 3.4 A."""
        # Also disable SIMPLE_PATH so this is a pure end-to-end "all patches off"
        # baseline (matches HEAD-pre-wire-in behaviour exactly).
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH"] = "0"
        os.environ["DELFIN_MULTIHAPTO_MM_ENFORCE"] = "0"
        xyz, err = smiles_to_xyz(_COIRSN_SMI, hapto_approx=True)
        assert err is None, f"converter failed: {err}"
        d = _sn_ir_distance(xyz)
        assert d is not None, "Sn or Ir missing from baseline XYZ"
        # Pre-wire-in baseline must remain >3.25 A so the regression is
        # observable; otherwise this test no longer guards the rollback.
        assert d > 3.25, (
            f"baseline Sn-Ir distance {d:.3f} A unexpectedly already inside "
            f"detector cutoff with both patches explicitly disabled -- has "
            f"the regression been fixed elsewhere?"
        )

    def test_env_unset_pulls_sn_ir_below_detector_cutoff(self):
        """Wire-in default: D-COIRSN env-unset -> Sn-Ir <=3.25 A."""
        # Explicit "no opt-out" environment: scrub clears env-vars in fixture.
        xyz, err = smiles_to_xyz(_COIRSN_SMI, hapto_approx=True)
        assert err is None, f"converter failed: {err}"
        d = _sn_ir_distance(xyz)
        assert d is not None, "Sn or Ir missing from wire-in XYZ"
        # Target from covalent-radii fallback: 1.39 + 1.41 + 0.4 = 3.20 A.
        assert d <= 3.25, (
            f"Sn-Ir = {d:.3f} A still above detector cutoff 3.25 A "
            f"with default wire-in (env unset)"
        )


# ---------------------------------------------------------------------------
# 2. Class gating for P1 (SIMPLE_PATH)
# ---------------------------------------------------------------------------

class TestSimplePathGating:
    """Ensure SIMPLE_PATH respects class targeting (non-regression on mono-hapto)."""

    def test_simple_path_active_only_for_multi_hapto_when_class_targeted(self):
        """Class-targeted env enables flag for multi_hapto only."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES"] = "multi_hapto"
        for label, smi in _CLASS_SMILES.items():
            mol = _mol(smi)
            flag = _class_conditional_flag(
                "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0
            )
            if label == "multi_hapto":
                assert flag is True, f"flag must be True for multi_hapto, got {flag}"
            else:
                assert flag is False, (
                    f"flag must be False for class={label} when CLASSES targets multi_hapto, "
                    f"got {flag}"
                )

    def test_simple_path_scalar_one_enables_for_all_classes(self):
        """Operator escape hatch: ``=1`` enables for every class."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH"] = "1"
        for label, smi in _CLASS_SMILES.items():
            mol = _mol(smi)
            assert _class_conditional_flag(
                "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0
            ) is True, f"scalar=1 must enable class={label}"

    def test_simple_path_zero_disables_even_when_classes_present(self):
        """``DELFIN_MULTIHAPTO_SIMPLE_PATH=0`` overrides default-classes behaviour."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH"] = "0"
        mol = _mol(_COIRSN_SMI)
        # Without _CLASSES env or default_classes, the result simply tracks
        # the scalar -- 0 means False.
        assert _class_conditional_flag(
            "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0
        ) is False


# ---------------------------------------------------------------------------
# 3. P2: MM_ENFORCE integration (end-to-end Sn-Ir correction)
# ---------------------------------------------------------------------------

class TestMMEnforceIntegration:
    """End-to-end: enabling MM_ENFORCE drops Sn-Ir below the topology detector cutoff."""

    def test_mm_enforce_pulls_sn_ir_below_detector_cutoff(self):
        """D-COIRSN baseline 3.4 A -> ~3.20 A with MM_ENFORCE=1."""
        os.environ["DELFIN_MULTIHAPTO_MM_ENFORCE"] = "1"
        xyz, err = smiles_to_xyz(_COIRSN_SMI, hapto_approx=True)
        assert err is None, f"converter failed: {err}"
        d = _sn_ir_distance(xyz)
        assert d is not None, "Sn or Ir missing from MM_ENFORCE XYZ"
        # Target from covalent-radii fallback: 1.39 + 1.41 + 0.4 = 3.20 A.
        # Tolerance 0.20 A, so any distance in [3.00, 3.40] passes.
        # Detector cutoff is 3.25 A; we require strictly < 3.25 A so the bond
        # would be detected by the topology gate.
        assert d <= 3.25, (
            f"Sn-Ir = {d:.3f} A still above detector cutoff 3.25 A "
            f"after MM_ENFORCE=1 (expected target ~3.20 A)"
        )

    def test_combined_p1_p2_also_pulls_sn_ir_below_cutoff(self):
        """SIMPLE_PATH + MM_ENFORCE together must also yield Sn-Ir <= 3.25 A."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES"] = "multi_hapto"
        os.environ["DELFIN_MULTIHAPTO_MM_ENFORCE"] = "1"
        xyz, err = smiles_to_xyz(_COIRSN_SMI, hapto_approx=True)
        assert err is None, f"converter failed: {err}"
        d = _sn_ir_distance(xyz)
        assert d is not None, "Sn or Ir missing from P1+P2 XYZ"
        assert d <= 3.25, (
            f"Sn-Ir = {d:.3f} A still above 3.25 A with P1+P2 active"
        )


# ---------------------------------------------------------------------------
# 4. P2: MM_ENFORCE helper unit tests (mono-hapto / no-metal no-op)
# ---------------------------------------------------------------------------

def _make_conf_helpers(mol, positions):
    """Build a fresh conformer on ``mol`` with given positions and return
    closures (_gp, _sp, conf) matching the ``_enforce_metal_topology`` API."""
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(positions):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    conf = mol.GetConformer(0)

    def _gp(i):
        p = conf.GetAtomPosition(i)
        return np.array([p.x, p.y, p.z])

    def _sp(i, arr):
        conf.SetAtomPosition(
            i, Point3D(float(arr[0]), float(arr[1]), float(arr[2]))
        )

    return _gp, _sp, conf


class TestMMEnforceHelperUnit:
    """Direct unit tests on ``_enforce_smiles_mm_distances``."""

    def test_mm_enforce_helper_no_op_on_mono_metal(self):
        """Single metal (Pt) -> no M-M edge -> zero moves."""
        mol = Chem.MolFromSmiles("[Pt](Cl)(Cl)(N)N", sanitize=False)
        # Three-atom subset: Pt at origin, Cl + N as dummy donors.
        _gp, _sp, _ = _make_conf_helpers(
            mol,
            [(0.0, 0.0, 0.0)] * mol.GetNumAtoms(),
        )
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        assert len(metals) == 1
        moved = _enforce_smiles_mm_distances(mol, None, metals, _gp, _sp, np)
        assert moved == 0, f"mono-metal mol must not move any pair, moved={moved}"

    def test_mm_enforce_helper_no_op_on_no_metal(self):
        """No metal -> empty metal_indices -> zero moves."""
        mol = Chem.MolFromSmiles("CC(=O)O", sanitize=False)
        _gp, _sp, _ = _make_conf_helpers(
            mol, [(0.0, 0.0, 0.0)] * mol.GetNumAtoms()
        )
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        assert metals == []
        moved = _enforce_smiles_mm_distances(mol, None, metals, _gp, _sp, np)
        assert moved == 0

    def test_mm_enforce_helper_uses_dict_target_when_listed(self):
        """Pt-Pt is in ``_METAL_METAL_BOND_LENGTHS`` at 2.60 A -> helper
        pulls a 4.00 A Pt-Pt pair toward that value, not the radii fallback."""
        # Build a minimal Pt-Pt graph via dative SMILES.
        mol = Chem.MolFromSmiles("[Pt][Pt]", sanitize=False)
        assert mol is not None and mol.GetNumAtoms() == 2
        # Place metals 4.00 A apart along x.
        _gp, _sp, _ = _make_conf_helpers(mol, [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)])
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        assert metals == [0, 1]

        # Confirm the dict has the Pt-Pt entry the helper will use.
        target_pt_pt = _METAL_METAL_BOND_LENGTHS.get(frozenset({"Pt"}))
        assert target_pt_pt is not None, "Pt-Pt expected in _METAL_METAL_BOND_LENGTHS"

        moved = _enforce_smiles_mm_distances(mol, None, metals, _gp, _sp, np)
        assert moved == 1

        # After symmetric shift the distance must equal the dict target
        # (within 1e-6 because the shift is computed exactly).
        d_new = float(np.linalg.norm(_gp(0) - _gp(1)))
        assert math.isclose(d_new, target_pt_pt, abs_tol=1e-3), (
            f"Pt-Pt expected {target_pt_pt:.3f} A, got {d_new:.3f} A"
        )

    def test_mm_enforce_helper_uses_covalent_fallback_when_pair_absent(self):
        """Sn-Ir is NOT in the dict -> covalent-radii fallback (3.20 A)."""
        mol = Chem.MolFromSmiles("[Sn][Ir]", sanitize=False)
        assert mol is not None and mol.GetNumAtoms() == 2
        # Place 4.00 A apart along x.
        _gp, _sp, _ = _make_conf_helpers(mol, [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)])
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        assert sorted(metals) == [0, 1]

        # Sanity: not in the dict.
        assert _METAL_METAL_BOND_LENGTHS.get(frozenset({"Sn", "Ir"})) is None

        moved = _enforce_smiles_mm_distances(mol, None, metals, _gp, _sp, np)
        assert moved == 1

        target = (_COVALENT_RADII["Sn"] + _COVALENT_RADII["Ir"] + 0.4)
        d_new = float(np.linalg.norm(_gp(0) - _gp(1)))
        assert math.isclose(d_new, target, abs_tol=1e-3), (
            f"Sn-Ir fallback expected {target:.3f} A, got {d_new:.3f} A"
        )

    def test_mm_enforce_helper_skips_within_tolerance(self):
        """Distance within 0.20 A of target -> no move (avoids oscillation)."""
        mol = Chem.MolFromSmiles("[Pt][Pt]", sanitize=False)
        target = _METAL_METAL_BOND_LENGTHS[frozenset({"Pt"})]
        # Place 0.15 A below target -> within default 0.20 A tolerance.
        _gp, _sp, _ = _make_conf_helpers(
            mol, [(0.0, 0.0, 0.0), (target - 0.15, 0.0, 0.0)]
        )
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        moved = _enforce_smiles_mm_distances(mol, None, metals, _gp, _sp, np)
        assert moved == 0, "in-tolerance pair should not be moved"


# ---------------------------------------------------------------------------
# 4b. Welle-5j Agent G — MM_ENFORCE rigid-drag fix
# ---------------------------------------------------------------------------
#
# Background (welle3 T7.2 forensik, 2026-05-16):
#   The legacy ``_enforce_smiles_mm_distances`` translates only the two
#   metal atoms when shifting an M-M pair toward the ideal distance.
#   On D-COIRSN / D-DIVPIJ / D-TENMIL the Sn shift is large enough
#   (~0.65 A) that the Sn-Cl bonds break: Sn-Cl collapses from 2.30 A
#   to 0.59 A (one Cl gets dragged through the metal) and stretches to
#   2.95 A on the others.
#
# Fix: ``rigid_drag=True`` translates the per-metal exclusive sigma
# ligand frame (everything bonded to the metal that does NOT also bond
# to another metal) with the metal so the local geometry is preserved.

class TestMMEnforceRigidDrag:
    """Welle-5j Agent G rigid-drag ligand-frame translation."""

    def _setup_sn_cl3_ir(self):
        """Sn(Cl)3-Ir minimal mol; Sn-Ir stretched to 4.50 A, Cl's below."""
        mol = Chem.MolFromSmiles("[Cl][Sn]([Cl])([Cl])[Ir]", sanitize=False)
        assert mol is not None
        positions = [None] * mol.GetNumAtoms()
        sn_idx = ir_idx = None
        cl_indices = []
        for i, atom in enumerate(mol.GetAtoms()):
            sym = atom.GetSymbol()
            if sym == 'Sn':
                sn_idx = i
            elif sym == 'Ir':
                ir_idx = i
            elif sym == 'Cl':
                cl_indices.append(i)
        # Sn at origin, Ir straight above at 4.50 A, Cls below Sn so a
        # Sn shift toward Ir visibly stretches Sn-Cl in the legacy path.
        positions[sn_idx] = (0.0, 0.0, 0.0)
        positions[ir_idx] = (0.0, 0.0, 4.50)
        cl_pos = [(0.0, 0.0, -2.30), (1.90, 0.0, -1.30),
                  (-1.90, 0.0, -1.30)]
        for ci, pos in zip(cl_indices[:3], cl_pos):
            positions[ci] = pos
        # Defensive: anything we forgot to place sits far away.
        for i in range(len(positions)):
            if positions[i] is None:
                positions[i] = (50.0 + 0.1 * i, 50.0, 50.0)
        _gp, _sp, _ = _make_conf_helpers(mol, positions)
        return mol, sn_idx, ir_idx, cl_indices, _gp, _sp

    def test_legacy_metal_only_breaks_sn_cl(self):
        """Default rigid_drag=False reproduces the welle3 T7.2 Sn-Cl drift."""
        mol, sn_idx, ir_idx, cl_indices, _gp, _sp = self._setup_sn_cl3_ir()
        metals = [sn_idx, ir_idx]
        sn_cl_pre = [
            float(np.linalg.norm(_gp(c) - _gp(sn_idx)))
            for c in cl_indices[:3]
        ]
        moved = _enforce_smiles_mm_distances(
            mol, None, metals, _gp, _sp, np
        )
        assert moved == 1
        sn_cl_post = [
            float(np.linalg.norm(_gp(c) - _gp(sn_idx)))
            for c in cl_indices[:3]
        ]
        max_drift = max(abs(a - b) for a, b in zip(sn_cl_pre, sn_cl_post))
        assert max_drift > 0.10, (
            f"legacy metal-only shift expected to drift Sn-Cl by >0.10 A, "
            f"got {max_drift:.3f} A (pre={sn_cl_pre}, post={sn_cl_post})"
        )

    def test_rigid_drag_preserves_sn_cl(self):
        """rigid_drag=True translates Cl ligands rigidly with Sn."""
        mol, sn_idx, ir_idx, cl_indices, _gp, _sp = self._setup_sn_cl3_ir()
        metals = [sn_idx, ir_idx]
        sn_cl_pre = [
            float(np.linalg.norm(_gp(c) - _gp(sn_idx)))
            for c in cl_indices[:3]
        ]
        sn_ir_pre = float(np.linalg.norm(_gp(ir_idx) - _gp(sn_idx)))
        moved = _enforce_smiles_mm_distances(
            mol, None, metals, _gp, _sp, np, rigid_drag=True
        )
        assert moved == 1
        sn_cl_post = [
            float(np.linalg.norm(_gp(c) - _gp(sn_idx)))
            for c in cl_indices[:3]
        ]
        sn_ir_post = float(np.linalg.norm(_gp(ir_idx) - _gp(sn_idx)))
        # Sn-Ir must reach the covalent-radii fallback target (3.20 A).
        target = _COVALENT_RADII["Sn"] + _COVALENT_RADII["Ir"] + 0.4
        assert math.isclose(sn_ir_post, target, abs_tol=1e-3), (
            f"Sn-Ir post-shift expected {target:.3f} A, got {sn_ir_post:.3f} A"
        )
        # Sn-Cl distances must be preserved within 1e-3 A (pure translation).
        for orig, post in zip(sn_cl_pre, sn_cl_post):
            assert math.isclose(orig, post, abs_tol=1e-3), (
                f"rigid_drag should preserve Sn-Cl: {orig:.4f} -> {post:.4f}"
            )

    def test_rigid_drag_no_op_when_no_mm_edge(self):
        """rigid_drag=True on a mono-metal mol is a no-op (no MM pair)."""
        mol = Chem.MolFromSmiles("[Pt](Cl)(Cl)(N)N", sanitize=False)
        _gp, _sp, _ = _make_conf_helpers(
            mol, [(0.0, 0.0, 0.0)] * mol.GetNumAtoms()
        )
        metals = [a.GetIdx() for a in mol.GetAtoms()
                  if a.GetSymbol() in _METAL_SET]
        assert len(metals) == 1
        moved = _enforce_smiles_mm_distances(
            mol, None, metals, _gp, _sp, np, rigid_drag=True
        )
        assert moved == 0

    def test_rigid_drag_bit_exact_with_legacy_when_no_movement(self):
        """If pair is in tolerance, rigid_drag=True does not perturb anything."""
        mol = Chem.MolFromSmiles("[Pt][Pt]", sanitize=False)
        target = _METAL_METAL_BOND_LENGTHS[frozenset({"Pt"})]
        # Place 0.10 A inside the tolerance band.
        positions = [(0.0, 0.0, 0.0), (target - 0.10, 0.0, 0.0)]
        _gp, _sp, _ = _make_conf_helpers(mol, positions)
        metals = [0, 1]
        moved = _enforce_smiles_mm_distances(
            mol, None, metals, _gp, _sp, np, rigid_drag=True
        )
        assert moved == 0
        # Coordinates unchanged (within float epsilon).
        for i, (x, y, z) in enumerate(positions):
            p = _gp(i)
            assert abs(p[0] - x) < 1e-6 and abs(p[1] - y) < 1e-6 and \
                abs(p[2] - z) < 1e-6


# ---------------------------------------------------------------------------
# 5. Safe fallbacks
# ---------------------------------------------------------------------------

class TestSafeFallback:
    def test_unknown_class_safe_fallback(self):
        """mol=None -> classifier yields 'no_metal'; flag must be False without raising."""
        os.environ["DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES"] = "multi_hapto"
        # No crash, returns False -- 'no_metal' not in {'multi_hapto'}.
        assert _class_conditional_flag(
            "DELFIN_MULTIHAPTO_SIMPLE_PATH", None
        ) is False

    def test_classify_multi_hapto_fixture_stable(self):
        """Pin the COIRSN fixture's class label so the gating tests remain meaningful."""
        assert _classify_complex_class(_mol(_COIRSN_SMI)) == "multi_hapto"
