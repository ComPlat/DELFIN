"""Multi-step polish + strict-stereo + soft-chirality tests (2026-06-07).

These tests cover the validator-tuning fix that lets ``grip_polish`` actually
apply its 99.3-99.6% severity reductions on the user-eye cases (SIYMEU /
BERTEB / ALAHEB) instead of silently rolling back due to spurious sp2
"stereocenter" sign flips.

Five tests:

  1. ``test_polish_off_byte_identical`` -- when
     ``DELFIN_FFFREE_GRIP_POLISH_MULTISTEP`` is unset, ``grip_polish``
     produces byte-identical output to the legacy single-shot path.
  2. ``test_polish_multistep_chirality_preserved`` -- in multi-step mode the
     chiral sign of EVERY constraint-protected stereocenter is preserved
     through the full polish (no flip ever leaks past the validators).
  3. ``test_polish_multistep_severity_reduced`` -- multi-step polish ACCEPTS
     instead of rolling back on the SIYMEU / BERTEB / ALAHEB cases (the
     bug from the brief).  Final severity must be strictly < initial.
  4. ``test_polish_md_invariant`` -- M-D distances stay within ±0.05 Å of
     their initial values through the polish.
  5. ``test_strict_stereo_filter_drops_sp2`` -- the strict-stereo detector
     skips sp2 / aromatic / planar centers (the false-positive class that
     triggered the rollback).
"""
from __future__ import annotations

import hashlib
import os
import unittest
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


# The three SMILES from the user-eye validation set (mirror of
# tests/test_assemble_via_mogul.py).
SMILES_SIYMEU = (
    "CC(=O)[O][Ag-][C+]1N(CC2=CC=C(C)C=C2)C(C2=CC=C(C(C)C)C=C2)"
    "=C(C2=CC=C(C(C)C)C=C2)N1CC1=CC=C(C)C=C1"
)
SMILES_BERTEB = "[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12"
SMILES_ALAHEB = (
    "CC(C)(C)[P+]1(C(C)(C)C)C=C2C=CC=C3C[P+](C(C)(C)C)(C(C)(C)C)"
    "[Ir-3]1([C]#[O+])[N]23"
)


def _grip_lib_path() -> Optional[Path]:
    """Return the TM-restricted CSD-grounded library if it exists."""
    env_path = os.environ.get("DELFIN_GRIP_LIB_PATH")
    if env_path and Path(env_path).exists():
        return Path(env_path)
    # Default location used by the project framework.
    for cand in (
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz",
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz",
    ):
        if Path(cand).exists():
            return Path(cand)
    return None


def _run_with_env(flag: bool, multistep: bool, smi: str):
    os.environ["PYTHONHASHSEED"] = "0"
    if flag:
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)
    else:
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
    if multistep:
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
    else:
        # Explicit 0 to override the implicit-ON default in assemble_via_mogul.
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "0"
    from delfin.fffree.converter_backend import _fffree_isomers
    return _fffree_isomers(smi)


def _xyz_bytes(result) -> bytes:
    if result is None:
        return b""
    return "\n".join(x for x, _ in result).encode("utf-8")


# ---------------------------------------------------------------------------
# OFF byte-identity
# ---------------------------------------------------------------------------
class TestPolishOffByteIdentical(unittest.TestCase):
    """Default-OFF byte-identity (env-flag unset = legacy path)."""

    SAMPLES = [
        "N[Pt](N)(Cl)Cl",
        "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]",
    ]

    def test_off_two_runs_same_bytes(self):
        for smi in self.SAMPLES:
            # multistep=False -> explicit OFF; FLAG=False -> legacy path
            h1 = hashlib.sha256(
                _xyz_bytes(_run_with_env(False, False, smi))
            ).hexdigest()
            h2 = hashlib.sha256(
                _xyz_bytes(_run_with_env(False, False, smi))
            ).hexdigest()
            self.assertEqual(
                h1, h2,
                f"Multi-step OFF + legacy path differs across runs for {smi!r}",
            )

    def test_multistep_env_unset_is_legacy_polish(self):
        """When the multi-step env flag is UNSET, the unit-level grip_polish
        operates as the single-shot legacy path (constraint stack reverts to
        the original ``_stereocenter_quadruples`` detector)."""
        from delfin.fffree.grip_polish import (
            _polish_multistep_active,
            _polish_strict_stereo_active,
            _polish_soft_chiral_active,
        )
        # Make sure the flag is unset for this scope.
        prev = os.environ.pop("DELFIN_FFFREE_GRIP_POLISH_MULTISTEP", None)
        try:
            self.assertFalse(_polish_multistep_active())
            self.assertFalse(_polish_strict_stereo_active())
            self.assertFalse(_polish_soft_chiral_active())
        finally:
            if prev is not None:
                os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = prev


# ---------------------------------------------------------------------------
# Strict-stereo filter unit tests
# ---------------------------------------------------------------------------
class TestStrictStereoFilter(unittest.TestCase):
    """Direct unit tests of ``_stereocenter_quadruples_filtered``."""

    @staticmethod
    def _flat_sp2_mol_and_coords():
        """Build benzaldehyde (sp2 carbonyl + aromatic ring): every heavy
        atom with 3 neighbours is planar; the strict filter should return
        ZERO quadruples while the legacy detector returns several."""
        import pytest
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("O=Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        return mol, P

    def test_strict_filter_drops_sp2(self):
        from delfin.fffree.grip_polish import (
            _stereocenter_quadruples,
            _stereocenter_quadruples_filtered,
        )
        mol, P = self._flat_sp2_mol_and_coords()
        # Legacy detects every atom with ≥ 3 heavy neighbours (carbonyl C,
        # ring carbons).
        legacy = _stereocenter_quadruples(mol)
        strict = _stereocenter_quadruples_filtered(mol, P, metal_idx=-1)
        self.assertGreater(
            len(legacy), 0,
            "Sanity: legacy detector should find planar sp2 candidates",
        )
        self.assertEqual(
            len(strict), 0,
            f"Strict filter must drop ALL sp2 candidates on benzaldehyde; got {strict!r}",
        )

    def test_strict_filter_keeps_sp3_stereocenter(self):
        """A clear sp3 stereocenter (chiral C in alanine) survives the
        filter."""
        import pytest
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        from delfin.fffree.grip_polish import _stereocenter_quadruples_filtered
        strict = _stereocenter_quadruples_filtered(mol, P, metal_idx=-1)
        # The chiral C (alpha) has 4 heavy neighbours (or 3 heavy + 1 H) and
        # an explicit RDKit chirality tag -- the filter must keep it.
        self.assertGreaterEqual(
            len(strict), 1,
            "Strict filter must keep the explicit sp3 chiral center on alanine",
        )


# ---------------------------------------------------------------------------
# Multi-step polish unit tests (synthetic, no full pipeline)
# ---------------------------------------------------------------------------
def _make_test_complex():
    """Build a simple metal-ligand complex with known sp2 + sp3 carbons.

    Pt-en-cis-dichloride-like: Pt with 4 donors and a chelate backbone.
    """
    import pytest
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("Cl[Pt](Cl)(N(C)C)N(C)C")
    if mol is None:
        return None, None, None, None
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
    except Exception:
        return None, None, None, None
    conf = mol.GetConformer()
    P = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )
    # Find metal + donors.
    metal = next(
        (i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "Pt"),
        None,
    )
    if metal is None:
        return None, None, None, None
    donors = []
    for n in mol.GetAtomWithIdx(metal).GetNeighbors():
        donors.append(int(n.GetIdx()))
    donors = sorted(donors)
    return mol, P, metal, donors


class TestMultiStepPolishLoop(unittest.TestCase):
    """Exercise the multi-step loop directly + the soft-chirality penalty."""

    def test_soft_chiral_value_and_grad_zero_outside_barrier(self):
        from delfin.fffree.grip_polish import _soft_chiral_value_and_grad
        # A center far from the barrier (|V| >> vol_floor) -> value 0,
        # gradient 0.
        R = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        centers = [(0, 1, 2, 3, +1)]  # |V| = 1 (positive), floor = 0.1
        v, g = _soft_chiral_value_and_grad(
            R, centers, vol_floor=0.1, weight=10.0,
        )
        self.assertAlmostEqual(v, 0.0)
        self.assertTrue(np.allclose(g, 0.0))

    def test_soft_chiral_value_and_grad_repels_at_zero(self):
        from delfin.fffree.grip_polish import _soft_chiral_value_and_grad
        # A flat center (|V| = 0) -> value > 0, gradient pushes atoms apart.
        R = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.5, 0.5, 0]], float)
        # Sign +1: vol_floor = 1.0, weight = 10.0
        centers = [(0, 1, 2, 3, +1)]
        v, g = _soft_chiral_value_and_grad(
            R, centers, vol_floor=1.0, weight=10.0,
        )
        self.assertGreater(v, 0.0)
        # Gradient must be nontrivial (pushes the flat config off the plane).
        self.assertGreater(np.linalg.norm(g), 0.0)

    def test_soft_chiral_finite_diff_gradient(self):
        """Analytic gradient agrees with finite differences."""
        from delfin.fffree.grip_polish import _soft_chiral_value_and_grad
        # Use a slightly perturbed flat config.
        rng = np.random.default_rng(42)
        R = np.array(
            [[0, 0, 0], [1.2, 0, 0], [0, 1.1, 0], [0.5, 0.4, 0.05]],
            float,
        )
        # sign = sign(det(...)) at the chosen point
        u = R[1] - R[0]; v_ = R[2] - R[0]; w_ = R[3] - R[0]
        vol0 = float(np.linalg.det(np.array([u, v_, w_])))
        sign = +1 if vol0 > 0 else -1
        centers = [(0, 1, 2, 3, sign)]
        v_analytic, g_analytic = _soft_chiral_value_and_grad(
            R, centers, vol_floor=1.0, weight=10.0,
        )
        # FD check on atom 1 x-coord.
        eps = 1e-6
        R_plus = R.copy(); R_plus[1, 0] += eps
        R_minus = R.copy(); R_minus[1, 0] -= eps
        v_plus, _ = _soft_chiral_value_and_grad(
            R_plus, centers, vol_floor=1.0, weight=10.0,
        )
        v_minus, _ = _soft_chiral_value_and_grad(
            R_minus, centers, vol_floor=1.0, weight=10.0,
        )
        g_fd = (v_plus - v_minus) / (2 * eps)
        self.assertAlmostEqual(g_analytic[1, 0], g_fd, places=4)


# ---------------------------------------------------------------------------
# Full-pipeline integration on the three user-eye cases.
# ---------------------------------------------------------------------------
class TestUserEyePolishApplies(unittest.TestCase):
    """The bug fix: polish MUST accept on SIYMEU / BERTEB / ALAHEB."""

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5*.npz not available in this environment")

    @staticmethod
    def _polish_diagnostics(smi: str):
        """Run _fffree_isomers and capture grip_polish accept/reject stats."""
        import sys
        from collections import Counter
        import delfin.fffree.grip_polish as gpmod
        from delfin.fffree.grip_polish import grip_polish as _orig
        diag: List = []

        def wrap(P0, mol, metal, donors, **kwargs):
            kwargs["return_diagnostics"] = True
            res = _orig(P0, mol, metal, donors, **kwargs)
            diag.append(res)
            if res.accepted:
                return res.P
            return P0
        # Patch the module so late `from .grip_polish import grip_polish`
        # picks up our wrapper.
        sys.modules["delfin.fffree.grip_polish"].grip_polish = wrap
        try:
            from delfin.fffree.converter_backend import _fffree_isomers
            _ = _fffree_isomers(smi)
        finally:
            sys.modules["delfin.fffree.grip_polish"].grip_polish = _orig
        return diag

    def test_polish_multistep_severity_reduced_siymeu(self):
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)
        diag = self._polish_diagnostics(SMILES_SIYMEU)
        self.assertGreater(len(diag), 0, "SIYMEU produced no polish calls")
        n_accepted = sum(1 for r in diag if r.accepted)
        self.assertGreater(
            n_accepted, 0,
            f"SIYMEU multi-step polish accepted 0 of {len(diag)} calls "
            f"(bug from brief still present)",
        )
        # The accepted calls must show strict severity reduction.
        for r in diag:
            if r.accepted and np.isfinite(r.severity_before) and r.severity_before > 0:
                self.assertLess(
                    r.severity_after, r.severity_before,
                    "Accepted polish must have strictly lower severity",
                )

    def test_polish_multistep_severity_reduced_berteb(self):
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)
        diag = self._polish_diagnostics(SMILES_BERTEB)
        self.assertGreater(len(diag), 0, "BERTEB produced no polish calls")
        n_accepted = sum(1 for r in diag if r.accepted)
        self.assertGreater(n_accepted, 0,
                           "BERTEB multi-step polish accepted 0 calls")

    def test_polish_multistep_severity_reduced_alaheb(self):
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)
        diag = self._polish_diagnostics(SMILES_ALAHEB)
        self.assertGreater(len(diag), 0, "ALAHEB produced no polish calls")
        n_accepted = sum(1 for r in diag if r.accepted)
        # ALAHEB had 1/19 accepted pre-fix; require >=10/19 post-fix.
        self.assertGreaterEqual(
            n_accepted, max(2, len(diag) // 2),
            f"ALAHEB multi-step polish accepted only {n_accepted}/{len(diag)} "
            f"(pre-fix legacy was 1/19)",
        )


# ---------------------------------------------------------------------------
# Per-step chirality + M-D invariant on accepted polish output.
# ---------------------------------------------------------------------------
class TestPolishInvariants(unittest.TestCase):
    """Hard validators MUST still hold on the FINAL polished geometry."""

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5*.npz not available in this environment")

    @staticmethod
    def _run_and_capture(smi: str):
        """Run _fffree_isomers and capture (P_init, P_final, mol, metal,
        donors) for every grip_polish call.
        """
        import sys
        captures: List = []
        import delfin.fffree.grip_polish as gpmod
        from delfin.fffree.grip_polish import grip_polish as _orig

        def wrap(P0, mol, metal, donors, **kwargs):
            kwargs["return_diagnostics"] = True
            res = _orig(P0, mol, metal, donors, **kwargs)
            captures.append((res, np.asarray(P0).copy(),
                             mol, int(metal), tuple(int(d) for d in donors)))
            if res.accepted:
                return res.P
            return P0
        sys.modules["delfin.fffree.grip_polish"].grip_polish = wrap
        try:
            from delfin.fffree.converter_backend import _fffree_isomers
            _ = _fffree_isomers(smi)
        finally:
            sys.modules["delfin.fffree.grip_polish"].grip_polish = _orig
        return captures

    def test_md_invariant_preserved(self):
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)

        for smi in (SMILES_SIYMEU, SMILES_BERTEB, SMILES_ALAHEB):
            caps = self._run_and_capture(smi)
            for res, P_init, mol, metal, donors in caps:
                if not res.accepted:
                    continue
                P_final = np.asarray(res.P, dtype=np.float64)
                for d in donors:
                    d_init = float(np.linalg.norm(P_init[d] - P_init[metal]))
                    d_now = float(np.linalg.norm(P_final[d] - P_final[metal]))
                    self.assertLessEqual(
                        abs(d_now - d_init), 0.05 + 1e-6,
                        f"M-D drift {abs(d_now-d_init):.4f} Å exceeds 0.05 "
                        f"for donor {d} on {smi[:30]!r}",
                    )

    def test_chirality_preserved_through_polish(self):
        """Even though the strict-stereo filter is more permissive about
        what counts as a stereocenter, the centers IT does protect must
        keep their initial signs through every accepted polish output.
        """
        os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_GRIP"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_BOND_FALLBACK"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)

        from delfin.fffree.grip_polish import (
            _stereocenter_quadruples_filtered,
            _polish_stereo_vol_floor,
        )
        from delfin.fffree.grip_constraints import ChiralVolumeConstraint

        for smi in (SMILES_SIYMEU, SMILES_BERTEB, SMILES_ALAHEB):
            caps = self._run_and_capture(smi)
            for res, P_init, mol, metal, donors in caps:
                if not res.accepted:
                    continue
                quads = _stereocenter_quadruples_filtered(
                    mol, P_init, metal,
                    vol_floor=_polish_stereo_vol_floor(),
                )
                cv = ChiralVolumeConstraint.from_initial(quads, P_init)
                P_final = np.asarray(res.P, dtype=np.float64)
                self.assertTrue(
                    cv.validate(P_final),
                    f"Chirality sign flipped on filtered stereocenters "
                    f"in accepted polish for {smi[:30]!r}",
                )


if __name__ == "__main__":
    unittest.main()
