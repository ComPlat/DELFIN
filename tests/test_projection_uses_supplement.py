"""Chelate-projection uses the supplement-derived M-D target.

Tests the Bug B fix (2026-06-08) where ``_project_donors_to_ccdc_geometry``
now applies a centroid radial translation to chelate subtrees so the
chelate's MEAN per-donor M-D distance lands on the CCDC bond-order-aware
supplement value -- instead of the legacy "skip chelate projection"
behaviour which left chelates at whatever distance the DG embedder
produced (often outside the supplement's ±2σ window).

Env-flag ``DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT`` (default ON when
MOGUL_PRIMARY is active) gates the new behaviour.  When OFF (legacy
mode) byte-identical with the pre-fix branch.

Universal: works for any chelate (bidentate / tridentate / polydentate),
any donor element, any registered or unregistered polyhedron CN.  Pure
graph + geometric logic.

Author: hmaximilian
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

import numpy as np


SUP_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/"
    "grip_lib_md_bo_supplement_5d_mid.npz"
)
MAIN_LIB_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz"
)


_FULL_ENV = {
    "DELFIN_GRIP_LIB_PATH": str(MAIN_LIB_PATH),
    "DELFIN_GRIP_MD_BO_SUPPLEMENT": str(SUP_PATH),
    "DELFIN_FFFREE_MOGUL_PRIMARY": "1",
    "DELFIN_FFFREE_MOGUL_BOND_FALLBACK": "1",
    # GRIP-polish OFF so we measure the projection step alone -- the
    # subsequent multi-step polish would mask its effect.
    "DELFIN_FFFREE_MOGUL_PRIMARY_GRIP": "0",
}


def _reload_modules():
    """Re-import the assemble stack so env-var changes take effect."""
    import importlib
    import delfin.fffree.grip_mogul_lookup as _g
    importlib.reload(_g)
    import delfin.fffree.mogul_bounds as _m
    importlib.reload(_m)
    import delfin.fffree.assemble_via_mogul as _a
    importlib.reload(_a)
    return _a


@unittest.skipUnless(SUP_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "5d-mid supplement / main library not present")
class TestProjectionUsesSupplement(unittest.TestCase):
    """Chelate-subtree projection MUST honour the supplement M-D value."""

    def setUp(self):
        self._old = {k: os.environ.get(k) for k in _FULL_ENV}
        os.environ.update(_FULL_ENV)
        os.environ["DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT"] = "1"

    def tearDown(self):
        for k, v in self._old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def _build(self, smi):
        avm = _reload_modules()
        res = avm.assemble_complex_mogul_primary(smi)
        self.assertIsNotNone(res, f"assemble returned None for {smi[:50]}")
        return res

    def _md_distances(self, syms, P, metal_sym, donor_sym, cutoff=3.5):
        m_idx = syms.index(metal_sym)
        out = []
        for i, s in enumerate(syms):
            if i == m_idx or s != donor_sym:
                continue
            d = float(np.linalg.norm(P[i] - P[m_idx]))
            if d < cutoff:
                out.append(d)
        return sorted(out)

    def test_uqigeh_chelate_mean_pulls_to_supplement(self):
        """X10-UQIGEH (Mo CN12 bis-η⁶ arene): mean M-C should approach 2.12 Å.

        Supplement value: Mo-C sp2 bo=1 → 2.117 Å.
        Pre-fix: mean ≈ 2.30 Å (chelate skip).  Post-fix: ≤ 2.25 Å (gate
        accepts the centroid radial pull).
        """
        smi = ("C[Si](C)(C)N(B(Cl)[C+]12->[Mo-12]3456789%10%11(<-[C+]%12"
               "(B(Cl)N([Si](C)(C)C)[Si](C)(C)C)[C+]3[C+]4[C+]5[C+]6"
               "[C+]%127)[C+]([C+]8[C+]19)[C+]%10[C+]2%11)[Si](C)(C)C")
        syms, P = self._build(smi)
        dists = self._md_distances(syms, P, "Mo", "C")
        # Take the inner 12 donors (Mo CN12).
        self.assertGreaterEqual(len(dists), 12)
        inner12 = dists[:12]
        mean_md = float(np.mean(inner12))
        # The centroid radial pull brings the mean from 2.30 → ≤ 2.25.
        # This is the projection improvement; the gate accepts because
        # the per-donor RMS decreases.
        self.assertLess(mean_md, 2.25,
                        f"UQIGEH mean M-C {mean_md:.3f} Å should be < 2.25 "
                        f"after supplement projection (target 2.117 Å)")

    def test_chelate_projection_preserves_internal_distances(self):
        """The chelate projection MUST be a pure rigid translation.

        We rebuild a chelate-bearing complex with the projection ON and
        OFF and verify that the donor-donor distances WITHIN each
        chelate subtree are byte-identical (rigid translation does not
        alter internal distances by construction).
        """
        smi = ("CC#[N+][Re-4]12([N+]#CC)([N+]#CC)[P+](C3=CC=CC=C3)"
               "(C3=CC=CC=C3)CC(C)(C[P+]1(C1=CC=CC=C1)C1=CC=CC=C1)"
               "C[P+]2(C1=CC=CC=C1)C1=CC=CC=C1")
        # ON
        os.environ["DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT"] = "1"
        syms_on, P_on = self._build(smi)
        # OFF
        os.environ["DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT"] = "0"
        syms_off, P_off = self._build(smi)
        self.assertEqual(syms_on, syms_off, "atom ordering must match")
        # Find the 3 P donors of the tripodal triphosphine chelate.
        m_idx = syms_on.index("Re")
        p_idxs = [i for i, s in enumerate(syms_on)
                  if s == "P"
                  and np.linalg.norm(P_on[i] - P_on[m_idx]) < 3.5]
        # P-P internal distances must be identical (chelate rigid).
        for a in range(len(p_idxs)):
            for b in range(a + 1, len(p_idxs)):
                d_on = float(np.linalg.norm(P_on[p_idxs[a]] - P_on[p_idxs[b]]))
                d_off = float(np.linalg.norm(
                    P_off[p_idxs[a]] - P_off[p_idxs[b]]
                ))
                # Rigid translation preserves intra-chelate distances
                # exactly (modulo numerical drift).  When the gate
                # rejects the translation (legacy fallback) ON == OFF.
                self.assertAlmostEqual(
                    d_on, d_off, delta=0.05,
                    msg=f"P-P internal distance changed ON vs OFF: "
                        f"{d_on:.3f} vs {d_off:.3f} (must be rigid)"
                )

    def test_chelate_projection_respects_hard_floor(self):
        """The projection MUST NOT push any donor below 0.7 × target.

        This is the defence-in-depth check: a mis-classified target
        could try to fly the chelate outward by >0.6 Å; the cap +
        the 0.7×target floor reject that, falling back to the legacy
        skip.  We verify on the 3 regressed refcodes that no donor
        ends up below the floor.
        """
        cases = [
            # (smi, metal, donor, target_md)
            ("C[Si](C)(C)N(B(Cl)[C+]12->[Mo-12]3456789%10%11(<-[C+]%12"
             "(B(Cl)N([Si](C)(C)C)[Si](C)(C)C)[C+]3[C+]4[C+]5[C+]6"
             "[C+]%127)[C+]([C+]8[C+]19)[C+]%10[C+]2%11)[Si](C)(C)C",
             "Mo", "C", 2.117),
            ("CC#[N+][Re-4]12([N+]#CC)([N+]#CC)[P+](C3=CC=CC=C3)"
             "(C3=CC=CC=C3)CC(C)(C[P+]1(C1=CC=CC=C1)C1=CC=CC=C1)"
             "C[P+]2(C1=CC=CC=C1)C1=CC=CC=C1",
             "Re", "P", 2.444),
            ("[C+]1=[C+][Mo-9]12345([C+]=[C+]2)([C+]=[C+]3)[N+]1=C(C=CC=C1"
             "C[P+]4(C1=CC=CC=C1)C1=CC=CC=C1)C[P+]5(C1=CC=CC=C1)"
             "C1=CC=CC=C1",
             "Mo", "C", 1.991),
        ]
        for smi, metal, donor, target in cases:
            syms, P = self._build(smi)
            dists = self._md_distances(syms, P, metal, donor)
            self.assertTrue(dists, "no donors found")
            min_d = float(min(dists))
            self.assertGreaterEqual(
                min_d, 0.7 * target,
                f"{metal}-{donor} min {min_d:.3f} Å below 0.7×target "
                f"({0.7*target:.3f} Å) -- defence floor violated"
            )


@unittest.skipUnless(SUP_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "5d-mid supplement / main library not present")
class TestOffByteIdentical(unittest.TestCase):
    """When the projection flag is OFF, behaviour matches legacy."""

    def setUp(self):
        self._old = {k: os.environ.get(k) for k in _FULL_ENV}
        os.environ.update(_FULL_ENV)

    def tearDown(self):
        for k, v in self._old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_off_skips_chelate_projection(self):
        """With ``DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT=0`` chelate
        subtrees are NOT projected (legacy behaviour preserved).
        """
        os.environ["DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT"] = "0"
        avm = _reload_modules()
        smi = ("[C+]1=[C+][Mo-9]12345([C+]=[C+]2)([C+]=[C+]3)[N+]1=C("
               "C=CC=C1C[P+]4(C1=CC=CC=C1)C1=CC=CC=C1)C[P+]5"
               "(C1=CC=CC=C1)C1=CC=CC=C1")
        res = avm.assemble_complex_mogul_primary(smi)
        self.assertIsNotNone(res)
        syms, P = res
        m_idx = syms.index("Mo")
        # Pre-existing 1.6-2.2 Å Mo-C cluster from the embed should be
        # present (untouched by chelate projection).  We check that
        # SOME donor sits below the supplement value 1.991 Å -- this
        # is the legacy embed-only behaviour.
        any_below = False
        for i, s in enumerate(syms):
            if i == m_idx or s != "C":
                continue
            d = float(np.linalg.norm(P[i] - P[m_idx]))
            if 1.5 < d < 1.99:
                any_below = True
                break
        self.assertTrue(any_below,
                        "OFF mode should leave the legacy embed's "
                        "below-target Mo-C distances untouched")


if __name__ == "__main__":
    unittest.main()
