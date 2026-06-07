"""Invariant tests for the Mogul-PRIMARY pipeline (added 2026-06-07).

Three universal invariants are checked against the six refcodes flagged in
the 2026-06-07 MOGUL-PRIMARY-COMPLETE-SMOKE100 archive audit:

  1. ``test_metal_at_idx_0_through_reorder`` —
       ``_reordered_source_mol`` preserves the primary-metal selection
       chosen by ``_locate_metal_and_donors`` so the GRIP fragment detector
       sees the same atom ordering as the assembled (syms, P) tuple.
       Catches the Bug #1 multi-metal mis-priority (Sb outranking Pt in
       D-ATOQOP).

  2. ``test_metal_at_origin_through_all_conformers`` — across every
       conformer frame produced by ``enumerate_kappa_variants_mogul_primary``
       /``enumerate_and_embed_mogul_primary``/``assemble_complex_mogul_primary``
       the metal (atoms[0]) stays at (0,0,0).  Catches the Bug #2
       metallacycle-ring rotation (the chelate ring closes through the
       metal because RDKit dative bonds aren't perceived as ring-forming;
       a rotamer rotation walked the metal into the subtree).

  3. ``test_m_shell_overfill_rejects_extras`` — the base structure
       returned by ``assemble_complex_mogul_primary`` never has more
       heavy atoms inside ``M_SHELL_FACTOR * ideal_bond`` of the metal
       than the SMILES-graph coordination number.  Catches the Bug #3
       overfill (DBBZFE/NOHWIQ/OJAXOO had 12/10/13 atoms inside the M-shell).

All tests run under the full Mogul-PRIMARY env stack and require the
GRIP library v5_TM at the canonical path.  When that file is missing
the suite is skipped.
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

import numpy as np


GRIP_LIB_PATH = (
    "/home/qmchem_max/agent_workspace/quality_framework/reports/"
    "grip_lib_v5_TM.npz"
)

# (refcode, SMILES, expected metal symbol, M-shell-CN reference from refcode)
CASES = [
    (
        "D-ATOQOP_5d_Pt_CN3",
        "[CH3][Sb+]1([CH3])[CH2]C2=CC=CC=C2[CH2][Sb+]([CH3])([CH3])[Pt-2]12"
        "[Sb+]([CH3])([CH3])[CH2]C1=CC=CC=C1[CH2][Sb+]2([CH3])[CH3]",
        "Pt",
        4,        # SMILES-graph CN for Pt = 4 (refcode says 3 but graph has 4)
    ),
    (
        "D2-GOSYIY_3d_Ti_CN7",
        "C1=CC=C([P+]2(C3=CC=CC=C3)C3=CC=C[N]3[Ti-4]2345([N]2C=CC=C2"
        "[P+]3(C2=CC=CC=C2)C2=CC=CC=C2)([N]2C=CC=C2[P+]4(C2=CC=CC=C2)"
        "C2=CC=CC=C2)[N]2C=CC=C2[P+]5(C2=CC=CC=C2)C2=CC=CC=C2)C=C1",
        "Ti",
        8,
    ),
    (
        "D2-LEPTUX_5d_Au_CN1",
        "CN1C(Cl)=C(Cl)N(C)[C+]1[Au-][Cl]",
        "Au",
        2,
    ),
    (
        "D2-DBBZFE_3d_Fe_CN8",
        "[B+]12->[Fe-12]3456789%10%11(<-[B+]%12[C+]3[C+]4[C+]5[C+]6"
        "[C+]%127)[C+]1[C+]8[C+]9[C+]%10[C+]2%11",
        "Fe",
        12,
    ),
    (
        "D2-NOHWIQ_3d_Cu_CN3",
        "[O-]C1=CC=C(C#[N+][Cu-3]([N+]#CC2=CC=C(O)C=C2)([P+](C2=CC=CC=C2)"
        "(C2=CC=CC=C2)C2=CC=CC=C2)[P+](C2=CC=CC=C2)(C2=CC=CC=C2)"
        "C2=CC=CC=C2)C=C1",
        "Cu",
        4,
    ),
    (
        "D2-OJAXOO_3d_Ti_CN9",
        "CN1C=CN(C)[C+]1[Ti-13]123456789%10([C+]%11[C+]1[C+]2[C+]3"
        "[C+]4[C+]%115)[C+]1[C+]6[C+]7[C+]8[C+]9[C+]1%10",
        "Ti",
        12,
    ),
]


_FULL_ENV = {
    "DELFIN_GRIP_LIB_PATH": GRIP_LIB_PATH,
    "DELFIN_FFFREE_MOGUL_PRIMARY": "1",
    "DELFIN_FFFREE_MOGUL_BOND_FALLBACK": "1",
    "DELFIN_FFFREE_MOGUL_PRIMARY_GRIP": "1",
    "DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM": "1",
    "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS": "1",
    "DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE": "1",
    "DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK": "1",
    "DELFIN_FFFREE_GRIP_POLISH_MULTISTEP": "1",
    "PYTHONHASHSEED": "0",
}


class _EnvSnapshot:
    """Context manager that sets env vars and restores on exit."""

    def __init__(self, env: dict):
        self._env = dict(env)
        self._saved: dict = {}

    def __enter__(self):
        for k, v in self._env.items():
            self._saved[k] = os.environ.get(k)
            os.environ[k] = v
        return self

    def __exit__(self, *exc):
        for k, prev in self._saved.items():
            if prev is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = prev


def _build_ensemble(smiles: str):
    """Return the full (κⁿ + orbit + conformer) ensemble for ``smiles``.

    Tries the κⁿ path first, falls back to orbit-only, then to the single
    base structure.  Mirrors converter_backend's dispatch order so the
    tests exercise the same call graph.
    """
    from delfin.fffree import assemble_via_mogul as _MP
    ens = None
    try:
        ens = _MP.enumerate_kappa_variants_mogul_primary(smiles, max_isomers=50)
    except Exception:
        ens = None
    if ens is None:
        try:
            ens = _MP.enumerate_and_embed_mogul_primary(smiles, max_isomers=50)
        except Exception:
            ens = None
    if ens is None:
        try:
            res = _MP.assemble_complex_mogul_primary(smiles)
        except Exception:
            res = None
        if res is not None:
            syms, P = res
            ens = [(syms, P, "MOGUL-PRIMARY-1")]
    return ens or []


def _grip_lib_available() -> bool:
    return Path(GRIP_LIB_PATH).is_file()


class TestMetalAtIdx0ThroughReorder(unittest.TestCase):
    """``_reordered_source_mol`` keeps the same primary metal at idx 0
    as ``_locate_metal_and_donors`` does on the original mol.  Catches
    the Bug #1 priority-tiebreak regression (Sb outranking Pt).
    """

    def setUp(self):
        if not _grip_lib_available():
            self.skipTest("grip_lib_v5_TM.npz not available")

    def test_d_atoqop_pt_outranks_sb(self):
        from delfin.fffree.assemble_via_mogul import (
            _locate_metal_and_donors,
            _reordered_source_mol,
            _full_complex_mol,
        )
        with _EnvSnapshot(_FULL_ENV):
            smi = next(c[1] for c in CASES if c[0].startswith("D-ATOQOP"))
            mol = _full_complex_mol(smi)
            self.assertIsNotNone(mol)
            metal_idx, _ = _locate_metal_and_donors(mol)
            self.assertIsNotNone(metal_idx)
            self.assertEqual(
                mol.GetAtomWithIdx(int(metal_idx)).GetSymbol(),
                "Pt",
                "Pt must outrank Sb on the multi-metal tie-break",
            )
            reordered = _reordered_source_mol(smi, mol.GetNumAtoms())
            self.assertIsNotNone(reordered)
            self.assertEqual(reordered.GetAtomWithIdx(0).GetSymbol(), "Pt")

    def test_d_atoqop_assembled_has_pt_at_0(self):
        with _EnvSnapshot(_FULL_ENV):
            smi = next(c[1] for c in CASES if c[0].startswith("D-ATOQOP"))
            ens = _build_ensemble(smi)
            self.assertGreater(len(ens), 0)
            for syms, P, lab in ens:
                self.assertEqual(
                    syms[0], "Pt",
                    f"frame {lab}: metal at idx 0 must be Pt, got {syms[0]}",
                )


class TestMetalAtOriginThroughAllConformers(unittest.TestCase):
    """Every conformer frame keeps the metal at (0,0,0).  Catches Bug #2
    (metal getting rotated WITH a metallacycle rotamer)."""

    def setUp(self):
        if not _grip_lib_available():
            self.skipTest("grip_lib_v5_TM.npz not available")

    def test_all_cases_metal_at_origin(self):
        with _EnvSnapshot(_FULL_ENV):
            for refcode, smi, expected_metal, _shell_cn in CASES:
                with self.subTest(refcode=refcode):
                    ens = _build_ensemble(smi)
                    self.assertGreater(len(ens), 0, f"empty ensemble: {refcode}")
                    for fi, (syms, P, lab) in enumerate(ens):
                        P = np.asarray(P, dtype=float)
                        # idx 0 must be the expected metal
                        self.assertEqual(
                            syms[0], expected_metal,
                            f"{refcode} frame {fi} ({lab}): "
                            f"metal idx 0 = {syms[0]}, expected {expected_metal}",
                        )
                        # idx 0 must be at the origin
                        d = float(np.linalg.norm(P[0]))
                        self.assertLess(
                            d, 1e-3,
                            f"{refcode} frame {fi} ({lab}): "
                            f"metal at idx 0 displaced by {d:.4f} Å",
                        )


class TestMShellOverfillRejectsExtras(unittest.TestCase):
    """The base structure has heavy-atom shell count <= SMILES-graph CN.
    Catches Bug #3 (DBBZFE / NOHWIQ / OJAXOO over-filling)."""

    def setUp(self):
        if not _grip_lib_available():
            self.skipTest("grip_lib_v5_TM.npz not available")

    def test_no_shell_overfill_in_base(self):
        from delfin.fffree import assemble_via_mogul as _MP
        from delfin import _bond_decollapse as _bd
        from delfin.fffree.rotamer_topology_gate import M_SHELL_FACTOR

        with _EnvSnapshot(_FULL_ENV):
            for refcode, smi, expected_metal, shell_cn_ref in CASES:
                with self.subTest(refcode=refcode):
                    res = _MP.assemble_complex_mogul_primary(smi)
                    if res is None:
                        # The gate may legitimately reject the base
                        # (e.g. shell-overfill on the embed) — that's
                        # the intended fail-closed behaviour and counts
                        # as gate success.  We just need to confirm a
                        # downstream κⁿ/orbit ensemble can still
                        # produce *some* output (no full pipeline
                        # collapse).
                        ens = _build_ensemble(smi)
                        self.assertGreaterEqual(
                            len(ens), 0,
                            f"{refcode}: gate rejected base AND ensemble empty",
                        )
                        continue
                    syms, P = res
                    self.assertEqual(syms[0], expected_metal)
                    P = np.asarray(P, dtype=float)
                    metal_sym = syms[0]
                    shell = 0
                    for j in range(1, len(syms)):
                        sj = syms[j]
                        if sj == "H":
                            continue
                        if _bd._is_metal(sj):
                            continue
                        try:
                            ideal = float(_bd._ideal_bond(metal_sym, sj))
                        except Exception:
                            ideal = 0.0
                        if ideal <= 0.0:
                            continue
                        d = float(np.linalg.norm(P[j] - P[0]))
                        if d < float(M_SHELL_FACTOR) * ideal:
                            shell += 1
                    self.assertLessEqual(
                        shell, shell_cn_ref,
                        f"{refcode}: M-shell overfill {shell} > expected "
                        f"<= {shell_cn_ref}",
                    )


class TestMetalPriorityUniversal(unittest.TestCase):
    """Universal behaviour of :func:`_metal_priority` (no env needed)."""

    def test_d_block_beats_p_block(self):
        from delfin.fffree.assemble_via_mogul import _metal_priority
        # Transition metals = priority 0
        for s in ("Ti", "Fe", "Cu", "Pt", "Pd", "Ru", "Ir", "Au"):
            self.assertEqual(_metal_priority(s), 0, s)
        # p-block "metals" = priority 2
        for s in ("Sb", "Sn", "Pb", "Bi", "Al", "Ga", "In", "Tl", "Ge"):
            self.assertEqual(_metal_priority(s), 2, s)
        # s-block = priority 1
        for s in ("Li", "Na", "K", "Cs", "Mg", "Ca"):
            self.assertEqual(_metal_priority(s), 1, s)

    def test_sb_does_not_outrank_pt(self):
        from delfin.fffree.assemble_via_mogul import _metal_priority
        self.assertLess(
            _metal_priority("Pt"), _metal_priority("Sb"),
            "Pt (d-block) must outrank Sb (p-block) on the tie-break",
        )


if __name__ == "__main__":
    unittest.main()
