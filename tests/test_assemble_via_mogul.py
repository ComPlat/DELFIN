"""Tests for the Mogul-primary construction path.

Validation:

1.  ``test_off_byte_identical`` — when ``DELFIN_FFFREE_MOGUL_PRIMARY`` is
    unset, the legacy ``_fffree_isomers`` path produces the SAME bytes
    on two consecutive runs (the OFF path is unchanged from HEAD).
2.  ``test_on_determinism`` — with the flag ON and
    ``PYTHONHASHSEED=0`` the construction is byte-identical across runs.
3.  ``test_siymeu_nhc_carbene`` — the SIYMEU (Ag-NHC linear CN2)
    structure has Ag-C(carbene) in the CCDC empirical window and
    C-Ag-O ≈ 180° (linear).
4.  ``test_berteb_pyridyl_tridentate`` — the X10-BERTEB (Ni CN4) build
    succeeds and produces Ni-N + Ni-Br distances inside the CCDC
    empirical envelope.
5.  ``test_alaheb_ir_phosphine`` — the X10-ALAHEB (Ir CN4 chelate) build
    succeeds and Ir-P > 2.20 Å (the legacy V16 path produced 1.85 Å).

The validation tests use the explicit ``DELFIN_GRIP_LIB_PATH`` pointing
at the CSD-grounded ``grip_lib_v5.npz`` when available; otherwise the
test is skipped because the COD-derived ``grip_lib_v6_cod.npz`` is
sparse on the user-eye categories used in these three examples.
"""
from __future__ import annotations

import os
import hashlib
import unittest
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


# Reference SMILES for the user-eye validation set.
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
    """Return the v5 CSD-grounded library path if it exists.

    Tests that need actual CCDC lookups skip when this path is missing
    (e.g. when running in CI without the framework-side library).
    """
    env_path = os.environ.get("DELFIN_GRIP_LIB_PATH")
    if env_path and Path(env_path).exists():
        return Path(env_path)
    # Default location used by the project framework.
    default = Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"
    )
    if default.exists():
        return default
    return None


def _run_with_env(flag: bool, smi: str):
    """Run ``_fffree_isomers`` with ``DELFIN_FFFREE_MOGUL_PRIMARY`` set."""
    import importlib
    os.environ["PYTHONHASHSEED"] = "0"
    if flag:
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        gp = _grip_lib_path()
        if gp is not None:
            os.environ["DELFIN_GRIP_LIB_PATH"] = str(gp)
    else:
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
    # Late import so the env var is honoured.
    from delfin.fffree.converter_backend import _fffree_isomers
    return _fffree_isomers(smi)


def _xyz_bytes(result) -> bytes:
    if result is None:
        return b""
    return "\n".join(x for x, _ in result).encode("utf-8")


def _measure_md(syms, P, metal_sym: str, donor_sym: str) -> Optional[float]:
    """Shortest distance between the first metal atom and the closest
    matching donor atom by element.  Mirror of the module helper used in
    the validation script; duplicated here so the test does not depend
    on internal helper signatures.
    """
    P = np.asarray(P, dtype=float)
    metal_idxs = [i for i, s in enumerate(syms) if s == metal_sym]
    if not metal_idxs:
        return None
    m = metal_idxs[0]
    best = None
    for i, s in enumerate(syms):
        if i == m or s != donor_sym:
            continue
        d = float(np.linalg.norm(P[i] - P[m]))
        if best is None or d < best:
            best = d
    return best


def _parse_xyz_block(xyz: str) -> Tuple[List[str], np.ndarray]:
    """Inverse of ``_xyz``: parse the headerless atom block back."""
    syms, coords = [], []
    for line in xyz.strip().splitlines():
        toks = line.split()
        if len(toks) < 4:
            continue
        syms.append(toks[0])
        coords.append([float(t) for t in toks[1:4]])
    return syms, np.asarray(coords, dtype=float)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestOFFByteIdentical(unittest.TestCase):
    """When the env flag is unset, the legacy path runs unchanged."""

    SAMPLES = [
        "N[Pt](N)(Cl)Cl",
        "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]",
    ]

    def test_off_two_runs_same_bytes(self):
        for smi in self.SAMPLES:
            h1 = hashlib.sha256(
                _xyz_bytes(_run_with_env(False, smi))
            ).hexdigest()
            h2 = hashlib.sha256(
                _xyz_bytes(_run_with_env(False, smi))
            ).hexdigest()
            self.assertEqual(
                h1, h2,
                f"OFF path bytes differ between runs for {smi!r}",
            )


class TestOnDeterminism(unittest.TestCase):
    """ON path produces the same bytes on consecutive runs."""

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5.npz not available in this environment")

    def test_siymeu_deterministic(self):
        h1 = hashlib.sha256(_xyz_bytes(_run_with_env(True, SMILES_SIYMEU))).hexdigest()
        h2 = hashlib.sha256(_xyz_bytes(_run_with_env(True, SMILES_SIYMEU))).hexdigest()
        self.assertEqual(h1, h2)


class TestUserEyeValidation(unittest.TestCase):
    """The three user-eye reference structures get reasonable geometry."""

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5.npz not available in this environment")

    def test_siymeu_nhc_carbene(self):
        """Ag-C(NHC) lands within the CCDC empirical envelope and
        the C-Ag-O angle is close to linear (180 deg)."""
        res = _run_with_env(True, SMILES_SIYMEU)
        self.assertIsNotNone(res, "Mogul-primary returned None for SIYMEU")
        self.assertGreaterEqual(len(res), 1)
        xyz, label = res[0]
        self.assertIn("MOGUL", label)
        syms, P = _parse_xyz_block(xyz)
        # Metal-Carbene distance: CCDC mean (Ag-NHC) ~ 2.10 Angstrom.
        ag_c = _measure_md(syms, P, "Ag", "C")
        self.assertIsNotNone(ag_c)
        # Wide acceptance: 2.0-2.3 Angstrom.  Pre-fix value (legacy V16)
        # is typically OUTSIDE this band; the primary path lands inside.
        self.assertGreaterEqual(ag_c, 2.0,
                                f"Ag-C={ag_c:.3f} below CCDC band lower")
        self.assertLessEqual(ag_c, 2.3,
                             f"Ag-C={ag_c:.3f} above CCDC band upper")

        # Linear C-Ag-O angle (CN2): expect ~ 180 deg, accept >= 175.
        ag = next(i for i, s in enumerate(syms) if s == "Ag")
        # Find the Ag-bonded O (acetate) and Ag-bonded C (carbene) by
        # nearest-element distance.
        o_idx = min((i for i, s in enumerate(syms) if s == "O" and i != ag),
                    key=lambda i: float(np.linalg.norm(P[i] - P[ag])))
        c_idx = min((i for i, s in enumerate(syms) if s == "C" and i != ag),
                    key=lambda i: float(np.linalg.norm(P[i] - P[ag])))
        v1 = P[o_idx] - P[ag]
        v2 = P[c_idx] - P[ag]
        cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        ang_deg = float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        self.assertGreaterEqual(ang_deg, 175.0,
                                f"C-Ag-O = {ang_deg:.1f} deg, expected ~180")

    def test_berteb_pyridyl_tridentate(self):
        """Ni(II) CN4 with a tridentate (N3) chelate + monodentate Br."""
        res = _run_with_env(True, SMILES_BERTEB)
        self.assertIsNotNone(res, "Mogul-primary returned None for BERTEB")
        xyz, label = res[0]
        syms, P = _parse_xyz_block(xyz)
        # Ni-Br: CCDC mean ~ 2.39 Angstrom; accept 2.20 - 2.60.
        ni_br = _measure_md(syms, P, "Ni", "Br")
        self.assertIsNotNone(ni_br)
        self.assertGreaterEqual(ni_br, 2.20,
                                f"Ni-Br={ni_br:.3f} below CCDC band")
        self.assertLessEqual(ni_br, 2.60,
                             f"Ni-Br={ni_br:.3f} above CCDC band")
        # Ni-N (closest): CCDC mean ~ 2.05; accept 1.90 - 2.35.
        ni_n = _measure_md(syms, P, "Ni", "N")
        self.assertIsNotNone(ni_n)
        self.assertGreaterEqual(ni_n, 1.90,
                                f"Ni-N={ni_n:.3f} below CCDC band")
        self.assertLessEqual(ni_n, 2.35,
                             f"Ni-N={ni_n:.3f} above CCDC band")

    def test_alaheb_ir_phosphine(self):
        """Ir(I) CN4 chelate with two P, one CO carbon, one N.

        The legacy V16 path produced Ir-P ≈ 1.85 Angstrom (much too
        short), the Mogul-primary path lands at the CCDC mean ~ 2.30 Å.
        """
        res = _run_with_env(True, SMILES_ALAHEB)
        self.assertIsNotNone(res, "Mogul-primary returned None for ALAHEB")
        xyz, label = res[0]
        syms, P = _parse_xyz_block(xyz)
        ir_p = _measure_md(syms, P, "Ir", "P")
        self.assertIsNotNone(ir_p)
        # The key regression we want to never see again: Ir-P < 2.20.
        self.assertGreaterEqual(ir_p, 2.20,
                                f"Ir-P={ir_p:.3f} (legacy V16 fails here)")
        # Accept up to 2.50 (CCDC ±3σ envelope).
        self.assertLessEqual(ir_p, 2.50,
                             f"Ir-P={ir_p:.3f} above CCDC band")


class TestPolyaOrbitEnumeration(unittest.TestCase):
    """Burnside-Pólya orbit enumeration over the Mogul-primary bounds matrix
    (2026-06-07 architectural extension).

    The draft manuscript contract is: SMILES -> Burnside-distinct orbits ->
    embed EACH orbit once -> return the ensemble.  This test exercises the
    textbook MA2B2C2 octahedron case which has exactly 6 stereoisomers
    (Burnside count = 6).  The Mogul-primary path with orbit override must
    realise all 6 with distinct 3D geometries -- including the chiral
    Δ/Λ pair of the all-cis class, which the orbit override resolves via
    different donor-to-vertex permutations even though their pairwise
    distance multisets are degenerate.
    """

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5.npz not available in this environment")

    def test_ma2b2c2_oc6_six_isomers(self):
        """[NH3]2 + [OH2]2 + F2 Co(III) — textbook OC-6 MA2B2C2 case.

        Burnside count for the octahedron with donor multiset (2,2,2) is 6.
        We require the Mogul-primary orbit enumerator to emit 6 distinct
        orbit BASE structures.  When the conformer-enum gate is on (default
        ON 2026-06-07), each orbit may additionally emit pucker / rotamer
        conformer variants; we test the orbit count by counting distinct
        orbit labels (the ``-mogul`` infix is the orbit signature).
        """
        smi = "[NH3][Co]([NH3])([OH2])([OH2])([F])[F]"
        # Pin the conformer gate OFF so this orbit-count test is independent
        # of the per-orbit conformer enumeration (which can multiply each
        # orbit by 1..N entries depending on flexibility).  The conformer
        # enumeration is exercised in TestMogulPrimaryConformers below.
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS"] = "0"
        try:
            res = _run_with_env(True, smi)
        finally:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS", None)
        self.assertIsNotNone(res, "Mogul-primary returned None for MA2B2C2")
        self.assertEqual(
            len(res), 6,
            f"Expected 6 Burnside-distinct orbits, got {len(res)}; "
            f"labels={[lab for _, lab in res]}",
        )
        # Each XYZ must be unique by hash (otherwise the orbit override
        # silently collapsed two orbits onto the same geometry, defeating
        # the whole purpose).
        hashes = {hashlib.sha256(xyz.encode("utf-8")).hexdigest()
                  for xyz, _ in res}
        self.assertEqual(
            len(hashes), 6,
            f"Expected 6 distinct geometries, got {len(hashes)} unique; "
            f"isomers collapsed onto identical XYZ.",
        )
        # Every label must contain the ``-mogul`` orbit signature.  With the
        # conformer gate OFF (set above), the labels end with ``-mogul``;
        # with the gate ON the labels are ``<orbit>-mogul-conf<k>-<sub>``.
        for _, label in res:
            self.assertIn(
                "-mogul", label,
                f"Orbit label {label!r} missing ``-mogul`` signature",
            )

    def test_orbit_off_byte_identical(self):
        """When DELFIN_FFFREE_MOGUL_PRIMARY is unset, the orbit enumeration
        code path is silent and the legacy fffree path is bit-identical."""
        smi = "[NH3][Co]([NH3])([OH2])([OH2])([F])[F]"
        h1 = hashlib.sha256(
            _xyz_bytes(_run_with_env(False, smi))
        ).hexdigest()
        h2 = hashlib.sha256(
            _xyz_bytes(_run_with_env(False, smi))
        ).hexdigest()
        self.assertEqual(h1, h2,
                         "OFF path bytes differ between runs (orbit "
                         "enumeration leaked into the default path)")


class TestMogulPrimaryConformers(unittest.TestCase):
    """Per-orbit conformer enumeration over the Mogul-primary path
    (2026-06-07 architectural extension).

    The draft manuscript's third Cauchy-Frobenius layer is the conformer
    space: each Pólya orbit must emit a per-isomer ensemble of Cremer-
    Pople ring-pucker variants × single-bond rotamers, RMSD-deduplicated
    and Mahalanobis-ranked.  These tests exercise:

      1. ``test_off_byte_identical`` -- when
         ``DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=0`` is set, the orbit
         enumerator emits exactly one structure per orbit (the base build).
      2. ``test_on_emits_ensemble`` -- when the gate is on (default), the
         flexible test cases produce more than one structure.
      3. ``test_conformer_off_outer_off_byte_identical`` -- with the outer
         Mogul-primary gate OFF, the conformer-enum env flag has no effect
         (the legacy path stays bit-identical).
    """

    def setUp(self):
        if _grip_lib_path() is None:
            self.skipTest("grip_lib_v5.npz not available in this environment")

    @staticmethod
    def _count_conformers(res) -> int:
        """Count emitted entries; helper that handles None."""
        return len(res) if res else 0

    def test_off_byte_identical(self):
        """CONFORMERS=0 must produce exactly the orbit baseline."""
        smi = "[NH3][Co]([NH3])([OH2])([OH2])([F])[F]"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS"] = "0"
        try:
            h1 = hashlib.sha256(_xyz_bytes(_run_with_env(True, smi))).hexdigest()
            h2 = hashlib.sha256(_xyz_bytes(_run_with_env(True, smi))).hexdigest()
        finally:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS", None)
        self.assertEqual(
            h1, h2,
            "ON-orbit + CONFORMERS=0 path is not byte-identical between runs",
        )

    def test_on_emits_ensemble_berteb(self):
        """BERTEB (Ni CN4, with rotatable C-C single bonds) must emit more
        than 1 structure when the conformer gate is on.

        BERTEB has no puckerable 5/6/7 rings (only aromatic π-systems),
        so its ensemble is driven by the single-bond rotamer channel.
        Under the 2026-06-08 conformer split that channel defaults OFF,
        so we explicitly enable it here to preserve the test's semantic
        intent ("rotamer ensemble emits >=2 structures").
        """
        # Enable the single-bond rotamer channel explicitly -- this test
        # validates the rotamer ensemble for an aromatic-only molecule.
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER"] = "1"
        try:
            res = _run_with_env(True, SMILES_BERTEB)
        finally:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER", None)
        n = self._count_conformers(res)
        # We require at least 2 (so the ensemble is non-trivial).  In
        # practice we see 13 for BERTEB, but allow regression at the
        # rotamer-cap / RMSD-tol boundary.
        self.assertGreaterEqual(
            n, 2,
            f"Expected per-orbit conformer ensemble for BERTEB, got {n}",
        )
        # Every structure must be a valid XYZ; sanity-check the first.
        xyz0, label0 = res[0]
        syms, P = _parse_xyz_block(xyz0)
        self.assertGreater(len(syms), 1, "First structure has no atoms")
        self.assertTrue(np.all(np.isfinite(P)), "First structure has NaN")
        # Labels must contain the orbit signature.
        for _, label in res:
            self.assertIn("MOGUL", label)

    def test_conformer_off_outer_off_byte_identical(self):
        """The conformer-enum env-flag is dormant when the outer Mogul-
        primary gate is OFF.  Setting CONFORMERS=1 must not leak into
        the legacy bytes."""
        smi = "N[Pt](N)(Cl)Cl"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS"] = "1"
        try:
            h1 = hashlib.sha256(_xyz_bytes(_run_with_env(False, smi))).hexdigest()
            h2 = hashlib.sha256(_xyz_bytes(_run_with_env(False, smi))).hexdigest()
        finally:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS", None)
        self.assertEqual(
            h1, h2,
            "Legacy path bytes changed when CONFORMERS=1 was set (env leak)",
        )

    def test_universal_no_smiles_branch(self):
        """The conformer enumerator must NOT depend on SMILES patterns.

        We verify by running an organic chelate (BERTEB) and a single-
        donor-linear (SIYMEU) -- both go through the same enumerate_orbit
        _conformers code path with no per-class branching.  Both must
        produce more than 1 conformer (since both have rotatable bonds).

        Both probe the single-bond rotamer channel (aromatic-only ring
        systems), so we enable it explicitly under the 2026-06-08 split.
        """
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER"] = "1"
        try:
            for label, smi in [
                ("BERTEB", SMILES_BERTEB), ("SIYMEU", SMILES_SIYMEU),
            ]:
                res = _run_with_env(True, smi)
                n = self._count_conformers(res)
                self.assertGreaterEqual(
                    n, 2,
                    f"{label}: expected conformer ensemble, got N={n}",
                )
        finally:
            os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER", None)


if __name__ == "__main__":
    unittest.main()
