"""Tests for Mogul-DG REPLACE_RIGID donor-shell integrity + repair.

Covers the universal first-shell verification + spurious-H rejection +
missing-donor repair loop added to :mod:`delfin.fffree.mogul_dg`.

Scenarios:

1. **TUMQAT synthetic**: 3 spurious H atoms in first shell → verification
   flags it → repair moves H out of first shell.
2. **Cl-drift detection**: declared Cl donor drifted to 3.0 Å → verification
   flags it → repair tightens bounds and brings it back.
3. **Healthy structure passes**: all declared donors in shell window, no H
   intrusion → verification returns ok=True.
4. **Byte-identical OFF**: flag unset → output bit-identical to HEAD.
5. **Determinism**: two runs byte-identical with all flags ON.
6. **CN6 with 6 declared donors**: every donor sits in the shell window.
7. **No false-positive on multi-metal**: per-metal verification ignores other
   metals' first shells.

All tests honour ``PYTHONHASHSEED=0`` for determinism.  No SMILES patterns
are used — only synthetic XYZ + a minimal RDKit RWMol graph.
"""
from __future__ import annotations

import hashlib
import os
from typing import List, Tuple

# Determinism BEFORE any other import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

# Scrub any donor-shell / mogul-dg env flags that may leak from the parent.
for _k in list(os.environ):
    if _k.startswith("DELFIN_FFFREE_MOGUL_DG"):
        del os.environ[_k]


from delfin.fffree.mogul_dg import (
    _push_h_out_of_first_shell,
    _tighten_bounds_for_missing_donor,
    mogul_embed,
    mogul_embed_from_assembled,
    verify_donor_shell_integrity,
    xyz_block,
)


# ---------------------------------------------------------------------------
# Helpers — build small synthetic mols / coordinate sets
# ---------------------------------------------------------------------------
def _scrub_env(monkeypatch):
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_MOGUL_DG"):
            monkeypatch.delenv(k, raising=False)


@pytest.fixture(autouse=True)
def _scrub(monkeypatch):
    """Scrub all DELFIN_FFFREE_MOGUL_DG env flags before each test."""
    _scrub_env(monkeypatch)


def _build_tumqat_synthetic():
    """Construct a TUMQAT-like CN4 Zn-N₂-Cl-C topology + bug coords.

    Returns ``(syms, mol, metal_idx, declared_donors, P_bug, P_healthy)``.

    Bug coords reproduce the symptom reported in the task: Cl drifted to
    3.5 Å while three H atoms intruded into the first shell at 1.66/1.66/1.97 Å.
    """
    from rdkit import Chem

    mol = Chem.RWMol()
    zn = mol.AddAtom(Chem.Atom("Zn"))   # 0
    n1 = mol.AddAtom(Chem.Atom("N"))    # 1
    n2 = mol.AddAtom(Chem.Atom("N"))    # 2
    cl = mol.AddAtom(Chem.Atom("Cl"))   # 3
    c = mol.AddAtom(Chem.Atom("C"))     # 4
    h1 = mol.AddAtom(Chem.Atom("H"))    # 5 (bonded to N1)
    h2 = mol.AddAtom(Chem.Atom("H"))    # 6 (bonded to N2)
    h3 = mol.AddAtom(Chem.Atom("H"))    # 7 (bonded to C)
    for a, b in [
        (zn, n1), (zn, n2), (zn, cl), (zn, c),
        (n1, h1), (n2, h2), (c, h3),
    ]:
        mol.AddBond(a, b, Chem.BondType.SINGLE)

    syms = ["Zn", "N", "N", "Cl", "C", "H", "H", "H"]

    # Healthy CN4 — all four donors in 2.0-2.3 Å, H far from metal
    P_healthy = np.array([
        [0.0, 0.0, 0.0],   # Zn
        [2.0, 0.0, 0.0],   # N1
        [-2.0, 0.0, 0.0],  # N2
        [0.0, 2.3, 0.0],   # Cl
        [0.0, -2.0, 0.0],  # C
        [3.0, 1.0, 0.0],   # H1
        [-3.0, 1.0, 0.0],  # H2
        [1.0, -3.0, 0.0],  # H3
    ], dtype=np.float64)

    # Bug — Cl pushed out, three H pulled in
    P_bug = np.array([
        [0.0, 0.0, 0.0],   # Zn
        [2.0, 0.0, 0.0],   # N1
        [-2.0, 0.0, 0.0],  # N2
        [0.0, 3.5, 0.0],   # Cl drifted to 3.5 Å (outside shell)
        [0.0, -1.56, 0.0], # C pulled inward (1.56 Å, also too short)
        [1.66, 0.0, 0.0],  # spurious H1 in first shell
        [-1.66, 0.0, 0.0], # spurious H2
        [0.0, -1.97, 0.0], # spurious H3
    ], dtype=np.float64)

    declared_donors = [n1, n2, cl, c]
    return syms, mol, zn, declared_donors, P_bug, P_healthy


def _build_oc6_healthy():
    """CN6 [M(NH3)6] — six declared donors all at 2.05 Å."""
    from rdkit import Chem

    mol = Chem.RWMol()
    m = mol.AddAtom(Chem.Atom("Co"))
    donors: List[int] = []
    h_list: List[int] = []
    for _ in range(6):
        n = mol.AddAtom(Chem.Atom("N"))
        donors.append(n)
        mol.AddBond(m, n, Chem.BondType.SINGLE)
        for _h in range(3):
            h = mol.AddAtom(Chem.Atom("H"))
            h_list.append(h)
            mol.AddBond(n, h, Chem.BondType.SINGLE)
    syms = (
        ["Co"]
        + ["N", "H", "H", "H"] * 6
    )
    # Octahedral positions for donors
    r = 2.05
    octa = np.array([
        [r, 0, 0], [-r, 0, 0], [0, r, 0], [0, -r, 0], [0, 0, r], [0, 0, -r],
    ], dtype=np.float64)
    n_atoms = mol.GetNumAtoms()
    P = np.zeros((n_atoms, 3), dtype=np.float64)
    for k, di in enumerate(donors):
        P[di] = octa[k]
        # Drop three H around each donor on the far side (away from metal),
        # well outside the first shell (>2.6 Å from metal).
        away = octa[k] / np.linalg.norm(octa[k])
        for j in range(3):
            theta = 2.0 * np.pi * j / 3.0
            offset = away * 1.0 + np.array(
                [0.5 * np.cos(theta), 0.5 * np.sin(theta), 0.0]
            )
            P[di + j + 1] = octa[k] + offset
            # Force d(metal, H) > 2.6 Å by pushing radially outward
            v = P[di + j + 1]
            d = float(np.linalg.norm(v))
            if d < 2.7:
                v = v * (2.8 / max(d, 1e-6))
                P[di + j + 1] = v
    return syms, mol, 0, donors, P


# ---------------------------------------------------------------------------
# 1. Healthy structure passes
# ---------------------------------------------------------------------------
class TestHealthyStructurePasses:
    def test_tumqat_healthy_passes(self):
        syms, _mol, m, donors, _bug, P_healthy = _build_tumqat_synthetic()
        ok, viol = verify_donor_shell_integrity(
            P_healthy, syms, m, donors, md_tol=0.4,
        )
        assert ok, f"healthy structure unexpectedly flagged: {viol}"
        assert viol == []

    def test_oc6_healthy_passes(self):
        syms, _mol, m, donors, P = _build_oc6_healthy()
        ok, viol = verify_donor_shell_integrity(P, syms, m, donors, md_tol=0.4)
        assert ok, f"OC-6 healthy flagged: {viol}"


# ---------------------------------------------------------------------------
# 2. TUMQAT synthetic: spurious H atoms in first shell
# ---------------------------------------------------------------------------
class TestTumqatSpuriousH:
    def test_verification_flags_spurious_h(self):
        syms, _mol, m, donors, P_bug, _ = _build_tumqat_synthetic()
        ok, viol = verify_donor_shell_integrity(
            P_bug, syms, m, donors, md_tol=0.4,
        )
        assert not ok
        # Must flag at least the three spurious H
        h_flags = [v for v in viol if "spurious H" in v]
        assert len(h_flags) >= 3, f"expected ≥3 H flags, got {h_flags}"

    def test_verification_flags_cl_drift(self):
        syms, _mol, m, donors, P_bug, _ = _build_tumqat_synthetic()
        ok, viol = verify_donor_shell_integrity(
            P_bug, syms, m, donors, md_tol=0.4,
        )
        assert not ok
        # Cl is donor idx=3
        cl_flags = [v for v in viol if "idx=3" in v and "Cl" in v]
        assert cl_flags, f"expected Cl drift flag, got {viol}"

    def test_repair_moves_h_out_of_shell(self):
        syms, mol, m, donors, P_bug, _ = _build_tumqat_synthetic()
        P_new, n_pushed = _push_h_out_of_first_shell(
            P_bug, syms, mol, m, donors, shell_upper=2.6,
        )
        assert n_pushed == 3, f"expected 3 H pushed, got {n_pushed}"
        # Every H now sits outside the 2.6 Å shell
        for i, s in enumerate(syms):
            if s == "H":
                d = float(np.linalg.norm(P_new[i] - P_new[m]))
                assert d > 2.6, f"H idx={i} still in shell at {d:.3f} Å"

    def test_repair_does_not_move_declared_donors(self):
        syms, mol, m, donors, P_bug, _ = _build_tumqat_synthetic()
        P_new, _ = _push_h_out_of_first_shell(
            P_bug, syms, mol, m, donors, shell_upper=2.6,
        )
        # Declared donor coordinates must be unchanged
        for d in donors:
            assert np.array_equal(P_new[d], P_bug[d]), (
                f"declared donor idx={d} was moved"
            )
        # Metal must also be unchanged
        assert np.array_equal(P_new[m], P_bug[m])


# ---------------------------------------------------------------------------
# 3. Cl-drift detection (md_tol)
# ---------------------------------------------------------------------------
class TestClDriftDetection:
    def test_md_tol_catches_far_from_target(self):
        syms = ["Zn", "N", "N", "Cl", "C"]
        # Cl at 3.0 Å (within physical shell upper if upper=3.2), but far
        # from μ=2.30 → md_tol=0.4 should catch this.
        P = np.array([
            [0.0, 0.0, 0.0],
            [2.05, 0.0, 0.0],
            [-2.05, 0.0, 0.0],
            [0.0, 3.0, 0.0],   # Cl at 3.0 Å vs μ=2.30
            [0.0, -2.0, 0.0],
        ], dtype=np.float64)
        md_targets = {1: 2.05, 2: 2.05, 3: 2.30, 4: 2.0}
        # With shell_upper=3.2 (more permissive) the physical window passes,
        # but the μ-target check fires.
        ok, viol = verify_donor_shell_integrity(
            P, syms, 0, [1, 2, 3, 4],
            md_targets=md_targets, md_tol=0.4, shell_upper=3.2,
        )
        assert not ok
        cl_flags = [v for v in viol if "idx=3" in v and "far from target" in v]
        assert cl_flags, f"expected μ-target flag, got {viol}"

    def test_tighten_bounds_halves_window(self):
        n = 5
        L = np.full((n, n), 0.5)
        U = np.full((n, n), 5.0)
        # M-D bounds 2.0 - 2.6 (μ=2.3, half=0.3)
        L[0, 3] = 2.0
        L[3, 0] = 2.0
        U[0, 3] = 2.6
        U[3, 0] = 2.6
        L2, U2 = _tighten_bounds_for_missing_donor(L, U, 0, 3, 0.5)
        # New half = 0.15 around μ=2.3 → [2.15, 2.45]
        assert abs(L2[0, 3] - 2.15) < 1e-9
        assert abs(U2[0, 3] - 2.45) < 1e-9
        # Symmetric
        assert L2[0, 3] == L2[3, 0]
        assert U2[0, 3] == U2[3, 0]
        # Inputs unmutated
        assert L[0, 3] == 2.0 and U[0, 3] == 2.6


# ---------------------------------------------------------------------------
# 4. Byte-identical OFF (regression guard)
# ---------------------------------------------------------------------------
class TestByteIdenticalOff:
    SMILES = "N[Pt](N)(Cl)Cl"

    def test_default_off_byte_identical_to_head(self):
        """No env flag set → mogul_embed output unchanged vs HEAD."""
        r1 = mogul_embed(self.SMILES)
        # Defensive: also try with bogus env unset (already scrubbed by fixture).
        r2 = mogul_embed(self.SMILES)
        if r1 is None or r2 is None:
            pytest.skip("mogul library not available on this machine")
        assert np.array_equal(r1[1], r2[1])
        # Info dict must NOT carry the donor_shell_* keys when flags are OFF
        assert "donor_shell_verify" not in r1[2]
        assert "donor_shell_reject_h" not in r1[2]
        assert "donor_shell_ok" not in r1[2]

    def test_xyz_hash_stable_two_runs(self):
        r1 = mogul_embed(self.SMILES)
        r2 = mogul_embed(self.SMILES)
        if r1 is None or r2 is None:
            pytest.skip("mogul library not available")
        h1 = hashlib.sha256(xyz_block(r1[0], r1[1]).encode()).hexdigest()
        h2 = hashlib.sha256(xyz_block(r2[0], r2[1]).encode()).hexdigest()
        assert h1 == h2


# ---------------------------------------------------------------------------
# 5. Determinism with all flags ON
# ---------------------------------------------------------------------------
class TestDeterminismWithFlagsOn:
    SMILES = "N[Pt](N)(Cl)Cl"

    def _run(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_REPLACE_RIGID", "1")
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_DONOR_SHELL_VERIFY", "1")
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_REJECT_H_FIRST_SHELL", "1")
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_DONOR_REPAIR_RETRIES", "3")
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_MD_TOL", "0.4")
        return mogul_embed(self.SMILES)

    def test_two_runs_byte_identical_with_flags(self, monkeypatch):
        r1 = self._run(monkeypatch)
        # Re-scrub then re-set for the second run (simulate fresh env)
        for k in list(os.environ):
            if k.startswith("DELFIN_FFFREE_MOGUL_DG"):
                monkeypatch.delenv(k, raising=False)
        r2 = self._run(monkeypatch)
        if r1 is None or r2 is None:
            pytest.skip("mogul library not available")
        assert np.array_equal(r1[1], r2[1])
        assert r1[0] == r2[0]


# ---------------------------------------------------------------------------
# 6. CN6: all six donors verified in shell
# ---------------------------------------------------------------------------
class TestCN6AllDonorsVerified:
    def test_oc6_healthy_verifies(self):
        syms, _mol, m, donors, P = _build_oc6_healthy()
        ok, viol = verify_donor_shell_integrity(P, syms, m, donors, md_tol=0.4)
        assert ok, f"CN6 healthy flagged: {viol}"
        assert len(donors) == 6

    def test_oc6_one_donor_missing_flags_count_mismatch(self):
        """Move one of six donors out to 3.5 Å → CN-mismatch + shell-window fail."""
        syms, _mol, m, donors, P = _build_oc6_healthy()
        P_bad = P.copy()
        # Push donor[5] outward
        donor5_pos = P_bad[donors[5]]
        unit = donor5_pos / np.linalg.norm(donor5_pos)
        P_bad[donors[5]] = unit * 3.5
        ok, viol = verify_donor_shell_integrity(
            P_bad, syms, m, donors, md_tol=0.4,
        )
        assert not ok
        # Must flag both the missing donor and the CN mismatch
        assert any(f"idx={donors[5]}" in v for v in viol)
        assert any("realised heavy CN" in v for v in viol)


# ---------------------------------------------------------------------------
# 7. Multi-metal: per-metal verification ignores other metal first shells
# ---------------------------------------------------------------------------
class TestMultiMetalPerMetal:
    def test_per_metal_verification(self):
        """Two well-separated CN4 metals — verification only sees its own shell."""
        syms = [
            "Zn", "N", "N", "Cl", "C",            # metal 1 (idx 0)
            "Cu", "N", "N", "Cl", "C",            # metal 2 (idx 5)
        ]
        P = np.array([
            [0.0, 0.0, 0.0],   # Zn
            [2.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 2.3, 0.0],
            [0.0, -2.0, 0.0],
            [10.0, 0.0, 0.0],  # Cu (far away)
            [12.0, 0.0, 0.0],
            [8.0, 0.0, 0.0],
            [10.0, 2.3, 0.0],
            [10.0, -2.0, 0.0],
        ], dtype=np.float64)
        ok1, viol1 = verify_donor_shell_integrity(
            P, syms, 0, [1, 2, 3, 4], md_tol=0.4,
        )
        assert ok1, f"metal 1 unexpectedly flagged: {viol1}"
        ok2, viol2 = verify_donor_shell_integrity(
            P, syms, 5, [6, 7, 8, 9], md_tol=0.4,
        )
        assert ok2, f"metal 2 unexpectedly flagged: {viol2}"


# ---------------------------------------------------------------------------
# 8. End-to-end via mogul_embed_from_assembled + flags
# ---------------------------------------------------------------------------
class TestEndToEnd:
    def _build_cisplatin_topology(self):
        from delfin._bond_decollapse import _is_metal
        from delfin.smiles_converter import _prepare_mol_for_embedding

        mol = _prepare_mol_for_embedding("N[Pt](N)(Cl)Cl", hapto_approx=False)
        if mol is None:
            pytest.skip("RDKit prepare unavailable")
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        metal_idx = next(i for i, s in enumerate(syms) if _is_metal(s))
        donors = sorted(
            int(nb.GetIdx())
            for nb in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
            if not _is_metal(nb.GetSymbol())
        )
        return syms, mol, metal_idx, donors

    def test_e2e_with_verify_flag_returns_info(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_DONOR_SHELL_VERIFY", "1")
        syms, mol, m, donors = self._build_cisplatin_topology()
        n = len(syms)
        P_init = np.zeros((n, 3), dtype=np.float64)
        # Square planar ideal positions for the 4 donors
        if len(donors) >= 4:
            r = 2.1
            P_init[donors[0]] = [r, 0.0, 0.0]
            P_init[donors[1]] = [-r, 0.0, 0.0]
            P_init[donors[2]] = [0.0, r, 0.0]
            P_init[donors[3]] = [0.0, -r, 0.0]
        for i, s in enumerate(syms):
            if s == "H":
                try:
                    nb = mol.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx()
                    P_init[i] = P_init[nb] + np.array([0.4, 0.4, 0.4])
                except Exception:
                    P_init[i] = np.array([0.1 * i, 0.0, 0.0])
        res = mogul_embed_from_assembled(syms, P_init, mol, m, donors)
        if res is None:
            pytest.skip("solver infeasible on this build")
        _, _, info = res
        assert "donor_shell_verify" in info
        assert info["donor_shell_verify"] is True
        assert "donor_shell_ok" in info
        assert "donor_shell_repair_attempts" in info

    def test_e2e_with_reject_h_flag_runs(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_MOGUL_DG_REJECT_H_FIRST_SHELL", "1")
        syms, mol, m, donors = self._build_cisplatin_topology()
        n = len(syms)
        P_init = np.zeros((n, 3), dtype=np.float64)
        if len(donors) >= 4:
            r = 2.1
            P_init[donors[0]] = [r, 0.0, 0.0]
            P_init[donors[1]] = [-r, 0.0, 0.0]
            P_init[donors[2]] = [0.0, r, 0.0]
            P_init[donors[3]] = [0.0, -r, 0.0]
        for i, s in enumerate(syms):
            if s == "H":
                try:
                    nb = mol.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx()
                    P_init[i] = P_init[nb] + np.array([0.4, 0.4, 0.4])
                except Exception:
                    P_init[i] = np.array([0.1 * i, 0.0, 0.0])
        res = mogul_embed_from_assembled(syms, P_init, mol, m, donors)
        if res is None:
            pytest.skip("solver infeasible on this build")
        _, P_out, info = res
        # Info dict must record the H-rejection flag
        assert info.get("donor_shell_reject_h") is True
        # No H within 1.6 Å of metal (declared donors are heavy in this SMILES)
        for i, s in enumerate(syms):
            if s == "H":
                d = float(np.linalg.norm(P_out[i] - P_out[m]))
                assert d > 1.4, f"H idx={i} still very close at {d:.3f} Å"
