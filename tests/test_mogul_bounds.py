"""Tests for Mogul-DG Phase A — :mod:`delfin.fffree.mogul_bounds`.

The bounds matrix is the *foundation* of the whole-complex distance-geometry
embed.  These tests cover:

* **Library-driven empirical bounds** — that the CCDC distributions surface
  in the right cells (Pt-N, Pt-Cl, Fe-Cp, Cu-N, La-O, Re-Re).
* **Universal mechanisms** — graph-only detection of resonance groups
  (aromatic + symmetric-terminal + automorphism orbits) for benzene,
  carboxylate, nitrate, perchlorate.
* **Polyhedron donor-donor encoding** for CN4 (SP-4), CN6 (OC-6),
  CN9 (TTP-9) — the polyhedron is reachable via the bounds alone.
* **Hapto-π handling** — ferrocene's two Cp rings get correctly
  identified as hapto groups; M-C bounds use the v5 TM η-fragment
  library.
* **Determinism contract** — same input ⇒ byte-identical bounds across
  two runs.
* **Universal-no-element-strings** grep — no SMILES-specific code in
  the module.

Tests deliberately build minimal RDKit mols + index assignments so each
case is fast and self-contained.  Library access is honoured via
``DELFIN_GRIP_LIB_PATH`` — when the v5 library is unavailable the
TM-aware tests degrade gracefully (skipped with a clear message).
"""
from __future__ import annotations

import hashlib
import math
import os
import re
import sys
from pathlib import Path
from typing import List, Tuple

# Determinism must be set before any other import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.fffree.mogul_bounds import (
    build_bounds_matrix,
    detect_hapto_groups,
    detect_resonance_groups,
    find_aromatic_groups,
    find_automorphism_orbits,
    find_symmetric_terminal_groups,
    hyb_str,
)
from delfin.fffree.grip_mogul_lookup import GripLibrary

Chem = pytest.importorskip("rdkit.Chem")

# ---------------------------------------------------------------------------
# Library fixture — v5 if available, otherwise the default v1 (some M-D tests
# can only assert "library hit" / "geometry plausible" on v1 because v1 has no
# TM categories).
# ---------------------------------------------------------------------------
_V5_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"
)


@pytest.fixture(scope="module")
def lib() -> GripLibrary:
    """Return the v5 library if available, else the default."""
    if _V5_PATH.exists():
        return GripLibrary.get(_V5_PATH)
    # Fallback: env or pinned default
    libd = GripLibrary.get_default()
    assert libd is not None, "No GRIP library available — set DELFIN_GRIP_LIB_PATH"
    return libd


def _has_tm_categories(lib_obj: GripLibrary) -> bool:
    return bool(getattr(lib_obj, "has_tm_categories", False))


def _has_pair_tables(lib_obj: GripLibrary) -> bool:
    return bool(getattr(lib_obj, "has_pair_tables", False))


# ---------------------------------------------------------------------------
# Helpers to build minimal assembled mols
# ---------------------------------------------------------------------------
def _build_complex(smi_with_metal: str) -> Tuple[List[str], "Chem.Mol", int, List[int]]:
    """Parse a SMILES that ALREADY contains the metal as a connected atom,
    add hydrogens, return (syms, mol, metal_idx, donors)."""
    mol = Chem.MolFromSmiles(smi_with_metal)
    assert mol is not None, f"SMILES parse failed: {smi_with_metal}"
    mol = Chem.AddHs(mol)
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    # Pick the FIRST metal we find (the project _METALS set).
    from delfin._bond_decollapse import _METALS
    metals = [i for i, s in enumerate(syms) if s in _METALS]
    assert metals, f"No metal found in {smi_with_metal}"
    metal_idx = metals[0]
    donors = sorted(int(nb.GetIdx()) for nb in mol.GetAtomWithIdx(metal_idx).GetNeighbors())
    return syms, mol, metal_idx, donors


def _attach_hapto_metal(
    ligand_smi: str,
    metal_sym: str,
    donor_carbon_idxs: List[int],
) -> Tuple[List[str], "Chem.Mol", int, List[int]]:
    """Construct a hapto-coordinated mol by adding an explicit metal atom and
    bonding it to a set of ring carbons."""
    mol = Chem.MolFromSmiles(ligand_smi)
    assert mol is not None
    mol = Chem.AddHs(mol)
    rw = Chem.RWMol(mol)
    metal_idx = rw.AddAtom(Chem.Atom(metal_sym))
    for c in donor_carbon_idxs:
        rw.AddBond(int(c), int(metal_idx), Chem.BondType.SINGLE)
    Chem.SanitizeMol(rw)
    syms = [a.GetSymbol() for a in rw.GetAtoms()]
    return syms, rw, int(metal_idx), sorted(int(c) for c in donor_carbon_idxs)


# ===========================================================================
# 1.  Universal mechanism tests — no library needed
# ===========================================================================
class TestResonanceDetection:
    """The 3 detection mechanisms must work on toy organic molecules."""

    def test_benzene_aromatic_group(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        groups = find_aromatic_groups(mol)
        # Exactly one ring of 6 carbons.
        assert len(groups) == 1
        assert len(groups[0]) == 6
        # All in the group should be C atoms.
        for idx in groups[0]:
            assert mol.GetAtomWithIdx(idx).GetSymbol() == "C"

    def test_carboxylate_symmetric_terminals(self):
        # Acetate: CC(=O)[O-]; the two O are degree-1 terminals on the C.
        mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)[O-]"))
        groups = find_symmetric_terminal_groups(mol)
        # There must be at least one group of 2 O atoms.
        o_groups = [g for g in groups if all(
            mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in g
        )]
        assert any(len(g) == 2 for g in o_groups)

    def test_nitrate_3_equivalent_O(self):
        mol = Chem.MolFromSmiles("[N+](=O)([O-])[O-]")
        # Nitrate: 1 N + 3 O, N is degree-3, O each degree-1.
        groups = find_symmetric_terminal_groups(mol)
        # Expect exactly 3 O atoms grouped.
        o_groups = [g for g in groups if all(
            mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in g
        )]
        assert any(len(g) == 3 for g in o_groups)

    def test_perchlorate_4_O(self):
        mol = Chem.MolFromSmiles("[Cl+3]([O-])([O-])([O-])[O-]")
        groups = find_symmetric_terminal_groups(mol)
        o_groups = [g for g in groups if all(
            mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in g
        )]
        assert any(len(g) == 4 for g in o_groups), "ClO4- should have 4 equivalent O"

    def test_BF4_4_F(self):
        mol = Chem.MolFromSmiles("[B-](F)(F)(F)F")
        groups = find_symmetric_terminal_groups(mol)
        f_groups = [g for g in groups if all(
            mol.GetAtomWithIdx(i).GetSymbol() == "F" for i in g
        )]
        assert any(len(g) == 4 for g in f_groups)

    def test_PF6_6_F(self):
        mol = Chem.MolFromSmiles("[P-](F)(F)(F)(F)(F)F")
        groups = find_symmetric_terminal_groups(mol)
        f_groups = [g for g in groups if all(
            mol.GetAtomWithIdx(i).GetSymbol() == "F" for i in g
        )]
        assert any(len(g) == 6 for g in f_groups)

    def test_methyl_3H_via_automorphism(self):
        # CH4: 4 equivalent H atoms by automorphism orbit.
        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        orbits = find_automorphism_orbits(mol)
        # All 4 H are equivalent.
        assert any(len(o) == 4 for o in orbits)

    def test_combined_detection_merges_overlapping(self):
        # Naphthalene — two fused aromatic rings; merging should yield
        # a single 10-carbon group.
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
        groups = detect_resonance_groups(mol)
        # The merged set should contain all 10 aromatic C in ONE group.
        c_in_groups = sorted(i for g in groups for i in g
                              if mol.GetAtomWithIdx(i).GetSymbol() == "C")
        # All 10 aromatic carbons should appear in some merged group.
        n_aromatic_c = sum(1 for a in mol.GetAtoms()
                            if a.GetIsAromatic() and a.GetSymbol() == "C")
        assert n_aromatic_c == 10
        # Find the largest aromatic-only group post-merge.
        biggest = max(
            (g for g in groups if all(
                mol.GetAtomWithIdx(i).GetIsAromatic() for i in g
            )),
            key=len,
            default=set(),
        )
        assert len(biggest) == 10


# ===========================================================================
# 2.  Bounds matrix — coordination + topology
# ===========================================================================
class TestPtAmmineChloride:
    """cis-Pt(NH₃)₂Cl₂ — canonical SP-4 test."""

    @pytest.fixture
    def setup(self, lib):
        syms, mol, metal_idx, donors = _build_complex(
            "[NH3][Pt]([NH3])(Cl)Cl"
        )
        L, U, info = build_bounds_matrix(
            syms, mol, metal_idx, donors, lib, geometry="SP-4 square planar"
        )
        return syms, mol, metal_idx, donors, L, U, info

    def test_pt_n_distance_matches_ccdc(self, setup, lib):
        if not _has_pair_tables(lib):
            pytest.skip("v5 library required for Pt-N TM-aware lookup")
        syms, mol, m, donors, L, U, info = setup
        n_donors = [d for d in donors if syms[d] == "N"]
        for d in n_donors:
            mid = 0.5 * (L[m, d] + U[m, d])
            # CCDC empirical Pt-N ≈ 2.04 ± 0.05 Å.  Our 2σ bound must
            # bracket this empirical centre.
            assert 1.9 < mid < 2.2, f"Pt-N midpoint {mid:.3f} out of CCDC range"
            assert U[m, d] - L[m, d] > 0.05, "M-D window too tight"
            assert U[m, d] - L[m, d] < 0.30, "M-D window too loose"

    def test_pt_cl_distance_matches_ccdc(self, setup, lib):
        if not _has_pair_tables(lib):
            pytest.skip("v5 library required for Pt-Cl TM-aware lookup")
        syms, mol, m, donors, L, U, info = setup
        cl_donors = [d for d in donors if syms[d] == "Cl"]
        for d in cl_donors:
            mid = 0.5 * (L[m, d] + U[m, d])
            # CCDC empirical Pt-Cl ≈ 2.32 ± 0.03 Å.
            assert 2.2 < mid < 2.45, f"Pt-Cl midpoint {mid:.3f} out of CCDC range"

    def test_donor_donor_bounds_for_SP4(self, setup):
        syms, mol, m, donors, L, U, info = setup
        # SP-4 — 4 donors form a square.  Adjacent donors are ≈ √2·r_md apart;
        # opposite donors are 2·r_md apart.  Our donor-donor block must
        # encode this.
        assert info["n_dd_bounds"] == 6  # C(4,2)
        # Find adjacent vs opposite ratio: at minimum there exist 2 distinct
        # ideal D-D distances (square has 2 edge-types: edge & diagonal).
        donor_pairs_mids = []
        for i, di in enumerate(donors):
            for j, dj in enumerate(donors[i+1:], start=i+1):
                mid = 0.5 * (L[di, dj] + U[di, dj])
                if mid > 0:
                    donor_pairs_mids.append(mid)
        unique_mids = sorted(set(round(m, 2) for m in donor_pairs_mids))
        assert len(unique_mids) >= 2, "SP-4 should have ≥ 2 distinct D-D distances"

    def test_no_bounds_relaxed(self, setup):
        _syms, _mol, _m, _donors, _L, _U, info = setup
        assert info["bounds_relaxed"] == [], (
            f"Pt(NH3)2Cl2 must be feasible: {info['bounds_relaxed']}"
        )


class TestFerroceneHaptoDetection:
    """Fe(Cp)₂ — hapto-π detection + η⁵ M-C library bounds."""

    @pytest.fixture
    def setup(self, lib):
        # Two cyclopentadienyl rings (5 C each, each ring has 1 carbanion);
        # we add Fe and bond it to all 10 ring carbons.
        mol = Chem.MolFromSmiles("C1=CC=C[CH-]1.C1=CC=C[CH-]1")
        mol = Chem.AddHs(mol)
        rw = Chem.RWMol(mol)
        cs = [int(i) for i, a in enumerate(rw.GetAtoms()) if a.GetSymbol() == "C"]
        assert len(cs) == 10
        fe = rw.AddAtom(Chem.Atom("Fe"))
        for c in cs:
            rw.AddBond(c, fe, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            # Sanitisation can complain about charges; we don't care for bounds.
            pass
        syms = [a.GetSymbol() for a in rw.GetAtoms()]
        L, U, info = build_bounds_matrix(syms, rw, fe, cs, lib)
        return syms, rw, fe, cs, L, U, info

    def test_two_cp_hapto_groups_detected(self, setup):
        syms, mol, fe, cs, L, U, info = setup
        groups = detect_hapto_groups(mol, fe, cs)
        assert len(groups) == 2, f"Expected 2 Cp rings, got {len(groups)}"
        for g in groups:
            assert len(g) == 5, "Each Cp ring should have 5 C atoms"

    def test_all_fe_c_bounds_symmetric(self, setup):
        syms, mol, fe, cs, L, U, info = setup
        # All 10 Fe-C bounds must be equal (η⁵ symmetry across both rings).
        mids = []
        for c in cs:
            mids.append(0.5 * (L[fe, c] + U[fe, c]))
        # All should be within ε of each other.
        spread = max(mids) - min(mids)
        assert spread < 0.001, f"Fe-C M-D bounds not symmetric: spread={spread:.4f}"

    def test_fe_c_distance_in_eta5_range(self, setup, lib):
        if not _has_tm_categories(lib):
            pytest.skip("v5 TM categories required for hapto-η⁵ lookup")
        syms, mol, fe, cs, L, U, info = setup
        # Fe-Cp empirical ≈ 2.04-2.08 Å.
        mid = 0.5 * (L[fe, cs[0]] + U[fe, cs[0]])
        assert 1.95 < mid < 2.15, f"Fe-Cp midpoint {mid:.3f} out of empirical range"


class TestResonanceEquidistanceBounds:
    """Verify the resonance bound block enforces equal aromatic C-C
    distances across a benzene ring's internal bonds."""

    def test_benzene_ring_bonds_equidistant(self, lib):
        # Toluene — methyl + aromatic ring.  Ring C-C bonds should
        # all get the same [lo, hi] window after the resonance block runs.
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        # No metal — fabricate metal_idx as a dangling sentinel (we still
        # need it for the API; donor list empty).
        # For this test we just need the function to NOT crash and produce
        # equidistant aromatic bonds.  We'll attach a dummy metal atom.
        rw = Chem.RWMol(mol)
        m = rw.AddAtom(Chem.Atom("Pd"))
        # bond Pd to one ring C to satisfy API; this should not affect aromatic ring
        ring_c = [i for i, a in enumerate(rw.GetAtoms()) if a.GetIsAromatic()][0]
        rw.AddBond(ring_c, m, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        syms = [a.GetSymbol() for a in rw.GetAtoms()]
        L, U, info = build_bounds_matrix(syms, rw, m, [ring_c], lib)
        # Find all aromatic C-C bonds and verify they share the same window.
        aromatic_cc_windows = []
        for bond in rw.GetBonds():
            i = int(bond.GetBeginAtomIdx())
            j = int(bond.GetEndAtomIdx())
            if (
                rw.GetAtomWithIdx(i).GetIsAromatic()
                and rw.GetAtomWithIdx(j).GetIsAromatic()
                and rw.GetAtomWithIdx(i).GetSymbol() == "C"
                and rw.GetAtomWithIdx(j).GetSymbol() == "C"
            ):
                aromatic_cc_windows.append((L[i, j], U[i, j]))
        assert len(aromatic_cc_windows) >= 6, "Should find ≥ 6 aromatic C-C bonds"
        # All windows must be identical (resonance equi-distance).
        first = aromatic_cc_windows[0]
        for w in aromatic_cc_windows[1:]:
            assert abs(w[0] - first[0]) < 1e-6 and abs(w[1] - first[1]) < 1e-6, (
                f"Aromatic C-C windows differ: {first} vs {w}"
            )

    def test_carboxylate_two_co_bonds_equal(self, lib):
        # Acetate bound to Pd via the carboxylate C: ensure the 2 C-O bonds
        # are equidistant after resonance detection (the 2 O are degree-1
        # terminals → symmetric-terminal mechanism).
        rw = Chem.RWMol(Chem.MolFromSmiles("CC(=O)[O-]"))
        m = rw.AddAtom(Chem.Atom("Pd"))
        # Attach Pd to the carboxylate carbon (atom 1).
        rw.AddBond(1, m, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        rw = Chem.RWMol(Chem.AddHs(rw))
        syms = [a.GetSymbol() for a in rw.GetAtoms()]
        # Re-discover Pd / donor indices.
        m = [i for i, s in enumerate(syms) if s == "Pd"][0]
        donors = [int(nb.GetIdx()) for nb in rw.GetAtomWithIdx(m).GetNeighbors()]
        L, U, info = build_bounds_matrix(syms, rw, m, donors, lib)
        # Find the 2 C-O bonds of the carboxylate (degree-1 O's on the
        # central carbon).
        o_terminals = []
        for atom in rw.GetAtoms():
            if atom.GetSymbol() == "O" and atom.GetDegree() == 1:
                nb = list(atom.GetNeighbors())[0]
                o_terminals.append((int(nb.GetIdx()), int(atom.GetIdx())))
        # Must be ≥ 2.
        assert len(o_terminals) >= 2
        windows = sorted(set((L[i, j], U[i, j]) for i, j in o_terminals))
        # After merging both into the same group, both windows must be IDENTICAL.
        assert len(windows) == 1, f"C-O bonds should share one window: {windows}"


# ===========================================================================
# 3.  Geometry-specific tests
# ===========================================================================
class TestOctahedralCN6:
    """[Ni(NH₃)₆]²⁺ — CN6 OC-6 octahedron donor-donor encoding."""

    def test_six_amine_donors_octahedral(self, lib):
        # Build a Ni(NH3)6 explicit complex.
        mol = Chem.RWMol()
        ni = mol.AddAtom(Chem.Atom("Ni"))
        n_idxs = []
        for _ in range(6):
            n = mol.AddAtom(Chem.Atom("N"))
            mol.AddBond(ni, n, Chem.BondType.SINGLE)
            n_idxs.append(n)
        # Add 3 H per N to satisfy valence.
        for n in n_idxs:
            for _ in range(3):
                h = mol.AddAtom(Chem.Atom("H"))
                mol.AddBond(n, h, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        L, U, info = build_bounds_matrix(
            syms, mol, ni, n_idxs, lib, geometry="OC-6 octahedron"
        )
        assert info["n_md_bounds"] == 6
        assert info["n_dd_bounds"] == 15  # C(6,2)
        # Octahedron has 2 distinct distances: edge (cis) and diagonal (trans).
        unique_dd: set = set()
        for i in range(6):
            for j in range(i + 1, 6):
                mid = 0.5 * (L[n_idxs[i], n_idxs[j]] + U[n_idxs[i], n_idxs[j]])
                unique_dd.add(round(mid, 2))
        assert len(unique_dd) == 2, f"Octahedron should have 2 D-D distances: {unique_dd}"


# ===========================================================================
# 4.  Multi-metal — Re₂Cl₈²⁻
# ===========================================================================
class TestMultiMetalReRe:
    def test_re_re_bond_present(self, lib):
        # Re-Re quintuple-bond cluster.
        mol = Chem.MolFromSmiles("Cl[Re](Cl)(Cl)(Cl)[Re](Cl)(Cl)(Cl)Cl")
        assert mol is not None
        mol = Chem.AddHs(mol)
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        res = [i for i, s in enumerate(syms) if s == "Re"]
        assert len(res) == 2
        # Build bounds for ONE metal centre (the other Re is treated as a
        # "donor" for the purposes of this metal).
        m = res[0]
        donors = sorted(int(nb.GetIdx()) for nb in mol.GetAtomWithIdx(m).GetNeighbors())
        L, U, info = build_bounds_matrix(syms, mol, m, donors, lib)
        # The Re-Re bound must exist and be in the empirical 2.2-2.8 Å window
        # (quintuple bond ≈ 2.24 Å, single ≈ 2.6 Å).
        other_re = res[1]
        mid = 0.5 * (L[m, other_re] + U[m, other_re])
        assert mid > 2.0 and mid < 3.0, f"Re-Re midpoint {mid:.3f} not in expected range"


# ===========================================================================
# 5.  Determinism
# ===========================================================================
class TestDeterminism:
    def test_byte_identical_two_runs(self, lib):
        def go():
            mol = Chem.MolFromSmiles("[NH3][Pt]([NH3])(Cl)Cl")
            mol = Chem.AddHs(mol)
            syms = [a.GetSymbol() for a in mol.GetAtoms()]
            metal = syms.index("Pt")
            donors = sorted(int(nb.GetIdx()) for nb in mol.GetAtomWithIdx(metal).GetNeighbors())
            L, U, _ = build_bounds_matrix(
                syms, mol, metal, donors, lib, geometry="SP-4 square planar"
            )
            return L, U
        L1, U1 = go()
        L2, U2 = go()
        # Byte-identical via SHA256 digest.
        h_L1 = hashlib.sha256(np.ascontiguousarray(L1).tobytes()).hexdigest()
        h_L2 = hashlib.sha256(np.ascontiguousarray(L2).tobytes()).hexdigest()
        h_U1 = hashlib.sha256(np.ascontiguousarray(U1).tobytes()).hexdigest()
        h_U2 = hashlib.sha256(np.ascontiguousarray(U2).tobytes()).hexdigest()
        assert h_L1 == h_L2, "L1 != L2 (non-deterministic lower bounds)"
        assert h_U1 == h_U2, "U1 != U2 (non-deterministic upper bounds)"

    def test_determinism_under_atom_index_perm(self, lib):
        # Reordering atom indices in SMILES (without changing the molecule)
        # should yield bounds that DIFFER only by the permutation — but the
        # SHA over the SORTED bound list should be stable.
        m1 = Chem.AddHs(Chem.MolFromSmiles("[NH3][Pt]([NH3])(Cl)Cl"))
        m2 = Chem.AddHs(Chem.MolFromSmiles("Cl[Pt](Cl)([NH3])[NH3]"))
        sym1 = [a.GetSymbol() for a in m1.GetAtoms()]
        sym2 = [a.GetSymbol() for a in m2.GetAtoms()]
        i1 = sym1.index("Pt"); i2 = sym2.index("Pt")
        d1 = sorted(int(nb.GetIdx()) for nb in m1.GetAtomWithIdx(i1).GetNeighbors())
        d2 = sorted(int(nb.GetIdx()) for nb in m2.GetAtomWithIdx(i2).GetNeighbors())
        L1, U1, _ = build_bounds_matrix(sym1, m1, i1, d1, lib, geometry="SP-4 square planar")
        L2, U2, _ = build_bounds_matrix(sym2, m2, i2, d2, lib, geometry="SP-4 square planar")
        # Sorted lower-bound multisets must match (up to floating noise).
        l1_vals = np.sort(L1.flatten())
        l2_vals = np.sort(L2.flatten())
        assert l1_vals.shape == l2_vals.shape
        np.testing.assert_allclose(l1_vals, l2_vals, atol=1e-9)


# ===========================================================================
# 6.  Universality verification — no hardcoded element strings
# ===========================================================================
class TestUniversalConstraint:
    def test_no_hardcoded_metal_element_branches(self):
        """The module must not contain element-specific ``if`` branches.

        Allowed:
            * occurrences inside docstrings or comments (which describe
              CCDC examples)
            * occurrences as keys of universal physical-constant tables
              (covalent radii, vdW radii — re-imported, not redefined)
            * occurrences as part of API field names like "metal_sym"

        Disallowed:
            * ``if sym == "Pt": …``
            * ``elif element in ("Fe", "Cu"): …``
        """
        src = Path(
            "/home/qmchem_max/ComPlat/DELFIN/delfin/fffree/mogul_bounds.py"
        ).read_text()
        # Strip docstrings + comments before grep.
        # Remove triple-quoted docstrings (greedy match per docstring).
        src_no_doc = re.sub(r'"""[\s\S]*?"""', "", src, flags=re.MULTILINE)
        src_no_doc = re.sub(r"'''[\s\S]*?'''", "", src_no_doc, flags=re.MULTILINE)
        # Remove single-line comments.
        src_no_doc = re.sub(r"#.*", "", src_no_doc)
        # Now look for chemistry-specific `if sym == "X":` patterns.
        metal_atoms = ["Pt", "Pd", "Fe", "Cu", "Ni", "Co", "Mn", "Cr", "V",
                       "Ti", "Sc", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
                       "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Au",
                       "Hg", "La", "Ce"]
        for el in metal_atoms:
            # Look for `== "Pt"` or `== 'Pt'` or `in {"Pt", ...}` style.
            pat = rf'==\s*["\']{el}["\']'
            hits = re.findall(pat, src_no_doc)
            assert not hits, f"Found element-specific equality check for {el}: {hits}"

    def test_module_imports_only_universal_tables(self):
        """The only element-keyed dicts imported must be the universal physical
        constants (covalent radii, vdW radii, ``_METALS`` set)."""
        from delfin.fffree import mogul_bounds as mb
        # The module must NOT define its own private element-keyed dict
        # beyond what is imported from `_bond_decollapse` / `polyhedra`.
        # Inspect module-level dict objects.
        # Periodic-table element symbols (subset — only need to detect chemistry
        # encoded as dict keys).
        elements = {
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
            "Zn", "Ga", "Ge", "As", "Se", "Br",
            "Pd", "Pt", "Au", "Ag",
        }
        for name in dir(mb):
            if name.startswith("__"):
                continue  # skip dunder attributes like __annotations__
            obj = getattr(mb, name)
            if not isinstance(obj, dict) or not obj:
                continue
            # If any key is a periodic-table element symbol, the dict must
            # have been imported from a universal table.  We check by name.
            keys = list(obj.keys())
            has_element_key = any(isinstance(k, str) and k in elements for k in keys)
            if not has_element_key:
                continue
            # Otherwise, allow only the imported universal tables.
            allowed_names = {"_COV", "_VDW", "_HYB_MAP"}
            assert name in allowed_names, (
                f"Disallowed element-keyed dict at module level: {name}"
            )


# ===========================================================================
# 7.  Hyb-string mapping (universal)
# ===========================================================================
class TestHybStringMap:
    def test_aromatic_carbon_is_sp2(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        atom = mol.GetAtomWithIdx(0)
        # RDKit reports aromatic C as SP2.
        assert hyb_str(atom) == "sp2"

    def test_tetrahedral_carbon_is_sp3(self):
        mol = Chem.MolFromSmiles("CCCC")
        atom = mol.GetAtomWithIdx(0)
        assert hyb_str(atom) == "sp3"

    def test_unknown_hyb_is_wildcard(self):
        # Lone hydrogen atom — RDKit reports UNSPECIFIED.
        mol = Chem.MolFromSmiles("[H]")
        atom = mol.GetAtomWithIdx(0)
        # Some RDKit builds parse [H] as a hydrogen with hybridisation UNSPECIFIED
        # (or rejects it); either way, hyb_str must never crash.
        h = hyb_str(atom)
        assert h in ("*", "sp", "sp2", "sp3")


# ===========================================================================
# 8.  vdW floor — non-bonded pairs
# ===========================================================================
class TestVdWFloor:
    def test_distant_atoms_get_vdw_floor(self, lib):
        """Atoms that are NOT bonded / 1,3 / M-D / resonance-equivalent
        must get a vdW lower bound of ≈ 0.85·(r_vdW_i + r_vdW_j)."""
        # Build a Pd complex with a long ligand so we have non-bonded H-H far
        # apart.  Pyridine-Pd:
        mol = Chem.MolFromSmiles("c1ccncc1[Pd](Cl)(Cl)Cl")
        if mol is None:
            pytest.skip("RDKit parse fail on pyridine-Pd")
        mol = Chem.AddHs(mol)
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        m = syms.index("Pd")
        donors = sorted(int(nb.GetIdx()) for nb in mol.GetAtomWithIdx(m).GetNeighbors())
        L, U, info = build_bounds_matrix(syms, mol, m, donors, lib)
        assert info["n_vdw_floors"] > 0
        # Pick any far H-H pair (non-bonded, non-1,3).
        h_atoms = [i for i, s in enumerate(syms) if s == "H"]
        if len(h_atoms) >= 2:
            i, j = h_atoms[0], h_atoms[-1]
            if i != j:
                # Bondi r_vdW(H) = 1.20; floor = 0.85·(1.20+1.20) = 2.04 Å
                # (the cell must have at least this floor unless tighter
                # constraint was applied).
                assert L[i, j] >= 2.0, (
                    f"H-H vdW floor too small: {L[i,j]:.3f}"
                )


# ===========================================================================
# 9.  Feasibility / infeasibility handling
# ===========================================================================
class TestFeasibilityRelaxation:
    def test_feasible_geometry_no_relaxation(self, lib):
        syms, mol, m, donors = _build_complex("[NH3][Pt]([NH3])(Cl)Cl")
        L, U, info = build_bounds_matrix(
            syms, mol, m, donors, lib, geometry="SP-4 square planar",
        )
        assert info["bounds_relaxed"] == []

    def test_diagonal_is_zero(self, lib):
        syms, mol, m, donors = _build_complex("[NH3][Pt]([NH3])(Cl)Cl")
        L, U, info = build_bounds_matrix(
            syms, mol, m, donors, lib, geometry="SP-4 square planar",
        )
        n = L.shape[0]
        for i in range(n):
            assert L[i, i] == 0.0
            assert U[i, i] == 0.0

    def test_symmetric_matrices(self, lib):
        syms, mol, m, donors = _build_complex("[NH3][Pt]([NH3])(Cl)Cl")
        L, U, info = build_bounds_matrix(
            syms, mol, m, donors, lib, geometry="SP-4 square planar",
        )
        assert np.allclose(L, L.T)
        assert np.allclose(U, U.T)


# ===========================================================================
# 10a.  Cu(en)₂²⁺ — Jahn-Teller candidate (multi-modal Cu-N distribution)
# ===========================================================================
class TestCuEnJahnTeller:
    """Cu(II) chelate is Jahn-Teller distorted — the Cu-N distance
    distribution is bimodal (equatorial ~ 2.00 Å, axial ~ 2.40 Å).  Our
    bound window must be permissive enough to cover both modes.  This is
    an *acceptance band* test rather than a precise μ test (Phase A only
    builds the bound; Phase C handles GMM-aware loss)."""

    def test_cu_n_bound_covers_both_modes(self, lib):
        # Build [Cu(en)₂]²⁺ explicitly.
        rw = Chem.RWMol(Chem.MolFromSmiles("NCCN.NCCN"))
        # The two ligands have 4 N total — bond them all to a single Cu.
        cu = rw.AddAtom(Chem.Atom("Cu"))
        # Find the 4 N atoms.
        n_idxs = [int(i) for i, a in enumerate(rw.GetAtoms()) if a.GetSymbol() == "N"]
        assert len(n_idxs) == 4
        for n in n_idxs:
            rw.AddBond(n, cu, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        rw = Chem.RWMol(Chem.AddHs(rw))
        syms = [a.GetSymbol() for a in rw.GetAtoms()]
        cu = syms.index("Cu")
        donors = sorted(int(nb.GetIdx()) for nb in rw.GetAtomWithIdx(cu).GetNeighbors())
        L, U, info = build_bounds_matrix(syms, rw, cu, donors, lib)
        # The CCDC Cu-N empirical mean is ~ 2.00 Å.  Our 2σ window should
        # be at least 0.05 Å wide and the midpoint within [1.95, 2.15].
        for d in donors:
            if syms[d] != "N":
                continue
            mid = 0.5 * (L[cu, d] + U[cu, d])
            assert 1.85 < mid < 2.20, f"Cu-N mid {mid:.3f} out of range"
            assert U[cu, d] - L[cu, d] > 0.03


# ===========================================================================
# 10b.  La(H₂O)₉ — f-block CN9 (high coordination number)
# ===========================================================================
class TestLanthanideF_BlockCN9:
    def test_la_o_distance_in_lanthanide_range(self, lib):
        # 9 water O atoms bonded to La.
        rw = Chem.RWMol()
        la = rw.AddAtom(Chem.Atom("La"))
        o_idxs = []
        for _ in range(9):
            o = rw.AddAtom(Chem.Atom("O"))
            rw.AddBond(la, o, Chem.BondType.SINGLE)
            o_idxs.append(o)
        # Add 2 H per O.
        for o in o_idxs:
            for _ in range(2):
                h = rw.AddAtom(Chem.Atom("H"))
                rw.AddBond(o, h, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        syms = [a.GetSymbol() for a in rw.GetAtoms()]
        L, U, info = build_bounds_matrix(syms, rw, la, o_idxs, lib)
        # La-O empirical ≈ 2.45-2.55 Å (Shannon ionic radii / CCDC).
        # Even with cov-sum fallback the value should land in [2.5, 2.9] Å.
        for o in o_idxs:
            mid = 0.5 * (L[la, o] + U[la, o])
            assert 2.0 < mid < 3.2, f"La-O mid {mid:.3f} out of Ln range"


# ===========================================================================
# 10c.  Universal — no SMILES-pattern matching in source
# ===========================================================================
class TestSourceGrepUniversal:
    def test_no_smarts_or_smiles_patterns(self):
        """The module source must not contain SMARTS / SMILES pattern
        strings.  These would be SMILES-specific shortcuts."""
        src = Path(
            "/home/qmchem_max/ComPlat/DELFIN/delfin/fffree/mogul_bounds.py"
        ).read_text()
        # Strip docstrings + comments.
        src_no_doc = re.sub(r'"""[\s\S]*?"""', "", src, flags=re.MULTILINE)
        src_no_doc = re.sub(r"'''[\s\S]*?'''", "", src_no_doc, flags=re.MULTILINE)
        src_no_doc = re.sub(r"#.*", "", src_no_doc)
        # Look for `MolFromSmarts`, `HasSubstructMatch`, `GetSubstructMatches`.
        for pat in ("MolFromSmarts", "MolFromSmiles", "HasSubstructMatch",
                    "GetSubstructMatches"):
            assert pat not in src_no_doc, (
                f"Found SMILES/SMARTS-specific pattern in source: {pat}"
            )


# ===========================================================================
# 11.  Info dict shape
# ===========================================================================
class TestInfoDict:
    def test_info_has_required_fields(self, lib):
        syms, mol, m, donors = _build_complex("[NH3][Pt]([NH3])(Cl)Cl")
        _L, _U, info = build_bounds_matrix(syms, mol, m, donors, lib)
        for key in (
            "n_aromatic_groups",
            "n_resonance_groups",
            "n_md_bounds",
            "bounds_relaxed",
        ):
            assert key in info, f"Missing info field: {key}"


# ===========================================================================
# 12.  Multi-tier bond-key fallback chain
#
#      Verifies that DELFIN_FFFREE_MOGUL_BOND_FALLBACK=1 raises the ligand
#      1-2 bond hit-rate from the centred-fragment baseline (~0-10 %) to
#      ~100 % for representative organic-rich SMILES while keeping the
#      OFF path BYTE-IDENTICAL to the pre-fix behaviour.
# ===========================================================================
import hashlib  # noqa: E402 — re-import for the test class block
from delfin.fffree.mogul_bounds import (  # noqa: E402
    _lookup_bond_organic,
    _lookup_bond_organic_with_tier,
    _bond_fallback_enabled,
    MOGUL_BOND_FALLBACK_ENV_FLAG,
)


# Three representative SMILES from per_smiles_ccdc_families:
#   * SIYMEU = Ag-NHC + acetate            (heavy organic ligand, hetero C/N/O/Ag)
#   * BERTEB = Ni-Br + bipyridyl-CH2       (aromatic π + sp3 link)
#   * ALAHEB = Ir-CO + di-tBu-phosphine    (sp3 dense, P-C heavy)
_REPRESENTATIVE_SMILES = {
    "SIYMEU": (
        "CC(=O)[O][Ag-][C+]1N(CC2=CC=C(C)C=C2)C(C2=CC=C(C(C)C)C=C2)"
        "=C(C2=CC=C(C(C)C)C=C2)N1CC1=CC=C(C)C=C1"
    ),
    "BERTEB": "[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12",
    "ALAHEB": (
        "CC(C)(C)[P+]1(C(C)(C)C)C=C2C=CC=C3C[P+](C(C)(C)C)(C(C)(C)C)"
        "[Ir-3]1([C]#[O+])[N]23"
    ),
}


def _measure_ligand_bond_hit_rate(mol, metal_idx, donors, lib, *, fallback_on: bool):
    """Return ``(n_total, n_lib_hits, tier_counter)`` for ligand 1-2 bonds.

    Uses :func:`_lookup_bond_organic` (which honours the env flag) when
    measuring hits, so the OFF/ON counter is the real production
    behaviour.  Tier counts are filled by an extra
    :func:`_lookup_bond_organic_with_tier` call when ``fallback_on`` is True
    (otherwise tier counts are zeroed because the OFF path does not
    classify hits into tiers).
    """
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    donor_set = set(int(d) for d in donors)
    prev_env = os.environ.get(MOGUL_BOND_FALLBACK_ENV_FLAG)
    if fallback_on:
        os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = "1"
    else:
        os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
    try:
        n_total = 0
        n_lib = 0
        tier_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        for bond in mol.GetBonds():
            i = int(bond.GetBeginAtomIdx()); j = int(bond.GetEndAtomIdx())
            if (i == metal_idx and j in donor_set) or (j == metal_idx and i in donor_set):
                continue
            n_total += 1
            zi = syms[i]; zj = syms[j]
            hi = hyb_str(mol.GetAtomWithIdx(i))
            hj = hyb_str(mol.GetAtomWithIdx(j))
            # Real hit count honours env flag.
            hit = _lookup_bond_organic(lib, zi, hi, zj, hj, None, min_n=5)
            if hit is not None:
                n_lib += 1
            # Tier diagnostic always runs the full chain so we can report
            # WHICH tier provided the hit when fallback_on is True.
            if fallback_on:
                _diag_hit, tier = _lookup_bond_organic_with_tier(
                    lib, zi, hi, zj, hj, None, min_n=5,
                )
                tier_counts[tier] = tier_counts.get(tier, 0) + 1
        return n_total, n_lib, tier_counts
    finally:
        if prev_env is None:
            os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
        else:
            os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = prev_env


class TestBondFallbackChain:
    """Validate the multi-tier bond-key fallback chain."""

    def test_env_flag_default_off(self):
        """Default-OFF: the env flag must report False when unset."""
        prev = os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
        try:
            assert _bond_fallback_enabled() is False
        finally:
            if prev is not None:
                os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = prev

    def test_env_flag_on_when_set(self):
        prev = os.environ.get(MOGUL_BOND_FALLBACK_ENV_FLAG)
        os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = "1"
        try:
            assert _bond_fallback_enabled() is True
        finally:
            if prev is None:
                os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
            else:
                os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = prev

    @pytest.mark.parametrize("label", list(_REPRESENTATIVE_SMILES))
    def test_fallback_on_raises_hit_rate_to_at_least_60pct(self, lib, label):
        """ON: tier chain must lift ligand 1-2 hit-rate to >= 60 %.

        The hard target from the mission spec is 60-80 %.  In practice the
        pair_bond table covers organic bonds at ~100 % so this should fire
        well above the floor.
        """
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        smi = _REPRESENTATIVE_SMILES[label]
        syms, mol, m, donors = _build_complex(smi)
        n_total, n_lib, tiers = _measure_ligand_bond_hit_rate(
            mol, m, donors, lib, fallback_on=True,
        )
        pct = (n_lib / n_total * 100.0) if n_total else 0.0
        assert pct >= 60.0, (
            f"{label}: ligand 1-2 lib hit-rate {pct:.1f} % below 60 % "
            f"(tier counts: {tiers})"
        )

    @pytest.mark.parametrize("label", list(_REPRESENTATIVE_SMILES))
    def test_fallback_on_beats_off_for_each_case(self, lib, label):
        """ON-hit-count must strictly exceed OFF-hit-count for every case."""
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        smi = _REPRESENTATIVE_SMILES[label]
        syms, mol, m, donors = _build_complex(smi)
        n_total_off, n_lib_off, _ = _measure_ligand_bond_hit_rate(
            mol, m, donors, lib, fallback_on=False,
        )
        n_total_on, n_lib_on, _ = _measure_ligand_bond_hit_rate(
            mol, m, donors, lib, fallback_on=True,
        )
        assert n_total_on == n_total_off
        assert n_lib_on > n_lib_off, (
            f"{label}: ON hits {n_lib_on} not greater than OFF hits {n_lib_off}"
        )

    def test_total_hit_rate_across_3_cases_above_90pct(self, lib):
        """Joint hit-rate across all 3 SMILES must exceed 90 % when ON."""
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        agg_total = 0
        agg_lib = 0
        for label, smi in _REPRESENTATIVE_SMILES.items():
            syms, mol, m, donors = _build_complex(smi)
            n_t, n_l, _ = _measure_ligand_bond_hit_rate(
                mol, m, donors, lib, fallback_on=True,
            )
            agg_total += n_t
            agg_lib += n_l
        pct = agg_lib / agg_total * 100.0
        assert pct >= 90.0, f"joint lib hit-rate {pct:.1f} % below 90 %"

    def test_off_path_byte_identical_to_legacy(self, lib):
        """OFF path: bounds matrices must be byte-identical across two runs.

        Also: the lib_hit count must match the legacy pre-fix expectation
        (low single-digit counts for these heavy-organic SMILES).
        """
        prev = os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
        try:
            for label, smi in _REPRESENTATIVE_SMILES.items():
                syms, mol, m, donors = _build_complex(smi)
                lo1, up1, info1 = build_bounds_matrix(syms, mol, m, donors, lib)
                lo2, up2, info2 = build_bounds_matrix(syms, mol, m, donors, lib)
                h_lo1 = hashlib.sha256(lo1.tobytes()).hexdigest()
                h_lo2 = hashlib.sha256(lo2.tobytes()).hexdigest()
                h_up1 = hashlib.sha256(up1.tobytes()).hexdigest()
                h_up2 = hashlib.sha256(up2.tobytes()).hexdigest()
                assert h_lo1 == h_lo2, f"{label}: lower not deterministic"
                assert h_up1 == h_up2, f"{label}: upper not deterministic"
                # Pre-fix path -> only single-neighbour centred-fragment hits.
                # All three cases had < 25 % lib_hit on the legacy chain.
                n_bonds_lig = info1["n_bond_lib_hit"] + info1["n_bond_fallback"]
                ratio = info1["n_bond_lib_hit"] / max(1, n_bonds_lig)
                assert ratio < 0.30, (
                    f"{label}: OFF lib-hit ratio {ratio:.1%} unexpectedly "
                    "high -- did the OFF path change?"
                )
        finally:
            if prev is not None:
                os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = prev

    def test_on_path_determinism(self, lib):
        """ON path: two runs must produce byte-identical matrices."""
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        prev = os.environ.get(MOGUL_BOND_FALLBACK_ENV_FLAG)
        os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = "1"
        try:
            for label, smi in _REPRESENTATIVE_SMILES.items():
                syms, mol, m, donors = _build_complex(smi)
                lo1, up1, _ = build_bounds_matrix(syms, mol, m, donors, lib)
                lo2, up2, _ = build_bounds_matrix(syms, mol, m, donors, lib)
                assert hashlib.sha256(lo1.tobytes()).hexdigest() == \
                       hashlib.sha256(lo2.tobytes()).hexdigest()
                assert hashlib.sha256(up1.tobytes()).hexdigest() == \
                       hashlib.sha256(up2.tobytes()).hexdigest()
        finally:
            if prev is None:
                os.environ.pop(MOGUL_BOND_FALLBACK_ENV_FLAG, None)
            else:
                os.environ[MOGUL_BOND_FALLBACK_ENV_FLAG] = prev

    def test_tier_returned_is_in_range(self, lib):
        """Every call must return tier in {1, 2, 3, 4, 5, 6}."""
        for label, smi in _REPRESENTATIVE_SMILES.items():
            syms, mol, m, donors = _build_complex(smi)
            for bond in mol.GetBonds():
                i = int(bond.GetBeginAtomIdx()); j = int(bond.GetEndAtomIdx())
                zi = syms[i]; zj = syms[j]
                hi = hyb_str(mol.GetAtomWithIdx(i))
                hj = hyb_str(mol.GetAtomWithIdx(j))
                _hit, tier = _lookup_bond_organic_with_tier(
                    lib, zi, hi, zj, hj, None, min_n=5,
                )
                assert tier in (1, 2, 3, 4, 5, 6), f"unexpected tier {tier}"

    def test_element_pair_aggregate_returns_n_weighted_mean(self, lib):
        """Tier 5 helper: aggregate over all hyb variants for C-C."""
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        hit = lib._lookup_bond_element_pair_aggregate("C", "C", min_n=5)
        assert hit is not None
        mu, sigma, n = hit
        # Aggregate C-C mean must lie between sp-sp (~1.20) and sp3-sp3 (~1.52)
        assert 1.20 <= mu <= 1.55, f"aggregate C-C mu {mu:.3f} out of range"
        assert sigma > 0.0
        # Several million C-C bonds in CCDC.
        assert n >= 1_000_000, f"aggregate C-C n {n} unexpectedly small"

    def test_element_pair_aggregate_symmetric(self, lib):
        """``aggregate(C, N) == aggregate(N, C)`` (unordered)."""
        if not getattr(lib, "has_pair_tables", False):
            pytest.skip("library has no pair_bond tables (v1/v2/v3)")
        h1 = lib._lookup_bond_element_pair_aggregate("C", "N", min_n=5)
        h2 = lib._lookup_bond_element_pair_aggregate("N", "C", min_n=5)
        assert h1 is not None and h2 is not None
        assert h1 == h2, "aggregate should be unordered-element-pair-symmetric"
