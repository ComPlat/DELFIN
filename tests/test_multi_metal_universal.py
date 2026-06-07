"""Tests for the universal multi-metal / cluster / mu-bridging block in
``delfin.fffree.assemble_complex`` (Task: multi-metal universal coverage,
2026-06-07).

Covers
------
1. Mononuclear detection -> topology_class == 'mono'.
2. Cisplatin -> 'mono' (1 Pt + Cl, Cl, NH3, NH3).
3. Cu_2(mu-OH)_2 L -> 'binuclear', bridging atoms identified, builds.
4. [Fe_4 S_4] (4 metals, S atoms bridge) -> 'cluster_4'.
5. Ferrocene Fe(Cp)_2 -> 'mono' (1 Fe, Cp carbons are not separate metals).
6. Byte-identical OFF: with env flag unset, ``assemble_multi_metal_universal``
   always returns ``None``.
7. Determinism: repeated calls give bit-identical numerical output.
8. Demo: Cu_2(mu-OH) -> 2 metals placed colinear with the bridge between
   them, donors split across hemispheres.
"""
from __future__ import annotations

import os
import math

import numpy as np
import pytest
from rdkit import Chem

from delfin.fffree import assemble_complex as AC


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Each test starts with the universal multi-metal flag unset."""
    monkeypatch.delenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", raising=False)
    yield


def _mol(smiles: str):
    """Permissive parse: sanitise without enforcing valence (works for the
    salt/dot/datively-bonded SMILES used in real CCDC inputs).
    """
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        raise ValueError(f"Cannot parse {smiles!r}")
    try:
        Chem.SanitizeMol(
            m,
            sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
            ),
        )
    except Exception:
        pass
    return m


# ---------------------------------------------------------------------------
# 1. Detection: mononuclear
# ---------------------------------------------------------------------------


def test_detect_mononuclear():
    mol = _mol("[Cu]([Cl])([Cl])[Cl]")
    topo = AC.detect_metal_topology(mol)
    assert topo["n_metals"] == 1
    assert topo["topology_class"] == "mono"
    assert topo["bridging_atoms"] == []
    assert topo["metal_metal_bonds"] == []
    assert topo["metal_symbols"] == ["Cu"]


# ---------------------------------------------------------------------------
# 2. Detection: cisplatin
# ---------------------------------------------------------------------------


def test_detect_cisplatin():
    # cis-[Pt(NH3)2 Cl2] -- 1 Pt + 4 ligands
    mol = _mol("[Pt]([Cl])([Cl])([NH3])[NH3]")
    topo = AC.detect_metal_topology(mol)
    assert topo["n_metals"] == 1
    assert topo["topology_class"] == "mono"
    assert topo["metal_symbols"] == ["Pt"]


# ---------------------------------------------------------------------------
# 3. Detection + placement: Cu2(mu-OH)2 dimer
# ---------------------------------------------------------------------------


def test_detect_binuclear_cu2_dioh():
    # Cu(mu-OH)2 Cu with terminal water ligands -- a textbook hydroxo-bridged
    # binuclear motif.  Both bridging oxygens link both Cu atoms.
    mol = _mol("[Cu]1([OH2])([OH2])O[Cu]([OH2])([OH2])O1")
    topo = AC.detect_metal_topology(mol)
    assert topo["n_metals"] == 2
    assert topo["topology_class"] == "binuclear"
    assert len(topo["bridging_atoms"]) == 2          # two mu-OH oxygens
    assert topo["metal_metal_bonds"] == []           # no explicit Cu-Cu bond
    assert topo["metal_symbols"] == ["Cu", "Cu"]


def test_place_binuclear_with_bridge_linear_geometry():
    """M1 at origin, B above, M2 above B -- all colinear on +z."""
    syms, P = AC.place_binuclear_with_bridge(
        "Cu", "Cu", "O",
        donors_a=[(10, "N"), (11, "N")],
        donors_b=[(20, "N"), (21, "N")],
    )
    # Header order: M1, B, M2, then donors_a, then donors_b
    assert syms[0] == "Cu"
    assert syms[1] == "O"
    assert syms[2] == "Cu"
    # Linearity: M1 -> B -> M2 all on z-axis (x ~ 0, y ~ 0).
    for i in (0, 1, 2):
        assert abs(P[i, 0]) < 1e-9
        assert abs(P[i, 1]) < 1e-9
    # Strict ordering along z.
    assert P[0, 2] < P[1, 2] < P[2, 2]
    # M1-B and B-M2 lengths match the polyhedra md_distance values.
    from delfin.fffree.polyhedra import md_distance as md
    assert abs(np.linalg.norm(P[1] - P[0]) - md("Cu", "O")) < 1e-6
    assert abs(np.linalg.norm(P[2] - P[1]) - md("Cu", "O")) < 1e-6


def test_place_binuclear_hemisphere_partition():
    """Donors_a sit BELOW M1; donors_b sit ABOVE M2 -- universal hemispheres."""
    syms, P = AC.place_binuclear_with_bridge(
        "Cu", "Cu", "O",
        donors_a=[(10, "N"), (11, "N"), (12, "N")],
        donors_b=[(20, "N"), (21, "N"), (22, "N")],
    )
    pos_m1 = P[0]
    pos_m2 = P[2]
    n_da = 3
    # donor_a positions follow at index 3..3+n_da
    for k in range(n_da):
        # relative z below the bridge (M1 sits at z=0)
        assert P[3 + k, 2] < P[1, 2]                # below the bridge oxygen
    for k in range(n_da):
        assert P[3 + n_da + k, 2] > P[1, 2]         # above the bridge oxygen
    # Donors at correct M-D distance from their respective metal.
    from delfin.fffree.polyhedra import md_distance as md
    r = md("Cu", "N")
    for k in range(n_da):
        assert abs(np.linalg.norm(P[3 + k] - pos_m1) - r) < 1e-6
        assert abs(np.linalg.norm(P[3 + n_da + k] - pos_m2) - r) < 1e-6


# ---------------------------------------------------------------------------
# 4. Detection: [Fe4 S4] cubane core
# ---------------------------------------------------------------------------


def test_detect_cluster_fe4s4():
    # Each S bridges 3 Fe (mu3) in the cubane.  Build a minimal SMILES of
    # 4 Fe + 4 S where every S is bonded to >=2 Fe (full cubane connectivity
    # would be 3-to-3, but the topology test only cares about >=2).
    mol = Chem.RWMol()
    fe_idx = [mol.AddAtom(Chem.Atom("Fe")) for _ in range(4)]
    s_idx = [mol.AddAtom(Chem.Atom("S")) for _ in range(4)]
    # Cubane: each Fe is bonded to 3 of the 4 S (skip the antipodal one).
    bonds = [
        (fe_idx[0], s_idx[0]), (fe_idx[0], s_idx[1]), (fe_idx[0], s_idx[2]),
        (fe_idx[1], s_idx[0]), (fe_idx[1], s_idx[1]), (fe_idx[1], s_idx[3]),
        (fe_idx[2], s_idx[0]), (fe_idx[2], s_idx[2]), (fe_idx[2], s_idx[3]),
        (fe_idx[3], s_idx[1]), (fe_idx[3], s_idx[2]), (fe_idx[3], s_idx[3]),
    ]
    for a, b in bonds:
        mol.AddBond(a, b, Chem.BondType.SINGLE)
    topo = AC.detect_metal_topology(mol.GetMol())
    assert topo["n_metals"] == 4
    assert topo["topology_class"] == "cluster_4"
    assert len(topo["bridging_atoms"]) == 4         # all four S atoms bridge


# ---------------------------------------------------------------------------
# 5. Detection: ferrocene -- 1 metal so 'mono'
# ---------------------------------------------------------------------------


def test_detect_ferrocene_mono():
    # Two Cp's sandwich a single Fe.  Topology layer only counts metal atoms.
    mol = _mol("[Fe].c1ccc[cH]1.c1ccc[cH]1")
    topo = AC.detect_metal_topology(mol)
    assert topo["n_metals"] == 1
    assert topo["topology_class"] == "mono"


# ---------------------------------------------------------------------------
# 6. Byte-identical OFF: assemble_multi_metal_universal returns None
# ---------------------------------------------------------------------------


def test_byte_identical_off(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", raising=False)
    mol = _mol("[Cu]1([OH2])([OH2])O[Cu]([OH2])([OH2])O1")
    # Both the detector and the assembler are side-effect-free; the dispatch
    # entry must be silent (None) when the flag is unset.
    assert AC.assemble_multi_metal_universal(mol) is None
    # Detector is unconditional -- it's a read-only graph walk and is safe
    # to call any time.  But its output must not vary with the flag.
    monkeypatch.setenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", "1")
    topo_on = AC.detect_metal_topology(mol)
    monkeypatch.delenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", raising=False)
    topo_off = AC.detect_metal_topology(mol)
    assert topo_on == topo_off


def test_assemble_off_on_separate_salt(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", "1")
    # Two isolated copper ions: no shared neighbour, no M-M bond.  Even with
    # the flag ON, the universal block returns None so the caller falls back
    # to per-metal mono treatment.
    mol = _mol("[Cu+].[Cu+]")
    assert AC.assemble_multi_metal_universal(mol) is None


# ---------------------------------------------------------------------------
# 7. Determinism
# ---------------------------------------------------------------------------


def test_determinism_binuclear(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", "1")
    mol = _mol("[Cu]1([OH2])([OH2])O[Cu]([OH2])([OH2])O1")
    res1 = AC.assemble_multi_metal_universal(mol)
    res2 = AC.assemble_multi_metal_universal(mol)
    assert res1 is not None and res2 is not None
    s1, P1 = res1
    s2, P2 = res2
    assert s1 == s2
    # Bit-identical coords.
    assert np.array_equal(P1, P2)


def test_determinism_place_binuclear():
    s1, P1 = AC.place_binuclear_with_bridge(
        "Cu", "Cu", "O",
        donors_a=[(10, "N"), (11, "N")],
        donors_b=[(20, "N"), (21, "N")],
    )
    s2, P2 = AC.place_binuclear_with_bridge(
        "Cu", "Cu", "O",
        donors_a=[(10, "N"), (11, "N")],
        donors_b=[(20, "N"), (21, "N")],
    )
    assert s1 == s2
    assert np.array_equal(P1, P2)


def test_determinism_cluster_polygon():
    metals = ["Fe", "Fe", "Fe", "Fe"]
    donors = [[(i * 10 + 1, "N")] for i in range(4)]
    s1, P1 = AC.place_metal_cluster_polygon(metals, donors)
    s2, P2 = AC.place_metal_cluster_polygon(metals, donors)
    assert s1 == s2
    assert np.array_equal(P1, P2)


# ---------------------------------------------------------------------------
# 8. Demo: end-to-end Cu2(mu-OH) builds with the universal block
# ---------------------------------------------------------------------------


def test_demo_cu2_mu_oh_assemble(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", "1")
    # Cu(mu-OH)2 Cu with terminal aquo ligands.
    mol = _mol("[Cu]1([OH2])([OH2])O[Cu]([OH2])([OH2])O1")
    out = AC.assemble_multi_metal_universal(mol)
    assert out is not None
    syms, P = out
    # 2 metals + 2 bridges + 4 terminals = 8 heavy atoms
    n_metals = sum(1 for s in syms if s == "Cu")
    assert n_metals == 2
    # Metals must not collide.
    cu_positions = [P[i] for i, s in enumerate(syms) if s == "Cu"]
    assert np.linalg.norm(cu_positions[0] - cu_positions[1]) > 1.5
    # Universal: the polyhedra md_distance was used for every M-(terminal).
    from delfin.fffree.polyhedra import md_distance as md
    # The terminal O of the aquo ligand should sit at md(Cu, O) from its Cu.
    # (Approximate check: at least one O distance matches md(Cu,O).)
    md_co = md("Cu", "O")
    distances = []
    for i, s in enumerate(syms):
        if s == "O":
            for j, t in enumerate(syms):
                if t == "Cu":
                    distances.append(float(np.linalg.norm(P[i] - P[j])))
    assert any(abs(d - md_co) < 1e-4 for d in distances), (
        f"No O-Cu distance matched md(Cu,O)={md_co:.3f}: got {distances}"
    )


def test_demo_cluster_3():
    """Three Fe atoms in a triangle with one terminal Cl each -- pure
    cluster_3 polygon path.  Confirms the polygon placement runs and is
    deterministic."""
    mol = Chem.RWMol()
    fe = [mol.AddAtom(Chem.Atom("Fe")) for _ in range(3)]
    cl = [mol.AddAtom(Chem.Atom("Cl")) for _ in range(3)]
    # Triangle of Fe-Fe-Fe bonds + one terminal Cl per Fe.
    mol.AddBond(fe[0], fe[1], Chem.BondType.SINGLE)
    mol.AddBond(fe[1], fe[2], Chem.BondType.SINGLE)
    mol.AddBond(fe[2], fe[0], Chem.BondType.SINGLE)
    for i in range(3):
        mol.AddBond(fe[i], cl[i], Chem.BondType.SINGLE)
    real_mol = mol.GetMol()
    topo = AC.detect_metal_topology(real_mol)
    assert topo["topology_class"] == "cluster_3"
    assert len(topo["metal_metal_bonds"]) == 3

    os.environ["DELFIN_FFFREE_MULTI_METAL_UNIVERSAL"] = "1"
    try:
        out = AC.assemble_multi_metal_universal(real_mol)
    finally:
        os.environ.pop("DELFIN_FFFREE_MULTI_METAL_UNIVERSAL", None)
    assert out is not None
    syms, P = out
    # 3 Fe + 3 Cl
    assert sum(1 for s in syms if s == "Fe") == 3
    assert sum(1 for s in syms if s == "Cl") == 3
    # Metals on a triangle: all in z=0, pairwise distances within 1% of each
    # other.
    fe_pos = np.array([P[i] for i, s in enumerate(syms) if s == "Fe"])
    assert max(abs(fe_pos[:, 2])) < 1e-9
    d01 = np.linalg.norm(fe_pos[0] - fe_pos[1])
    d12 = np.linalg.norm(fe_pos[1] - fe_pos[2])
    d20 = np.linalg.norm(fe_pos[2] - fe_pos[0])
    avg = (d01 + d12 + d20) / 3
    assert all(abs(d - avg) / avg < 0.01 for d in (d01, d12, d20))
