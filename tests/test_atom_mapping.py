"""Tests for delfin.atom_mapping (universal atom-order mapper).

Covers:
  * pure renumbering round-trip (symmetric organic),
  * a metal complex (metal bonds excluded from the graph),
  * universality on a metal outside DELFIN's tuned radii table (lanthanide),
  * a reaction pair differing by exactly one formed bond,
  * formula-mismatch error handling.
"""

import numpy as np
import pytest

pytest.importorskip("networkx")
pytest.importorskip("scipy")
pytest.importorskip("rdkit")

from delfin import atom_mapping as am


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _roundtrip(syms, coords, seed=0):
    """Permute the target, map it back, and assert exact recovery."""
    coords = np.asarray(coords, dtype=float)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(syms))
    tsyms = [syms[i] for i in perm]
    tcoords = coords[perm]

    res = am.map_atoms(syms, coords, tsyms, tcoords)

    assert res["verified"] is True
    assert res["n_bond_edits"] == 0
    order = res["order"]
    # reordered target must reproduce the reference exactly
    assert [tsyms[j] for j in order] == list(syms)
    assert np.allclose(tcoords[order], coords, atol=1e-6)


def _mol_from_smiles(smiles, seed=0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    return syms, coords


# minimal metal complex: M + 2 equivalent Cl (isolated once metal bonds are
# excluded) + an NH2 group (2 equivalent H) + an OH group.  Exercises metal
# handling and degenerate-class (Cl/Cl, H/H) resolution.
def _metal_complex(metal="Fe"):
    syms = [metal, "Cl", "Cl", "N", "H", "H", "O", "H"]
    coords = np.array(
        [
            [0.00, 0.00, 0.00],   # M
            [2.20, 0.00, 0.00],   # Cl
            [-2.20, 0.00, 0.00],  # Cl
            [0.00, 2.00, 0.00],   # N
            [0.60, 2.60, 0.40],   # N-H
            [-0.60, 2.60, 0.40],  # N-H
            [0.00, -2.00, 0.00],  # O
            [0.20, -2.90, 0.30],  # O-H
        ],
        dtype=float,
    )
    return syms, coords


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
def test_element_classification_and_radii():
    assert am.is_metal("Pd") and am.is_metal("Fe") and am.is_metal("Dy")
    assert am.is_metal("U") and am.is_metal("Na")
    assert not am.is_metal("C") and not am.is_metal("N")
    assert not am.is_metal("B") and not am.is_metal("Si")  # covalent metalloids
    # every radius resolves to a finite positive number (universal)
    for el in ("H", "C", "Pd", "Dy", "U", "Xx"):
        assert am.cov_radius(el) > 0.0


def test_roundtrip_symmetric_organic():
    # toluene: symmetric methyl H (class of 3) + ring symmetry
    syms, coords = _mol_from_smiles("Cc1ccccc1")
    _roundtrip(syms, coords, seed=1)
    _roundtrip(syms, coords, seed=2)


def test_metal_complex_roundtrip():
    syms, coords = _metal_complex("Fe")
    _roundtrip(syms, coords, seed=3)


def test_universality_lanthanide_outside_tuned_table():
    # Dy is not in delfin.manta.polyhedra.COV -> must fall through to the RDKit
    # periodic table; category-based is_metal must still treat it as a metal.
    from delfin.manta.polyhedra import COV as tuned
    assert "Dy" not in tuned
    syms, coords = _metal_complex("Dy")
    _roundtrip(syms, coords, seed=4)


def test_reaction_one_formed_bond():
    # ethanol product; educt = same atoms with the C-C bond broken by pulling
    # the methyl group away.  Mapping must detect exactly one bond edit.
    syms, coords = _mol_from_smiles("CCO")
    prod_syms, prod_xyz = list(syms), coords.copy()

    G = am.connectivity(syms, coords)
    heavy_deg = {
        i: sum(1 for j in G[i] if syms[j] != "H")
        for i in range(len(syms))
        if syms[i] == "C"
    }
    methyl_c = min(heavy_deg, key=heavy_deg.get)      # 1 heavy neighbour
    other_c = max(heavy_deg, key=heavy_deg.get)       # 2 heavy neighbours
    assert heavy_deg[methyl_c] == 1 and heavy_deg[other_c] == 2

    group = [methyl_c] + [j for j in G[methyl_c] if syms[j] == "H"]
    axis = coords[methyl_c] - coords[other_c]
    axis /= np.linalg.norm(axis)

    educt_xyz = coords.copy()
    educt_xyz[group] += 1.8 * axis  # break the C-C bond, keep it a candidate pair

    # sanity: the educt really has no C-C bond, the product does
    assert not am.connectivity(syms, educt_xyz).has_edge(methyl_c, other_c)
    assert am.connectivity(prod_syms, prod_xyz).has_edge(methyl_c, other_c)

    res = am.map_atoms(syms, educt_xyz, prod_syms, prod_xyz)
    assert res["verified"] is True
    assert res["n_bond_edits"] == 1


def test_formula_mismatch_raises():
    a_syms, a_xyz = _metal_complex("Fe")
    b_syms, b_xyz = _metal_complex("Fe")
    b_syms = list(b_syms)
    b_syms[1] = "Br"  # Cl -> Br: composition no longer matches
    with pytest.raises(ValueError):
        am.map_atoms(a_syms, a_xyz, b_syms, b_xyz)
