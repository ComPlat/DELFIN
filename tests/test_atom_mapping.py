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


@pytest.mark.parametrize(
    "smiles",
    ["c1ccccc1", "c1ccc2ccccc2c1", "C1CCCCC1", "CC(C)(C)C", "C1C2CC3CC1CC(C2)C3"],
)
def test_symmetric_molecules_resolve_via_seeding(smiles):
    # benzene, naphthalene, cyclohexane, neopentane, adamantane: no unique
    # anchors -> must be recovered by the anchor-free multi-orientation seeding.
    syms, coords = _mol_from_smiles(smiles)
    rng = np.random.default_rng(5)
    perm = rng.permutation(len(syms))
    tsyms = [syms[i] for i in perm]
    tcoords = coords[perm]

    res = am.map_atoms(syms, coords, tsyms, tcoords)

    assert res["verified"] is True
    assert res["n_bond_edits"] == 0
    # a valid isomorphism onto an exact permuted copy must reach ~0 RMSD
    assert res["rmsd"] < 1e-3


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


def test_twin_quotient_drops_the_interchangeable_hydrogens():
    """The VF2 fallback searches the twin quotient, not the full graph.

    Neopentane has 31104 element-preserving automorphisms, and 1296 of them only
    permute hydrogens within a methyl -- mappings that differ in nothing a graph
    can see.  The quotient keeps the 24 that are genuinely different.  Scoring
    the full set is what made a numbering check on a real ligand take minutes.
    """
    import networkx as nx

    syms, coords = _mol_from_smiles("CC(C)(C)C")
    G = am.connectivity(syms, coords)

    groups = am.false_twin_groups(G, syms)
    assert sorted(len(g) for g in groups if len(g) > 1) == [3, 3, 3, 3]

    quotient = am.twin_quotient(G, syms, {i: syms[i] for i in G.nodes}, groups)
    assert quotient.number_of_nodes() == 9   # 5 C + one node per methyl H triple
    assert quotient.number_of_edges() == 8   # the skeleton is untouched

    node_match = lambda a, b: a["el"] == b["el"] and a["mult"] == b["mult"]
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        quotient, quotient, node_match=node_match
    )
    assert sum(1 for _ in matcher.isomorphisms_iter()) == 24

    full = nx.algorithms.isomorphism.GraphMatcher(G, G, node_match=am._elem_match)
    assert sum(1 for _ in zip(full.isomorphisms_iter(), range(2000))) == 2000


def test_twins_are_never_adjacent_so_the_quotient_keeps_every_bond():
    # Members of a twin group share a neighbour set, which no atom can do with
    # something it is bonded to -- so collapsing a group never swallows a bond.
    syms, coords = _mol_from_smiles("CC(C)(C)c1ccccc1")
    G = am.connectivity(syms, coords)
    groups = am.false_twin_groups(G, syms)
    for members in groups:
        for i in members:
            for j in members:
                assert i == j or not G.has_edge(i, j)
    assert sum(len(g) for g in groups) == G.number_of_nodes()


def test_formula_mismatch_raises():
    a_syms, a_xyz = _metal_complex("Fe")
    b_syms, b_xyz = _metal_complex("Fe")
    b_syms = list(b_syms)
    b_syms[1] = "Br"  # Cl -> Br: composition no longer matches
    with pytest.raises(ValueError):
        am.map_atoms(a_syms, a_xyz, b_syms, b_xyz)
