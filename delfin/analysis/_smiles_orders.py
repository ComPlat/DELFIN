"""Bond orders from the calc's SMILES, transferred to the XYZ graph.

The geometric heuristic ``_rdkit_bond_orders`` in ``results_overview.py``
fails on metal-coordinated heteroaromatics (pyrrole, pyrazole, pyridine
ligands of a Co complex, calc 428): the coordinating N has no valence
deficit the heuristic can see, so the ring double bonds never appear.
The calc's own SMILES (``DELFIN_Data.json`` ``input.smiles``) carries the
correct orders.  This module kekulizes that SMILES and transfers the
orders onto the XYZ heavy-atom graph by a backtracking isomorphism that
matches element AND neighbourhood, so the mapping is unambiguous even
for symmetric ligands.
"""

import itertools
import json
import os

DELFIN_DATA_NAME = "DELFIN_Data.json"


def _load_smiles(calc_dir):
    """SMILES string from DELFIN_Data.json ``input.smiles``, or None."""
    path = os.path.join(calc_dir, DELFIN_DATA_NAME)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    value = (data.get("input") or {}).get("smiles") or ""
    return str(value) or None


def _kekule_orders(smiles):
    """Kekulize the SMILES; return (atoms, orders) or None.

    atoms: list of element symbols.  orders: {frozenset({a, b}): int}.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append(atom.GetSymbol())
    orders = {}
    BOND_ORDERS = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
    }
    for bond in mol.GetBonds():
        order = BOND_ORDERS.get(bond.GetBondType())
        if order is None:
            continue
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        orders[frozenset((a, b))] = order
    return atoms, orders


def _neighbours(n_atoms, bonds):
    """Adjacency sets over heavy-atom indices."""
    adj = [set() for _ in range(n_atoms)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def smiles_bond_orders(sym, bonds, calc_dir):
    """Bond orders for the XYZ graph, taken from the calc's SMILES.

    Parameters
    ----------
    sym : list of element symbols (H included or not, metals included).
    bonds : iterable of (i, j) index pairs over ``sym``.
    calc_dir : directory holding ``DELFIN_Data.json``.

    Returns
    -------
    dict mapping (i, j) (both orderings, i < j) to order 2 or 3, or
    None when no SMILES is available, RDKit is missing, or no
    neighbourhood-preserving isomorphism between the XYZ heavy-atom
    graph and the SMILES molecule graph is found.
    """
    if not bonds:
        return None
    parsed = _kekule_orders(_load_smiles(calc_dir))
    if parsed is None:
        return None
    s_atoms, s_orders = parsed

    # Heavy atoms only: metals and H carry no orders to transfer and
    # are excluded from the isomorphism on both sides.
    _NONHEAVY = {"H"} | {
        "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Al", "Ga", "In", "Tl", "Pb", "Bi",
        "Sn", "Sb", "B", "Si", "P", "S", "Cl", "Br", "I", "F", "C",
        "N", "O", "Se", "Te",
    }
    _METAL = _NONHEAVY - {"H", "B", "C", "N", "O", "F", "P", "S", "Cl",
                          "Br", "I", "Se", "Si", "Te"}
    x_keep = [i for i, s in enumerate(sym) if s not in _METAL and s != "H"]
    s_keep = [k for k, s in enumerate(s_atoms)
              if s not in _METAL and s != "H"]
    if len(x_keep) != len(s_keep):
        return None

    x_adj = _neighbours(len(sym), bonds)
    x_adj = [[j for j in x_adj[i] if j in set(x_keep)] for i in range(len(sym))]
    s_adj_full = [[] for _ in s_atoms]
    for pair in s_orders:
        a, b = tuple(pair)
        s_adj_full[a].append(b)
        s_adj_full[b].append(a)
    s_adj = {k: [j for j in s_adj_full[k] if j in set(s_keep)]
             for k in s_keep}

    # Compatibility: element equal, degree equal, and the multiset of
    # neighbour elements equal.  This prunes the backtracking hard.
    s_elem = {k: s_atoms[k] for k in s_keep}
    x_elem = {i: sym[i] for i in x_keep}
    x_deg = {i: len(x_adj[i]) for i in x_keep}
    s_deg = {k: len(s_adj[k]) for k in s_keep}
    x_nbelem = {i: sorted(sym[j] for j in x_adj[i]) for i in x_keep}
    s_nbelem = {k: sorted(s_atoms[j] for j in s_adj[k]) for k in s_keep}

    def _consistent(xi, si, mapping):
        if x_elem[xi] != s_elem[si]:
            return False
        if x_deg[xi] != s_deg[si]:
            return False
        if x_nbelem[xi] != s_nbelem[si]:
            return False
        # already-mapped neighbours must agree too
        for j in x_adj[xi]:
            sj = mapping.get(j)
            if sj is not None and sj not in s_adj[si]:
                return False
        for sk in s_adj[si]:
            xj = None
            for j, mv in mapping.items():
                if mv == sk:
                    xj = j
                    break
            if xj is not None and xj not in x_adj[xi]:
                return False
        return True

    candidates = {xi: [si for si in s_keep
                       if x_elem[xi] == s_elem[si]
                       and x_deg[xi] == s_deg[si]
                       and x_nbelem[xi] == s_nbelem[si]]
                  for xi in x_keep}
    if any(not c for c in candidates.values()):
        return None

    # Order the XYZ atoms by fewest candidates first (fail fast).
    order = sorted(x_keep, key=lambda xi: len(candidates[xi]))
    mapping = {}
    used = set()

    def _backtrack(pos):
        if pos == len(order):
            return True
        xi = order[pos]
        for si in candidates[xi]:
            if si in used:
                continue
            if not _consistent(xi, si, mapping):
                continue
            mapping[xi] = si
            used.add(si)
            if _backtrack(pos + 1):
                return True
            del mapping[xi]
            used.discard(si)
        return False

    if not _backtrack(0):
        return None

    # Transfer the orders through the isomorphism.  Only bonds that
    # exist in the XYZ graph get an order; a SMILES edge with no XYZ
    # counterpart is skipped silently.
    inv = {si: xi for xi, si in mapping.items()}
    _bond_set = {tuple(sorted((i, j))) for i, j in bonds}
    result = {}
    for pair, order_val in s_orders.items():
        a, b = tuple(pair)
        xa, xb = inv.get(a), inv.get(b)
        if xa is None or xb is None:
            continue
        if tuple(sorted((xa, xb))) not in _bond_set:
            continue
        if order_val >= 2:
            result[(xa, xb)] = order_val
            result[(xb, xa)] = order_val
    return result or None
