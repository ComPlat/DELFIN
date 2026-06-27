"""Per-atom mass weights for Baustein 5 PBD post-optimizer.

Determines how much each atom may move during constraint resolution.
Mass-weighted projection: alpha_i = w_i / (w_i + w_j) — heavier (lower w)
atoms move less.

Anchor hierarchy:
    Metals                  -> w = 0.0   (anchor, immovable)
    sigma-donors            -> w = 0.3   (semi-fixed; 0.5 in hapto class)
    hapto-atoms             -> w = 0.1   (rigid as group)
    bridge atoms (mu-X)     -> w = 0.2   (semi-fixed)
    H atoms                 -> w = 1.5   (light, easy correction)
    ligand-internal heavy   -> w = 1.0   (free)

Detection is chemistry-graph driven (not class-label primary). The
`class_label` argument only modulates the sigma-donor weight in hapto
contexts, mirroring the Wave-7 V e6761e4 finding that hapto-class
complexes need softer sigma-donors.
"""

from __future__ import annotations

from typing import Dict, Optional

# Mass weights (constants — keep in sync with module docstring)
W_METAL = 0.0
W_SIGMA = 0.3
W_SIGMA_HAPTO = 0.5   # softer sigma in hapto-class context
W_HAPTO = 0.1
W_BRIDGE = 0.2
W_HYDROGEN = 1.5
W_LIGAND = 1.0

# Metal symbol set: prefer canonical _METAL_SET from smiles_converter, with
# inline fallback if smiles_converter cannot be imported (e.g., during unit
# tests with stripped-down environments).
try:
    from delfin.smiles_converter import _METAL_SET as METAL_SYMBOLS  # type: ignore
except Exception:  # pragma: no cover - defensive fallback
    METAL_SYMBOLS = frozenset({
        # 3d transition metals
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        # 4d
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        # 5d
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        # f-block
        "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu", "Th", "U",
        # main-group sometimes-metal
        "Sn", "Pb", "Bi", "Tl", "In",
        # alkali / alkaline-earth
        "Li", "Na", "K", "Rb", "Cs",
        "Be", "Mg", "Ca", "Sr", "Ba",
    })


# ---------------------------------------------------------------------------
# Detection helpers
# ---------------------------------------------------------------------------

def _is_sigma_donor(mol, atom_idx: int) -> bool:
    """True if atom is bonded to >=1 metal (and is not a metal itself)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetSymbol() in METAL_SYMBOLS:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() in METAL_SYMBOLS:
            return True
    return False


def _is_hapto_atom(mol, atom_idx: int) -> bool:
    """True if atom is part of an eta-fragment.

    Heuristic: atom is in a ring, AND that ring contains >=2 atoms bonded
    to the same metal (metal-pi coordination).
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    if not atom.IsInRing():
        return False

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if atom_idx not in ring:
            continue
        # Count, per metal, how many ring atoms attach to it.
        metal_bonded_in_ring: Dict[int, list] = {}
        for ring_atom_idx in ring:
            ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetSymbol() in METAL_SYMBOLS:
                    metal_bonded_in_ring.setdefault(
                        nbr.GetIdx(), []
                    ).append(ring_atom_idx)
        for ring_atoms in metal_bonded_in_ring.values():
            if len(ring_atoms) >= 2:
                return True
    return False


def _is_bridge_atom(mol, atom_idx: int) -> bool:
    """True if atom is bonded to >=2 metals (mu-X bridging)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetSymbol() in METAL_SYMBOLS:
        return False
    metal_count = sum(
        1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() in METAL_SYMBOLS
    )
    return metal_count >= 2


# ---------------------------------------------------------------------------
# Main weight function
# ---------------------------------------------------------------------------

def get_atom_weight(mol, atom_idx: int, class_label: str = "sigma") -> float:
    """Return mass weight for atom. Higher = atom can move more.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule with explicit M-L bonds (rdkit-style).
    atom_idx : int
        Atom index in `mol`.
    class_label : str, default "sigma"
        Coordination class label. Only modulates the sigma-donor weight:
        in "hapto" / "multi_hapto" classes sigma-donors are softer (0.5
        instead of 0.3), per Wave-7 V e6761e4 finding.

    Returns
    -------
    float
        Mass weight in [0.0, 1.5].
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    sym = atom.GetSymbol()

    # Metals = anchor, never move.
    if sym in METAL_SYMBOLS:
        return W_METAL

    # Bridge atoms (mu-X bonded to >=2 metals) take precedence over plain
    # sigma — they are part of the polyhedron skeleton.
    if _is_bridge_atom(mol, atom_idx):
        return W_BRIDGE

    # sigma-donor: direct neighbor of a metal, non-metal itself.
    if _is_sigma_donor(mol, atom_idx):
        if class_label in ("hapto", "multi_hapto"):
            return W_SIGMA_HAPTO
        return W_SIGMA

    # Hapto-atom: member of an eta-fragment. Checked after sigma so that
    # actual metal-bonded ring atoms get sigma weight; only non-bonded
    # ring members of an eta-fragment receive the rigid-group weight.
    if _is_hapto_atom(mol, atom_idx):
        return W_HAPTO

    # H atom: light.
    if sym == "H":
        return W_HYDROGEN

    # Generic ligand-internal heavy atom.
    return W_LIGAND


def get_all_atom_weights(
    mol, class_label: str = "sigma"
) -> Dict[int, float]:
    """Batch helper. Return {atom_idx: weight} for every atom in `mol`."""
    return {
        i: get_atom_weight(mol, i, class_label=class_label)
        for i in range(mol.GetNumAtoms())
    }


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":  # pragma: no cover
    # Minimal synthetic test using a tiny RDKit mol.
    from rdkit import Chem

    # Build a tiny Fe complex with two pyridine sigma-donors via RWMol
    # so we don't depend on SMILES dative-bond parsing quirks.
    rw = Chem.RWMol()
    fe = rw.AddAtom(Chem.Atom("Fe"))
    n1 = rw.AddAtom(Chem.Atom("N"))
    n2 = rw.AddAtom(Chem.Atom("N"))
    c_lig = rw.AddAtom(Chem.Atom("C"))  # ligand-internal C
    rw.AddBond(fe, n1, Chem.BondType.SINGLE)
    rw.AddBond(fe, n2, Chem.BondType.SINGLE)
    rw.AddBond(n1, c_lig, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL
            ^ Chem.SANITIZE_PROPERTIES
            ^ Chem.SANITIZE_KEKULIZE,
        )
    except Exception:
        pass
    if mol is None:
        print("FAIL: could not build synthetic mol")
    else:
        try:
            mol = Chem.AddHs(mol)
        except Exception:
            pass
        ws = get_all_atom_weights(mol, class_label="sigma")
        print(f"n_atoms = {mol.GetNumAtoms()}")
        for i, w in ws.items():
            a = mol.GetAtomWithIdx(i)
            tag = []
            if a.GetSymbol() in METAL_SYMBOLS:
                tag.append("M")
            if _is_sigma_donor(mol, i):
                tag.append("sigma")
            if _is_hapto_atom(mol, i):
                tag.append("hapto")
            if _is_bridge_atom(mol, i):
                tag.append("bridge")
            print(f"  [{i:2d}] {a.GetSymbol():2s}  w={w:.2f}  {','.join(tag)}")

        # Quick assertions.
        fe_idx = next(
            i for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetSymbol() == "Fe"
        )
        assert get_atom_weight(mol, fe_idx) == W_METAL, "metal must be 0.0"

        # First non-metal donor neighbor of Fe.
        donor_idx = next(
            n.GetIdx()
            for n in mol.GetAtomWithIdx(fe_idx).GetNeighbors()
            if n.GetSymbol() not in METAL_SYMBOLS
        )
        assert get_atom_weight(mol, donor_idx, "sigma") == W_SIGMA
        assert (
            get_atom_weight(mol, donor_idx, "hapto") == W_SIGMA_HAPTO
        ), "hapto class must soften sigma to 0.5"

        # Any H should be 1.5.
        h_idx = next(
            (i for i in range(mol.GetNumAtoms())
             if mol.GetAtomWithIdx(i).GetSymbol() == "H"),
            None,
        )
        if h_idx is not None:
            assert get_atom_weight(mol, h_idx) == W_HYDROGEN

        print("OK: self-test passed")
