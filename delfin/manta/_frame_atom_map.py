"""_frame_atom_map.py — Atom mapping FRAME -> MOL for the FF-free path.

THE PROBLEM, MEASURED ON 16.08.2026.  The builder's post-hoc correctors need
`mol` ATOM INDICES and apply them to XYZ COORDINATES.  On the legacy path that is correct,
because the XYZ comes from the same `mol`.  On the FF-FREE path it is not: there the
metal sits at 0, followed by one `AddHs` block per ligand in BUILD ORDER
(`converter_backend._heteroleptic_block_offsets`, `_config_block_offsets`) -- a completely
different ordering from RDKit's.

WHAT DAMAGE THIS DOES.  `_ffree_shared_tail` warns about it in so many words; the same break
has already turned the ring-pucker emitter there into a null lever (185 of 187 systems
byte-identical).  And the planarizers' guard checks only the NUMBER of atoms: with the same
count and a permuted order they would have flattened the WRONG atoms -- not silent
ineffectiveness, but silent destruction.  (Since 16.08. a latch there checks the order,
`864cba3f`.)

⚠ WHY NOT CARRY THE MAPPING ALONG.  The obvious approach would be to transport the offsets
known at BUILD TIME together with the frame.  Two reasons against it: the XYZ comment line
already carries the LABEL ("BEBGUL frame0 alt-bind-C"), and labels are parsed downstream
(`_arrangement_key`, arrangement families) -- appending anything there endangers them.  And
the tuple width `(xyz, label)` is assumed in many places.

THE WAY TAKEN HERE: RECONSTRUCT, DO NOT GUESS.  From the frame's geometry a graph is
built (the same adjacency the correctors use anyway), from that an RDKit molecule
with GENERIC bonds, and then RDKit searches for the subgraph match.  A
COMPLETE substructure match over all atoms IS the isomorphism -- it preserves elements
and connectivity by definition.  This is not a heuristic with a fallback chain.

⚠⚠ GENERIC BONDS ARE MANDATORY.  The frame graph only knows "bonded yes/no" (it comes
from distances), `mol` knows single/double/aromatic/dative.  A match with bond orders
would ALWAYS fail -- and silently, as "no mapping".  `makeBondsGeneric` lifts that.

WHEN IN DOUBT, NOTHING.  No match, ambiguous match, unequal atom count or missing
RDKit -> `None`.  The caller then falls back to its current behavior, and that is
byte-identical.  A WRONG mapping would be worse than none: it would let the correctors
touch exactly the wrong atoms -- the error that the latch from 16.08. catches.
"""
from __future__ import annotations

import logging
from typing import List, Optional, Sequence

_LOG = logging.getLogger(__name__)


def frame_to_mol_map(mol, syms: Sequence[str], nbrs: List[List[int]]) -> Optional[List[int]]:
    """Mapping FRAME index -> MOL index, or None.

    ``nbrs`` is the frame's adjacency (from `_build_geometric_adjacency`), i.e. exactly the
    graph the correctors work with anyway -- no second bond definition.

    Returns: ``m[frame_idx] = mol_idx``.  None EXPLICITLY means "not determinable".
    """
    if mol is None:
        return None
    try:
        from rdkit import Chem
    except Exception:
        return None
    n = len(syms)
    if n == 0 or mol.GetNumAtoms() != n:
        return None
    # Rebuild the frame graph as a molecule -- elements and edges, no bond orders.
    try:
        rw = Chem.RWMol()
        for s in syms:
            a = Chem.Atom(s)
            a.SetNoImplicit(True)
            rw.AddAtom(a)
        seen = set()
        for i in range(n):
            for j in nbrs[i]:
                if j <= i or (i, j) in seen:
                    continue
                seen.add((i, j))
                rw.AddBond(i, j, Chem.BondType.SINGLE)
        frame = rw.GetMol()
    except Exception as exc:
        _LOG.debug("frame-map: Frame-Graph nicht baubar: %s", exc)
        return None
    # Generic bonds on BOTH sides: the frame knows no bond orders, `mol` does.
    # Without this the match always fails -- and silently.
    try:
        params = Chem.AdjustQueryParameters.NoAdjustments()
        params.makeBondsGeneric = True
        params.makeDummiesQueries = False
        q = Chem.AdjustQueryProperties(frame, params)
        target = Chem.AdjustQueryProperties(Chem.Mol(mol), params)
    except Exception as exc:
        _LOG.debug("frame-map: Bindungen nicht generisch zu machen: %s", exc)
        return None
    try:
        matches = target.GetSubstructMatches(q, useChirality=False, uniquify=False,
                                             maxMatches=2)
    except Exception as exc:
        _LOG.debug("frame-map: Substruktursuche fehlgeschlagen: %s", exc)
        return None
    if not matches:
        return None
    m = list(matches[0])
    if len(m) != n:
        return None
    # Element equality is guaranteed by the search; it is checked here anyway.  A
    # mapping that is not recomputed is a claim -- and a measurement has already failed
    # on exactly that once today.
    try:
        for fi, mi in enumerate(m):
            if mol.GetAtomWithIdx(int(mi)).GetSymbol() != syms[fi]:
                return None
    except Exception:
        return None
    return m


if __name__ == "__main__":  # pragma: no cover -- self-test
    # TWO CASES, and the second is the real one:
    #   (a) same order      -> the mapping MUST be the identity
    #   (b) PERMUTED order  -> it must return the permutation
    # If (b) fails, the module is worthless: that is exactly what it exists for.
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from delfin.manta._coord_angle_corrector import _build_geometric_adjacency
    import numpy as _np

    smi = "CC(=O)C=C(C)[O][Cu][O]C(C)=CC(C)=O"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print("  SMILES nicht lesbar"); raise SystemExit(1)
    molH = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    if AllChem.EmbedMolecule(molH, randomSeed=1) != 0:
        print("  Einbettung fehlgeschlagen"); raise SystemExit(1)
    conf = molH.GetConformer()
    syms = [a.GetSymbol() for a in molH.GetAtoms()]
    pts = _np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                      conf.GetAtomPosition(i).z] for i in range(molH.GetNumAtoms())])
    nbrs, _ = _build_geometric_adjacency(syms, pts)

    m = frame_to_mol_map(molH, syms, nbrs)
    ident = (m == list(range(len(syms)))) if m else False
    print(f"  (a) gleiche Ordnung  -> {'IDENTITAET' if ident else m if m else 'None'}")

    # (b) Reorder the frame: metal to the front, as the FF-free builder does.
    mi = next(i for i, s in enumerate(syms) if s == "Cu")
    order = [mi] + [i for i in range(len(syms)) if i != mi]
    inv = {old: new for new, old in enumerate(order)}
    syms2 = [syms[o] for o in order]
    nbrs2 = [[inv[k] for k in nbrs[o]] for o in order]
    m2 = frame_to_mol_map(molH, syms2, nbrs2)
    ok2 = bool(m2) and all(molH.GetAtomWithIdx(m2[f]).GetSymbol() == syms2[f]
                           for f in range(len(syms2)))
    print(f"  (b) Metall-auf-0    -> {'ZUORDNUNG GEFUNDEN' if ok2 else (m2 if m2 else 'None')}")
    if m2:
        print(f"      Frame 0 ({syms2[0]}) -> mol {m2[0]} ({molH.GetAtomWithIdx(m2[0]).GetSymbol()})")
