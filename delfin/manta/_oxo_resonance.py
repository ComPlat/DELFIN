"""delfin.manta._oxo_resonance -- resonance lengths of degenerate terminal groups on the LEGACY path
(DELFIN_FFFREE_OXO_RESONANCE, 28.10.2026).

THE DEFECT.  A nitrate, nitrite, nitro group, carboxylate, carbonate, sulfate, phosphate, amidinate ... is written
in the SMILES as ONE Lewis structure (O=[N+]([O-])[O-]: one double, two single bonds).  Every length the legacy
path gives these bonds comes from that one structure: the whole-complex ETKDG template (measured: nitrate kappa2
N-O 1.17-2.14 A on OPEVEN / WEMKUX / AGOXIE) and Open Babel UFF, which relaxes each "single" N-O to its single
length (1.36-1.50 A in the built frames; the crystal has ~1.27 A for all three).  The eye reads ligspec_nitrate
bad_bond_len (dev300: 87 frames without the legacy fixes, 409 with the fold window + LP weight, where more frames
survive).  Nothing after UFF restores the resonance.

THE RULE is the one DELFIN_FFFREE_MESOMERY_SEAT (assemble_ligand_confs) applies on the FF-free path, which the
legacy path never reaches: for every heavy atom X its TERMINAL heavy neighbours (no other non-metal heavy
neighbour, no H -- a protonated O/N is not resonance-equivalent) are grouped by element; a group with >= 2
members AND different bond orders is degenerate, its bond order is the mean of the SMILES orders, and its target
length the Pyykko single/double radius interpolation at that order.  Each terminal atom moves along its X-Y
direction to the target.  A frame is replaced only if the moved frame still passes the legacy topology verifier;
the frame's atom order must be the SMILES mol's (legacy batch results), else the frame is left unchanged.
Off -> byte-identical.  Public radii only (Pyykko & Atsumi 2009 via delfin.manta.polyhedra); no CSD/CCDC data.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import List, Tuple

import numpy as np

from delfin.manta.ml_tables import _METAL_SET
from delfin.manta.polyhedra import PYYKKO_DOUBLE, PYYKKO_SINGLE


def oxo_resonance_enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_OXO_RESONANCE", "0") == "1"


def degenerate_terminal_sets(mol) -> List[Tuple[int, List[int], float]]:
    """(centre X, [terminal Y of one element], target X-Y length A) for every degenerate group."""
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in _METAL_SET or a.GetAtomicNum() <= 1:
            continue
        groups = defaultdict(list)
        for b in a.GetBonds():
            y = b.GetOtherAtom(a)
            if y.GetAtomicNum() <= 1 or y.GetSymbol() in _METAL_SET:
                continue
            if any(n.GetAtomicNum() == 1 for n in y.GetNeighbors()):
                continue
            if any(n.GetIdx() != a.GetIdx() and n.GetAtomicNum() > 1 and n.GetSymbol() not in _METAL_SET
                   for n in y.GetNeighbors()):
                continue
            bo = b.GetBondTypeAsDouble()
            groups[y.GetSymbol()].append((y.GetIdx(), bo if bo > 0 else 1.0))
        for ysym, mem in groups.items():
            orders = [o for _i, o in mem]
            if len(mem) < 2 or len(set(orders)) < 2:
                continue
            rs_x, rd_x = PYYKKO_SINGLE.get(a.GetSymbol()), PYYKKO_DOUBLE.get(a.GetSymbol())
            rs_y, rd_y = PYYKKO_SINGLE.get(ysym), PYYKKO_DOUBLE.get(ysym)
            if None in (rs_x, rd_x, rs_y, rd_y):
                continue
            f = min(1.0, max(0.0, float(np.mean(orders)) - 1.0))
            target = (rs_x + rs_y) - f * ((rs_x + rs_y) - (rd_x + rd_y))
            out.append((a.GetIdx(), [i for i, _o in mem], target))
    return out


def resonance_xyz(xyz: str, mol, sets) -> str:
    lines = [ln for ln in xyz.strip().splitlines() if ln.strip()]
    if len(lines) != mol.GetNumAtoms():
        return xyz
    syms = [ln.split()[0] for ln in lines]
    if any(s != mol.GetAtomWithIdx(i).GetSymbol() for i, s in enumerate(syms)):
        return xyz
    P = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines])
    moved = False
    for c, term, target in sets:
        for y in term:
            v = P[y] - P[c]
            d = float(np.linalg.norm(v))
            if d < 1e-6 or abs(d - target) < 0.02 or abs(d - target) > 0.9:
                continue
            P[y] = P[c] + v / d * target
            moved = True
    if not moved:
        return xyz
    return "\n".join(f"{s:4s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(syms, P)) + "\n"


def resonance_results(mol, results):
    """Apply to every (xyz, label) frame; a frame stays unchanged when the verifier rejects the moved one."""
    if not results or mol is None:
        return results
    try:
        sets = degenerate_terminal_sets(mol)
    except Exception:
        return results
    if not sets:
        return results
    from delfin.manta.topo_isomers import _verify_topology_from_graph
    out = []
    for item in results:
        try:
            xyz, label = item
            new = resonance_xyz(xyz, mol, sets)
            if new is not xyz and _verify_topology_from_graph(new, mol):
                out.append((new, label))
            else:
                out.append(item)
        except Exception:
            out.append(item)
    return out
