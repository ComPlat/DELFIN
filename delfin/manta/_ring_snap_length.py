"""_ring_snap_length -- a ring projected onto its plane keeps its bond lengths
(DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH, default OFF -> byte-identical).

THE DEFECT.  Three passes of the legacy post-processing flatten aromatic rings by a
"minimal-movement projection": every ring atom loses the component of its position along
the ring normal (``ligand_placement._snap_aromatic_rings_to_plane``,
``ligand_placement._snap_fused_aromatic_groups_to_plane``, ``_arom_planarize._flatten_system``).
Their docstrings say the bond lengths are preserved "to first order" -- true for a ring that
is flat to within a few hundredths of an Angstrom, false for the rings the legacy builder
hands over: a bond tilted by 25 degrees out of the fitted plane shrinks by cos 25 = 0.91,
an aromatic C-N of 1.34 A becomes 1.22 A.  Measured on the H1 projection step of NAMFUE
(Ru bipyridine/terpyridine, DELFIN_LEGCONF_TRACE census, 23.10.2026): frames with a ligand
bond below 0.88 x its bond-order reference go 25 -> 58 of 119 in this one step; the frames
handed in were puckered, not compressed.

THE FIX.  After the projection the ring (or fused group) is scaled IN ITS PLANE about a
pivot so that its mean intra-ring bond length is the one it had before the projection.
The pivot is the ring's one metal-bonded atom when it has exactly one (the M-D distance
is then untouched by the scaling), the centroid when it has none, and with two or more
metal-bonded ring atoms the scaling is skipped (it would move a donor).  Hydrogens bonded
to a ring atom move with that atom's displacement, so every ring X-H keeps its length.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def ring_snap_keep_length_enabled() -> bool:
    """THE one read site of DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH (default 0)."""
    return os.environ.get("DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH", "0") == "1"


def ring_scale_restoring_lengths(
    before: np.ndarray,
    after: np.ndarray,
    bonds: Sequence[Tuple[int, int]],
    pivot: Optional[int],
    skip: bool = False,
) -> np.ndarray:
    """Positions ``after`` (rows = ring atoms, already projected onto one plane) scaled about
    row ``pivot`` (None -> their centroid) so the mean length of ``bonds`` (row pairs) equals
    its value in ``before``.  Returns ``after`` unchanged when ``skip``, when there are no
    bonds, or when the projection did not shorten them."""
    if skip or not bonds:
        return after
    lb = float(np.mean([np.linalg.norm(before[i] - before[j]) for i, j in bonds]))
    la = float(np.mean([np.linalg.norm(after[i] - after[j]) for i, j in bonds]))
    if la < 1e-9 or lb <= la * (1.0 + 1e-6):
        return after
    c = after[pivot] if pivot is not None else after.mean(axis=0)
    return c + (lb / la) * (after - c)


def ride_hydrogens(
    displacement: Dict[int, np.ndarray],
    hydrogens_of: Dict[int, List[int]],
) -> Dict[int, np.ndarray]:
    """Per hydrogen the displacement of its ring parent (first parent wins)."""
    out: Dict[int, np.ndarray] = {}
    for parent, d in displacement.items():
        for h in hydrogens_of.get(parent, ()):
            out.setdefault(h, d)
    return out
