"""delfin.manta._me_bond_snap — terminal M=E multiple bonds on the FF-FREE path.

WHAT FOR (measured 17.08.2026).  The kind parameter ``kind="me"`` was "wired up" at three sites
on 16.08.  Recounted, the three are not three:

  * ``_clamp_metalloid_md_xyz``   -- METALLOID-only (Sb/As/Bi/Te/Se/Ge/Sn/Pb).  Sees not a
    single oxo metal; Re/W/V/Mo/Os are not in ``_METALLOID_MD_DONORS``.
  * ``_md_distance_in_tolerance`` -- by signature ``-> bool``.  A PREDICATE, it moves
    no atom.
  * ``_manual_metal_embed``       -- the CN4 path.  Terminal oxo complexes are CN5-7.

What remains is ``_snap_md_distances_to_ideal`` -- with EXACTLY ONE call site
(``smiles_converter.py:33066``), and that one lies BEHIND the FF-free return (``:32399``,
``return _ff, None``); both in the same body ``_smiles_to_xyz_isomers_impl:32156``.
=> The only real setter is LEGACY-ONLY, and FF-free is 98.8 % of the builds.
Measured consequence: ``me42b`` reached 2 of 42, ``byte_identical`` 40.

Fourth instance of the 10.08. pattern ("three additive modules only on LEGACY -> count CALL
SITES, never lines").  The mistake was counting OCCURRENCES of the kind parameter instead of
SETTING SITES on the live path.

WHY THIS MODULE TAKES NO ``mol``.  ``_ffree_shared_tail`` warns verbatim about the atom
order: a mol parsed from SMILES carries RDKit's order, FF-free frames carry metal-at-0 plus
AddHs blocks.  Precisely on this, the ring-pucker emitter at this spot once already turned
into a null test (185 of 187 byte-identical).  But the criterion for a terminal M=E does not
need the graph at all -- it is STRUCTURAL and readable from the frame itself.  Thus the trap
disappears instead of being worked around.

THE CRITERION is identical to ``smiles_converter._ml_bond_kind``, only read from the geometry
instead of from the graph: donor from {O, N, C}, NO bonded hydrogen, and its ONLY heavy
neighbour is a metal.  An OH/NH2 is thereby correctly not an oxo/imido; a bridging mu-oxo
(two metals) drops out via "exactly one heavy neighbour" -- likewise as in the original,
where a single translation could not satisfy two M-D ideals anyway.

WHY THE CORRECTION IS PARTICULARLY SAFE HERE.  A terminal donor by definition has no further
heavy neighbour and no hydrogen.  Hence exactly ONE atom is moved -- no BFS fragment, no
ligand body.  It cannot tear a bond, because the atom has none apart from the M-D bond.

Default OFF -> not called -> byte-identical.  If no calibrated band exists for a pair, the
caller returns ``None`` and this module touches nothing.  LICENSE: the values are
CCDC-derived, are NOT in this repo and are not read here either -- the caller passes them
in.  Deterministic (sorted order, no RNG), never a non-finite coordinate.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple

import numpy as np

from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

# Elements for which terminal M=E chemistry exists at all (oxo / nitrido / imido / carbyne).
_ME_ELEMENTS = frozenset({"O", "N", "C"})

# Hard lower bound for a heavy-heavy contact.  If the shortening brings the donor closer to
# a THIRD atom than this, it is rolled back -- the shortening does, after all, pull it into
# the coordination sphere.
_CLASH_FLOOR_A = 1.30

# Minimum relative change before anything is set at all.  Prevents noise on pairs whose
# band is practically 1.0 (62 % of the measured pairs lie between 0.95 and 1.05).
_MIN_REL_DELTA = 0.02


def terminal_me_pairs(syms: List[str], nbrs: List[List[int]]) -> List[Tuple[int, int]]:
    """(metal_idx, donor_idx) for every structurally terminal M=E donor, sorted."""
    out: List[Tuple[int, int]] = []
    for d, s in enumerate(syms):
        if s not in _ME_ELEMENTS:
            continue
        nb = nbrs[d] if d < len(nbrs) else []
        if any(syms[x] == "H" for x in nb):
            continue                                  # OH / NH2 / CH is not an oxo/imido
        heavy = [x for x in nb if syms[x] != "H"]
        if len(heavy) != 1:
            continue                                  # terminal: exactly ONE heavy neighbour
        m = heavy[0]
        if not _is_metal_sym(syms[m]):
            continue                                  # and that one must be the metal
        out.append((m, d))
    out.sort()
    return out


def snap_me_bonds(
    xyz_str: str,
    target_for: Callable[[str, str], Optional[float]],
) -> str:
    """Set every terminal M=E bond to its calibrated length.

    ``target_for(metal_symbol, donor_symbol)`` returns the target length in Angstrom, or
    ``None`` if NO calibrated band exists for the pair.  ``None`` explicitly means
    "do nothing" -- without a table the pass is byte-identical.

    Returns the unchanged input on every failure.
    """
    if not xyz_str:
        return xyz_str
    try:
        syms, pts, lines = _parse_xyz(xyz_str)
    except Exception:
        return xyz_str
    if not syms or pts is None or len(syms) < 2:
        return xyz_str
    try:
        nbrs, _blen = _build_geometric_adjacency(syms, pts)
    except Exception:
        return xyz_str

    pairs = terminal_me_pairs(syms, nbrs)
    if not pairs:
        return xyz_str

    new_pts = np.array(pts, dtype=float, copy=True)
    moved = False
    for m, d in pairs:
        try:
            target = target_for(syms[m], syms[d])
        except Exception:
            target = None
        if target is None:
            continue
        try:
            target = float(target)
        except Exception:
            continue
        if not np.isfinite(target) or target <= 0.0:
            continue
        v = new_pts[d] - new_pts[m]
        cur = float(np.linalg.norm(v))
        if cur < 1e-8:
            continue
        if abs(cur - target) / target < _MIN_REL_DELTA:
            continue
        cand = new_pts[m] + v * (target / cur)
        if not np.all(np.isfinite(cand)):
            continue
        # ROLLBACK: the donor must not come closer to a THIRD atom than the hard
        # lower bound.  The metal is exempt -- for it, the new distance is the target.
        others = [k for k in range(len(syms)) if k != d and k != m]
        if others:
            before = float(np.min(np.linalg.norm(new_pts[others] - new_pts[d], axis=1)))
            after = float(np.min(np.linalg.norm(new_pts[others] - cand, axis=1)))
            if after < before and after < _CLASH_FLOOR_A:
                continue                              # rejected, this donor stays
        new_pts[d] = cand
        moved = True

    if not moved:
        return xyz_str
    try:
        return _format_xyz(lines, syms, new_pts)
    except Exception:
        return xyz_str


# ---------------------------------------------------------------------------
# Self-test:  python delfin/manta/_me_bond_snap.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    def _dist(xyz, i, j):
        s, p, _l = _parse_xyz(xyz)
        return float(np.linalg.norm(p[i] - p[j]))

    fails = 0

    # 1) TERMINAL OXO is set to the target length.
    base = _xyz([("W", 0.0, 0.0, 0.0),
                 ("O", 2.10, 0.0, 0.0),      # terminal -> should move
                 ("Cl", 0.0, 2.30, 0.0),
                 ("Cl", 0.0, -2.30, 0.0),
                 ("Cl", 0.0, 0.0, 2.30)])
    got = snap_me_bonds(base, lambda m, d: 1.905 if (m, d) == ("W", "O") else None)
    d1 = _dist(got, 0, 1)
    ok = abs(d1 - 1.905) < 1e-4
    print(f"1 terminales W=O 2.100 -> {d1:.4f} (Ziel 1.905)  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 2) WITHOUT A BAND (target None) NOTHING happens -- byte-identical.
    same = snap_me_bonds(base, lambda m, d: None)
    ok = (same == base)
    print(f"2 ohne Band byte-identisch: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 3) HYDROXO is NOT touched (O carries an H).
    oh = _xyz([("W", 0.0, 0.0, 0.0),
               ("O", 2.10, 0.0, 0.0),
               ("H", 2.70, 0.90, 0.0),
               ("Cl", 0.0, 2.30, 0.0),
               ("Cl", 0.0, -2.30, 0.0)])
    ok = (snap_me_bonds(oh, lambda m, d: 1.905) == oh)
    print(f"3 Hydroxo unangetastet: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 4) BRIDGING mu-oxo (two metals) is NOT touched.
    mu = _xyz([("W", 0.0, 0.0, 0.0),
               ("O", 1.95, 0.0, 0.0),
               ("W", 3.90, 0.0, 0.0),
               ("Cl", 0.0, 2.30, 0.0),
               ("Cl", 3.90, 2.30, 0.0)])
    ok = (snap_me_bonds(mu, lambda m, d: 1.60) == mu)
    print(f"4 mu-Oxo unangetastet: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 5) The SECOND call changes nothing any more (idempotent).
    twice = snap_me_bonds(got, lambda m, d: 1.905 if (m, d) == ("W", "O") else None)
    ok = (twice == got)
    print(f"5 idempotent: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 6) Only the donor moves -- all other atoms stay exactly still.
    s0, p0, _ = _parse_xyz(base)
    s1, p1, _ = _parse_xyz(got)
    moved = [i for i in range(len(s0)) if float(np.linalg.norm(p0[i] - p1[i])) > 1e-9]
    ok = (moved == [1])
    print(f"6 bewegte Atome = {moved} (erwartet [1])  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    print(f"\n{6 - fails}/6 bestanden")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
