"""Jahn-Teller frames for Cu(II) (DELFIN_FFFREE_ML_TARGETS_BV_JT, default OFF).

WHY.  Cu(II) is d9: one hole in the e_g set, and the octahedron distorts.  Among the CSD crystals of
Cu(II) CN6 with three trans pairs, 36 % are elongated (two long axial bonds, ~2.39 A against 2.02 A
for N), 10 % compressed, 12 % undistorted, the rest rhombic.  The bond-valence targets give every
donor of one class the same, JT-averaged length, so no Cu(II) frame can match an elongated crystal.

WHAT.  After the build, each CN6 Cu(II) frame gets an axis: the trans pair whose two donors can move
on their own -- the only donor of its ligand (the ligand is translated as a rigid body along its M-D
bond) or one donor of a bidentate whose other donor is equatorial (the chelate is rotated rigidly
about that other donor) -- and whose elements the table's J rows set furthest apart between axial and
equatorial; ties go to the lower atom index.  Its donors move to the axial length, every other
movable donor to the equatorial length.  A frame without a movable trans pair is left as it is.
``elong``: the frames are elongated in place (label suffix -jte).  ``all``: the frame as built is kept
and an elongated and a compressed copy are appended.  ``both``: elongated in place and
a compressed copy of each appended (-jtc), the way the spin pass appends its frames; the
finalisation chain judges every frame on its own.

Default OFF, and without a bond-valence context for a Cu(II) CN6 complex the frames are untouched.
"""
from __future__ import annotations

import os
from typing import List, Optional, Tuple

import numpy as np

JT_ENV = "DELFIN_FFFREE_ML_TARGETS_BV_JT"


def jt_mode() -> str:
    v = os.environ.get(JT_ENV, "").strip().lower()
    return v if v in ("elong", "both", "all") else ""


def _jt_lengths(table, element: str, form: str):
    ax = table.get(("J", "Cu", "2", "6", element, f"{form}_ax"))
    eq = table.get(("J", "Cu", "2", "6", element, f"{form}_eq"))
    return (ax[1] if ax else None), (eq[1] if eq else None)


def _rotate_to_length(X, atoms, metal, a, pivot, target):
    """Rigidly rotate ``atoms`` about the axis through ``pivot`` perpendicular to the plane
    (metal, a, pivot) so that |metal - a| becomes ``target``; None if unreachable."""
    M, A, P = X[metal], X[a], X[pivot]
    n = np.cross(A - P, M - P)
    if np.linalg.norm(n) < 1e-6:
        return None
    n /= np.linalg.norm(n)
    rr = A - P
    best = None
    for k in range(1, 721):                      # scan both senses, 0.25 deg steps up to 90 deg
        th = np.radians(0.25 * ((k + 1) // 2) * (1 if k % 2 else -1))
        c, s = np.cos(th), np.sin(th)
        Ar = P + rr * c + np.cross(n, rr) * s + n * (n @ rr) * (1 - c)
        if abs(np.linalg.norm(Ar - M) - target) < 0.005:
            best = th
            break
    if best is None:
        return None
    c, s = np.cos(best), np.sin(best)
    Y = X.copy()
    for j in atoms:
        v = X[j] - P
        Y[j] = P + v * c + np.cross(n, v) * s + n * (n @ v) * (1 - c)
    return Y


def distort_frame(syms, X, table, form: str) -> Optional[np.ndarray]:
    """New coordinates of one Cu(II) CN6 frame in the given JT form, or None (left as it is).

    A donor can take its JT length when it is the only donor of its ligand (rigid translation along
    its M-D bond) or one of the two donors of a bidentate ligand whose other donor is not on the axis
    (rigid rotation of the chelate about that other donor).  The axis is the trans pair of such donors
    whose elements differ most between the axial and equatorial lengths; ties to the lower index."""
    from delfin.manta import _elements as EL
    from delfin.manta._d8_square_planar import _fragments
    mi = [i for i, s in enumerate(syms) if s == "Cu"]
    if len(mi) != 1 or sum(1 for s in syms if EL.is_metal(s)) != 1:
        return None
    mi = mi[0]
    d = np.linalg.norm(X - X[mi], axis=1)
    cand = [i for i in range(len(syms)) if i != mi and syms[i] != "H"
            and d[i] < 1.25 * (EL.covalent_radius("Cu") + EL.covalent_radius(syms[i]))]
    if len(cand) != 6:
        return None
    frags, frag_of = _fragments(syms, X, mi, cand)
    dons = {i: [j for j in cand if frag_of[j] == frag_of[i]] for i in cand}
    u = {i: (X[i] - X[mi]) / d[i] for i in cand}

    def handle(i, other_axial):
        if len(dons[i]) == 1:
            return "mono"
        if len(dons[i]) == 2 and other_axial not in dons[i]:
            return "bi"
        return None

    best, best_key = None, None
    for a in cand:
        b = min((x for x in cand if x != a), key=lambda x: float(u[a] @ u[x]))
        if not (a < b and float(u[a] @ u[b]) < -0.8):
            continue
        if handle(a, b) is None or handle(b, a) is None:
            continue
        gap, ok = 0.0, True
        for i in (a, b):
            ax, eq = _jt_lengths(table, syms[i], form)
            if ax is None or eq is None:
                ok = False
                break
            gap += abs(ax - eq)
        if ok and (best_key is None or (gap, -a) > best_key):
            best, best_key = (a, b), (gap, -a)
    if best is None:
        return None
    Y = X.copy()
    for i in best:
        ax, _eq = _jt_lengths(table, syms[i], form)
        if handle(i, best[0] if i == best[1] else best[1]) == "mono":
            shift = (ax - d[i]) * u[i]
            for j in frags[frag_of[i]]:
                Y[j] = Y[j] + shift
        else:
            pivot = [j for j in dons[i] if j != i][0]
            Z = _rotate_to_length(Y, frags[frag_of[i]], mi, i, pivot, ax)
            if Z is None:
                return None
            Y = Z
    for i in cand:                               # equatorial single-donor ligands: their length
        if i in best or len(dons[i]) != 1:
            continue
        _ax, eq = _jt_lengths(table, syms[i], form)
        if eq is not None:
            dd = float(np.linalg.norm(Y[i] - Y[mi]))
            shift = (eq - dd) * (Y[i] - Y[mi]) / dd
            for j in frags[frag_of[i]]:
                Y[j] = Y[j] + shift
    return Y


def apply_jt(results, ctx, table) -> List[Tuple]:
    """Elongate (and with ``both`` add compressed copies of) the frames of a Cu(II) CN6 complex."""
    mode = jt_mode()
    if not mode or not results or not ctx or ctx.get("metal") != "Cu" or ctx.get("ox") != 2 \
            or ctx.get("cn") != 6:
        return results
    from delfin.manta._fix_sp2n_planarize import _format_xyz, _parse_xyz
    out, extra = [], []
    for e in results:
        try:
            xyz, lbl = e[0], e[1]
            syms, X, lines = _parse_xyz(xyz)
            Ye = distort_frame(syms, X, table, "elong")
            if mode == "all":
                out.append(e)
                if Ye is not None:
                    extra.append((_format_xyz(lines, syms, Ye), (f"{lbl}-jte" if lbl else "jte")) + tuple(e[2:]))
            else:
                out.append(((_format_xyz(lines, syms, Ye), (f"{lbl}-jte" if lbl else "jte")) + tuple(e[2:]))
                           if Ye is not None else e)
            if mode in ("both", "all"):
                Yc = distort_frame(syms, X, table, "comp")
                if Yc is not None:
                    extra.append((_format_xyz(lines, syms, Yc), (f"{lbl}-jtc" if lbl else "jtc"))
                                 + tuple(e[2:]))
        except Exception:
            if not out or out[-1] is not e:
                out.append(e)
    return out + extra
