"""Jahn-Teller frames for Cu(II) CN6 and square-pyramidal CN5 (DELFIN_FFFREE_ML_TARGETS_BV_JT, default OFF).

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


def _jt_len(table, cn: int, element: str, label: str):
    row = table.get(("J", "Cu", "2", str(cn), element, label))
    return row[1] if row else None


def _rotate_about(X, atoms, metal, a, p1, p2, target):
    """Rigidly rotate ``atoms`` about the line through p1 (and p2; with p2 None the line through p1
    perpendicular to the plane metal-a-p1) so that |metal - a| becomes ``target``; None if not
    reachable within 90 degrees."""
    M, A, P = X[metal], X[a], X[p1]
    n = (X[p2] - P) if p2 is not None else np.cross(A - P, M - P)
    if np.linalg.norm(n) < 1e-6:
        return None
    n = n / np.linalg.norm(n)

    def rot(v, th):
        c, s = np.cos(th), np.sin(th)
        return v * c + np.cross(n, v) * s + n * (n @ v) * (1 - c)

    best = None
    for k in range(1, 721):                      # both senses, 0.25 deg steps up to 90 deg
        th = np.radians(0.25 * ((k + 1) // 2) * (1 if k % 2 else -1))
        if abs(np.linalg.norm(P + rot(A - P, th) - M) - target) < 0.005:
            best = th
            break
    if best is None:
        return None
    Y = X.copy()
    for j in atoms:
        Y[j] = P + rot(X[j] - P, best)
    return Y


def _move_donor(Y, frags, frag_of, dons, mi, i, target, fixed):
    """Bring donor i to ``target`` A from the metal: translate a single-donor ligand; rotate a
    bi- or tridentate about its other donor(s) when none of them is in ``fixed`` motion."""
    others = [j for j in dons[i] if j != i]
    if any(j in fixed for j in others):
        return None
    if not others:
        dd = float(np.linalg.norm(Y[i] - Y[mi]))
        Z = Y.copy()
        shift = (target - dd) * (Y[i] - Y[mi]) / dd
        for j in frags[frag_of[i]]:
            Z[j] = Z[j] + shift
        return Z
    if len(others) == 1:
        return _rotate_about(Y, frags[frag_of[i]], mi, i, others[0], None, target)
    if len(others) == 2:
        return _rotate_about(Y, frags[frag_of[i]], mi, i, others[0], others[1], target)
    return None


def distort_frame(syms, X, table, form: str) -> Optional[np.ndarray]:
    """New coordinates of one Cu(II) frame in the given JT form, or None (left as it is).

    CN6 (``elong`` / ``comp``): the axis is a trans pair; each axial donor reaches the J-row axial
    length -- a single-donor ligand by rigid translation along its M-D bond, a bi- or tridentate by
    rigid rotation about its other donor(s), which stay in place.  The pair whose elements differ
    most between axial and equatorial length wins, ties to the lower atom index.  Equatorial
    single-donor ligands then take the equatorial length.
    CN5 square pyramid (``elong`` only): the apical donor (the one without a trans partner) takes the
    apical length, basal single-donor ligands the basal length."""
    from delfin.manta import _elements as EL
    from delfin.manta._d8_square_planar import _fragments
    mi = [i for i, s in enumerate(syms) if s == "Cu"]
    if len(mi) != 1 or sum(1 for s in syms if EL.is_metal(s)) != 1:
        return None
    mi = mi[0]
    d = np.linalg.norm(X - X[mi], axis=1)
    cand = [i for i in range(len(syms)) if i != mi and syms[i] != "H"
            and d[i] < 1.25 * (EL.covalent_radius("Cu") + EL.covalent_radius(syms[i]))]
    if len(cand) not in (5, 6):
        return None
    cn = len(cand)
    frags, frag_of = _fragments(syms, X, mi, cand)
    dons = {i: [j for j in cand if frag_of[j] == frag_of[i]] for i in cand}
    u = {i: (X[i] - X[mi]) / d[i] for i in cand}
    trans = {}
    for a in cand:
        b = min((x for x in cand if x != a), key=lambda x: float(u[a] @ u[x]))
        trans[a] = b if float(u[a] @ u[b]) < -0.8 else None
    if cn == 5:
        if form != "elong":
            return None
        lone = [i for i in cand if trans[i] is None or trans.get(trans[i]) != i]
        if len(lone) != 1:
            return None
        ax_label, eq_label = "sp5_ap", "sp5_bas"
        axis = (lone[0],)
    else:
        ax_label, eq_label = f"{form}_ax", f"{form}_eq"
        axis, best_key = None, None
        for a in cand:
            b = trans[a]
            if b is None or not a < b or trans.get(b) != a:
                continue
            if b in dons[a]:
                continue                                   # both axial donors in one ligand
            gap, ok = 0.0, True
            for i in (a, b):
                ax, eq = _jt_len(table, 6, syms[i], ax_label), _jt_len(table, 6, syms[i], eq_label)
                if ax is None or eq is None or len(dons[i]) > 3:
                    ok = False
                    break
                gap += abs(ax - eq)
            if ok and (best_key is None or (gap, -a) > best_key):
                axis, best_key = (a, b), (gap, -a)
        if axis is None:
            return None
    Y = X.copy()
    for i in axis:
        t = _jt_len(table, cn, syms[i], ax_label)
        if t is None:
            return None
        Z = _move_donor(Y, frags, frag_of, dons, mi, i, t, set(axis) - {i})
        if Z is None:
            return None
        Y = Z
    for i in cand:
        if i in axis or len(dons[i]) != 1:
            continue
        t = _jt_len(table, cn, syms[i], eq_label)
        if t is not None:
            Y = _move_donor(Y, frags, frag_of, dons, mi, i, t, set())
    return Y


def apply_jt(results, ctx, table) -> List[Tuple]:
    """Elongate (and with ``both`` add compressed copies of) the frames of a Cu(II) CN6 complex."""
    mode = jt_mode()
    if not mode or not results or not ctx or ctx.get("metal") != "Cu" or ctx.get("ox") != 2 \
            or ctx.get("cn") not in (5, 6):
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
