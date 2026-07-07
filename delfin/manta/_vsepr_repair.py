"""Universal VSEPR terminal-group repair (license-clean, pure-numpy, no FF).

The metal-complex embed/seating sometimes distorts a LIGAND's local geometry
while the coordination sphere is fine -- e.g. a CF3 left with an F-C-F of 157deg
(LUXWIL), a distorted SO3, a splayed CCl3.  A terminal EX3 group is RIGID: it is
tetrahedral (E has a 4th bond) or trigonal (BF3), with no conformational freedom,
so any distortion is unambiguously a construction defect and can be repaired
deterministically by placing the terminal atoms at their ideal VSEPR positions
around the central atom -- anchored by the central atom's other bond so the
group's attachment (and the rest of the molecule) is untouched.

Fundamental + universal: works for CF3 / CCl3 / CBr3 / CH3 / SO3 / PO3 / any
central atom carrying >=3 terminal same-element neighbours, metal or organic.
Only fires on a genuinely-distorted group (worst intra-group angle > tol);
byte-identical when every terminal group is already near-ideal.
"""
from typing import List, Optional, Tuple

try:
    import numpy as _np
except Exception:                                    # pragma: no cover
    _np = None

_METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Th", "U",
}
_COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05,
        "Cl": 1.02, "Br": 1.20, "I": 1.39, "P": 1.07, "B": 0.84, "Si": 1.11,
        "Se": 1.20}


def _parse(xyz: str):
    lines = [l for l in xyz.strip().splitlines() if l.strip()]
    # tolerate a leading count / comment header
    if lines and lines[0].split()[0].isdigit():
        try:
            n = int(lines[0].split()[0])
            lines = lines[2:2 + n] if len(lines) >= 2 + n else lines[1:]
        except Exception:
            pass
    syms, P = [], []
    for l in lines:
        t = l.split()
        if len(t) >= 4:
            syms.append(t[0])
            P.append([float(t[1]), float(t[2]), float(t[3])])
    return syms, _np.asarray(P, dtype=float)


def _format(syms, P) -> str:
    return "".join(f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n"
                   for s, p in zip(syms, P))


def _adjacency(syms, P):
    n = len(syms)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(_np.linalg.norm(P[i] - P[j]))
            cut = (_COV.get(syms[i], 0.9) + _COV.get(syms[j], 0.9)) * 1.3
            if d < cut:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _worst_angle_dev(P, c, members, ideal):
    worst = 0.0
    for a in range(len(members)):
        for b in range(a + 1, len(members)):
            v1 = P[members[a]] - P[c]
            v2 = P[members[b]] - P[c]
            nn = float(_np.linalg.norm(v1) * _np.linalg.norm(v2))
            if nn < 1e-9:
                continue
            ang = _np.degrees(_np.arccos(max(-1.0, min(1.0, float(_np.dot(v1, v2) / nn)))))
            worst = max(worst, abs(ang - ideal))
    return worst


def repair_terminal_groups(xyz: str, tol: float = 20.0) -> str:
    """Repair distorted terminal EX3 groups to their ideal VSEPR geometry.
    Returns the repaired XYZ (same atom order) or the input unchanged when every
    terminal group is already within ``tol`` degrees of ideal."""
    if _np is None:
        return xyz
    try:
        syms, P = _parse(xyz)
        if len(syms) < 4:
            return xyz
        adj = _adjacency(syms, P)
        n = len(syms)
        changed = False
        Pn = P.copy()
        for c in range(n):
            if syms[c] == "H" or syms[c] in _METALS:
                continue
            # terminal neighbours (no heavy neighbour other than c), grouped by element
            groups = {}
            others = []
            for k in adj[c]:
                heavy_other = [x for x in adj[k] if x != c and syms[x] != "H"]
                # terminal neighbour (nothing else heavy attached) -> a group
                # member.  H IS included so a distorted methyl CH3 / ammonium NH3
                # is repaired too (user: the CF3 problem also hits CH3 in some
                # systems); heavy terminals (F/Cl/O of CF3/SO3) as before.
                if not heavy_other and syms[k] not in _METALS:
                    groups.setdefault(syms[k], []).append(k)
                elif syms[k] != "H":
                    others.append(k)
            for el, members in groups.items():
                if len(members) != 3:
                    continue
                # sp3 tetrahedral case needs a heavy anchor bond (CF3-C, SO3-S,
                # ...); a bare 3-terminal centre (BF3, NX3) is trigonal/pyramidal
                # -> handled separately, skip.  A SEVERELY distorted group can
                # collide with a 4th atom so distance-perception sees a spurious
                # extra "bond" (LUXWIL C37: 3F + real O anchor + 1 close contact
                # = 2 others) -> pick the anchor as the most BOND-LIKE neighbour
                # (distance closest to the ideal single bond); the repair then
                # tetrahedralises the group and relieves the spurious contact.
                if not others:
                    continue
                ideal = 109.5
                if _worst_angle_dev(Pn, c, members, ideal) <= tol:
                    continue
                anchor = min(others, key=lambda o: abs(
                    float(_np.linalg.norm(Pn[o] - Pn[c]))
                    - (_COV.get(syms[c], 0.76) + _COV.get(syms[o], 0.76))))
                u = Pn[anchor] - Pn[c]
                nu = float(_np.linalg.norm(u))
                if nu < 1e-6:
                    continue
                u = u / nu
                # orthonormal basis perpendicular to the anchor axis
                tmp = _np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else _np.array([0.0, 1.0, 0.0])
                e1 = tmp - _np.dot(tmp, u) * u
                e1 /= float(_np.linalg.norm(e1))
                e2 = _np.cross(u, e1)
                # preserve the group's azimuthal orientation (rotamer) from member 0
                v0 = Pn[members[0]] - Pn[c]
                v0p = v0 - _np.dot(v0, u) * u
                phi0 = _np.arctan2(float(_np.dot(v0p, e2)), float(_np.dot(v0p, e1))) if float(_np.linalg.norm(v0p)) > 1e-6 else 0.0
                bl = float(_np.mean([_np.linalg.norm(Pn[m] - Pn[c]) for m in members]))
                # tetrahedral: 109.47deg from the anchor -> cos = -1/3, sin = sqrt(8)/3
                for k, m in enumerate(members):
                    phi = phi0 + 2.0 * _np.pi * k / 3.0
                    d = (-1.0 / 3.0) * u + (_np.sqrt(8.0) / 3.0) * (_np.cos(phi) * e1 + _np.sin(phi) * e2)
                    Pn[m] = Pn[c] + bl * d
                changed = True
        if not changed:
            return xyz
        return _format(syms, Pn)
    except Exception:
        return xyz
