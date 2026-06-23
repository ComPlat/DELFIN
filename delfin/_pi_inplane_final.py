"""Final post-conformer MONODENTATE aromatic σ-donor in-plane re-assertion.

WHY (measured 2026-06-23): a σ aromatic donor binds through an IN-PLANE sp² lone
pair, so the metal must lie IN the donor ring's mean plane.  ``assemble_complex``
places it correctly (metal out-of-plane = 0.00 Å for a fresh pyridine), but the
conformer / reembed passes (``_conf``) run AFTER the mid-pipeline coplanar-M pass
and re-tilt the ring, so 43-48 % of σ-aromatic-donor complexes end with the metal
0.6-2.8 Å OUT of the donor plane.  97.7 % of these are MONODENTATE (one donor in
the conjugated system) — exactly the case a per-donor rotation fixes with nothing
to break (no second M-D bond).  The existing ``_pi_coplanar_m`` does not perceive/
correct them; this is a self-contained, geometry-only re-assertion run on the
FINAL ensemble.

METHOD (universal, deterministic, FF-free, no SMILES/refcode knowledge):
For each ligand whose ONLY metal-coordinating atom is a single heteroatom D that
sits in a flat 5/6-membered aromatic ring (η / face-on π donors excluded), rotate
the WHOLE rigid ligand about the fixed donor D so the donor's in-plane lone-pair
(exterior bisector of D's two ring neighbours) aligns with the D→M direction.  M
then lies in the ring plane by construction; M-D length is preserved exactly
(rotation is about D, M is fixed).  The remaining spin about the M-D axis is free
and is used to keep the inter-ligand clash never-worse.

SAFETY (asymmetric — never make a frame worse): apply only if the metal
out-of-plane STRICTLY decreases AND the minimum inter-ligand heavy-heavy distance
does not drop below a hard floor; otherwise keep the frame unchanged.  Never
raises; returns the input on any error.  Default-OFF / byte-identical unless
``DELFIN_FFFREE_PI_RIGID_PLACE=1`` (gated by the caller).
"""

import math
import numpy as np

_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U",
}
_COV = {"C": 0.76, "N": 0.71, "O": 0.66, "H": 0.31, "P": 1.07, "S": 1.05,
        "Cl": 1.02, "F": 0.57, "B": 0.84, "Se": 1.20, "Br": 1.20, "Si": 1.11,
        "As": 1.19, "Te": 1.38, "I": 1.39}
_HET = {"N", "O", "S", "P", "Se", "As", "Te"}

_MD_SHELL = 2.7      # heavy donor within this of the metal = coordinated
_RING_FLAT = 0.35    # ring's own max out-of-plane to be judged "flat"
_MOOP_MIN = 0.30     # only act if metal is at least this far out of plane
_CLASH_FLOOR = 1.95  # heavy-heavy inter-ligand floor (Å); below = clash
_SPINS = [0.0, 30.0, -30.0, 60.0, -60.0, 90.0, -90.0, 120.0, -120.0, 150.0,
          -150.0, 180.0]  # deterministic spin grid about M-D for clash relief


def _cov(s):
    return _COV.get(s, 0.80)


def _parse(block):
    """Return (lines, atoms) where atoms = list of (line_idx, sym, xyz np.array).
    Non-atom lines (count/comment) are kept by index for verbatim re-emit."""
    lines = block.splitlines()
    atoms = []
    for li, ln in enumerate(lines):
        p = ln.split()
        if len(p) >= 4:
            try:
                x, y, z = float(p[1]), float(p[2]), float(p[3])
            except ValueError:
                continue
            if math.isfinite(x) and math.isfinite(y) and math.isfinite(z):
                atoms.append((li, p[0], np.array([x, y, z])))
    return lines, atoms


def _emit(lines, atoms, coords, endswith_nl):
    """Rebuild the block: atom lines get new coords (6dp), others verbatim."""
    out = list(lines)
    for (li, sym, _), c in zip(atoms, coords):
        out[li] = f"{sym} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}"
    return "\n".join(out) + ("\n" if endswith_nl else "")


def _rot_align(a, b):
    """Rotation matrix taking unit vector a onto unit vector b (Rodrigues)."""
    a = a / (np.linalg.norm(a) + 1e-12)
    b = b / (np.linalg.norm(b) + 1e-12)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        # antiparallel: 180° about any axis ⟂ a
        ax = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            ax = np.array([0.0, 1.0, 0.0])
        ax = ax - np.dot(ax, a) * a
        ax /= (np.linalg.norm(ax) + 1e-12)
        return _axis_rot(ax, math.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def _axis_rot(axis, theta):
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    x, y, z = axis
    c, s = math.cos(theta), math.sin(theta)
    C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def _components(n, adj, metal_set):
    """Connected components of the non-metal atoms (ligands)."""
    seen = [False] * n
    comps = []
    for s in range(n):
        if s in metal_set or seen[s]:
            continue
        stack = [s]
        seen[s] = True
        comp = []
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in adj[a]:
                if not seen[b]:
                    seen[b] = True
                    stack.append(b)
        comps.append(comp)
    return comps


def _rings(comp_set, adj, syms):
    """5/6-membered rings fully inside comp_set (heavy only)."""
    R = []
    seen = set()

    def dfs(s, cur, sn):
        last = cur[-1]
        for nb in adj[last]:
            if nb not in comp_set or syms[nb] == "H":
                continue
            if nb == s and 5 <= len(cur) <= 6:
                k = frozenset(cur)
                if k not in seen:
                    seen.add(k)
                    R.append(list(cur))
            elif nb not in sn and len(cur) < 6:
                dfs(s, cur + [nb], sn | {nb})

    for s in comp_set:
        if syms[s] != "H":
            dfs(s, [s], {s})
    return R


def correct_xyz(block):
    """Re-assert monodentate aromatic σ-donor metal-in-plane on one xyz frame.
    Returns a corrected block string, or the input unchanged on no-op/error."""
    endswith_nl = block.endswith("\n")
    try:
        lines, atoms = _parse(block)
        n = len(atoms)
        if n < 5:
            return block
        syms = [a[1] for a in atoms]
        P = np.array([a[2] for a in atoms])
        metal_idx = [i for i in range(n) if syms[i] in _METALS]
        if len(metal_idx) != 1:
            return block               # multi-metal/no-metal: out of scope
        mi = metal_idx[0]
        metal_set = {mi}
        M = P[mi]
        D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)

        # non-metal heavy+H bond graph
        adj = [[] for _ in range(n)]
        for i in range(n):
            if i == mi:
                continue
            for j in range(i + 1, n):
                if j == mi:
                    continue
                if D[i, j] < 1.25 * (_cov(syms[i]) + _cov(syms[j])):
                    adj[i].append(j)
                    adj[j].append(i)

        comps = _components(n, adj, metal_set)
        coords = P.copy()
        changed = False

        for comp in comps:
            comp_set = set(comp)
            # the ligand's metal-coordinating heavy atoms
            donors = [k for k in comp if syms[k] != "H" and D[mi, k] < _MD_SHELL]
            if len(donors) != 1:
                continue                # only clean MONODENTATE here
            d = donors[0]
            if syms[d] not in _HET:
                continue
            # donor must sit in a flat aromatic ring; pick the flattest containing d
            rings = [rg for rg in _rings(comp_set, adj, syms) if d in rg]
            best_ring = None
            best_flat = 1e9
            for rg in rings:
                pts = coords[rg]
                c = pts.mean(0)
                _, _, vt = np.linalg.svd(pts - c)
                flat = float(np.abs((pts - c) @ vt[2]).max())
                if flat < best_flat:
                    best_flat = flat
                    best_ring = rg
            if best_ring is None or best_flat > _RING_FLAT:
                continue
            # η / face-on exclusion: no OTHER ring atom may be near the metal
            if any(k != d and D[mi, k] < 2.8 for k in best_ring):
                continue

            # current metal out-of-plane
            pts = coords[best_ring]
            cen = pts.mean(0)
            _, _, vt = np.linalg.svd(pts - cen)
            nrm = vt[2]
            m_oop = abs(float(np.dot(M - cen, nrm)))
            if m_oop < _MOOP_MIN:
                continue                # already in-plane

            # in-plane lone pair at d = exterior bisector (away from ring nbrs)
            ring_nbrs = [k for k in best_ring if k in adj[d]]
            if len(ring_nbrs) < 2:
                continue
            lp = coords[d] - np.mean([coords[k] for k in ring_nbrs], axis=0)
            lpn = np.linalg.norm(lp)
            if lpn < 1e-6:
                continue
            lp = lp / lpn
            target = M - coords[d]
            tn = np.linalg.norm(target)
            if tn < 1e-6:
                continue
            target = target / tn

            # rigid rotation about d aligning lp -> target (M into ring plane)
            R0 = _rot_align(lp, target)
            atoms_to_move = sorted(comp_set)
            base = coords[d].copy()

            def apply(rot):
                Q = coords.copy()
                for k in atoms_to_move:
                    Q[k] = base + rot @ (coords[k] - base)
                return Q

            # candidate orientations: align, then spin about M-D for clash relief
            mdaxis = target
            best = None
            best_clash = -1.0
            for sp in _SPINS:
                rot = _axis_rot(mdaxis, math.radians(sp)) @ R0
                Q = apply(rot)
                # verify metal now ~in-plane
                rp = Q[best_ring]
                rc = rp.mean(0)
                _, _, rvt = np.linalg.svd(rp - rc)
                noop = abs(float(np.dot(M - rc, rvt[2])))
                if noop > m_oop - 0.1:      # must strictly improve
                    continue
                clash = _min_interlig(Q, syms, comps, metal_set)
                if clash > best_clash:
                    best_clash = clash
                    best = Q
            if best is None:
                continue
            # asymmetric safety: do not introduce a clash below the hard floor
            cur_clash = _min_interlig(coords, syms, comps, metal_set)
            if best_clash < _CLASH_FLOOR and best_clash < cur_clash - 1e-6:
                continue
            coords = best
            changed = True

        if not changed:
            return block
        if not np.all(np.isfinite(coords)):
            return block
        return _emit(lines, atoms, coords, endswith_nl)
    except Exception:
        return block


def _min_interlig(coords, syms, comps, metal_set):
    """Min heavy-heavy distance between atoms of DIFFERENT ligand components."""
    mind = 1e9
    heavy = [[k for k in comp if syms[k] != "H"] for comp in comps]
    for a in range(len(heavy)):
        for b in range(a + 1, len(heavy)):
            for i in heavy[a]:
                ci = coords[i]
                for j in heavy[b]:
                    dd = float(np.linalg.norm(ci - coords[j]))
                    if dd < mind:
                        mind = dd
    return mind


def correct_results(isomers):
    """Map :func:`correct_xyz` over an emitted isomer ensemble.  Each item is an
    xyz string or a tuple whose first element is the xyz string."""
    if not isomers:
        return isomers
    out = []
    for item in isomers:
        try:
            if isinstance(item, (tuple, list)) and item:
                new = correct_xyz(str(item[0]))
                out.append((new,) + tuple(item[1:]))
            else:
                out.append(correct_xyz(str(item)))
        except Exception:
            out.append(item)
    return out
