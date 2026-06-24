"""Final per-arm π-coplanarity polish — metal INTO every σ-aromatic donor plane.

WHY (measured, pipeline-anchored 2026-06-24): a σ aromatic donor binds through an
IN-PLANE sp² lone pair, so the metal must lie IN the donor ring's mean plane.
``assemble_from_config`` places it correctly mid-pipeline (metal out-of-plane ≈ 0
for a fresh pyridine), but the downstream backbone reembed / conformer passes
(``backbone_reembed.reembed_complex``) run AFTER and re-tilt the rings, so 43-48 %
of σ-aromatic-donor complexes EMIT with the metal 0.6-2.8 Å out of plane.

The existing final pass ``_pi_inplane_final`` fixes only a MONODENTATE LIGAND (the
whole ligand is one rigid body rotated about its single donor).  That is just 8.6 %
of the defect; **91 % are MULTI-ARM polydentates** (e.g. a Zr tripodal: 3 pyridyl
arms on an sp³ backbone) where rotating the whole ligand about one donor would
break the OTHER M-D bonds.  This module closes that gap.

METHOD (universal, geometry-only, deterministic, FF-free — no SMILES/refcode):
Operate on the FINAL emitted frame (after everything → nothing downstream can undo
it).  For each conjugated π-system that has EXACTLY ONE coordinating donor D:
  * the rigid body = the π-system + its donor-less attached substituents (H, alkyl,
    …) — i.e. the arm's *exclusive* atoms.  Any fragment that leads on to ANOTHER
    coordinating donor (the shared backbone) STAYS FIXED.
  * rotate that rigid body ABOUT D (D is the pivot → M-D length is preserved
    exactly; M is fixed) so the donor's in-plane lone-pair direction (exterior
    bisector of D's two ring neighbours) aligns with the D→M vector ⇒ M lies in the
    π-system plane by construction.  The residual spin about the M-D axis keeps M
    in-plane and is used to minimise the arm→backbone linker stretch + inter-ligand
    clash.
For a true monodentate ligand the rigid body is the whole ligand (reduces to the
``_pi_inplane_final`` behaviour); for a multi-arm only that one arm moves.

SAFETY (asymmetric — never make a frame worse):  per arm, apply ONLY if the metal
out-of-plane STRICTLY decreases AND the arm→backbone linker bond(s) do not stretch
past a tolerance AND the inter-ligand heavy-heavy minimum does not drop below a
hard floor.  Otherwise the arm is left unchanged.  Never raises; returns the input
on any error.  Default-OFF / byte-identical unless ``DELFIN_FFFREE_PI_COPLANAR_FINAL=1``
(read by the caller).
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

_MD_SHELL = 2.7       # heavy donor within this of the metal = coordinated (Å)
_RING_FLAT = 0.35     # a ring's own max out-of-plane to be judged "flat" (Å)
_MOOP_MIN = 0.30      # only act if the metal is at least this far out of plane (Å)
_MOOP_EPS = 0.10      # required strict improvement in metal out-of-plane (Å)
_CLASH_FLOOR = 1.95   # heavy-heavy inter-ligand floor (Å); below = clash
_LINKER_STRETCH = 0.45  # max allowed arm→backbone linker bond growth (Å)
# deterministic spin grid about the M-D axis (keeps M in-plane; relieves linker+clash)
_SPINS = [0.0, 20.0, -20.0, 40.0, -40.0, 60.0, -60.0, 90.0, -90.0,
          120.0, -120.0, 150.0, -150.0, 180.0]


def _cov(s):
    return _COV.get(s, 0.80)


def _parse(block):
    """Return (lines, atoms) where atoms = [(line_idx, sym, xyz np.array), ...].
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
                if b not in metal_set and not seen[b]:
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


def _flat(coords, ring):
    """(centroid, unit-normal, max-out-of-plane) of a ring atom set."""
    pts = coords[ring]
    cen = pts.mean(0)
    _, _, vt = np.linalg.svd(pts - cen)
    nrm = vt[2]
    oop = float(np.abs((pts - cen) @ nrm).max())
    return cen, nrm, oop


def _pi_systems(comp_set, adj, syms, coords):
    """Conjugated π-systems in a ligand component = 5/6-membered flat rings fused
    by a shared edge (≥2 common atoms), unioned.  Returns list of atom-index sets."""
    rings = [rg for rg in _rings(comp_set, adj, syms)
             if _flat(coords, rg)[2] <= _RING_FLAT]
    if not rings:
        return []
    parent = list(range(len(rings)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        parent[find(i)] = find(j)

    for i in range(len(rings)):
        si = set(rings[i])
        for j in range(i + 1, len(rings)):
            if len(si & set(rings[j])) >= 2:   # shared edge -> fused
                union(i, j)
    groups = {}
    for i in range(len(rings)):
        groups.setdefault(find(i), set()).update(rings[i])
    return list(groups.values())


def _global_nb_heavy_min(coords, syms, metal_set):
    """Min heavy-heavy distance over all non-metal atom pairs that are NOT 1-2 or
    1-3 graph neighbours (true non-bonded clash measure).  A tight 1.15·cov bond
    graph is built from ``coords`` so genuine bonds — even unusually short build
    artefacts (nitro N-O, nitrile C-N) — are excluded, while a real through-space
    overlap is caught.  Self-contained; immune to component bookkeeping."""
    n = len(syms)
    P = np.asarray(coords, float)
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 1.15 * (_cov(syms[i]) + _cov(syms[j])):
                adj[i].add(j); adj[j].add(i)
    near = [set(adj[i]) for i in range(n)]    # 1-2 and 1-3 graph neighbours
    for i in range(n):
        for j in list(adj[i]):
            near[i] |= adj[j]
    heavy = [i for i in range(n) if syms[i] != "H" and i not in metal_set]
    mn = 1e9
    for a in range(len(heavy)):
        i = heavy[a]
        for b in range(a + 1, len(heavy)):
            j = heavy[b]
            if j in near[i]:
                continue
            if D[i, j] < mn:
                mn = float(D[i, j])
    return mn


def _arm_clash(coords, mov, syms, adj, metal_set):
    """Min heavy-heavy distance between the MOVABLE arm atoms and every other
    non-metal heavy atom they are NOT directly bonded to.  This is the exact
    contact the arm rotation could create — immune to component bookkeeping."""
    movh = [a for a in mov if syms[a] != "H"]
    oth = [a for a in range(len(syms))
           if a not in mov and a not in metal_set and syms[a] != "H"]
    if not movh or not oth:
        return 1e9
    A = coords[movh]; B = coords[oth]
    d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    movset = set(mov)
    for i, a in enumerate(movh):
        bonded = set(adj[a])
        for j, b in enumerate(oth):
            if b in bonded:           # 1-2 linker bond: not a clash
                d[i, j] = 1e9
    return float(d.min())


def _min_interlig(coords, syms, comps, metal_set):
    """Min heavy-heavy distance between atoms of DIFFERENT ligand components."""
    mind = 1e9
    heavy = [[k for k in comp if syms[k] != "H"] for comp in comps]
    for a in range(len(heavy)):
        ca = coords[heavy[a]] if heavy[a] else np.empty((0, 3))
        for b in range(a + 1, len(heavy)):
            cb = coords[heavy[b]] if heavy[b] else np.empty((0, 3))
            if len(ca) == 0 or len(cb) == 0:
                continue
            d = np.linalg.norm(ca[:, None, :] - cb[None, :, :], axis=-1)
            m = float(d.min())
            if m < mind:
                mind = m
    return mind


def _movable_set(pi, comp_set, adj, coord_donors):
    """Rigid body to rotate about the π-system's donor = the π-system + every
    attached fragment that does NOT contain another coordinating donor (the arm's
    exclusive substituents).  Fragments leading to other donors (the shared
    backbone) stay fixed.  For a monodentate ligand all fragments qualify ⇒ the
    whole ligand moves (= _pi_inplane_final behaviour)."""
    rest = [a for a in comp_set if a not in pi]
    restset = set(rest)
    seen = set()
    movable = set(pi)
    for s in rest:
        if s in seen:
            continue
        # grow the fragment (within rest)
        stack = [s]
        seen.add(s)
        frag = [s]
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b in restset and b not in seen:
                    seen.add(b)
                    stack.append(b)
                    frag.append(b)
        attaches = any(nb in pi for a in frag for nb in adj[a])
        if not attaches:
            continue                       # not bonded to this π-system: irrelevant
        if any(a in coord_donors for a in frag):
            continue                       # backbone toward another donor: stays
        movable.update(frag)               # terminal substituent: moves with the ring
    return movable


def correct_xyz(block):
    """Final per-arm metal-in-plane polish on one xyz frame.  Returns a corrected
    block string, or the input unchanged on no-op/error."""
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
            return block                   # multi-metal / no-metal: out of scope
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
        # all coordinating donors (heavy atoms inside the M-D shell), complex-wide
        coord_donors = {k for k in range(n)
                        if k != mi and syms[k] != "H" and D[mi, k] < _MD_SHELL}
        coords = P.copy()
        changed = False

        for comp in comps:
            comp_set = set(comp)
            systems = _pi_systems(comp_set, adj, syms, coords)
            for pi in systems:
                # the π-system's coordinating donors
                pdon = [k for k in pi if k in coord_donors and syms[k] in _HET]
                if len(pdon) != 1:
                    continue               # 0 = non-coordinating; ≥2 = rigid chelate (mid-pipeline)
                d = pdon[0]
                # η / face-on exclusion: no OTHER π atom may be near the metal
                if any(k != d and D[mi, k] < 2.8 for k in pi):
                    continue

                cen, nrm, _ = _flat(coords, list(pi))
                m_oop = abs(float(np.dot(M - cen, nrm)))
                if m_oop < _MOOP_MIN:
                    continue               # already in-plane

                # in-plane lone pair at d = exterior bisector of d's two π-neighbours
                ring_nbrs = [k for k in adj[d] if k in pi]
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

                movable = _movable_set(pi, comp_set, adj, coord_donors)
                # linker bonds = movable→fixed bonds (the arm's hinges to the backbone)
                linkers = []
                for a in movable:
                    for b in adj[a]:
                        if b not in movable and b not in metal_set:
                            d0 = float(D[a, b])
                            linkers.append((a, b, d0))

                R0 = _rot_align(lp, target)     # bring M into the π-system plane
                base = coords[d].copy()
                idx = sorted(movable)
                rel = coords[idx] - base
                cur_clash = _arm_clash(coords, movable, syms, adj, metal_set)
                # never push the arm's contact with the rest below the worse of
                # (current, floor): an arm move may not create or deepen a clash.
                clash_min = min(cur_clash, _CLASH_FLOOR)

                best = None
                best_key = None
                for sp in _SPINS:
                    rot = _axis_rot(target, math.radians(sp)) @ R0
                    Q = coords.copy()
                    Q[idx] = base + rel @ rot.T
                    # metal must now be ~in the π-system plane (strict improvement)
                    cen2, nrm2, ringoop = _flat(Q, list(pi))
                    if ringoop > _RING_FLAT:
                        continue
                    noop = abs(float(np.dot(M - cen2, nrm2)))
                    if noop > m_oop - _MOOP_EPS:
                        continue
                    # linker stretch guard
                    max_link = 0.0
                    for a, b, d0 in linkers:
                        dl = float(np.linalg.norm(Q[a] - Q[b])) - d0
                        if dl > max_link:
                            max_link = dl
                    if max_link > _LINKER_STRETCH:
                        continue
                    # asymmetric clash safety: never worsen the arm's contact past the floor
                    clash = _arm_clash(Q, movable, syms, adj, metal_set)
                    if clash < clash_min:
                        continue
                    # rank: lowest metal-oop, then least linker stretch, then most clash
                    key = (round(noop, 4), round(max_link, 4), -round(clash, 4))
                    if best_key is None or key < best_key:
                        best_key = key
                        best = Q
                if best is None:
                    continue
                Q = best
                if not np.all(np.isfinite(Q)):
                    continue
                # per-arm global clash gate: applying THIS arm may not create or
                # deepen any sub-floor non-bonded contact anywhere in the complex
                # (greedy arms can interact); if it would, skip this arm and keep
                # the others.  Bulletproof, order-independent never-worse.
                g_before = _global_nb_heavy_min(coords, syms, metal_set)
                g_after = _global_nb_heavy_min(Q, syms, metal_set)
                if g_after < _CLASH_FLOOR and g_after < g_before - 1e-3:
                    continue
                coords = Q
                changed = True

        if not changed:
            return block
        if not np.all(np.isfinite(coords)):
            return block
        # FINAL frame-level clash guard: greedy per-arm moves can interact
        # (arm A moved toward arm B's original pose); if the whole-frame minimum
        # non-bonded heavy contact dropped below the floor AND below the original,
        # revert the ENTIRE frame (never emit a frame worse than the input).
        g0 = _global_nb_heavy_min(P, syms, metal_set)
        g1 = _global_nb_heavy_min(coords, syms, metal_set)
        if g1 < _CLASH_FLOOR and g1 < g0 - 1e-3:
            return block
        return _emit(lines, atoms, coords, endswith_nl)
    except Exception:
        return block


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
