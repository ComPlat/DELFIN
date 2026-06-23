"""Phase 2 — per-RING coplanar metal-centered conformer search for FLEXIBLE
multi-arm polydentate ligands (the 91 % of metal-out-of-plane defects).

Why this module: a σ aromatic donor binds via an in-plane sp² lone pair, so the
metal must lie IN that donor ring's plane.  ~44 % of σ-aromatic-donor complexes
end with the metal out of plane; 91 % of those are MULTI-ARM polydentates (e.g. a
Zr tripodal with 3 pyridyl arms on an sp³ backbone).  No post-hoc rotation can fix
them (any ring re-orientation moves its donor → breaks M-D).  The existing
`_coplanar_metal_centered_conformer` only handles RIGID planar polydentates (all
donors in ONE plane).  This generalises it to "metal coplanar with EACH arm's ring"
— no common donor plane required.

Core solve (FF-free, deterministic, closed-form Gauss-Newton):
  Given, per arm i: donor position D_i, the donor's aromatic-ring plane (centroid
  c_i, unit normal n_i), and ideal bond md_i, find the metal position M minimising
      Σ_i (‖M − D_i‖ − md_i)²   +   w · Σ_i ((M − c_i)·n_i)²
  i.e. M sits at its ideal M-D distance from every donor AND in every ring's plane.
  For a rigid planar terpy (all n_i parallel) this reduces to the single-plane
  case; for a tripodal the planes differ and M is pulled to their common near-
  intersection.  The conformer with the smallest coplanarity residual wins.

This module is standalone-testable (``solve_metal`` is pure numpy) and is wired
into the polydentate placement path of ``assemble_complex`` (default-OFF, byte-id;
falls back to the current placement on any infeasibility — in-doubt-keep).
"""

import numpy as np

SEED = 42


def solve_metal(D, centroids, normals, md, has_plane=None, w=1.0, iters=200):
    """Least-squares metal position (see module docstring).

    D, centroids, normals : (k,3) arrays (normals unit).  md : (k,) ideals.
    has_plane : optional (k,) bool — include the coplanarity term only for those
    donors (e.g. aromatic σ-donors); non-aromatic donors (amine bridgehead, …)
    contribute only the M-D distance term.  Returns (M, dist_residual,
    coplanar_residual).  Deterministic.
    """
    D = np.asarray(D, float)
    C = np.asarray(centroids, float)
    N = np.asarray(normals, float)
    md = np.asarray(md, float)
    k = len(D)
    if k == 0:
        return None, 1e9, 1e9
    if has_plane is None:
        has_plane = [True] * k
    N = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-12)
    sw = float(np.sqrt(max(w, 0.0)))
    M = D.mean(axis=0) + np.array([1e-3, 1e-3, 1e-3])
    for _ in range(iters):
        rows = []
        res = []
        for i in range(k):
            diff = M - D[i]
            nn = float(np.linalg.norm(diff))
            if nn < 1e-9:
                nn = 1e-9
            res.append(nn - md[i])
            rows.append(diff / nn)
        for i in range(k):
            if not has_plane[i]:
                continue
            res.append(sw * float(np.dot(M - C[i], N[i])))
            rows.append(sw * N[i])
        J = np.array(rows)
        r = np.array(res)
        try:
            step = np.linalg.lstsq(J, -r, rcond=None)[0]
        except Exception:
            break
        M = M + step
        if float(np.linalg.norm(step)) < 1e-10:
            break
    dist_res = float(np.sqrt(np.mean([(np.linalg.norm(M - D[i]) - md[i]) ** 2
                                      for i in range(k)])))
    pl = [i for i in range(k) if has_plane[i]]
    cop_res = float(np.sqrt(np.mean([(np.dot(M - C[i], N[i])) ** 2
                                     for i in pl]))) if pl else 1e9
    return M, dist_res, cop_res


def _ring_of_donor(mol, donor, ring_info):
    """Smallest aromatic-ish ring (5/6) containing the donor atom, as an index
    list; else None."""
    cand = [list(r) for r in ring_info.AtomRings()
            if donor in r and 5 <= len(r) <= 6]
    if not cand:
        return None
    # prefer an aromatic ring; tie-break smallest
    arom = [r for r in cand
            if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)]
    pool = arom or cand
    pool.sort(key=len)
    return pool[0]


def best_conformer(mol, donor_idxs, mds, k=40, w=2.0, dist_tol=0.25):
    """Embed K free-ligand conformers, solve the metal into every arm's ring plane
    per conformer, and return the SAME contract as
    ``_coplanar_metal_centered_conformer``: ``(lsyms, [coords])`` with the
    placeholder metal EXCLUDED and coords recentred so the metal is at the origin
    (donor positions are then the M->donor vectors).  ``None`` on failure.

    Universal, geometry/graph-only.  Deterministic (fixed seed, single thread).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return None
    try:
        dons = [int(x) for x in donor_idxs]
        if len(dons) < 2 or len(mds) != len(dons):
            return None
        mh = Chem.AddHs(mol)
        params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
        params.randomSeed = SEED
        params.numThreads = 1
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=max(int(k), 1),
                                               params=params))
        if not cids:
            if AllChem.EmbedMolecule(mh, randomSeed=SEED, useRandomCoords=True) != 0:
                return None
            cids = [0]
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mh, numThreads=1)
        except Exception:
            pass
        ri = mol.GetRingInfo()
        # each donor's aromatic ring (None for non-aromatic donors, e.g. an amine
        # bridgehead — they get only the M-D distance term, no plane term).  AddHs
        # preserves heavy-atom indices, so the same indices apply to mh heavy atoms.
        rings = [_ring_of_donor(mol, d, ri) for d in dons]
        if not any(rg is not None for rg in rings):
            return None                       # no aromatic arm to coplanarise
        keep = list(range(mh.GetNumAtoms()))
        best = None  # (cop_res, coords)
        for cid in cids:
            P = np.array(mh.GetConformer(cid).GetPositions(), float)
            Dpos = np.array([P[d] for d in dons])
            C = np.zeros((len(dons), 3))
            Nn = np.tile(np.array([0.0, 0.0, 1.0]), (len(dons), 1))
            has_plane = [False] * len(dons)
            ok = True
            for i, rg in enumerate(rings):
                if rg is None:
                    continue                  # non-aromatic donor: distance only
                pts = P[rg]
                cen = pts.mean(axis=0)
                try:
                    _, sv, Vt = np.linalg.svd(pts - cen)
                except Exception:
                    ok = False
                    break
                if sv[0] < 1e-9 or (sv[2] / sv[0]) > 0.25:  # ring not flat -> skip arm
                    continue
                C[i] = cen
                Nn[i] = Vt[2]
                has_plane[i] = True
            if not ok or not any(has_plane):
                continue
            M, dres, cres = solve_metal(Dpos, C, Nn, mds, has_plane=has_plane, w=w)
            if M is None or not np.all(np.isfinite(M)):
                continue
            if dres > dist_tol:
                continue
            coords = P[keep] - M
            if not np.all(np.isfinite(coords)):
                continue
            if best is None or cres < best[0]:
                best = (cres, coords)
        if best is None:
            return None
        lsyms = [mh.GetAtomWithIdx(i).GetSymbol() for i in keep]
        return lsyms, [best[1]]
    except Exception:
        return None
