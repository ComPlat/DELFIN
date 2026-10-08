"""Monodentate and chelate assembly of the FF-free constructor: clash relief, bite-aware and trilaterated donor targets, lone-pair orientation, bite closure and multichelate assembly.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import numpy as np
import os
from typing import List, Tuple

from delfin.manta import metal_sphere_builder as MSB
from delfin.manta.assemble_donor_plane import (
    _bite_aware_targets,
)
from delfin.manta.assemble_ligand_embed import (
    _ligand_3d,
    _rot_align,
)
from delfin.manta.assemble_orient import (
    _axis_rot,
    _donor_and_lp,
)


def assemble_monodentate(metal: str, ligand_smiles: str, donor_idx: int,
                         geometry: str, thetas=None) -> Tuple[List[str], np.ndarray]:
    """Assemble; thetas[i] = rotation of ligand i about its own M-D axis (free DOF
    used for clash relief — geometric, metal-FF-free)."""
    ref = MSB._ref_vectors(geometry)
    n = len(ref)
    lsyms, lP, lmol = _ligand_3d(ligand_smiles)
    donor_elem = lsyms[donor_idx]
    md = MSB.md_distance(metal, donor_elem,
                         atom=lmol.GetAtomWithIdx(donor_idx), mol=lmol)
    if thetas is None:
        thetas = [0.0] * n
    out_syms = [metal]
    blocks = [np.zeros((1, 3))]
    lp = _donor_and_lp(lsyms, lP, lmol, donor_idx)
    for i in range(n):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        vertex = Vunit * md
        R = _rot_align(lp, -Vunit)
        Q = (lP - lP[donor_idx]) @ R.T            # donor at origin, lp along -Vunit
        Q = (Q) @ _axis_rot(-Vunit, thetas[i]).T  # spin about the M-D axis
        Q = Q + vertex
        out_syms += lsyms
        blocks.append(Q)
    return out_syms, np.vstack(blocks)


def _min_interligand(syms, P, n_lig, lig_natoms):
    """Min heavy-heavy distance between atoms of DIFFERENT ligands (clash proxy)."""
    blocks = []  # (start,end) per ligand in P (atom 0 is metal)
    s = 1
    for _ in range(n_lig):
        blocks.append((s, s + lig_natoms)); s += lig_natoms
    mind = 1e9
    for a in range(n_lig):
        for b in range(a + 1, n_lig):
            for i in range(*blocks[a]):
                if syms[i] == "H": continue
                for j in range(*blocks[b]):
                    if syms[j] == "H": continue
                    d = float(np.linalg.norm(P[i] - P[j]))
                    if d < mind: mind = d
    return mind


def clash_relief(metal, ligand_smiles, donor_idx, geometry, grid=12, passes=3):
    """Coordinate-descent over per-ligand M-D-axis rotations to maximize the min
    inter-ligand heavy distance. Deterministic (fixed grid + order)."""
    ref = MSB._ref_vectors(geometry); n = len(ref)
    lsyms, _, _ = _ligand_3d(ligand_smiles)
    lig_natoms = len(lsyms)
    angles = [2 * np.pi * k / grid for k in range(grid)]
    thetas = [0.0] * n
    for _ in range(passes):
        improved = False
        for i in range(n):
            best_t, best_score = thetas[i], -1.0
            for t in angles:
                trial = list(thetas); trial[i] = t
                syms, P = assemble_monodentate(metal, ligand_smiles, donor_idx, geometry, trial)
                sc = _min_interligand(syms, P, n, lig_natoms)
                if sc > best_score:
                    best_score, best_t = sc, t
            if abs(best_t - thetas[i]) > 1e-9:
                thetas[i] = best_t; improved = True
        if not improved:
            break
    return assemble_monodentate(metal, ligand_smiles, donor_idx, geometry, thetas), thetas


def _lone_pair_dir(P, mol, idx):
    """Direction of the donor's lone pair, from its own bonds only.

    Minus the sum of the unit vectors to ALL its neighbours -- hydrogens included, because
    for an sp3 amine the two N-H bonds are what fix where the lone pair points; dropping
    them would leave only the C-N axis and give the wrong direction entirely.

    One expression covers every case, which is the test this project applies to any rule:
    an aromatic pyridine N (two ring neighbours) gets the in-plane outward bisector, an sp3
    amine (three neighbours) gets the fourth tetrahedral direction, a linear nitrile N (one
    neighbour) gets the bond axis reversed.  No hybridisation label is read, no functional
    group is named, and it reads the same for a donor nobody has classified yet.
    """
    v = np.zeros(3)
    try:
        for nb in mol.GetAtomWithIdx(int(idx)).GetNeighbors():
            d = P[int(nb.GetIdx())] - P[int(idx)]
            n = float(np.linalg.norm(d))
            if n > 1.0e-9:
                v += d / n
    except Exception:
        return None
    n = float(np.linalg.norm(v))
    if n < 1.0e-6:
        return None                      # symmetric surroundings: no defined direction
    return -v / n


def _lp_orient_seated_bidentate(Q, mol, d1, d2, lsyms=None):
    """Set the ONE rotational freedom a SEATED bidentate still has, from the lone pairs.

    WHERE THE TILT ACTUALLY COMES FROM.  The live path seats a chelate with
    _orient_chelate_to_vertices, a Kabsch fit of the donor directions onto the target vertex
    directions.  With TWO donors that fit is UNDERDETERMINED: rotating the ligand about the
    donor-donor axis leaves both donors exactly where they are, so it changes nothing the fit
    can see, and whatever the SVD happens to return decides it.  That leftover rotation is
    precisely the one that decides whether the metal lies in a conjugated donor's pi plane --
    i.e. beta, the largest single realism gap this project has (18-19 % of frames against
    0.98 % of crystals).  It was being set by a numerical tie-break.

    The law that fixes it already existed and was already measured -- over 995 built systems
    the tilt tracks co-ligand contact (median 12.1 deg where vdW shells overlap, 0.9 deg where
    they do not) -- but it sat in _place_chelate_block, which assemble_from_config reaches ONLY
    as a fallback after the embed has already failed.  Three independent flags in that function
    (LP_ORIENT, BITE_LAW, BITE_MEASURED) all measured ZERO reach on a 187-system sweep, while
    flags elsewhere in the same file measured 123 and 7 -- the function is live code on a path
    the champion does not take.

    SAFE BY CONSTRUCTION for the seating: both donors lie ON the rotation axis, so they do not
    move at all.  The M-D distances, the bite angle and the arm-to-vertex correspondence are
    all preserved EXACTLY; only the backbone turns.  Returns None when no lone-pair direction
    is defined, or when the turn would introduce a collapsed bond the seated frame did not
    have -- the caller then keeps what it had.
    """
    Q = np.asarray(Q, float)
    axis = Q[d2] - Q[d1]
    if float(np.linalg.norm(axis)) < 1.0e-9:
        return None
    mid = 0.5 * (Q[d1] + Q[d2])
    th = _lp_aligned_angle(Q, mol, d1, d2, axis, mid)
    if th is None:
        return None
    Qr = (Q - mid) @ _axis_rot(axis, th).T + mid
    if not np.all(np.isfinite(Qr)):
        return None
    # NO STERIC ESCAPE HERE, AND THAT IS THE POINT.  A first version of this rejected the
    # derived orientation whenever it brought a backbone atom near the metal, i.e. it fell
    # back to the old sweep exactly where the sweep and the law disagree -- and measured
    # affected=0 on 187 systems, changing nothing at all.  The docstring of _lp_aligned_angle
    # records the same experiment and the same outcome one layer down: "it is why the first
    # version changed almost nothing".  Twice is enough.  The M-D sigma bond goes THROUGH the
    # lone pair, so the plane is a BONDING constraint and a clash is not a licence to break
    # it; a real molecule relieves the clash by turning substituents, opening the bite or
    # twisting the polyhedron.  If a collision remains it now STAYS VISIBLE for the self-gate
    # and the clash axis to report, instead of being hidden inside a tilt.
    return Qr


def _lp_aligned_angle(Q, mol, d1, d2, axis, mid):
    """The rotation about the donor-donor axis that points BOTH lone pairs at the metal.

    THE DEFECT THIS REPLACES.  The orientation used to be picked by a 36-step sweep that
    maximised the minimum metal-distance of the backbone -- collision avoidance, i.e. a
    STERIC heuristic where the constraint is a BONDING one.  Three consequences, and all
    three look like the flapping the manifold is full of: the criterion is simply the wrong
    one, so a planar chelate's plane tilts out; 36 steps quantise to 10 deg; and scoring by
    the MINIMUM over backbone atoms lets a single atom decide, so a small change in the
    ligand flips the chosen step and the orientation jumps between similar systems.

    A conjugated donor's lone pair lies IN its pi plane, so the metal lies in that plane
    too -- there is no such thing as a tilted aromatic chelate.  That is a hard geometric
    condition, not a preference, and it has a closed form.  Both donors sit ON the rotation
    axis, so they do not move and their target directions t_i are constant; each lone pair
    rotates by Rodrigues, and the total alignment is

        sum_i lp_i(theta) . t_i  =  A cos(theta) + B sin(theta) + C

    maximised at ``theta = atan2(B, A)``.  One atan2 instead of 36 samples: exact, no
    quantisation, and no discontinuity to jump across.

    Returns ``None`` when no lone-pair direction is defined at either donor, so the caller
    keeps the old sweep rather than inventing an orientation.
    """
    k = axis / float(np.linalg.norm(axis))
    A = B = C = 0.0
    used = 0
    for d in (d1, d2):
        v = _lone_pair_dir(Q, mol, d)
        if v is None:
            continue
        t = mid * 0.0 - Q[d]             # metal sits at the origin
        nt = float(np.linalg.norm(t))
        if nt < 1.0e-9:
            continue
        t = t / nt
        kv, kt = float(np.dot(k, v)), float(np.dot(k, t))
        A += float(np.dot(v, t)) - kv * kt
        B += float(np.dot(np.cross(k, v), t))
        C += kv * kt
        used += 1
    if used == 0 or (abs(A) < 1.0e-12 and abs(B) < 1.0e-12):
        return None
    return math.atan2(B, A)


def _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=None):
    """Rigid-fit a chelating ligand's donors onto targets T1,T2; return placed
    coords. (donor-donor vector -> target vector, midpoints matched, backbone
    rotated away from metal at origin.)

    Targets are first contracted to the chelate's natural bite (see
    _bite_aware_targets) so a tight ring is not over-stretched onto the ideal
    vertices.  The axial sweep then maximises the MINIMUM metal-distance of the
    non-donor heavy atoms (not the centroid distance): donors sit ON the rotation
    axis so only the backbone moves, and pushing its CLOSEST atom as far from the
    metal as possible keeps backbone atoms out of the coordination shell (the
    self-gate's cshm picks the cn closest atoms; a single intruder -> shape
    outlier / over-coordination).  Universal across all chelate geometries,
    deterministic."""
    T1, T2 = _bite_aware_targets(lP, d1, d2, T1, T2)
    Q = lP.copy()
    R1 = _rot_align(Q[d2] - Q[d1], T2 - T1)
    Q = (Q - Q[d1]) @ R1.T + Q[d1]
    Q = Q + (0.5 * (T1 + T2) - 0.5 * (Q[d1] + Q[d2]))
    axis = T2 - T1
    mid = 0.5 * (T1 + T2)
    body_idx = [i for i in range(len(lsyms))
                if i not in (d1, d2) and lsyms[i] != "H"]
    if not body_idx:
        return Q                              # diatomic chelate: no backbone to rotate

    # LONE-PAIR ORIENTATION (DELFIN_FFREE_LP_ORIENT, default OFF -> byte-identical).
    #
    # THE PLANE WINS, WITHOUT A STERIC ESCAPE.  A first version kept a floor that fell back
    # to the sweep when the derived orientation pushed a backbone atom near the metal.  That
    # is the wrong hierarchy and it is why the first version changed almost nothing: the
    # M-D sigma bond goes THROUGH the lone pair, so a conjugated donor holds the metal in
    # its pi plane, and a real molecule relieves a clash by turning substituents, opening
    # the bite or twisting the polyhedron -- never by tipping the metal out of the plane,
    # which would destroy the overlap that binds it.  Measured over 995 built systems: the
    # tilt tracks the co-ligand contact (Pearson -0.20; median tilt 12.1 deg where the vdW
    # shells already overlap, 0.9 deg where they do not), i.e. the seat is spending the
    # ONE degree of freedom it has on steric relief and paying for it with the bond.
    #
    # So the plane is set, unconditionally.  If a collision remains it now STAYS VISIBLE
    # for the clash axis to report, instead of being hidden inside a tilt -- and it has to
    # be resolved where reality resolves it (cis-edge choice, conformer), not here.
    if mol is not None and os.environ.get("DELFIN_FFREE_LP_ORIENT", "0") == "1":
        _th = _lp_aligned_angle(Q, mol, d1, d2, axis, mid)
        if _th is not None:
            return (Q - mid) @ _axis_rot(axis, _th).T + mid

    best, bestQ = -1e9, Q
    for k in range(36):
        Qr = (Q - mid) @ _axis_rot(axis, 2 * np.pi * k / 36).T + mid
        score = min(float(np.linalg.norm(Qr[i])) for i in body_idx)   # closest backbone atom to metal
        if score > best:
            best, bestQ = score, Qr
    return bestQ


_BITE_MAX_RING = 7      # above this a "shared ring" is a macrocycle, not a bite


def bite_close(u, v, theta_deg):
    """Rotate two unit directions symmetrically in their own plane until the angle between
    them is ``theta_deg``, keeping their bisector fixed.

    WHY.  An ideal octahedron puts two cis vertices 90 deg apart, and a five-membered
    chelate ring measured over real crystals bites at 70.6 deg (p10 69.4, p90 72.6 -- a
    3 deg wide band, one of the best defined quantities we have).  So every five-ring
    chelate is seated 19.4 deg too open, and chelates are most of our systems.  This is the
    quantified form of the earthquake root, "rigid polydentates forced onto ideal vertices",
    and of poly_cshm_vs_ccdc being the worst axis in all seven UFF-constraint A/Bs.

    Keeping the bisector fixed matters: the chelate must still occupy the SAME edge of the
    polyhedron that the isomer enumeration assigned to it.  Only the opening changes, so
    the isomer identity is untouched and this stays a pure seating correction.
    """
    u = np.asarray(u, float); v = np.asarray(v, float)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-9 or nv < 1e-9:
        return u, v
    u = u / nu; v = v / nv
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    cur = math.degrees(math.acos(c))
    if abs(cur - theta_deg) < 1e-6:
        return u, v
    bis = u + v
    nb = np.linalg.norm(bis)
    if nb < 1e-9:                       # exactly opposed: no bisector, nothing to close on
        return u, v
    bis = bis / nb
    perp = u - float(np.dot(u, bis)) * bis
    np_ = np.linalg.norm(perp)
    if np_ < 1e-9:
        return u, v
    perp = perp / np_
    half = math.radians(theta_deg) / 2.0
    return (math.cos(half) * bis + math.sin(half) * perp,
            math.cos(half) * bis - math.sin(half) * perp)


def _chelate_ring_size(lmol, d1, d2):
    """Size of the chelate ring the two donors would close WITH the metal.

    The ligand mol has no metal, so the donors are not yet in a common ring: the chelate
    ring is the shortest donor-to-donor path through the ligand, plus the metal itself.
    Path of k bonds -> ring of k+2 atoms (k+1 ligand atoms and the metal).
    """
    try:
        from rdkit import Chem
        path = Chem.GetShortestPath(lmol, int(d1), int(d2))
        return (len(path) + 1) if path else 0
    except Exception:
        return 0


def _measured_bite(metal, cn, ring_size, e1, e2):
    """Measured p50 bite angle (deg) for this chelate, or None.

    Only for SMALL rings: shared rings of 11-14 come out 68-82 deg WIDE, because two
    donors of a macrocycle can sit cis or trans while still sharing a ring.  Beyond
    _BITE_MAX_RING there is no single bite to aim at, and forcing one would be the
    bimodal-p50 mistake this table was designed to avoid.
    """
    if not ring_size or ring_size > _BITE_MAX_RING:
        return None
    path = os.environ.get("DELFIN_FFREE_DMD_BANDS", "")
    if not path:
        return None
    try:
        from delfin.manta.polyhedra import _load_band_tsv
        tbl = _load_band_tsv(path)
    except Exception:
        return None
    ea, eb = sorted((e1, e2))
    for k in (f"{metal},{cn}|{ea},{eb}|{ring_size}",
              f"{metal},{cn}|*|{ring_size}",
              f"{metal}|*|{ring_size}",
              f"*,{cn}|*|{ring_size}"):
        row = tbl.get(k)
        if row is not None:
            return float(row[1])                     # p50
    return None


def assemble_multichelate(metal: str, geometry: str, chelate_specs):
    """chelate_specs: list of (smiles, [d1,d2], [v1,v2]).  Place several chelating
    ligands on their assigned cis vertex-pairs (e.g. tris-chelate [M(en)3])."""
    ref = MSB._ref_vectors(geometry)
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    for smi, dons, verts in chelate_specs:
        lsyms, lP, lmol = _ligand_3d(smi)
        d1, d2 = dons
        u1 = ref[verts[0]] / np.linalg.norm(ref[verts[0]])
        u2 = ref[verts[1]] / np.linalg.norm(ref[verts[1]])
        # MEASURED BITE (DELFIN_FFREE_BITE_MEASURED, default OFF -> byte-identical).
        # The ideal polyhedron separates two cis vertices by 90 deg; a five-membered
        # chelate measured over real crystals bites at 70.6 deg with a 3 deg wide band.
        # Closing the pair onto the measured angle -- bisector fixed, so the ligand keeps
        # the very edge the isomer enumeration assigned it -- is a pure SETTING: no weight,
        # no gate, no minimiser.
        _r1 = MSB.md_distance(metal, lsyms[d1],
                              atom=lmol.GetAtomWithIdx(d1), mol=lmol)
        _r2 = MSB.md_distance(metal, lsyms[d2],
                              atom=lmol.GetAtomWithIdx(d2), mol=lmol)
        # DERIVED BITE (DELFIN_FFREE_BITE_LAW, default OFF -> byte-identical).
        #
        # A chelate closes a ring THROUGH the metal, so the bite is not an independent
        # quantity at all -- it is ring closure:
        #
        #     cos(bite) = (r1^2 + r2^2 - d_DD^2) / (2 r1 r2)
        #
        # with d_DD the donor-donor distance the LIGAND ALREADY HAS, read straight off its
        # own FF-free 3D construction.  Verified against 131769 crystals without building
        # anything: inverting this on the measured D-M-D table returns chemically exact
        # donor-donor distances (1.402 A for an eta2 C=C, 2.614 A for an en/bipy 5-ring)
        # and reproduces the measured bite to 0.24 deg across 18 metals.  Nine of eleven
        # ligand types sit at or below the precision of the M-D input itself.
        #
        # WHY THIS AND NOT THE TABLE.  Not precision -- the error in r enters either way.
        # COUPLING.  mdAB measured what happens when they are set independently: the M-D
        # length was corrected while the vertex angle stayed at the ideal 90 deg, which
        # forces d_DD = |u1*r1 - u2*r2| to a value the ligand DOES NOT HAVE, so the ligand
        # absorbs the difference (valid 155 -> 146).  Deriving the bite from d_DD and r
        # makes that contradiction impossible by construction: the ligand keeps its own
        # donor separation, whatever r turns out to be.
        #
        # It also extrapolates where the table has holes (bins below n=200 -- precisely the
        # exotic systems carrying our defects), separates cis from trans by itself in a
        # macrocycle (they simply ARE two different d_DD), and is a formula rather than a
        # dataset.
        if os.environ.get("DELFIN_FFREE_BITE_LAW", "0") == "1":
            _dd = float(np.linalg.norm(lP[d1] - lP[d2]))
            if _dd > 1.0e-6 and _r1 > 1.0e-6 and _r2 > 1.0e-6:
                _c = (_r1 * _r1 + _r2 * _r2 - _dd * _dd) / (2.0 * _r1 * _r2)
                # Outside [-1, 1] the triangle does not close: the ligand's own donor
                # separation cannot be reached at these M-D lengths.  Forcing a clamped
                # angle there would silently invent a geometry, so leave the vertex
                # direction alone and let the existing path handle it.
                if -1.0 <= _c <= 1.0:
                    u1, u2 = bite_close(u1, u2, math.degrees(math.acos(_c)))
        elif os.environ.get("DELFIN_FFREE_BITE_MEASURED", "0") == "1":
            # MEASURED BITE.  The ideal polyhedron separates two cis vertices by 90 deg; a
            # five-membered chelate measured over real crystals bites at 70.6 deg with a
            # 3 deg wide band.  Closing the pair onto the measured angle -- bisector fixed,
            # so the ligand keeps the very edge the isomer enumeration assigned it -- is a
            # pure SETTING: no weight, no gate, no minimiser.
            _ring = _chelate_ring_size(lmol, d1, d2)
            _bite = _measured_bite(metal, len(ref), _ring, lsyms[d1], lsyms[d2])
            if _bite is not None:
                u1, u2 = bite_close(u1, u2, _bite)
        T1 = u1 * _r1
        T2 = u2 * _r2
        Q = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=lmol)
        out_syms += lsyms; blocks.append(Q)
    return out_syms, np.vstack(blocks)


def assemble_chelate(metal: str, ligand_smiles: str, donor_indices: List[int],
                     vertex_indices: List[int], geometry: str):
    """Place a chelating ligand: rigid-fit its donor atoms onto the assigned
    polyhedron vertices (align donor-donor vector to vertex-vertex vector, match
    midpoints), then rotate about that axis so the backbone points away from the
    metal.  Bidentate; donors land near (not exactly on) vertices if the
    ligand bite != vertex spacing (a constrained relax closes that)."""
    ref = MSB._ref_vectors(geometry)
    lsyms, lP, lmol = _ligand_3d(ligand_smiles)
    d1, d2 = donor_indices[0], donor_indices[1]
    md1 = MSB.md_distance(metal, lsyms[d1], atom=lmol.GetAtomWithIdx(d1), mol=lmol)
    md2 = MSB.md_distance(metal, lsyms[d2], atom=lmol.GetAtomWithIdx(d2), mol=lmol)
    T1 = ref[vertex_indices[0]] / np.linalg.norm(ref[vertex_indices[0]]) * md1
    T2 = ref[vertex_indices[1]] / np.linalg.norm(ref[vertex_indices[1]]) * md2
    # rigid-fit + backbone-away-from-metal axial sweep (single source of truth)
    bestQ = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=lmol)
    syms = [metal] + lsyms
    P = np.vstack([np.zeros((1, 3)), bestQ])
    md_act = (float(np.linalg.norm(bestQ[d1])), float(np.linalg.norm(bestQ[d2])))
    return syms, P, md_act
