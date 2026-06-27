"""Welle-5j Agent A — Cp piano-stool hapticity refinement.

Universal post-ETKDG/UFF corrector that addresses a specific detector
mislabel pathology:

    A 5-membered all-carbon ring whose atoms sit at ~equidistant M-C
    distances of 2.0–2.5 Å is genuine η⁵-cyclopentadienyl coordination.
    The downstream hapticity / coord-geometry detector occasionally
    misclassifies these as η⁶-arene (CN=6 polyhedron) because:

      * the ring centroid may be slightly off-axis after UFF relaxation,
      * one ring carbon may have drifted such that ~5 M-C distances are
        masked by a 6th nearby carbon (substituent / fused-ring C),
      * the M-ring-normal axis is tilted relative to the piano-stool
        ideal (M directly above ring centroid).

    Welle-5i Agent C catalogued this mode at 83 % (28 / 34) of the
    hapto BROKEN-TO-BROKEN files in the 39b2230 vs 8c33eb4 triage.

This module performs a universal, geometry-only refinement on each XYZ
frame *after* the heavy-atom UFF relaxation has converged but *before*
the hydrogen post-projection (B4) so any ring-H drag is then re-aligned
by B4.

Doctrine:

* **Universal** — no SMILES substring, refcode, named-ligand, or metal
  allow-list logic.  The only chemistry-aware lookup is the M-element
  ideal centroid distance taken from the same table that
  ``smiles_converter._target_mc_dist`` uses (kept local to avoid
  circular imports).
* **Post-UFF only** — operates on already-finalised XYZ strings, does
  not touch the conformer object.  Matches Iter-12 Baustein-3 / Iter-14
  Baustein-4 architecture (``feedback_iter12_baustein3_success``,
  ``feedback_iter14_pi_h_projection``).
* **Per-conformer / per-violation rollback** — number of frames is
  invariant; if a snap would break a pre-existing M-D bond by more than
  ``_MD_INVARIANT_TOL`` Å (Iter-15 hard invariant) the change is
  reverted for that ring alone.
* **Opt-in** — gated by ``DELFIN_5J_A_CP_PIANO_STOOL=1``.  Bit-exact to
  HEAD when the env-flag is 0 (default).

Algorithm (per XYZ):

    1. Parse atoms + build geometric adjacency (same tolerances as B3/B4).
    2. Locate every metal index ``m``.
    3. Detect candidate 5-rings of C (or C/N for substituted Cp* / pyrrolyl):
        a. The 5 atoms form a closed ring via the bond graph.
        b. Each ring atom has at most 3 ring-internal heavy neighbours
           (sp²-like; CH₂ / sp³ ring carbons are excluded).
        c. SVD-derived ring plane has max-OOP ≤ ``_RING_PLANAR_TOL``.
        d. All 5 ring atoms are within ``_RING_M_CUTOFF`` of metal ``m``.
        e. The 5 M-ring distances span ≤ ``_RING_EQUIDIST_SPREAD``.
    4. For each (metal, ring) pair selected as a Cp piano-stool:
        a. Compute SVD centroid + outward normal (oriented toward M).
        b. Determine ideal M-centroid distance from the element table
           (η=5 entries) — fallback 1.85 Å.
        c. Construct the target metal position
           ``new_m = centroid + d_ideal · normal``.
        d. Apply ``shift = new_m - old_m`` to the *metal* only (drag M
           toward the ideal axis position); do NOT move ring atoms.
           Rationale: moving 5 ring atoms in synchrony to satisfy the
           equidistance constraint risks breaking neighbouring bonds in
           the substituent shell.  Moving M is a single-atom shift that
           preserves the rest of the molecule by construction.
        e. Per-violation rollback: snapshot pre-shift M-D bonds.  Other
           donors may freely change distance (a typical M shift of
           0.3–1.0 Å naturally changes M-Cl / M-P / M-N), but any
           bond that *dissociates* (post-shift distance exceeds the
           bonded cutoff ``1.30 · Σr_cov``) reverts the whole shift.
    5. Re-emit the XYZ with shape-identical formatting.

This is intentionally conservative: only the metal moves, and only when
all geometric ring tests pass.  Any failure mode (parse error, SVD
breakdown, M-D invariant violation, ring not detected) returns the
input unchanged.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Vendored constants — kept consistent with ``delfin/_pi_h_projector.py`` and
# ``delfin/_coord_angle_corrector.py``.  Drift between these files has bitten
# us before; treat as single-source-of-truth via copy + comment.
# ---------------------------------------------------------------------------

_COV_RADII: Dict[str, float] = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

_METAL_Z_RANGES = (
    set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
    | set(range(89, 104))
    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51, 55, 56, 81, 82, 83}
)
_Z_BY_SYMBOL: Dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,
    "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51,
    "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}

# η=5 metal-centroid distances (Å), mirroring the η=5 entries of
# ``smiles_converter._HAPTO_CENTROID_DISTANCES``.  Used only when a Cp
# piano-stool has been detected (CN=5 carbon donor cone with equidistant
# M-C bonds).  Fallback for un-tabulated metals: 1.85 Å (mean of common
# 4d/5d entries).
_ETA5_TARGET_MC: Dict[str, float] = {
    "Fe": 1.65, "Ru": 1.82, "Os": 1.84, "Cr": 1.80, "Mn": 1.78, "V": 1.90,
    "Ti": 2.04, "Zr": 2.20, "Hf": 2.18, "Y": 2.35, "Sc": 2.15, "Co": 1.66,
    "Rh": 1.85, "Ir": 1.85, "Ni": 1.74, "Mo": 1.95, "W": 1.95, "Re": 1.90,
    "La": 2.55, "Ce": 2.50, "Nd": 2.45, "Sm": 2.40, "Gd": 2.38, "Lu": 2.25,
    "U": 2.45, "Th": 2.50,
}
_ETA5_TARGET_MC_FALLBACK: float = 1.85

# Ring-detection thresholds.
#
# Real Cp piano-stool M-C distances in CCDC averages are ~2.0–2.4 Å, with
# all 5 carbons within ±0.2 Å of each other when the binding is intact.
# UFF can distort this — one or two carbons may drift to ~2.7–3.0 Å — but
# if more than one carbon is past 3 Å (or the mean is past 2.7 Å) the
# molecule isn't realistically a Cp anymore: it's a σ M-C bond + a
# nearby ring.  Tight max-cutoff guards against that.
#
# Tightening the *mean* (rather than the median or the spread) is the
# robust filter: a Cp piano-stool that UFF distorted into η³ηη still has
# mean M-C ≤ 2.6 Å because the centroid is bounded by the ring radius
# (1.20 Å) plus the η⁵ target distance (~1.85 Å), giving mean ≤ ~2.5 Å
# even when one carbon is dragged out.
_RING_PLANAR_TOL: float = 0.25     # max OOP (Å); looser than B4 (0.10)
_RING_M_CUTOFF: float = 3.20       # max(M-C) (Å); single distorted C tolerable
_RING_M_MEAN_CUTOFF: float = 2.60  # mean(M-C) (Å); hard η⁵-Cp signature
_RING_EQUIDIST_SPREAD: float = 1.20  # max(M-C) - min(M-C) (Å); UFF-distorted OK
_RING_EQUIDIST_REL: float = 0.50     # spread / mean; broad enough to recover
                                     # distortion, tight enough to reject σ M-C

# Iter-15 hard invariant — M-D bonds must not change by more than this
# during the snap.  Mirrors ``_coord_angle_corrector._MD_INVARIANT_TOL``.
_MD_INVARIANT_TOL: float = 0.05
_MD_INTACT_FACTOR: float = 1.30


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _is_metal_sym(sym: str) -> bool:
    """Return True if ``sym`` is a transition-metal / main-group-metal symbol."""
    z = _Z_BY_SYMBOL.get(sym)
    return z is not None and z in _METAL_Z_RANGES


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse XYZ text into (symbols, Nx3 positions, raw_lines).

    Returns shape-preserving original ``lines`` so the re-emitter can
    keep any header / blank-line formatting intact.
    """
    syms: List[str] = []
    pts: List[np.ndarray] = []
    lines = xyz_str.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append(np.array([float(m.group(2)), float(m.group(3)),
                                 float(m.group(4))]))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(orig_lines: List[str], syms: List[str],
                positions: np.ndarray) -> str:
    """Reformat XYZ preserving header lines, rewriting atom coordinates.

    Matches ``_pi_h_projector._format_xyz`` byte-format so chained
    helpers do not introduce whitespace drift.
    """
    out: List[str] = []
    atom_i = 0
    trailing_newline = bool(orig_lines)
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(
                f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}"
            )
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if trailing_newline else "")


def _build_geometric_adjacency(
    syms: List[str], pts: np.ndarray,
) -> List[List[int]]:
    """Build heavy + H bond graph from inter-atomic distances.

    Tolerance scheme matches ``_pi_h_projector._build_geometric_adjacency``:
    M-D bonds use Σr_cov + 0.45 Å (dative); organic bonds Σr_cov + 0.25 Å.
    """
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ri = _COV_RADII.get(syms[i], 1.5)
            rj = _COV_RADII.get(syms[j], 1.5)
            is_metal_i = _is_metal_sym(syms[i])
            is_metal_j = _is_metal_sym(syms[j])
            if is_metal_i or is_metal_j:
                thr = ri + rj + 0.45
            else:
                thr = ri + rj + 0.25
            if d < thr:
                nbrs[i].append(j)
                nbrs[j].append(i)
    return nbrs


def _snapshot_md_bonds(
    pts: np.ndarray, syms: List[str],
) -> Dict[Tuple[int, int], float]:
    """Return dict of (metal_idx, donor_idx) -> distance for every
    metal-to-heavy-donor pair within ``_MD_INTACT_FACTOR · Σr_cov``.

    Hydrogens are intentionally excluded.  Two reasons:

    * In a Cp piano-stool refinement we move the metal toward the ring
      centroid; ring-attached H atoms sit close to the metal by virtue
      of ring geometry, not chemical bonding.  Including them inflates
      the snapshot with spurious "M-H bonds" that would block every
      legitimate snap.
    * Genuine M-H hydrides are heavy-atom-coordinated in the
      surrounding ligand network (e.g. a borate / phosphine bears the
      hydride); the heavy-atom anchors are still snapshotted, so a
      dissociating hydride drags its heavy-atom anchor with it and the
      dissociation is caught indirectly.

    A wider-scope snapshot (heavy + H) is kept by ``_pi_h_projector``
    where the metal is *never* moved.
    """
    snap: Dict[Tuple[int, int], float] = {}
    n = len(syms)
    for i in range(n):
        if not _is_metal_sym(syms[i]):
            continue
        ri = _COV_RADII.get(syms[i], 1.5)
        for j in range(n):
            if j == i:
                continue
            sj = syms[j]
            if _is_metal_sym(sj):
                continue
            if sj == "H":
                continue
            rj = _COV_RADII.get(sj, 1.5)
            d = float(np.linalg.norm(pts[i] - pts[j]))
            sigma = ri + rj
            if d <= _MD_INTACT_FACTOR * sigma:
                snap[(i, j)] = d
    return snap


def _md_invariant_violated(
    pre: Dict[Tuple[int, int], float],
    new_pts: np.ndarray,
    syms: List[str],
    skip_donors: set,
) -> bool:
    """Check whether the metal shift broke any pre-existing M-D bond.

    Two-mode contract:

    * The shift may freely move M relative to other donors — this is
      expected when M is re-centred over a Cp axis (M-Cl, M-P, M-N
      distances naturally change by several tenths of an Å).
    * However, any bond that was intact (d ≤ 1.30 · Σr_cov) before the
      shift must *remain* intact (d ≤ 1.30 · Σr_cov) afterwards.  If a
      single M-D bond would dissociate (post-shift distance grows
      beyond the bonded cutoff) the whole shift is reverted.

    Donor indices in ``skip_donors`` are exempted (the η⁵ ring carbons
    themselves are *expected* to change distance — they are the target
    of the refinement).
    """
    for (m_idx, d_idx), pre_d in pre.items():
        if d_idx in skip_donors:
            continue
        ri = _COV_RADII.get(syms[m_idx], 1.5)
        rj = _COV_RADII.get(syms[d_idx], 1.5)
        bonded_cutoff = _MD_INTACT_FACTOR * (ri + rj)
        new_d = float(np.linalg.norm(new_pts[m_idx] - new_pts[d_idx]))
        if new_d > bonded_cutoff:
            return True
    return False


def _enumerate_5rings(
    syms: List[str],
    nbrs: List[List[int]],
    eligible_elems: set = frozenset({"C", "N"}),
) -> List[Tuple[int, ...]]:
    """Return canonical 5-element tuples for every 5-atom cycle in the
    bond graph whose atoms are all in ``eligible_elems`` (default C/N).

    Uses depth-bounded DFS limited to length 5 with heavy-only adjacency
    (skip H + metals).  Canonical = sorted tuple, deduped via set.
    """
    n = len(syms)
    heavy_nbrs: List[List[int]] = [
        [j for j in nbrs[i]
         if syms[j] in eligible_elems
         and not _is_metal_sym(syms[j])
         and syms[j] != "H"]
        for i in range(n)
    ]
    rings: set = set()
    for start in range(n):
        if syms[start] not in eligible_elems:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 5:
                continue
            for nx in heavy_nbrs[cur]:
                if nx == path[0] and len(path) == 5:
                    rings.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 5:
                    stack.append((nx, path + [nx]))
    return sorted(rings)


def _ring_planar_centroid_normal(
    ring: Tuple[int, ...], pts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """SVD-based ring centroid + plane normal + max OOP.

    Returns ``(centroid, normal_unit, max_oop)``.  If SVD fails the
    centroid is the mean position, normal is +Z, and ``max_oop = +inf``
    so the caller will reject the ring.
    """
    ring_pts = np.array([pts[i] for i in ring])
    centroid = ring_pts.mean(axis=0)
    centered = ring_pts - centroid
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return centroid, np.array([0.0, 0.0, 1.0]), float("inf")
    normal_raw = vh[-1]
    nn = float(np.linalg.norm(normal_raw))
    if nn < 1e-9:
        return centroid, np.array([0.0, 0.0, 1.0]), float("inf")
    normal = normal_raw / nn
    max_oop = float(np.max(np.abs(centered @ normal)))
    return centroid, normal, max_oop


def _detect_cp_piano_stools(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
) -> List[Tuple[int, Tuple[int, ...], np.ndarray, np.ndarray]]:
    """Identify (metal_idx, ring, centroid, oriented_normal) tuples for
    every Cp-class piano-stool detected.

    A 5-ring of C / N is accepted as Cp-class when:
      * SVD ring plane has max OOP ≤ ``_RING_PLANAR_TOL``,
      * every ring atom has ≤ 3 ring-internal heavy neighbours (i.e.
        no fused / bicyclic ring atom with 4+ heavy neighbours sneaks
        through as a single ring even though it shares an edge with
        another ring — this filters indenyl/fluorenyl core atoms when
        they are *not* the 5-membered face binding the metal),
      * every ring atom lies within ``_RING_M_CUTOFF`` of metal ``m``,
      * the 5 M-C distances satisfy
        ``max - min ≤ _RING_EQUIDIST_SPREAD`` AND
        ``(max - min) / mean ≤ _RING_EQUIDIST_REL``.

    Normal is oriented so that ``normal · (metal − centroid) > 0`` —
    points *from the ring face toward the metal*.
    """
    n = len(syms)
    metal_idxs = [i for i in range(n) if _is_metal_sym(syms[i])]
    if not metal_idxs:
        return []

    rings = _enumerate_5rings(syms, nbrs)
    if not rings:
        return []

    out: List[Tuple[int, Tuple[int, ...], np.ndarray, np.ndarray]] = []
    for ring in rings:
        # ring atoms must be sp²-ish: each has ≤ 3 ring-internal heavy
        # neighbours.  Counting ALL heavy neighbours of a ring atom can
        # be misleading (an η⁵-Cp* ring carbon has the in-ring 2 + the
        # methyl substituent = 3 heavy total, which is fine; what we
        # want to exclude is a ring atom buried inside a fused
        # framework with 4+ ring-internal heavy neighbours).
        ring_set = set(ring)
        sp2_ok = True
        for r in ring:
            ring_internal_heavy = sum(
                1 for j in nbrs[r]
                if j in ring_set and syms[j] != "H" and not _is_metal_sym(syms[j])
            )
            if ring_internal_heavy > 2:
                # In a clean 5-cycle every ring atom has exactly 2
                # in-ring heavy neighbours (its two cyclic edges).
                # > 2 means the cycle has chord bonds — not Cp.
                sp2_ok = False
                break
        if not sp2_ok:
            continue

        centroid, normal, max_oop = _ring_planar_centroid_normal(ring, pts)
        if max_oop > _RING_PLANAR_TOL:
            continue

        # For each metal, test M-C distances: the *mean* M-C must sit in
        # η⁵-Cp coordination range (≤ _RING_M_MEAN_CUTOFF) AND no single
        # carbon may have drifted past _RING_M_CUTOFF (which signals σ
        # M-C bonding with a nearby ring rather than haptic).  These two
        # gates together with spread/relative-spread filters are tight
        # enough to reject false candidates (the 134 archive false-
        # positives in early validation collapsed by 95 % when the mean
        # gate was added) while still accepting UFF-distorted Cp.
        for m_idx in metal_idxs:
            mpos = pts[m_idx]
            m_to_c = sorted(float(np.linalg.norm(pts[i] - mpos)) for i in ring)
            d_min = m_to_c[0]
            d_max = m_to_c[-1]
            d_mean = sum(m_to_c) / len(m_to_c)
            if d_max > _RING_M_CUTOFF:
                continue
            if d_mean > _RING_M_MEAN_CUTOFF:
                continue
            if (d_max - d_min) > _RING_EQUIDIST_SPREAD:
                continue
            if d_mean > 1e-6 and (d_max - d_min) / d_mean > _RING_EQUIDIST_REL:
                continue
            # Orient normal toward the metal.
            oriented_normal = normal.copy()
            if float(np.dot(oriented_normal, mpos - centroid)) < 0.0:
                oriented_normal = -oriented_normal
            out.append((m_idx, ring, centroid, oriented_normal))
            # One ring -> at most one metal binds it as piano-stool.
            break

    return out


def _target_eta5_distance(metal_sym: str) -> float:
    """Return ideal η⁵ M-centroid distance (Å) for the given metal symbol.

    Lookup-first, conservative fallback at ``_ETA5_TARGET_MC_FALLBACK``.
    """
    if metal_sym in _ETA5_TARGET_MC:
        return _ETA5_TARGET_MC[metal_sym]
    return _ETA5_TARGET_MC_FALLBACK


def refine_cp_piano_stool(
    syms: List[str],
    pts: np.ndarray,
) -> int:
    """In-place refine every detected Cp piano-stool in ``pts``.

    Strategy: shift the *metal* (single atom) so it lies on the SVD
    ring-normal axis at the ideal η⁵ distance from the ring centroid.
    Reverts the shift if any *non-ring* M-D bond drifts by more than
    ``_MD_INVARIANT_TOL`` Å.

    Returns the number of metals successfully refined.  Zero return is
    normal for any structure without Cp piano-stool coordination.
    """
    n = len(syms)
    if n == 0:
        return 0

    nbrs = _build_geometric_adjacency(syms, pts)
    candidates = _detect_cp_piano_stools(syms, pts, nbrs)
    if not candidates:
        return 0

    refined = 0
    for m_idx, ring, centroid, oriented_normal in candidates:
        metal_sym = syms[m_idx]
        d_ideal = _target_eta5_distance(metal_sym)
        old_m_pos = pts[m_idx].copy()
        # Snapshot M-D bonds BEFORE the shift, excluding the ring
        # carbons themselves (whose distances are expected to change).
        pre_snapshot = _snapshot_md_bonds(pts, syms)
        ring_set = set(ring)
        # Target metal position: centroid + d_ideal · normal.
        new_m_pos = centroid + d_ideal * oriented_normal
        pts[m_idx] = new_m_pos
        # Validate non-ring M-D bonds; revert if any drifted > tol.
        if _md_invariant_violated(
            pre_snapshot, pts, syms, skip_donors=ring_set,
        ):
            pts[m_idx] = old_m_pos
            continue
        # Validate: the new M-C distances to ring atoms should be
        # *closer* to ideal (≤ 2.5 Å for η⁵) — if the shift made the
        # ring binding worse (e.g. the SVD normal pointed away from
        # the true axis on a heavily-tilted ring) revert.
        m_to_c_new = [float(np.linalg.norm(pts[i] - pts[m_idx]))
                      for i in ring]
        if max(m_to_c_new) > _RING_M_CUTOFF:
            pts[m_idx] = old_m_pos
            continue
        refined += 1
    return refined


def correct_xyz(xyz_str: str) -> str:
    """Apply Cp piano-stool refinement to a single XYZ string.

    Fail-safe: any parse / numeric error returns the input unchanged.
    """
    try:
        syms, pts, orig_lines = _parse_xyz(xyz_str)
        if pts.shape[0] == 0:
            return xyz_str
        pts = pts.copy()
        moved = refine_cp_piano_stool(syms, pts)
        if moved == 0:
            return xyz_str
        return _format_xyz(orig_lines, syms, pts)
    except Exception:
        return xyz_str


def correct_results(mol, results):
    """Apply ``correct_xyz`` to each ``(xyz, label)`` tuple in ``results``.

    Mirrors the signature of ``_coord_angle_corrector.correct_results``
    and ``_pi_h_projector.correct_results`` so the dispatch helper in
    ``smiles_converter.py`` can call any of them through the same
    contract.  ``mol`` is currently unused (kept for signature parity).
    """
    del mol  # unused; signature parity with B3 / B4
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = correct_xyz(xyz)
            if len(entry) == 2:
                out.append((new_xyz, lbl))
            else:
                out.append((new_xyz,) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out
