"""PI-PLANE-FULL — the WHOLE conjugated plane, not just the ring.

WHY (measured 19.08.2026 on 930 systems from a 9600-system run):

    conjugated axis NOT planar        395 of 732    54.0 %   (654 axes)
    metal out of the pi plane > 0.2 A  89           9.6 %    (mean 0.65 A)
    sp2 donor out of the plane         47           5.1 %    (81 donors)
    biaryl OVER-planarized             20           --       (OPPOSITE DIRECTION)

Today DELFIN planarizes the *ring* (``_arom_planarize``, ``_aromatic_ring_flattener``)
and individual sp2 centers (``_fix_sp2n_planarize``, ``_fix_sp2c_planarize``).  A fully
conjugated system reaches further: the exocyclic C=O/C=N, the amide nitrogen, the
bridging C=C between two rings lie in ONE plane together with the ring, not in several
ring planes.  Exactly this plane is what this module builds.

NAMING RULE IN THIS FILE.  Everything that starts with ``_axis_type``, ``_pi_lobe_normal``,
``_fold_deg``, ``_rings_5_6_aromatic``, ``_conj_bond_cut``, ``_biaryl_dihedrals_indexed``
is a VERBATIM copy of a detector definition of the eye -- deliberately renamed so that the
name keeps visible what is a measurement term and what is our own rule.
The eye itself is only READ, never changed.

WHICH DEFINITION THIS IS BUILT AGAINST
------------------------------------
``MANTA2/weddell/detectors/conjugated_torsion.py``
  * ``_plane_normal``            :234  pi lobe = smallest-variance axis (SVD) of the
                                       bond vectors, normalized to length 1, to the
                                       HEAVY NON-METAL neighbors (>= 2 required).
  * ``_twist`` / ``_twist_h``    :439/:523  fold = angle between the two
                                       pi-lobe normals, folded into [0, 90].
  * ``classify_axis``            :382  amide / ester / amidinate / aryl_nitro /
                                       aryl_carboxyl / enone / biaryl.
  * ``_donor_axis`` :533, ``_generic_conj`` :565  arylamine / aryl_ether / "conj".
  * ``_FOLD_MIN = 12.0``         :166  below that: thermal noise, NEVER reported.
  * ``_PLANAR_DELTA_TOL``        :171  |fold_build - fold_crystal| tolerance,
                                       28-38 degrees, TWO-SIDED (over-planarized too!).
  * ``_JAG_TOL``                 :190  reference-free window; TIGHT only for
                                       amide 47.6 / ester 22.9 / aryl_nitro 47.8,
                                       WIDE (84-89) for biaryl/conj/enone/arylamine.
``MANTA2/weddell/detectors/find_ligand_quality.py``
  * ``sp2_donor_metal_oop``      :313  metal distance from the plane that the donor
                                       (N/O, 0.5 < d(M,D) < 2.75 A) spans with the FIRST TWO
                                       neighbors in the aromatic window 1.24-1.44 A.
  * ``sp2_donor_planarity``      :348  ``oop_cut = 0.40`` A -> n_sp2_donor_out_of_plane;
                                       BEST-OF over the frames.
  * ``biaryl_dihedrals``         :38   biaryl bond = C-C 1.44-1.53 A, both C with
                                       >= 2 neighbors in the aromatic window; dihedral
                                       ra[0]-a-b-rb[0], folded into [0, 90].
  * ``_OVER_TOL = 8.0``          :25   from 8 degrees difference to the crystal = over-planarized.

THE DECISIVE CONSEQUENCE OF THE TOLERANCE.  The crystal-anchored arm of the eye is
TWO-SIDED with a 28-38 degree tolerance.  Whoever blindly flattens an axis whose crystal
is twisted by 40 degrees produces exactly the defect they believe they are healing.  That is
why every hinge rotation here is CAPPED at ``MAXFOLD`` (default 25 degrees): 25 < 28 means
the correction cannot make a hitherto SILENT axis fire -- it can only heal an axis whose
delta lies between 28 and 53 degrees.  That is not a calibration, that is a bound.

WHAT IS BUILT -- three stages, by mass, each with its own fallback
----------------------------------------------------------------------
A1  SHEET PROJECTION (the 395).  A "sheet" is the maximal set of atoms connected by
    CHEMICALLY RIGID bonds: ring bonds of an aromatic ring, genuine double bonds
    (exocyclic C=O / C=N / N=O and the bridging C=C between two rings) and the tight
    conjugated single bonds amide / ester / aryl_nitro / amidinate.  All atoms of a sheet
    MUST be coplanar -- that is chemistry, not statistics.  The deviation is INTERNAL
    (pyramidalization, residual kink); no rigid-body motion can remove it, so it is
    PROJECTED.  Cheapest form under the cost law (ordering ~ 0 < isometry +0.98 pp <
    rigid rotation with new conformation +6.57 pp < re-embedding +11.9 pp).
    Capped at ``MAXPROJ`` (0.25 A): a more strongly folded sheet is no longer a
    planarity defect but a torsion -- that belongs to A2, and a projection there would
    compress the bond lengths.

A2  HINGE ROTATION (the rest of the 395).  Between two sheets sits a conjugated
    SINGLE bond.  Here the smaller side is rotated RIGIDLY about the bond axis --
    bond-length-exact, capped at ``MAXFOLD``.
    ⚠ BIARYL IS EXCLUDED.  A biaryl is twisted about its axis and must not be planar;
    ``n_biaryl_overplanar`` is the opposite direction.  In addition a hard latch keeps
    watch: if ANY biaryl dihedral (definition of the eye, ``biaryl_dihedrals``) changes
    by more than ``BIARYL_TOL`` (default 5 degrees), the move is rolled back.

B   METAL INTO THE PLANE (the 89).  ⚠ HERE THE CHEMISTRY DECIDES, NOT THE GEOMETRY
    ALONE, and a mechanism that treats both cases alike destroys the one while it
    repairs the other:
      * sigma-bound sp2 donor (pyridine N, imine N, carboxylate O): the lone pair
        points outward IN the plane -> the metal lies IN the plane.
      * eta coordination (eta5-Cp, eta6-arene): the metal sits on the axis PERPENDICULAR
        through the ring centroid -> it is precisely NOT in the plane.
    The distinction is TRIPLY safeguarded, see ``_is_sigma_donor``.

WHAT THIS MODULE DOES NOT REINVENT (rule: before every new build, check the dark switches)
------------------------------------------------------------------------------------------
Stage B overlaps with ``_pi_coplanar_final`` (DELFIN_FFFREE_PI_COPLANAR_FINAL),
``_pi_inplane_final`` (DELFIN_FFFREE_PI_RIGID_PLACE) and ``_pi_coplanar_m``
(DELFIN_FFFREE_PI_COPLANAR_M) -- all three are built, wired and SWITCHED OFF.
The eta discriminator ``_face_on`` is IMPORTED from there, not rebuilt, so that both
mechanisms draw the same boundary.  The only thing new in stage B is that the plane is the
WHOLE sheet (ring + exocyclic conjugation) instead of only the fused 5/6-rings.
That is why stage B has its OWN default OFF -- it is the second build, not the first.

MEASURED 19.08.2026 (250 champion systems from ``results/archive_builderAB_champ``,
all numbers collected with the DETECTORS OF THE EYE themselves, ``--measure``)
------------------------------------------------------------------------------------------
    changed                               94 of 250    (37.6 %)
    sheets seen 20905, projected             239       (rejected: 347 too strongly folded,
                                                        162 sp3 latch, 2196 metal contact
                                                        -- the rest is already planar)
    hinges seen 6888, rotated               1626       (1345 at the metal on both sides)
    [1] axes with fold > 12 degrees    6301 -> 5989    (-312, -5.0 %)
        systems with >= 1 such axis     154 ->  152
    [2] metal-oop best-of mean        0.116 -> 0.116   (NO effect, see below)
    [3] sp2 donors oop > 0.40 A        5307 -> 5280
    [4] biaryl moved > 8 degrees                0      (the latch held 47 times)

TWO FINDINGS THAT THE BUILD FORCED
---------------------------------------
1. THE 395 ARE A TORSION PROBLEM, NOT A SURFACE PROBLEM.  Of 20905 sheets only 347 are
   folded within themselves by more than 0.25 A; the rest are already planar.  The
   non-planar conjugated axis almost always sits in the HINGE between two planar sheets,
   not in the surface.  Whoever says "planarize better" here really means "set the torsion
   right" -- and that is a rotation, not a projection.
2. STAGE B IS A SECOND BUILD WITH ZERO EFFECT; THE FIRST BUILD IS MERELY SWITCHED OFF.
   Measured with ``--dark`` on THE SAME 250 systems, ``_pi_coplanar_final``
   (DELFIN_FFFREE_PI_COPLANAR_FINAL, default 0, NOT in the champion) gives:
       changed 75 of 250 (30.0 %)
       metal-oop best-of mean      0.116 -> 0.095 A
       systems best-of > 0.20 A       28 -> 22
       systems frame-0  > 0.20 A       90 -> 77
   Stage B here, by contrast, moves NOTHING (154 hits, best-of unchanged): ``_movable_arm``
   refuses as soon as a SECOND metal contact lies in the movable part, and precisely that
   is the normal case for a chelate (5085 ``second_metal``, 4130 ``eta_multi_contact``).
   Whoever wants the 89 systems switches on DELFIN_FFFREE_PI_COPLANAR_FINAL -- and builds nothing.

DEFAULT OFF.  ``DELFIN_FFFREE_PI_PLANE_FULL`` is 0; then every entry function is the
identity on the input text (byte-identical, not "roughly equal").  Proven on
27903 real champion frames from 995 archives: 0 deviations (``--identity``).

⚠ REACH TODAY: ZERO.  This module has NO call site.  The line that is needed belongs in
``delfin/smiles_converter.py`` (foreign territory, hence only named here): in the
outermost public chain at :32325/:32327, right next to ``_apply_pi_coplanar_final``
(whose dispatcher is at :31600) -- this chain runs through BOTH paths (the FF-free
branch leaves ``_impl`` at :32703 and arrives back there), it needs no ``mol``
and therefore does not fall into the atom-order trap that ``_ffree_shared_tail`` warns about.
What is needed is a dispatcher following the pattern of ``_apply_pi_coplanar_final``:

    def _apply_pi_plane_full(isomers):
        if not isomers or os.environ.get("DELFIN_FFFREE_PI_PLANE_FULL", "0") != "1":
            return isomers
        try:
            from delfin.manta._pi_plane_full import correct_isomer_results
            return correct_isomer_results(isomers)
        except Exception:
            return isomers

and its call in both chains, immediately around ``_apply_pi_coplanar_final(...)``.

Self-test / measurement (no inline python):
    PYTHONPATH=/home/qmchem_max/DELFIN_dev python -m delfin.manta._pi_plane_full
    PYTHONPATH=... python -m delfin.manta._pi_plane_full --measure <xyz-dir> [N]
    PYTHONPATH=... python -m delfin.manta._pi_plane_full --identity <xyz-dir> [N]
"""
from __future__ import annotations

import math
import os
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

from delfin.manta._fix_sp2n_planarize import _format_xyz, _parse_xyz, _is_metal_sym

# ---------------------------------------------------------------------------------------
# SWITCHES.  Default 0 -> identity -> byte-identical.
# ---------------------------------------------------------------------------------------
ENV_MAIN = "DELFIN_FFFREE_PI_PLANE_FULL"
ENV_STAGE_A1 = "DELFIN_FFFREE_PI_PLANE_FULL_SHEET"     # sheet projection    (default on)
ENV_STAGE_A2 = "DELFIN_FFFREE_PI_PLANE_FULL_HINGE"     # hinge rotation      (default on)
ENV_STAGE_B = "DELFIN_FFFREE_PI_PLANE_FULL_METAL"      # metal into the plane (default OFF,
#                                                        _pi_coplanar_final can already do that)
ENV_MAXPROJ = "DELFIN_FFFREE_PI_PLANE_FULL_MAXPROJ_A"  # 0.25 A
ENV_MAXFOLD = "DELFIN_FFFREE_PI_PLANE_FULL_MAXFOLD"    # 25.0 degrees (< 28 = smallest eye tolerance)
ENV_MAXROT = "DELFIN_FFFREE_PI_PLANE_FULL_MAXROT"      # 30.0 degrees
ENV_BIARYL_TOL = "DELFIN_FFFREE_PI_PLANE_FULL_BIARYL_TOL"   # 5.0 degrees

# ---------------------------------------------------------------------------------------
# CONSTANTS -- taken VERBATIM from the eye (conjugated_torsion.py :67-:97, :151,
# :166, :371; find_ligand_quality.py :22-:25, :348), so that the same definition is measured
# and repaired.  If the eye changes its numbers, these must move with them.
# ---------------------------------------------------------------------------------------
_COV_R: Dict[str, float] = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}
_DOUBLE_CUT = {
    ("C", "O"): 1.29, ("C", "N"): 1.35, ("C", "C"): 1.40, ("C", "S"): 1.68,
    ("N", "O"): 1.32, ("N", "N"): 1.32, ("C", "Se"): 1.80,
}
_CONJ_CUT = {
    ("C", "C"): 1.49, ("C", "N"): 1.45, ("C", "O"): 1.40, ("C", "S"): 1.78,
    ("N", "N"): 1.42, ("N", "O"): 1.42, ("B", "N"): 1.50,
}
_CONJ_AXIS_MAX = {("C", "N"): 1.42, ("C", "O"): 1.40, ("C", "C"): 1.49, ("N", "N"): 1.40}
_SP2_SUM_MIN = 348.0                     # conjugated_torsion.py :151
_FOLD_MIN = 12.0                         # conjugated_torsion.py :166
_AROM_MIN, _AROM_MAX = 1.24, 1.44        # find_ligand_quality.py :23
_BIARYL_MIN, _BIARYL_MAX = 1.44, 1.53    # find_ligand_quality.py :22
_SP2_OOP_CUT = 0.40                      # find_ligand_quality.py :348 (oop_cut)

# Bonds that chemically FORCE a conjugated plane (sheet formation, stage A1).
# amidinate is included: the NCN is delocalized and planar; the WIDE eye window (86.1)
# comes from the metal coordination, not from genuine rotatability.
_STIFF_AXES = frozenset(("amide", "ester", "aryl_nitro", "amidinate"))
# Axes with GENUINE rotatability -> never connect rigidly, and do not rotate biaryl either.
_NEVER_TOUCH_AXES = frozenset(("biaryl",))

# Guard thresholds
_BOND_TOL_A = 0.05          # max. change of a 1-2 bond length
_MD_TOL_A = 0.02            # max. change of a metal-donor length
_CLASH_SLACK_A = 0.02       # non-bonded minimum distance may drop by this much ...
_CLASH_SAFE_A = 2.20        # ... or it lies above this value anyway
_CHIR_KEEP = 0.25           # signed volume must not fall below this fraction


# ---------------------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------------------
def _env_flag(name: str, default: str = "0") -> bool:
    return str(os.environ.get(name, default)).strip().lower() in ("1", "true", "yes", "on")


def _env_float(name: str, default: float) -> float:
    try:
        return float(str(os.environ.get(name, "")).strip())
    except Exception:
        return float(default)


def pi_plane_full_enabled() -> bool:
    """The ONE read point of the main switch."""
    return _env_flag(ENV_MAIN, "0")


# ---------------------------------------------------------------------------------------
# Graph + geometry -- mirror of conjugated_torsion._Mol
# ---------------------------------------------------------------------------------------
def _pair_key(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def _double_bond_cut(a: str, b: str) -> Optional[float]:
    """Verbatim conjugated_torsion._double_cut (:202)."""
    k = _pair_key(a, b)
    if k in _DOUBLE_CUT:
        return _DOUBLE_CUT[k]
    ra, rb = _COV_R.get(a), _COV_R.get(b)
    return 0.90 * (ra + rb) if (ra and rb) else None


def _conj_bond_cut(a: str, b: str) -> Optional[float]:
    """Verbatim conjugated_torsion._conj_cut (:210)."""
    k = _pair_key(a, b)
    if k in _CONJ_CUT:
        return _CONJ_CUT[k]
    ra, rb = _COV_R.get(a), _COV_R.get(b)
    return 0.965 * (ra + rb) if (ra and rb) else None


def _conj_axis_len_ok(a: str, b: str, d: float) -> bool:
    """Verbatim conjugated_torsion._conj_axis_ok (:374)."""
    lim = _CONJ_AXIS_MAX.get(_pair_key(a, b))
    if lim is not None:
        return d <= lim
    cc = _conj_bond_cut(a, b)
    return cc is not None and d <= cc


def _pi_lobe_normal(P: np.ndarray, center: int,
                    neighbors: Sequence[int]) -> Optional[np.ndarray]:
    """pi lobe = smallest-variance axis of the normalized bond vectors.
    Verbatim conjugated_torsion._plane_normal (:234) -- the same quantity the eye
    measures; any other definition would repair something other than what is measured."""
    if len(neighbors) < 2:
        return None
    V = P[list(neighbors)] - P[center]
    norms = np.linalg.norm(V, axis=1)
    if np.any(norms < 1e-6):
        return None
    V = V / norms[:, None]
    try:
        _u, _s, vt = np.linalg.svd(V)
    except np.linalg.LinAlgError:
        return None
    return vt[-1]


def _angle_sum_deg(P: np.ndarray, center: int, neighbors: Sequence[int]) -> float:
    """Verbatim conjugated_torsion._sum_neighbor_angles (:256)."""
    if len(neighbors) < 3:
        return 0.0
    s = 0.0
    for a in range(len(neighbors)):
        for b in range(a + 1, len(neighbors)):
            v1 = P[neighbors[a]] - P[center]
            v2 = P[neighbors[b]] - P[center]
            nn = float(np.linalg.norm(v1) * np.linalg.norm(v2))
            if nn < 1e-9:
                continue
            s += math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2)) / nn))))
    return s


def _rings_5_6_aromatic(nbr: List[List[int]], syms: Sequence[str],
                        P: np.ndarray) -> List[FrozenSet[int]]:
    """Aromatic-like 5/6-rings -- verbatim conjugated_torsion._rings (:271)."""
    arom_el = {"C", "N", "O", "S", "Se"}
    adjA = {i: [k for k in nbr[i] if syms[k] in arom_el]
            for i in range(len(syms)) if syms[i] in arom_el}
    found: Set[FrozenSet[int]] = set()
    for start in adjA:
        stack = [(start, (start,))]
        while stack:
            node, path = stack.pop()
            for nb in adjA.get(node, ()):
                if nb == start and 5 <= len(path) <= 6:
                    found.add(frozenset(path))
                elif nb not in path and len(path) < 6:
                    stack.append((nb, path + (nb,)))
    rings: List[FrozenSet[int]] = []
    for rg in found:
        idx = list(rg)
        bl = [float(np.linalg.norm(P[a] - P[b])) for ai, a in enumerate(idx)
              for b in idx[ai + 1:] if b in nbr[a]]
        if not bl or (sum(bl) / len(bl)) > 1.46:
            continue
        Q = P[idx]
        c = Q.mean(0)
        try:
            _u, _s, vt = np.linalg.svd(Q - c)
        except np.linalg.LinAlgError:
            continue
        if float(np.sqrt(np.mean(np.dot(Q - c, vt[2]) ** 2))) > 0.30:
            continue
        rings.append(rg)
    return rings


class _Frame:
    """Graph + geometry of ONE frame.  Mirror of conjugated_torsion._Mol (:304)."""

    def __init__(self, syms: Sequence[str], P: np.ndarray):
        self.syms = list(syms)
        self.P = np.asarray(P, float)
        n = len(self.syms)
        self.n = n
        r = np.array([_COV_R.get(s, 1.5) for s in self.syms])
        D = np.sqrt(((self.P[:, None, :] - self.P[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(D, 1e9)
        self.D = D
        cut = 1.30 * (r[:, None] + r[None, :])
        nbr: List[List[int]] = [[] for _ in range(n)]
        ii, jj = np.where(D <= cut)
        for i, j in zip(ii.tolist(), jj.tolist()):
            if i < j:
                nbr[i].append(j)
                nbr[j].append(i)
        self.nbr = nbr
        self.metal = [bool(_is_metal_sym(s)) for s in self.syms]
        self.heavy = [[k for k in nbr[i] if self.syms[k] != "H"] for i in range(n)]
        self.hnm = [[k for k in self.heavy[i] if not self.metal[k]] for i in range(n)]
        self.rings = _rings_5_6_aromatic(self.nbr, self.syms, self.P)
        self.ring_of: List[List[FrozenSet[int]]] = [[] for _ in range(n)]
        for rg in self.rings:
            for a in rg:
                self.ring_of[a].append(rg)

    # -- Predicates, verbatim from the eye ----------------------------------------------
    def is_sp2(self, i: int) -> bool:
        s = self.syms[i]
        if self.metal[i] or s == "H":
            return False
        sig = self.nbr[i]
        if len(sig) == 3:
            if _angle_sum_deg(self.P, i, sig) >= _SP2_SUM_MIN:
                return True
            for k in self.heavy[i]:
                cc = _conj_bond_cut(s, self.syms[k])
                if cc is not None and float(self.D[i, k]) <= cc:
                    return True
            return bool(self.ring_of[i])
        if len(sig) == 2 and s in ("O", "N", "S"):
            for k in self.heavy[i]:
                cc = _conj_bond_cut(s, self.syms[k])
                if cc is not None and float(self.D[i, k]) <= cc:
                    return True
        return False

    def has_double(self, c: int, elem: str, exclude: Optional[int] = None) -> Optional[int]:
        for k in self.heavy[c]:
            if k == exclude or self.syms[k] != elem:
                continue
            dc = _double_bond_cut(self.syms[c], elem)
            if dc is not None and float(self.D[c, k]) <= dc:
                return k
        return None

    def count_elem(self, c: int, elem: str) -> int:
        return sum(1 for k in self.heavy[c] if self.syms[k] == elem)

    def share_ring(self, i: int, j: int) -> bool:
        return any(i in rg and j in rg for rg in self.rings)


def _axis_type(fr: _Frame, i: int, j: int) -> Optional[str]:
    """Conjugation type of the bond i-j.  Verbatim copy of
    conjugated_torsion.classify_axis (:382)."""
    si, sj = fr.syms[i], fr.syms[j]
    d_ij = float(fr.D[i, j])
    dc = _double_bond_cut(si, sj)
    if dc is not None and d_ij <= dc:
        if not (_pair_key(si, sj) == ("C", "C") and d_ij > 1.355):
            return None
    for a, b in ((i, j), (j, i)):
        sa, sb = fr.syms[a], fr.syms[b]
        if (sa == "C" and sb == "N" and fr.count_elem(b, "H") <= 2
                and _conj_axis_len_ok("C", "N", d_ij)):
            if fr.has_double(a, "O") is not None or fr.has_double(a, "S") is not None:
                return "amide"
            if fr.count_elem(a, "N") >= 2 and fr.has_double(a, "N") is not None:
                return "amidinate"
        if sa == "C" and sb == "O" and _conj_axis_len_ok("C", "O", d_ij):
            if fr.has_double(a, "O", exclude=b) is not None and len(fr.hnm[b]) >= 2:
                return "ester"
        if sa == "N" and sb == "C":
            if fr.count_elem(a, "O") >= 2 and fr.is_sp2(b):
                return "aryl_nitro"
        if sa == "C" and sb == "C":
            if fr.count_elem(a, "O") >= 2 and fr.is_sp2(b):
                return "aryl_carboxyl"
            if (fr.has_double(a, "O") is not None
                    and fr.has_double(b, "C", exclude=a) is not None):
                return "enone"
    if fr.ring_of[i] and fr.ring_of[j] and not fr.share_ring(i, j):
        return "biaryl"
    return None


def _axis_type_generic(fr: _Frame, i: int, j: int) -> Optional[str]:
    """Verbatim copy of conjugated_torsion._generic_conj (:565)."""
    si, sj = fr.syms[i], fr.syms[j]
    if fr.share_ring(i, j):
        return None
    if not (fr.is_sp2(i) and fr.is_sp2(j)):
        return None
    d = float(fr.D[i, j])
    dc = _double_bond_cut(si, sj)
    if dc is not None and d <= dc:
        if not (_pair_key(si, sj) == ("C", "C") and d > 1.355):
            return None
    cc = _conj_bond_cut(si, sj)
    if cc is None or d > cc:
        return None
    return "conj"


def _fold_deg(fr: _Frame, P: np.ndarray, i: int, j: int) -> Optional[float]:
    """Fold in degrees, [0, 90] -- verbatim conjugated_torsion._twist (:439), but on
    ONE passed-in coordinate matrix, so that before/after are comparable."""
    ni = _pi_lobe_normal(P, i, fr.hnm[i])
    nj = _pi_lobe_normal(P, j, fr.hnm[j])
    if ni is None or nj is None:
        return None
    c = max(0.0, min(1.0, abs(float(np.dot(ni, nj)))))
    return math.degrees(math.acos(c))


# ---------------------------------------------------------------------------------------
# SHEETS: the maximal set of atoms connected by chemically RIGID bonds
# ---------------------------------------------------------------------------------------
def _is_rigid_bond(fr: _Frame, i: int, j: int) -> bool:
    """Does the bond i-j force a common plane?

    Three sources, each chemical (not statistical):
      1. ring bond of an aromatic/conjugated 5/6-ring -> the ring plane;
      2. genuine DOUBLE bond between two sp2 centers -> exocyclic C=O / C=N / N=O
         and the bridging C=C between two rings (the pi bond has no rotation
         axis -- E/Z is a configuration, not a conformation);
      3. one of the TIGHT conjugated single bonds (amide/ester/aryl_nitro/
         amidinate), which the eye itself measures with a narrow window.
    Everything else -- biaryl, enone, aryl_carboxyl, arylamine, generically conjugated --
    rotates in real crystals and is NOT connected rigidly."""
    if fr.metal[i] or fr.metal[j]:
        return False
    if fr.syms[i] == "H" or fr.syms[j] == "H":
        return False
    if fr.share_ring(i, j):
        return True
    d = float(fr.D[i, j])
    dc = _double_bond_cut(fr.syms[i], fr.syms[j])
    if dc is not None and d <= dc:
        # ⚠ A TERMINAL =O IS NOT "sp2" IN THE SENSE OF THE EYE.  ``is_sp2`` (and with it
        # ``_Mol.is_sp2`` in the detector) demands THREE or TWO neighbors -- a carbonyl or
        # nitro oxygen has ONE.  That is consistent within the eye (its ``_twist`` needs
        # >= 2 neighbors anyway to span a plane), but here it would be a hole: precisely
        # this O is the third substituent of the sp2 center and by definition lies in
        # its plane.  So: the PARTNER must be sp2, the terminal atom may be.
        ok_i = fr.is_sp2(i) or len(fr.hnm[i]) <= 1
        ok_j = fr.is_sp2(j) or len(fr.hnm[j]) <= 1
        if ok_i and ok_j and (fr.is_sp2(i) or fr.is_sp2(j)):
            return True
    return _axis_type(fr, i, j) in _STIFF_AXES


def _pi_sheets(fr: _Frame, min_atoms: int = 4) -> List[List[int]]:
    """Connected components over the rigid bonds (union-find).
    Only sheets with >= min_atoms heavy atoms -- three atoms are trivially planar."""
    parent = list(range(fr.n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(fr.n):
        if fr.metal[i] or fr.syms[i] == "H":
            continue
        for j in fr.heavy[i]:
            if j <= i or fr.metal[j]:
                continue
            if _is_rigid_bond(fr, i, j):
                union(i, j)
    groups: Dict[int, List[int]] = {}
    for i in range(fr.n):
        if fr.metal[i] or fr.syms[i] == "H":
            continue
        groups.setdefault(find(i), []).append(i)
    out = [sorted(v) for v in groups.values() if len(v) >= min_atoms]
    out.sort(key=lambda g: (-len(g), g[0]))     # deterministic
    return out


def _sheet_with_h(fr: _Frame, sheet: Sequence[int]) -> List[int]:
    """Sheet + the H attached to it (an H on an sp2 center lies in its plane)."""
    s = set(sheet)
    out = list(sheet)
    for a in sheet:
        for k in fr.nbr[a]:
            if fr.syms[k] == "H" and k not in s:
                s.add(k)
                out.append(k)
    return sorted(out)


def _fit_plane(pts: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    cen = pts.mean(0)
    _u, _s, vt = np.linalg.svd(pts - cen)
    return cen, vt[2]


# ---------------------------------------------------------------------------------------
# GUARDS -- every move is rolled back if one of them trips
# ---------------------------------------------------------------------------------------
def _bond_pairs(fr: _Frame) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(fr.n) for j in fr.nbr[i] if j > i]


def _min_nonbonded(P: np.ndarray, syms: Sequence[str], nbr: List[List[int]],
                   metal: Sequence[bool]) -> float:
    """Smallest non-bonded heavy-heavy distance (1-2 and 1-3 excluded).
    Same construction as _pi_coplanar_final._global_nb_heavy_min (:229)."""
    n = len(syms)
    near = [set(nbr[i]) for i in range(n)]
    for i in range(n):
        for j in list(nbr[i]):
            near[i] |= set(nbr[j])
    heavy = [i for i in range(n) if syms[i] != "H" and not metal[i]]
    if len(heavy) < 2:
        return 1e9
    Q = P[heavy]
    D = np.sqrt(((Q[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
    mn = 1e9
    for a in range(len(heavy)):
        i = heavy[a]
        for b in range(a + 1, len(heavy)):
            if heavy[b] in near[i]:
                continue
            if D[a, b] < mn:
                mn = float(D[a, b])
    return mn


def _chirality_volumes(fr: _Frame, P: np.ndarray) -> Dict[int, float]:
    """Signed volume of every 4-coordinate center.  THIS is the check that was missing
    in the aromatic case -- a projection may neither mirror an sp3 stereocenter nor
    press it flat."""
    out: Dict[int, float] = {}
    for i in range(fr.n):
        if fr.metal[i] or fr.syms[i] == "H":
            continue
        nb = fr.nbr[i]
        if len(nb) != 4:
            continue
        v = [P[k] - P[i] for k in nb]
        out[i] = float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))
    return out


def _torsion_folded(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray,
                    p4: np.ndarray) -> float:
    """Dihedral angle, folded into [0, 90] -- verbatim find_ligand_quality._dihedral (:29)."""
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2n = b2 / (np.linalg.norm(b2) + 1e-12)
    m1 = np.cross(n1, b2n)
    d = abs(math.degrees(math.atan2(float(np.dot(m1, n2)), float(np.dot(n1, n2)))))
    return 180.0 - d if d > 90.0 else d


def _biaryl_dihedrals_indexed(fr: _Frame, P: np.ndarray) -> Dict[Tuple[int, int], float]:
    """Biaryl dihedrals ACCORDING TO THE DEFINITION OF THE EYE (find_ligand_quality.py :38),
    but index-stable: the key is the C-C pair, so that before/after are comparable.
    This is the latch against the OPPOSITE DIRECTION (n_biaryl_overplanar)."""
    S = fr.syms
    n = fr.n
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    arom = [set(np.where((D[i] >= _AROM_MIN) & (D[i] <= _AROM_MAX))[0].tolist()) - {i}
            for i in range(n)]
    out: Dict[Tuple[int, int], float] = {}
    for a in range(n):
        if S[a] != "C":
            continue
        for b in range(a + 1, n):
            if S[b] != "C" or not (_BIARYL_MIN <= D[a, b] <= _BIARYL_MAX):
                continue
            ra = sorted(arom[a] - {b})
            rb = sorted(arom[b] - {a})
            if len(ra) < 2 or len(rb) < 2:
                continue
            out[(a, b)] = _torsion_folded(P[ra[0]], P[a], P[b], P[rb[0]])
    return out


class _Guard:
    """The state BEFORE the move; ``verdict(new)`` says whether the move may stay."""

    def __init__(self, fr: _Frame):
        self.fr = fr
        self.P0 = fr.P.copy()
        self.bond_len0 = {(i, j): float(fr.D[i, j]) for i, j in _bond_pairs(fr)}
        self.md0 = {(m, k): float(fr.D[m, k])
                    for m in range(fr.n) if fr.metal[m]
                    for k in fr.nbr[m] if not fr.metal[k]}
        self.clash0 = _min_nonbonded(self.P0, fr.syms, fr.nbr, fr.metal)
        self.chir0 = _chirality_volumes(fr, self.P0)
        self.biaryl0 = _biaryl_dihedrals_indexed(fr, self.P0)
        self.biaryl_tol = _env_float(ENV_BIARYL_TOL, 5.0)

    def verdict(self, P: np.ndarray) -> Tuple[bool, str]:
        if not np.all(np.isfinite(P)):
            return False, "nonfinite"
        for (i, j), l0 in self.bond_len0.items():
            if abs(float(np.linalg.norm(P[i] - P[j])) - l0) > _BOND_TOL_A:
                return False, "bond_length"
        for (m, k), l0 in self.md0.items():
            if abs(float(np.linalg.norm(P[m] - P[k])) - l0) > _MD_TOL_A:
                return False, "metal_donor"
        cl = _min_nonbonded(P, self.fr.syms, self.fr.nbr, self.fr.metal)
        if cl < self.clash0 - _CLASH_SLACK_A and cl < _CLASH_SAFE_A:
            return False, "clash"
        ch = _chirality_volumes(self.fr, P)
        for i, v0 in self.chir0.items():
            v1 = ch.get(i)
            if v1 is None:
                return False, "chirality_lost"
            if v0 * v1 < 0.0:
                return False, "chirality_flip"
            if abs(v1) < _CHIR_KEEP * abs(v0):
                return False, "chirality_flat"
        ba = _biaryl_dihedrals_indexed(self.fr, P)
        for kk, d0 in self.biaryl0.items():
            d1 = ba.get(kk)
            if d1 is None:
                return False, "biaryl_lost"
            if abs(d1 - d0) > self.biaryl_tol:
                return False, "biaryl_moved"
        return True, "ok"


# ---------------------------------------------------------------------------------------
# STAGE A1 -- SHEET PROJECTION
# ---------------------------------------------------------------------------------------
def _stage_sheet_projection(fr: _Frame, P: np.ndarray, guard: _Guard,
                            report: Dict) -> np.ndarray:
    maxproj = _env_float(ENV_MAXPROJ, 0.25)
    for sheet in _pi_sheets(fr):
        report["n_sheets"] += 1
        pts = P[sheet]
        cen, nrm = _fit_plane(pts)
        worst = float(np.abs((pts - cen) @ nrm).max())
        if worst <= 0.05:
            continue                                    # already planar
        if worst > maxproj:
            report["n_sheet_too_folded"] += 1
            continue                                    # that is a torsion -> stage A2
        move = _sheet_with_h(fr, sheet)
        # sp3 LATCH before every projection: a 4-coordinate center is never projected.
        # ⚠ METALS DO NOT COUNT.  Measured 19.08. on 250 systems: with ``fr.nbr``
        # (which contains the metal contact) the latch tripped 976 times -- almost all
        # eta ring carbons (2 ring + 1 H + 1 metal = 4), i.e. NO sp3 centers.
        # A latch that trips on the wrong case hides the case it is meant for.
        if any(sum(1 for k in fr.nbr[a] if not fr.metal[k]) >= 4
               and fr.syms[a] != "H" and not fr.metal[a] for a in move):
            report["n_sheet_sp3_block"] += 1
            continue
        # a sheet atom bound to a metal must not move (M-D stays exact)
        if any(any(fr.metal[k] for k in fr.nbr[a]) for a in move):
            report["n_sheet_metal_block"] += 1
            continue
        cand = P.copy()
        for a in move:
            cand[a] = P[a] - float(np.dot(P[a] - cen, nrm)) * nrm
        okk, why = guard.verdict(cand)
        if not okk:
            report["rollback"][why] = report["rollback"].get(why, 0) + 1
            continue
        new_worst = float(np.abs((cand[sheet] - cen) @ nrm).max())
        if new_worst >= worst - 1e-9:
            continue
        P = cand
        report["n_sheet_fixed"] += 1
        report["sheet_gain_A"] += worst - new_worst
    return P


# ---------------------------------------------------------------------------------------
# STAGE A2 -- HINGE ROTATION
# ---------------------------------------------------------------------------------------
def _rot_matrix(axis: np.ndarray, theta_deg: float) -> np.ndarray:
    a = np.asarray(axis, float)
    a = a / (np.linalg.norm(a) + 1e-12)
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    x, y, z = a
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def _bond_sides(fr: _Frame, i: int, j: int) -> Optional[Tuple[List[int], List[int]]]:
    """Both sides of the bond i-j in the ligand graph (metals are NOT a bridge).
    None if the bond lies in a cycle -- then no rotation is possible."""
    def bfs(start: int) -> Set[int]:
        seen = {start}
        stack = [start]
        while stack:
            a = stack.pop()
            for k in fr.nbr[a]:
                if fr.metal[k]:
                    continue
                if (a == i and k == j) or (a == j and k == i):
                    continue
                if k not in seen:
                    seen.add(k)
                    stack.append(k)
        return seen

    si = bfs(i)
    if j in si:
        return None                                     # ring -> not rotatable
    return sorted(si), sorted(bfs(j))


def _hinge_bonds(fr: _Frame) -> List[Tuple[int, int, str]]:
    """Conjugated SINGLE bonds between two sheets that the eye lists as
    planarity-preferring -- WITHOUT biaryl (the opposite direction)."""
    out: List[Tuple[int, int, str]] = []
    for i in range(fr.n):
        if fr.metal[i] or fr.syms[i] == "H":
            continue
        for j in fr.heavy[i]:
            if j <= i or fr.metal[j]:
                continue
            if _is_rigid_bond(fr, i, j):
                continue                                # inside a sheet
            t = _axis_type(fr, i, j) or _axis_type_generic(fr, i, j)
            if t is None or t in _NEVER_TOUCH_AXES:
                continue
            if len(fr.hnm[i]) < 2 or len(fr.hnm[j]) < 2:
                continue                                # terminal donor: plane depends on the H
            out.append((i, j, t))
    out.sort()
    return out


def _stage_hinge(fr: _Frame, P: np.ndarray, guard: _Guard, report: Dict) -> np.ndarray:
    maxfold = _env_float(ENV_MAXFOLD, 25.0)
    for i, j, _t in _hinge_bonds(fr):
        report["n_hinges"] += 1
        fold = _fold_deg(fr, P, i, j)
        if fold is None or fold <= _FOLD_MIN:
            continue
        sides = _bond_sides(fr, i, j)
        if sides is None:
            continue
        si, sj = sides
        mi = any(any(fr.metal[k] for k in fr.nbr[a]) for a in si)
        mj = any(any(fr.metal[k] for k in fr.nbr[a]) for a in sj)
        if mi and mj:
            report["n_hinge_metal_block"] += 1
            continue                                    # both sides at the metal -> strained
        if mi:
            mov, piv_a, piv_b = sj, i, j
        elif mj:
            mov, piv_a, piv_b = si, j, i
        else:
            mov, piv_a, piv_b = (sj, i, j) if len(sj) <= len(si) else (si, j, i)
        mov = [a for a in mov if a != piv_a]
        if not mov:
            continue
        axis = P[piv_b] - P[piv_a]
        if float(np.linalg.norm(axis)) < 1e-6:
            continue
        step = min(fold, maxfold)
        org = P[piv_a]
        best = None
        for sgn in (+1.0, -1.0):
            R = _rot_matrix(axis, sgn * step)
            cand = P.copy()
            cand[mov] = (P[mov] - org) @ R.T + org
            f2 = _fold_deg(fr, cand, i, j)
            if f2 is None:
                continue
            if best is None or f2 < best[0]:
                best = (f2, cand)
        if best is None or best[0] >= fold - 1e-9:
            continue
        okk, why = guard.verdict(best[1])
        if not okk:
            report["rollback"][why] = report["rollback"].get(why, 0) + 1
            continue
        P = best[1]
        report["n_hinge_fixed"] += 1
        report["hinge_gain_deg"] += fold - best[0]
    return P


# ---------------------------------------------------------------------------------------
# STAGE B -- THE METAL INTO THE PLANE.  Here the chemistry decides.
# ---------------------------------------------------------------------------------------
try:                                    # the SAME eta discriminator as _pi_coplanar_final
    from delfin.manta._pi_coplanar_final import _face_on as _face_on
except Exception:                       # pragma: no cover
    def _face_on(M, pts):               # type: ignore[misc]
        pts = np.asarray(pts, float)
        cen = pts.mean(0)
        _u, _s, vt = np.linalg.svd(pts - cen)
        nrm = vt[2]
        in_plane = (M - cen) - np.dot(M - cen, nrm) * nrm
        radius = float(np.mean(np.linalg.norm(pts - cen, axis=1)))
        return float(np.linalg.norm(in_plane)) < 0.7 * radius


def _is_sigma_donor(fr: _Frame, P: np.ndarray, sheet: Sequence[int],
                    m: int, donors: Sequence[int]) -> Tuple[bool, str]:
    """sigma donor or eta coordination?  THREE independent criteria, ALL must speak for
    sigma -- a mechanism that confuses the two cases destroys the one while it
    repairs the other.

    (1) COUNTING (the hapticity marker, as far as it can be read off a frame at
        all -- via the contacts, not via a label).  A sigma donor is EXACTLY ONE
        contact atom of the sheet to the metal; eta2 has two, eta5-Cp five,
        eta6-arene six.
    (2) FOOT OF THE PERPENDICULAR (``_face_on``, imported VERBATIM from _pi_coplanar_final :174).
        If the projection of the metal onto the ring plane falls INSIDE the ring (closer
        than 0.7 ring radii to the centroid), the metal sits ABOVE the face -> eta.
        A sigma donor has the metal OUTSIDE the ring on the extension of the lone
        pair, i.e. about one ring radius away from the centroid.
    (3) DIRECTION OF THE LONE PAIR.  For the sigma donor D->M points OUTWARD, away from
        the ring centroid: (D->M) . (centroid->D) > 0.20.  In the eta case D->M points
        essentially along the normal; the radial component is small or negative.

    Returns (is_sigma, reason)."""
    if len(donors) != 1:
        return False, "eta_multi_contact"
    d = int(donors[0])
    ring = [a for a in sheet if any(d in rg and a in rg for rg in fr.ring_of[d])]
    if len(ring) < 4:
        ring = [a for a in sheet if fr.syms[a] != "H"]
    if len(ring) < 3:
        return False, "no_plane"
    pts = P[list(ring)]
    if _face_on(P[m], pts):
        return False, "eta_face_on"
    cen = pts.mean(0)
    u = P[m] - P[d]
    nu = float(np.linalg.norm(u))
    if nu < 1e-6:
        return False, "degenerate"
    u = u / nu
    radial = P[d] - cen
    nr = float(np.linalg.norm(radial))
    if nr < 1e-6:
        return False, "donor_is_centroid"
    if float(np.dot(u, radial / nr)) <= 0.20:
        return False, "eta_not_outward"
    return True, "sigma"


def _movable_arm(fr: _Frame, sheet: Sequence[int], donor: int) -> Optional[List[int]]:
    """The piece that may be rotated along with the sheet: the sheet plus all fragments
    hanging off it.  None if a SECOND metal contact would be included -- the rotation
    would then tear apart another M-D bond."""
    s = set(sheet)
    stack = list(sheet)
    while stack:
        a = stack.pop()
        for k in fr.nbr[a]:
            if fr.metal[k] or k in s:
                continue
            s.add(k)
            stack.append(k)
    for a in s:
        if a == donor:
            continue
        if any(fr.metal[k] for k in fr.nbr[a]):
            return None
    return sorted(s)


def _metal_oop(P: np.ndarray, m: int, d: int, fr: _Frame) -> Optional[float]:
    """Metal distance from the donor plane -- VERBATIM the quantity that
    find_ligand_quality.sp2_donor_metal_oop (:313) measures: plane from the donor and the
    FIRST TWO neighbors in the aromatic window 1.24-1.44 A, in index order."""
    D = np.sqrt(((P[d] - P) ** 2).sum(-1))
    ring = [k for k in range(fr.n)
            if k not in (m, d) and fr.syms[k] != "H" and _AROM_MIN <= D[k] <= _AROM_MAX]
    if len(ring) < 2:
        return None
    a, b = ring[0], ring[1]
    nn = np.cross(P[a] - P[d], P[b] - P[d])
    nl = float(np.linalg.norm(nn))
    if nl < 1e-6:
        return None
    return abs(float(np.dot(P[m] - P[d], nn)) / nl)


def _stage_metal_in_plane(fr: _Frame, P: np.ndarray, guard: _Guard,
                          report: Dict) -> np.ndarray:
    maxrot = _env_float(ENV_MAXROT, 30.0)
    metals = [i for i in range(fr.n) if fr.metal[i]]
    if not metals:
        return P
    for sheet in _pi_sheets(fr):
        for m in metals:
            donors = [a for a in sheet if m in fr.nbr[a]]
            if not donors:
                continue
            report["n_metal_sheets"] += 1
            sig, why = _is_sigma_donor(fr, P, sheet, m, donors)
            if not sig:
                report["eta_skip"][why] = report["eta_skip"].get(why, 0) + 1
                continue
            d = int(donors[0])
            cen, nrm = _fit_plane(P[list(sheet)])
            u = P[m] - P[d]
            nu = float(np.linalg.norm(u))
            if nu < 1e-6:
                continue
            u = u / nu
            theta = math.degrees(math.acos(max(-1.0, min(1.0, abs(float(np.dot(nrm, u)))))))
            delta = 90.0 - theta                        # this much is missing to the in-plane position
            if abs(delta) < 3.0:
                continue
            arm = _movable_arm(fr, sheet, d)
            if arm is None:
                report["eta_skip"]["second_metal"] = \
                    report["eta_skip"].get("second_metal", 0) + 1
                continue
            axis = np.cross(u, nrm)
            if float(np.linalg.norm(axis)) < 1e-6:
                continue
            mov = [a for a in arm if a != d]
            if not mov:
                continue
            step = min(abs(delta), maxrot)
            oop0 = _metal_oop(P, m, d, fr)
            if oop0 is None:
                continue
            org = P[d]
            best = None
            for sgn in (+1.0, -1.0):
                R = _rot_matrix(axis, sgn * step)
                cand = P.copy()
                cand[mov] = (P[mov] - org) @ R.T + org
                o2 = _metal_oop(cand, m, d, fr)
                if o2 is None:
                    continue
                if best is None or o2 < best[0]:
                    best = (o2, cand)
            if best is None or best[0] >= oop0 - 1e-9:
                continue
            okk, why2 = guard.verdict(best[1])
            if not okk:
                report["rollback"][why2] = report["rollback"].get(why2, 0) + 1
                continue
            P = best[1]
            report["n_metal_fixed"] += 1
            report["metal_gain_A"] += oop0 - best[0]
    return P


# ---------------------------------------------------------------------------------------
# ENTRY
# ---------------------------------------------------------------------------------------
def _new_report() -> Dict:
    return {"n_sheets": 0, "n_sheet_fixed": 0, "n_sheet_too_folded": 0,
            "n_sheet_sp3_block": 0, "n_sheet_metal_block": 0, "sheet_gain_A": 0.0,
            "n_hinges": 0, "n_hinge_fixed": 0, "n_hinge_metal_block": 0,
            "hinge_gain_deg": 0.0,
            "n_metal_sheets": 0, "n_metal_fixed": 0, "metal_gain_A": 0.0,
            "rollback": {}, "eta_skip": {}}


def correct_xyz(block: str, force: bool = False) -> Tuple[str, Dict]:
    """ONE frame.  ``(new_text, report)``.  Without the switch -> input unchanged.
    Same signature family as ``_pi_coplanar_final.correct_xyz`` / ``_pi_coplanar_m``,
    so that the wiring has the same shape as for the siblings."""
    report = _new_report()
    if not block or not (force or pi_plane_full_enabled()):
        return block, report
    try:
        syms, pts, lines = _parse_xyz(block)
    except Exception:
        return block, report
    if pts.shape[0] < 4:
        return block, report
    try:
        fr = _Frame(syms, pts)
        guard = _Guard(fr)
        P = pts.copy()
        if _env_flag(ENV_STAGE_A1, "1"):
            P = _stage_sheet_projection(fr, P, guard, report)
        if _env_flag(ENV_STAGE_A2, "1"):
            P = _stage_hinge(fr, P, guard, report)
        if _env_flag(ENV_STAGE_B, "0"):
            P = _stage_metal_in_plane(fr, P, guard, report)
    except Exception:
        return block, report
    if report["n_sheet_fixed"] + report["n_hinge_fixed"] + report["n_metal_fixed"] == 0:
        return block, report
    if not np.all(np.isfinite(P)):
        return block, report
    okk, why = guard.verdict(P)             # overall latch over ALL moves together
    if not okk:
        report["rollback"]["final_" + why] = report["rollback"].get("final_" + why, 0) + 1
        return block, report
    try:
        return _format_xyz(lines, syms, P), report
    except Exception:
        return block, report


def correct_isomer_results(isomers, force: bool = False):
    """``[(xyz, label), ...]`` -> the same list, count and order unchanged.

    Without the switch the IDENTITY (the same object), so that a call is byte-identical.
    Deliberately NOT named ``correct_results`` like the siblings (_pi_coplanar_final,
    _isolated_reseat, _aromatic_ring_flattener): the name should show, when wiring,
    that a DIFFERENT corrector hangs here and not accidentally the one of the same name.
    THIS is the function that needs a call site (see module header)."""
    if not isomers or not (force or pi_plane_full_enabled()):
        return isomers
    out = []
    changed = False
    for item in isomers:
        try:
            xyz, label = item[0], item[1]
        except Exception:
            out.append(item)
            continue
        try:
            new_xyz, _rep = correct_xyz(xyz, force=force)
        except Exception:
            new_xyz = xyz
        if new_xyz != xyz:
            changed = True
        out.append((new_xyz, label) if len(item) == 2 else item)
    return out if changed else isomers


# =======================================================================================
# SELF-TEST / MEASUREMENT -- no inline python, everything in the module
# =======================================================================================
def _ring_coords(n: int, r: float, z: float = 0.0) -> np.ndarray:
    return np.array([[r * math.cos(2 * math.pi * k / n), r * math.sin(2 * math.pi * k / n), z]
                     for k in range(n)])


def _mk_xyz(syms: Sequence[str], P: np.ndarray, comment: str = "test") -> str:
    out = [str(len(syms)), comment]
    for s, p in zip(syms, P):
        out.append(f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
    return "\n".join(out) + "\n"


def _mk_benzamide(twist_deg: float) -> Tuple[List[str], np.ndarray]:
    """Benzene + exocyclic C(=O)NH2, the amide twisted about the ring-C bond."""
    ring = _ring_coords(6, 1.39)
    syms = ["C"] * 6
    P = [ring[k] for k in range(6)]
    dirn = ring[0] / np.linalg.norm(ring[0])
    perp = np.array([-dirn[1], dirn[0], 0.0])
    R = _rot_matrix(dirn, twist_deg)
    c_carb = ring[0] + dirn * 1.48
    o = c_carb + R @ (perp * 1.22)
    nn = c_carb + R @ (-perp * 0.68 + dirn * 1.15)
    h1 = nn + R @ (-perp * 0.80 + dirn * 0.72)
    h2 = nn + R @ (perp * 0.10 + dirn * 1.00)
    syms += ["C", "O", "N", "H", "H"]
    P += [c_carb, o, nn, h1, h2]
    for k in range(1, 6):
        d = ring[k] / np.linalg.norm(ring[k])
        syms.append("H")
        P.append(ring[k] + d * 1.08)
    return syms, np.array(P)


def _mk_nitrobenzene(tilt_deg: float) -> Tuple[List[str], np.ndarray]:
    """Nitrobenzene -- the textbook case for "the plane extends beyond the ring":
    C(ring)-N is ``aryl_nitro`` (tight axis), N=O are double bonds, so ring + N + O + O
    belong in ONE sheet.  ``tilt_deg`` tilts the NO2 group about the C-N axis."""
    ring = _ring_coords(6, 1.39)
    syms = ["C"] * 6
    P = [ring[k] for k in range(6)]
    outward = ring[0] / np.linalg.norm(ring[0])
    perp = np.array([-outward[1], outward[0], 0.0])
    nn = ring[0] + outward * 1.47
    R = _rot_matrix(outward, tilt_deg)
    o1 = nn + R @ (outward * 0.61 + perp * 1.057)
    o2 = nn + R @ (outward * 0.61 - perp * 1.057)
    syms += ["N", "O", "O"]
    P += [nn, o1, o2]
    for k in range(1, 6):
        d = ring[k] / np.linalg.norm(ring[k])
        syms.append("H")
        P.append(ring[k] + d * 1.08)
    return syms, np.array(P)


def _mk_biphenyl(twist_deg: float) -> Tuple[List[str], np.ndarray]:
    A = _ring_coords(6, 1.39)
    R = _rot_matrix(np.array([1.0, 0.0, 0.0]), twist_deg)
    B = (A @ R.T) + np.array([2.87, 0.0, 0.0])   # C1-C1' = 1.48 A along x
    return ["C"] * 12, np.vstack([A, B])


def _mk_cp_metal() -> Tuple[List[str], np.ndarray]:
    A = _ring_coords(5, 1.19)
    syms = ["C"] * 5 + ["Fe"]
    P = np.vstack([A, np.array([0.0, 0.0, 1.70])])
    for k in range(5):
        d = A[k] / np.linalg.norm(A[k])
        syms.append("H")
        P = np.vstack([P, A[k] + d * 1.08])
    return syms, P


def _mk_pyridine_metal(tilt_deg: float) -> Tuple[List[str], np.ndarray]:
    """Pyridine with N at index 0, metal in the outward direction, tilted out of the
    ring plane by ``tilt`` -- that is the 89-system defect in its purest form."""
    A = _ring_coords(6, 1.39)
    syms = ["N"] + ["C"] * 5
    outward = A[0] / np.linalg.norm(A[0])
    axis = np.cross(outward, np.array([0.0, 0.0, 1.0]))
    m = A[0] + _rot_matrix(axis, tilt_deg) @ (outward * 2.05)
    P = np.vstack([A, m])
    syms.append("Ni")
    for k in range(1, 6):
        d = A[k] / np.linalg.norm(A[k])
        syms.append("H")
        P = np.vstack([P, A[k] + d * 1.08])
    return syms, P


def _selftest() -> int:
    bad = 0

    def _ck(name: str, cond: bool, extra: str = "") -> None:
        nonlocal bad
        print(f"  [{'OK ' if cond else 'FAIL'}] {name} {extra}")
        if not cond:
            bad += 1

    print("== 0. Vorgabe AUS = Identitaet ==")
    os.environ.pop(ENV_MAIN, None)
    syms, P = _mk_benzamide(55.0)
    xyz = _mk_xyz(syms, P)
    out, rep = correct_xyz(xyz)
    _ck("Schalter aus -> Text identisch", out == xyz)
    _ck("Schalter aus -> nichts angefasst",
        rep["n_sheet_fixed"] + rep["n_hinge_fixed"] + rep["n_metal_fixed"] == 0)
    probe = [(xyz, "a"), (xyz, "b")]
    _ck("correct_isomer_results aus -> dasselbe Objekt",
        correct_isomer_results(probe) is probe)

    print("== 1. Blatt reicht UEBER den Ring hinaus (Nitrobenzol) ==")
    syms, P = _mk_nitrobenzene(0.0)
    fr = _Frame(syms, P)
    sh = _pi_sheets(fr)
    big = set(max(sh, key=len)) if sh else set()
    _ck("Blatt = Ring + N + O + O", set(range(9)) <= big, f"(Blatt = {sorted(big)})")

    print("== 1b. gekippte NO2 wird in die Ringebene projiziert ==")
    syms, P = _mk_nitrobenzene(10.0)
    fr = _Frame(syms, P)
    dev0 = float(np.abs((P[:9] - P[:9].mean(0)) @ _fit_plane(P[:9])[1]).max())
    out, rep = correct_xyz(_mk_xyz(syms, P), force=True)
    _s2, P2, _l = _parse_xyz(out)
    dev1 = float(np.abs((P2[:9] - P2[:9].mean(0)) @ _fit_plane(P2[:9])[1]).max())
    _ck("Blatt-Projektion feuert", rep["n_sheet_fixed"] >= 1, f"({rep['n_sheet_fixed']})")
    _ck("Rest-Auslenkung sinkt", dev1 < dev0 - 0.02, f"({dev0:.3f} -> {dev1:.3f} A)")
    _ck("N-O-Laenge bleibt in Toleranz",
        abs(float(np.linalg.norm(P2[7] - P2[6])) - float(np.linalg.norm(P[7] - P[6])))
        <= _BOND_TOL_A)

    print("== 2. verdrehtes Amid: Faltung sinkt, Deckel haelt ==")
    syms, P = _mk_benzamide(55.0)
    fr = _Frame(syms, P)
    hs = _hinge_bonds(fr)
    f0 = max((_fold_deg(fr, P, i, j) or 0.0) for i, j, _t in hs) if hs else 0.0
    out, rep = correct_xyz(_mk_xyz(syms, P), force=True)
    _s2, P2, _l = _parse_xyz(out)
    f1 = max((_fold_deg(fr, P2, i, j) or 0.0) for i, j, _t in hs) if hs else 0.0
    _ck("Scharnier gefunden", bool(hs), f"({[(i, j, t) for i, j, t in hs]})")
    _ck("Faltung sinkt", f1 < f0 - 1.0, f"({f0:.1f} -> {f1:.1f} Grad)")
    _ck("Deckel eingehalten (<= MAXFOLD)",
        (f0 - f1) <= _env_float(ENV_MAXFOLD, 25.0) + 1e-6, f"(Delta {f0 - f1:.1f})")

    print("== 3. Biaryl bleibt verdrillt (die GEGENRICHTUNG) ==")
    syms, P = _mk_biphenyl(38.0)
    fr = _Frame(syms, P)
    b0 = _biaryl_dihedrals_indexed(fr, P)
    src = _mk_xyz(syms, P)
    out, rep = correct_xyz(src, force=True)
    _s2, P2, _l = _parse_xyz(out)
    b1 = _biaryl_dihedrals_indexed(fr, P2)
    worst = max((abs(b1.get(k, v) - v) for k, v in b0.items()), default=0.0)
    _ck("Biaryl ist kein Scharnier", all(t != "biaryl" for _i, _j, t in _hinge_bonds(fr)))
    _ck("Biaryl-Diederwinkel unveraendert", worst <= 1e-6, f"(max Delta {worst:.4f} Grad)")

    print("== 4. eta5-Cp: das Metall wird NICHT in die Ebene gezogen ==")
    os.environ[ENV_STAGE_B] = "1"
    syms, P = _mk_cp_metal()
    fr = _Frame(syms, P)
    verdicts = [_is_sigma_donor(fr, P, shh, 5, [a for a in shh if 5 in fr.nbr[a]])
                for shh in _pi_sheets(fr) if any(5 in fr.nbr[a] for a in shh)]
    _ck("Cp wird als eta erkannt", bool(verdicts) and all(not v[0] for v in verdicts),
        f"({verdicts})")
    src = _mk_xyz(syms, P)
    out, rep = correct_xyz(src, force=True)
    _ck("Cp-Frame unveraendert", out == src,
        f"(metal_fixed={rep['n_metal_fixed']}, eta_skip={rep['eta_skip']})")

    print("== 5. sigma-Pyridin: das Metall WIRD in die Ebene gezogen ==")
    syms, P = _mk_pyridine_metal(28.0)
    fr = _Frame(syms, P)
    sig = [_is_sigma_donor(fr, P, shh, 6, [a for a in shh if 6 in fr.nbr[a]])
           for shh in _pi_sheets(fr) if any(6 in fr.nbr[a] for a in shh)]
    _ck("Pyridin wird als sigma erkannt", bool(sig) and all(v[0] for v in sig), f"({sig})")
    o0 = _metal_oop(P, 6, 0, fr)
    out, rep = correct_xyz(_mk_xyz(syms, P), force=True)
    _s2, P2, _l = _parse_xyz(out)
    o1 = _metal_oop(P2, 6, 0, fr)
    _ck("Metall-oop sinkt", o0 is not None and o1 is not None and o1 < o0 - 0.05,
        f"({(o0 if o0 is not None else -1):.3f} -> {(o1 if o1 is not None else -1):.3f} A)")
    md0 = float(np.linalg.norm(P[6] - P[0]))
    md1 = float(np.linalg.norm(P2[6] - P2[0]))
    _ck("M-D-Laenge erhalten", abs(md1 - md0) <= _MD_TOL_A, f"({md0:.4f} -> {md1:.4f})")
    os.environ.pop(ENV_STAGE_B, None)

    print("== 6. sp3-Stereozentrum: der Wachhund schlaegt an ==")
    syms = ["C", "H", "F", "Cl", "Br"]
    dirs = np.array([[1.0, 1.0, 1.0], [1.0, -1.0, -1.0],
                     [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]]) / math.sqrt(3.0)
    lens = [1.09, 1.35, 1.77, 1.94]              # real bond lengths -> all four bond
    P = np.vstack([np.zeros(3)] + [dirs[k] * lens[k] for k in range(4)])
    frx = _Frame(syms, P)
    _ck("Zentrum hat vier Nachbarn", len(frx.nbr[0]) == 4, f"({len(frx.nbr[0])})")
    g = _Guard(frx)
    flip = P.copy()
    flip[:, 2] *= -1.0                                  # mirroring = sign change
    okk, why = g.verdict(flip)
    _ck("Spiegelung wird abgelehnt", (not okk) and "chirality" in why, f"({why})")

    print(f"\n{'ALLE TESTS BESTANDEN' if bad == 0 else str(bad) + ' TEST(S) FEHLGESCHLAGEN'}")
    return bad


# --- Measurement against the detectors of the eye (READ-ONLY) ---------------------------
def _load_eye_detectors():
    import sys
    p = "/home/qmchem_max/agent_workspace/MANTA2/weddell/detectors"
    if p not in sys.path:
        sys.path.insert(0, p)
    import conjugated_torsion as CT           # type: ignore
    import find_ligand_quality as FLQ         # type: ignore
    return CT, FLQ


def _split_frames(text: str) -> List[str]:
    """Multi-frame XYZ -> list of single-frame blocks."""
    lines = text.splitlines()
    out: List[str] = []
    k = 0
    while k < len(lines):
        try:
            n = int(lines[k].strip())
        except Exception:
            k += 1
            continue
        blk = lines[k:k + n + 2]
        if len(blk) == n + 2:
            out.append("\n".join(blk) + "\n")
        k += n + 2
    return out


def _atom_tuples(block: str):
    syms, P, _l = _parse_xyz(block)
    return [(syms[i], float(P[i][0]), float(P[i][1]), float(P[i][2]))
            for i in range(len(syms))]


def _measure_against_eye(dirpath: str, limit: int) -> None:
    import glob
    CT, FLQ = _load_eye_detectors()
    files = sorted(glob.glob(os.path.join(dirpath, "*.xyz")))[:limit]
    t = {"sys": 0, "axes": 0, "off0": 0, "off1": 0, "soff0": 0, "soff1": 0,
         "oop0": 0.0, "oop1": 0.0, "soop0": 0, "soop1": 0,
         "nd0": 0, "nd1": 0, "bmv": 0, "chg": 0,
         "sheets": 0, "sheet_fix": 0, "hinges": 0, "hinge_fix": 0, "metal_fix": 0,
         "sheet_folded": 0, "sheet_sp3": 0, "sheet_metal": 0, "hinge_metal": 0,
         "metal_sheets": 0}
    eta: Dict[str, int] = {}
    rb: Dict[str, int] = {}
    for fp in files:
        try:
            with open(fp) as fh:
                text = fh.read()
        except Exception:
            continue
        frames = _split_frames(text)
        if not frames:
            continue
        t["sys"] += 1
        a0 = a1 = nd0 = nd1 = bmv = 0
        oops0: List[float] = []
        oops1: List[float] = []
        chg = False
        for blk in frames:
            try:
                new, rep = correct_xyz(blk, force=True)
            except Exception:
                new, rep = blk, _new_report()
            t["sheets"] += rep["n_sheets"]
            t["sheet_fix"] += rep["n_sheet_fixed"]
            t["hinges"] += rep["n_hinges"]
            t["hinge_fix"] += rep["n_hinge_fixed"]
            t["metal_fix"] += rep["n_metal_fixed"]
            t["sheet_folded"] += rep["n_sheet_too_folded"]
            t["sheet_sp3"] += rep["n_sheet_sp3_block"]
            t["sheet_metal"] += rep["n_sheet_metal_block"]
            t["hinge_metal"] += rep["n_hinge_metal_block"]
            t["metal_sheets"] += rep["n_metal_sheets"]
            for k2, v2 in rep["rollback"].items():
                rb[k2] = rb.get(k2, 0) + v2
            for k4, v4 in rep["eta_skip"].items():
                eta[k4] = eta.get(k4, 0) + v4
            if new != blk:
                chg = True
            at0, at1 = _atom_tuples(blk), _atom_tuples(new)
            try:
                x0, _m0 = CT._axes_typed(at0)
                x1, _m1 = CT._axes_typed(at1)
            except Exception:
                x0 = x1 = []
            t["axes"] += len(x0)
            a0 += sum(1 for _i, _j, _tt, f in x0 if f > CT._FOLD_MIN)
            a1 += sum(1 for _i, _j, _tt, f in x1 if f > CT._FOLD_MIN)
            try:
                o0 = FLQ.sp2_donor_metal_oop(at0)
                o1 = FLQ.sp2_donor_metal_oop(at1)
            except Exception:
                o0 = o1 = []
            oops0.append(max((x[1] for x in o0), default=0.0))
            oops1.append(max((x[1] for x in o1), default=0.0))
            nd0 += sum(1 for x in o0 if x[1] > _SP2_OOP_CUT)
            nd1 += sum(1 for x in o1 if x[1] > _SP2_OOP_CUT)
            try:
                b0 = sorted(FLQ.biaryl_dihedrals(at0), reverse=True)
                b1 = sorted(FLQ.biaryl_dihedrals(at1), reverse=True)
                bmv += sum(1 for k3 in range(min(len(b0), len(b1)))
                           if abs(b0[k3] - b1[k3]) > FLQ._OVER_TOL)
            except Exception:
                pass
        t["off0"] += a0
        t["off1"] += a1
        t["soff0"] += 1 if a0 else 0
        t["soff1"] += 1 if a1 else 0
        m0 = min(oops0) if oops0 else 0.0          # best-of, like sp2_donor_planarity
        m1 = min(oops1) if oops1 else 0.0
        t["oop0"] += m0
        t["oop1"] += m1
        t["soop0"] += 1 if m0 > 0.20 else 0
        t["soop1"] += 1 if m1 > 0.20 else 0
        t["nd0"] += nd0
        t["nd1"] += nd1
        t["bmv"] += bmv
        t["chg"] += 1 if chg else 0
    n = max(1, t["sys"])
    print(f"\nSysteme: {t['sys']}   veraendert: {t['chg']} ({100.0 * t['chg'] / n:.1f} %)")
    print(f"Blaetter {t['sheets']} (projiziert {t['sheet_fix']}) | "
          f"Scharniere {t['hinges']} (gedreht {t['hinge_fix']}) | "
          f"Metall {t['metal_fix']}")
    print(f"  Blatt abgelehnt: zu stark gefaltet {t['sheet_folded']}, sp3-Riegel "
          f"{t['sheet_sp3']}, Metallkontakt {t['sheet_metal']}")
    print(f"  Scharnier abgelehnt: beide Seiten am Metall {t['hinge_metal']}")
    print(f"  Metall-Blaetter gesehen {t['metal_sheets']}; eta/Abbruch: {eta}")
    print(f"Rueckrollungen: {rb}")
    print(f"[1] konjugierte Achsen gesamt      : {t['axes']}")
    print(f"    Achsen mit Faltung > {CT._FOLD_MIN:.0f} Grad  : {t['off0']} -> {t['off1']}")
    print(f"    Systeme mit >= 1 solcher Achse : {t['soff0']} -> {t['soff1']}")
    print(f"[2] Metall-oop best-of Mittel (A)  : {t['oop0'] / n:.3f} -> {t['oop1'] / n:.3f}")
    print(f"    Systeme mit oop > 0,20 A       : {t['soop0']} -> {t['soop1']}")
    print(f"[3] sp2-Donoren oop > 0,40 A       : {t['nd0']} -> {t['nd1']}")
    print(f"[4] Biaryl > 8 Grad bewegt         : {t['bmv']}   (MUSS 0 sein)")


def _measure_dark_sibling(dirpath: str, limit: int) -> None:
    """BEFORE EVERY NEW BUILD, CHECK THE DARK SWITCHES.  Measures what the already built
    but SWITCHED-OFF ``_pi_coplanar_final`` (DELFIN_FFFREE_PI_COPLANAR_FINAL) moves on
    the same frames in metal-oop -- so that the report carries a NUMBER instead of a
    guess as to whether stage B is needed at all."""
    import glob
    _CT, FLQ = _load_eye_detectors()
    from delfin.manta._pi_coplanar_final import correct_xyz as _pcf
    files = sorted(glob.glob(os.path.join(dirpath, "*.xyz")))[:limit]
    n_sys = n_chg = s0 = s1 = f0_bad = f1_bad = 0
    sum0 = sum1 = 0.0
    for fp in files:
        try:
            with open(fp) as fh:
                text = fh.read()
        except Exception:
            continue
        frames = _split_frames(text)
        if not frames:
            continue
        n_sys += 1
        best0: List[float] = []
        best1: List[float] = []
        chg = False
        for k, blk in enumerate(frames):
            try:
                new = _pcf(blk)
            except Exception:
                new = blk
            if new != blk:
                chg = True
            w0 = max((x[1] for x in FLQ.sp2_donor_metal_oop(_atom_tuples(blk))), default=0.0)
            w1 = max((x[1] for x in FLQ.sp2_donor_metal_oop(_atom_tuples(new))), default=0.0)
            best0.append(w0)
            best1.append(w1)
            if k == 0:                       # the DELIVERED frame, not best-of
                f0_bad += 1 if w0 > 0.20 else 0
                f1_bad += 1 if w1 > 0.20 else 0
        m0 = min(best0) if best0 else 0.0
        m1 = min(best1) if best1 else 0.0
        sum0 += m0
        sum1 += m1
        s0 += 1 if m0 > 0.20 else 0
        s1 += 1 if m1 > 0.20 else 0
        n_chg += 1 if chg else 0
    n = max(1, n_sys)
    print(f"\n_pi_coplanar_final (DUNKLER Schalter DELFIN_FFFREE_PI_COPLANAR_FINAL)")
    print(f"Systeme: {n_sys}   veraendert: {n_chg} ({100.0 * n_chg / n:.1f} %)")
    print(f"Metall-oop best-of Mittel (A)   : {sum0 / n:.3f} -> {sum1 / n:.3f}")
    print(f"Systeme best-of oop > 0,20 A    : {s0} -> {s1}")
    print(f"Systeme Frame-0  oop > 0,20 A   : {f0_bad} -> {f1_bad}")


def _check_identity_off(dirpath: str, limit: int) -> int:
    """With the switch OFF, EVERY frame must remain byte-identical."""
    import glob
    os.environ.pop(ENV_MAIN, None)
    files = sorted(glob.glob(os.path.join(dirpath, "*.xyz")))[:limit]
    nf = bad = 0
    for fp in files:
        try:
            with open(fp) as fh:
                text = fh.read()
        except Exception:
            continue
        blocks = _split_frames(text)
        for blk in blocks:
            nf += 1
            out, _rep = correct_xyz(blk)
            if out != blk:
                bad += 1
        res = [(blk, "l") for blk in blocks]
        if correct_isomer_results(res) is not res:
            bad += 1
    print(f"Byte-Identitaet AUS: {nf} Frames aus {len(files)} Dateien, "
          f"{bad} Abweichung(en)")
    return bad


if __name__ == "__main__":
    import sys
    argv = sys.argv[1:]
    if argv and argv[0] == "--measure":
        _measure_against_eye(argv[1], int(argv[2]) if len(argv) > 2 else 100)
        raise SystemExit(0)
    if argv and argv[0] == "--dark":
        _measure_dark_sibling(argv[1], int(argv[2]) if len(argv) > 2 else 100)
        raise SystemExit(0)
    if argv and argv[0] == "--identity":
        raise SystemExit(
            1 if _check_identity_off(argv[1], int(argv[2]) if len(argv) > 2 else 100) else 0)
    raise SystemExit(1 if _selftest() else 0)
