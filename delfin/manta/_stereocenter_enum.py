"""Stereocentre Enumeration — coordination-created X-H R/S fold completeness.

THE GAP (measured by the eye, find_isomer_coverage `_nh_sign` + the `ccdc_isomer_realized`
hard floor with ``DELFIN_EYE_NH_STEREO=1``): a donor that becomes a tetrahedral stereocentre
ONLY once coordinated — a secondary amine ``[NH+]`` bonded to the metal + two backbone carbons —
has an up/down N-H fold (R/S) that the coordination fingerprint alone is blind to.  ETKDG samples
folds only by accident (USEMOW: the crystal's alternating ``N1+N1-N1+N1-`` fold is NEVER embedded ->
``ccdc_isomer_realized = FALSE`` even though the OC-6 polyhedron + coordination isomers are perfect).

COMPLETENESS LAW (User 2026-07-13): "we always need all that are possible.  not always all work, but
always + and - must be built for completeness."  So every stereocentre must have BOTH folds ATTEMPTED;
the multidentate rigidity (a small macrocycle cannot point every N-H the same way) prunes the
infeasible ones by a hard steric clash — that is physics deciding, not us pre-selecting.

THE FIX (this module, env-gated ``DELFIN_STEREOCENTER_ENUM``, default-OFF):
  * X-H fold: reflect the donor H across the (M, donor, backbone-centroid) plane.  That plane contains
    M and D, so the reflection is an isometry that EXACTLY negates the eye's signed triple product
    ``dot(M-D, cross(H-D, B-D))`` -> flips the R/S sign, while preserving |H-M| and angle(M-D-H)
    (so the amine-H realism gate cannot false-reject it).
  * per-fold light constrained UFF relax (reuse ``assemble_complex._constrained_uff_relax``) with EVERY
    heavy atom frozen — only the flipped H's move — so no aromatic/backbone distortion, just local
    N-H settling.  Fall back to the pure reflection if the relax would un-flip the sign.
  * clash-prune: a fold whose flipped H is driven into a hard steric overlap (the macrocycle forbids it)
    is dropped — "not all work".

ARCHITECTURE: ADDITIVE (originals are never modified/removed -> never-worse-safe by construction) and
DETERMINISTIC (no RNG; reflection + fixed-conformer UFF minimize are byte-stable).  Dispatched AFTER the
final dedup (next to the rotamer/conformer expansions) so the folds — which are heavy-atom-close to their
base — are never collapsed.  Bit-exact no-op when the env-flag is unset OR the structure has no
coordination-created X-H stereocentre.
"""
from __future__ import annotations

import itertools
import logging
import os
from typing import List, Optional, Tuple

import numpy as np

# Single source of truth for the geometry primitives (no drift, no runtime dep on the
# validation framework): the same helpers Baustein 3 uses.
from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

_LOG = logging.getLogger(__name__)

_SIGN_EPS = 1e-6

# Elements that form a CONFIGURATIONALLY STABLE coordination-created stereocentre: the group-15
# pnictogens.  Coordination donates their lone pair to the metal, which LOCKS the pyramidal
# configuration (an uncoordinated secondary amine inverts freely; the metal-bound one does not).
# Oxygen / chalcogen "X-H stereocentres" invert essentially barrierlessly (not stable), and a
# coordinated C-H is not this axis at all -- flipping those does NOT create a real stereoisomer, it
# only perturbs the eye's manifold-derived spec and spuriously breaks the CCDC-isomer floor (measured:
# WIPHOW [OOOO] + JAJMUG [CC] went ccdc true->false).  Whole-group scope -> universal, not element-
# specific.  (Chalcogen / other centres are a separate axis needing their own eye + validation.)
_STEREO_DONOR_ELEMENTS = frozenset({"N", "P", "As", "Sb", "Bi"})

# ===== CHALCOGENS: THE SECOND HALF OF THE SAME AXIS (20.08.2026) =====
# The exclusion above says "Oxygen / chalcogen X-H stereocentres invert essentially
# barrierlessly" and backs it with TWO cases: WIPHOW [OOOO] and JAJMUG [CC].  Those are
# OXYGEN and CARBON.  **Sulfur and selenium were never tested** -- they were swept along
# under "chalcogen" by analogy, and the analogy is chemically weak: the inversion barrier
# at a COORDINATED thioether lies far above the one at an ether, which is why S-configured
# complexes are isolable as separate diastereomers.  The line above admits it itself:
# "Chalcogen / other centres are a separate axis needing their own eye + validation".
#
# 🔑 THE EYE ALREADY DEMANDS THEM.  find_isomer_coverage._stereo_kind is by default
# ELEMENT-FREE (:759) -- it counts every donor as a centre and demands both folds.
# It says there: "on 38 of the 83 systems that carry a sign it demands a fold ... the
# builder deliberately refuses to make.  Two instruments, never held against each other."
# Measured consequence: systems with such a centre realise the natural isomer in 52 % of
# cases against 90 % otherwise.  This axis closes the half of the dispute for which there
# is NO counterexample -- O and C stay excluded.
#
# DENOMINATOR, measured before building (harness/prochiral_zensus.py on pool_full_10000):
#   272 of 10000 systems = 2.72 % carry >=1 distinguishable chalcogen centre
#   457 centres, of which 410 = 89.7 % have at most ONE arm in the chelate ring
#   (for comparison: the pnictogens already covered are only 34.9 % that favourable)
# Calibration of the census: pnictogens 3.42 % against the known census of 3.03 %.
_STEREO_CHALCOGEN_ELEMENTS = frozenset({"S", "Se"})


# ---------------------------------------------------------------------------
# Env gate + tunables (all default to the calibrated USEMOW-validated values)
# ---------------------------------------------------------------------------

def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    s = str(raw).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return 1
    if s in ("0", "false", "no", "off", ""):
        return 0
    try:
        return int(s)
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return float(raw)
    except Exception:
        return default


def _is_enabled() -> bool:
    """Master switch.  Accepts BOTH ``DELFIN_STEREOCENTER_ENUM`` (the validated A/B name) and the
    champion-list form ``DELFIN_FFFREE_STEREOCENTER_ENUM`` (``_CHAMPION_FLAGS`` prefixes its entries
    with ``DELFIN_FFFREE_``).  Same code path either way; default OFF (byte-exact) when neither is set."""
    return (_env_int("DELFIN_STEREOCENTER_ENUM", 0) == 1
            or _env_int("DELFIN_FFFREE_STEREOCENTER_ENUM", 0) == 1)


# ---------------------------------------------------------------------------
# Stereocentre sign — MUST match the eye (find_isomer_coverage._nh_sign)
# ---------------------------------------------------------------------------

def _subtree(nbrs: List[List[int]], start: int, blocked: set) -> set:
    """BFS component reachable from ``start`` without crossing any atom in ``blocked`` (inclusive
    of ``start``).  Used to collect a substituent's whole atom set for an X-R inversion."""
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for y in nbrs[x]:
            if y in blocked or y in seen:
                continue
            seen.add(y)
            stack.append(y)
    return seen


def _reaches_metal(nbrs: List[List[int]], start: int, block_d: int, metal_set: set) -> bool:
    """True if a BFS from ``start`` (blocking only the donor ``block_d``) reaches a metal — i.e. the
    substituent is a chelate-ring arm, so the X-R centre cannot be inverted by a clean reflection."""
    seen = {start, block_d}
    stack = [start]
    while stack:
        x = stack.pop()
        for y in nbrs[x]:
            if y in seen:
                continue
            if y in metal_set:
                return True
            seen.add(y)
            stack.append(y)
    return False


def _center_sign(pts: np.ndarray, c: dict) -> str:
    """Signed handedness of a coordination-created tetrahedral stereocentre — the scalar triple
    product of the metal anchor and the two highest-information substituent directions.  IDENTICAL
    formula to the eye (find_isomer_coverage._nh_sign), for BOTH kinds:
      * X-H: directions = the single H and the non-metal heavy-neighbour centroid.
      * X-R: directions = the two highest-priority of the 3 distinguishable heavy substituents.
    So a built fold maps to the eye's stereocentre fingerprint."""
    D = pts[c["d"]]
    M = pts[c["m"]]
    if c["type"] == "XH":
        a = pts[c["h"]] - D
        b = np.mean(pts[c["heavy"]], axis=0) - D
    else:                                              # XR
        a = pts[c["keyed"][0]] - D
        b = pts[c["keyed"][1]] - D
    s = float(np.dot(M - D, np.cross(a, b)))
    if s > _SIGN_EPS:
        return "+"
    if s < -_SIGN_EPS:
        return "-"
    return "0"


def _center_plane_normal(pts: np.ndarray, c: dict):
    """Unit normal of the reflection plane that INVERTS this centre's sign: the plane through D that
    contains M and one in-plane reference direction (the heavy centroid for X-H, the highest-priority
    substituent for X-R).  Reflecting the centre's substituent atoms across it is an isometry fixing M
    and D, so it exactly negates the triple product.  Returns None if the plane is degenerate."""
    D = pts[c["d"]]
    M = pts[c["m"]]
    # `plane_ref` selects the IN-PLANE reference separately from `keyed`.  Only the chalcogen
    # branch with ONE chelate arm needs it: there the plane must contain the CHELATE ARM so
    # that the ring stays exactly in place and only the free arm tips over.  `keyed` is left
    # untouched by this, because it carries the eye's sign convention.  If the key is absent,
    # the expression is literally the old one -> byte-identical.
    ref = (np.mean(pts[c["heavy"]], axis=0) if c["type"] == "XH"
           else pts[c.get("plane_ref", c["keyed"][0])])
    n = np.cross(M - D, ref - D)
    nn = float(np.linalg.norm(n))
    if nn < 1e-9:
        return None
    return n / nn


def _find_centers(syms: List[str], pts: np.ndarray, nbrs: List[List[int]],
                  metal_set: set) -> List[dict]:
    """Coordination-created tetrahedral stereocentres — UNIVERSAL (not only N-H): a non-metal donor
    bonded to >=1 metal that becomes chiral once the metal occupies the 4th tetrahedral position.
      * X-H  : exactly one H + >=1 non-metal heavy neighbour (secondary amine / P-H / ...).
      * X-R  : no H + exactly 3 DISTINGUISHABLE (radius-1 env) heavy substituents (tertiary amine / P
               with three different R groups).  Only INVERTIBLE X-R centres are returned — the three
               substituents must be independent (no chelate ring back to the metal, pairwise-disjoint
               arms) so a reflection cleanly inverts them; ring-embedded X-R needs a re-embed (future).
    Returns center dicts sorted by donor index (deterministic)."""
    def _pkey(x):
        return (syms[x], tuple(sorted(syms[y] for y in nbrs[x] if y not in metal_set)))
    centers: List[dict] = []
    # ===== SB AND BI WERE IN THE DONOR LIST AND WERE SKIPPED AS METAL =====
    # Measured 18.08.2026 on 620 diastereomer cases: _STEREO_DONOR_ELEMENTS (:61) lists
    # N, P, As, Sb, Bi -- but ``metal_set`` comes from _coord_angle_corrector._is_metal_sym,
    # and its _METAL_Z_RANGES (:69-73) contains Z 51 (Sb) and Z 83 (Bi).  The metal check
    # stands TWO LINES BEFORE the element check, so Sb and Bi were unreachable: they are in
    # the list, but the path never arrives at them.  What was in effect was {N, P, As},
    # not {N, P, As, Sb, Bi}.
    # Affected: 33 centres in 26 systems, of which 19 are free SbR3 donors without a ring --
    # exactly the XR case the branch below was built for.
    #
    # DELFIN_STEREOCENTER_PNICTOGEN_METALLOID (default 0 -> byte-identical) lets the element
    # check take precedence.  An Sb that is ITSELF the central metal still drops out: the
    # test ``ms`` further down demands a metal neighbour, and the neighbours of a central
    # metal are donors.
    _pnict_metalloid = (os.environ.get(
        "DELFIN_STEREOCENTER_PNICTOGEN_METALLOID", "0") == "1")
    # Default 0 -> `_donor_elems` IS `_STEREO_DONOR_ELEMENTS`, the chalcogen branch below is
    # unreachable, and the run stays byte-identical.
    _chalc = (os.environ.get("DELFIN_STEREOCENTER_CHALCOGEN", "0") == "1")
    _donor_elems = (_STEREO_DONOR_ELEMENTS | _STEREO_CHALCOGEN_ELEMENTS) if _chalc \
        else _STEREO_DONOR_ELEMENTS
    for d in range(len(syms)):
        if d in metal_set and not (_pnict_metalloid
                                   and syms[d] in _donor_elems):
            continue
        if syms[d] not in _donor_elems:
            continue                                   # stable donor centres only
        nb = nbrs[d]
        ms = [x for x in nb if x in metal_set]
        if not ms:
            continue                                   # not coordinated -> no metal-created centre
        hs = [x for x in nb if syms[x] == "H"]
        heavy = [x for x in nb if (x not in metal_set) and syms[x] != "H"]
        if len(hs) == 1 and len(heavy) >= 1:
            centers.append({"type": "XH", "d": d, "m": ms[0], "h": hs[0],
                            "heavy": heavy, "flip": [hs[0]]})
        elif len(hs) == 0 and len(heavy) == 3 and len({_pkey(x) for x in heavy}) == 3:
            # X-R: invertible only if the three arms are independent (no chelate ring, disjoint).
            if any(_reaches_metal(nbrs, x, d, metal_set) for x in heavy):
                continue
            blocked = {d} | metal_set
            subs = [_subtree(nbrs, x, blocked) for x in heavy]
            if len(subs[0] | subs[1] | subs[2]) != sum(len(s) for s in subs):
                continue                               # arms interconnect (ring) -> not cleanly invertible
            keyed = sorted(heavy, key=_pkey)           # priority order == the eye's
            flip = sorted(subs[0] | subs[1] | subs[2])
            centers.append({"type": "XR", "d": d, "m": ms[0], "heavy": heavy,
                            "keyed": keyed, "flip": flip})
        elif (_chalc and syms[d] in _STEREO_CHALCOGEN_ELEMENTS
              and len(hs) == 0 and len(heavy) == 2
              and len({_pkey(x) for x in heavy}) == 2):
            # ===== CHALCOGEN: TWO ARMS + METAL + LONE PAIR =====
            # Pyramidal like a pnictogen, but with one arm fewer.  The inversion pushes
            # the lone pair to the other side -- the metal switches faces.  Three cases,
            # and only the chelate ring decides, not "in a ring":
            chel = [x for x in heavy if _reaches_metal(nbrs, x, d, metal_set)]
            if len(chel) == 2:
                # BOTH arms lead back to the metal.  Then the "inversion" is nothing
                # other than the RING FOLD of the chelate -- it belongs on the fold
                # axis (DELFIN_STEREOCENTER_FAMILY_PARTITION), not here.  Building
                # both would be double-counted completeness.
                # Measured: only 47 of 457 centres, i.e. 10.3 %.
                continue
            blocked = {d} | metal_set
            subs = [_subtree(nbrs, x, blocked) for x in heavy]
            if len(subs[0] | subs[1]) != sum(len(s) for s in subs):
                continue                               # arms interconnect (carbocycle)
            keyed = sorted(heavy, key=_pkey)           # order == the eye's
            if len(chel) == 1:
                # ONE arm in the chelate: the ring STAYS in place, only the free arm is
                # mirrored across the plane through M, D and the chelate arm.
                # _center_plane_normal holds M and D fixed, so this is a true isometry --
                # inner ligand geometry untouched, price +0.98 pp.  `plane_ref` selects the
                # plane, `keyed` stays the eye's order: re-sorting for this would flip the
                # sign in the cross product and the build would contradict the fingerprint.
                frei = [x for x in heavy if x not in chel][0]
                centers.append({"type": "XR", "d": d, "m": ms[0], "heavy": heavy,
                                "keyed": keyed, "plane_ref": chel[0],
                                "flip": sorted(_subtree(nbrs, frei, blocked))})
            else:
                centers.append({"type": "XR", "d": d, "m": ms[0], "heavy": heavy,
                                "keyed": keyed, "flip": sorted(subs[0] | subs[1])})
    centers.sort(key=lambda c: c["d"])
    return centers


def _relax_flipped(syms: List[str], pts: np.ndarray, nbrs: List[List[int]],
                   metal_set: set, free_hs: List[int]) -> np.ndarray:
    """Light constrained UFF relax of ONLY the flipped H atoms (every other atom frozen), on a
    metal-free ligand mol built from the geometric bond graph.  Reuses
    ``assemble_complex._constrained_uff_relax`` (the _hapto_rigid_v2 pattern).  Fail-safe:
    returns ``pts`` unchanged on any error (missing rdkit type, sanitize failure, ff None)."""
    try:
        from rdkit import Chem  # local import (rdkit is a hard dep; keeps import cost off no-op path)
        from delfin.manta.assemble_complex import _constrained_uff_relax
    except Exception:
        return pts
    try:
        keep = [i for i in range(len(syms)) if i not in metal_set]
        old2new = {o: n for n, o in enumerate(keep)}
        rw = Chem.RWMol()
        for o in keep:
            rw.AddAtom(Chem.Atom(str(syms[o])))
        seen = set()
        for o in keep:
            for j in nbrs[o]:
                if j in metal_set:
                    continue
                a, b = (o, j) if o < j else (j, o)
                if (a, b) in seen:
                    continue
                seen.add((a, b))
                rw.AddBond(old2new[a], old2new[b], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, catchErrors=True)
        except Exception:
            pass
        conf = Chem.Conformer(mol.GetNumAtoms())
        for o in keep:
            x, y, z = pts[o]
            conf.SetAtomPosition(old2new[o], (float(x), float(y), float(z)))
        mol.AddConformer(conf, assignId=True)
        free_new = {old2new[h] for h in free_hs if h in old2new}
        fixed = [old2new[o] for o in keep if old2new[o] not in free_new]
        if not _constrained_uff_relax(mol, fixed, max_its=200):
            return pts
        outc = mol.GetConformer()
        newpts = pts.copy()
        for h in free_hs:
            if h in old2new:
                p = outc.GetAtomPosition(old2new[h])
                newpts[h] = np.array([p.x, p.y, p.z], float)
        return newpts
    except Exception:
        return pts


def _flip_clash(syms: List[str], pts: np.ndarray, nbrs: List[List[int]],
                flipped_hs: List[int], h_heavy_min: float, h_h_min: float) -> bool:
    """True if any flipped H sits in a hard steric overlap with a non-bonded atom (the
    multidentate rigidity forbids this fold).  Absolute floors well inside vdW so only genuinely
    impossible folds are pruned — the crystal's fold is real and never trips these."""
    n = len(syms)
    for h in flipped_hs:
        bonded = set(nbrs[h])
        bonded.add(h)
        for j in range(n):
            if j in bonded:
                continue
            dist = float(np.linalg.norm(pts[h] - pts[j]))
            if syms[j] == "H":
                if dist < h_h_min:
                    return True
            elif dist < h_heavy_min:
                return True
    return False


def _coord_iso_key(syms: List[str], pts: np.ndarray, nbrs: List[List[int]],
                   metal_set: set) -> tuple:
    """A cheap, conformer-INVARIANT, enantiomer-SENSITIVE signature of the coordination isomer, so
    the fold axis (which is orthogonal to the coordination + conformer axes) is enumerated ONCE per
    coordination isomer — not once per conformer.  Per metal: the sorted donor elements, the sorted
    multiset of rounded pairwise donor-M-donor cosines (the cis/trans/fac/mer pattern — reflection-
    invariant), AND a global chirality sign (sum of donor-triple scalar-triple-products — flips under
    reflection) so Δ/Λ enantiomers stay DISTINCT groups (their eye fingerprints differ, so each needs
    its own folds).  Two conformers of one isomer -> same key; two isomers -> different key."""
    keys = []
    for m in sorted(metal_set):
        donors = sorted(j for j in nbrs[m] if j not in metal_set)
        if not donors:
            continue
        M = pts[m]
        u = {}
        for d in donors:
            v = pts[d] - M
            nv = float(np.linalg.norm(v))
            u[d] = v / nv if nv > 1e-9 else v
        # cis/trans BUCKET (not a rounded cosine): robust to the heavy angular distortion real builds
        # carry (USEMOW angle_dev ~26deg -> a rounded cosine fragments one isomer across conformers).
        # cos < -0.5  (angle > 120deg) = trans ; else = cis.  Captures cis/trans/fac/mer for mixed donors.
        cosangs = []
        for a in range(len(donors)):
            for b in range(a + 1, len(donors)):
                da, db = donors[a], donors[b]
                ea, eb = sorted((syms[da], syms[db]))
                rel = "t" if float(np.dot(u[da], u[db])) < -0.5 else "c"
                cosangs.append((ea, eb, rel))
        cosangs.sort()
        chir = 0.0
        for a in range(len(donors)):
            for b in range(a + 1, len(donors)):
                for c in range(b + 1, len(donors)):
                    chir += float(np.dot(u[donors[a]], np.cross(u[donors[b]], u[donors[c]])))
        chir_sign = "+" if chir > 1e-3 else ("-" if chir < -1e-3 else "0")
        keys.append((tuple(sorted(syms[d] for d in donors)), tuple(cosangs), chir_sign))
    return tuple(keys)


def _analyze_frame(xyz: str):
    """Parse one frame + locate its coordination-created X-H stereocentres.  Returns a dict (or None
    when the frame has no such stereocentre) carrying everything the fold builder needs, plus the
    coordination-isomer key and this frame's own base sign-vector."""
    syms, pts, lines = _parse_xyz(xyz)
    if len(syms) < 3:
        return None
    nbrs, _bd = _build_geometric_adjacency(syms, pts)
    metal_set = {i for i, s in enumerate(syms) if _is_metal_sym(s)}
    if not metal_set:
        return None
    centers: List[dict] = []
    base_signs: List[str] = []
    for c in _find_centers(syms, pts, nbrs, metal_set):
        s = _center_sign(pts, c)
        if s in ("+", "-"):                            # non-degenerate centre only
            centers.append(c)
            base_signs.append(s)
    if not centers:
        return None
    return {
        "syms": syms, "pts": pts, "lines": lines, "nbrs": nbrs, "metal_set": metal_set,
        "centers": centers, "base_signs": base_signs,
        "group": _coord_iso_key(syms, pts, nbrs, metal_set),
    }


def _build_fold(A: dict, target: List[str], h_heavy_min: float, h_h_min: float):
    """Realise one target fold from a representative frame ``A`` by reflecting each flipped centre's
    substituent atoms across its sign-inverting plane (a single H for X-H; the substituent subtrees for
    X-R), a light frozen-boundary relax of the moved atoms, and a clash-prune.  Returns (fold_xyz, tag)
    or None if the fold is degenerate / sterically infeasible (the multidentate rigidity forbids it)."""
    syms, pts, lines = A["syms"], A["pts"], A["lines"]
    nbrs, metal_set, centers, base = A["nbrs"], A["metal_set"], A["centers"], A["base_signs"]
    k = len(centers)
    flips = [i for i in range(k) if target[i] != base[i]]
    if not flips:
        return None
    reflected = pts.copy()
    flipped_atoms: List[int] = []
    for i in flips:
        c = centers[i]
        normal = _center_plane_normal(pts, c)
        if normal is None:
            return None                                # degenerate plane — cannot realise this fold
        D = pts[c["d"]]
        for a in c["flip"]:
            v = pts[a] - D
            reflected[a] = pts[a] - 2.0 * float(np.dot(v, normal)) * normal
            flipped_atoms.append(a)
    if [_center_sign(reflected, centers[i]) for i in range(k)] != target:
        return None                                    # reflection did not realise the target signs
    relaxed = _relax_flipped(syms, reflected, nbrs, metal_set, flipped_atoms)
    use = relaxed if [_center_sign(relaxed, centers[i]) for i in range(k)] == target else reflected
    if _flip_clash(syms, use, nbrs, flipped_atoms, h_heavy_min, h_h_min):
        return None
    tag = "".join("u" if s == "+" else "d" for s in target)
    return _format_xyz(lines, syms, use), tag


# ---------------------------------------------------------------------------
# MULTI-CENTRE FOLDS BEYOND KMAX (DELFIN_STEREOCENTER_MULTI_CENTRE)
# ---------------------------------------------------------------------------
#
# THE FINDING (18.08.2026).  The corpus carries NO stereochemistry: 9 of 129 314
# SMILES.  Handedness is therefore 100 % a GENERATION DUTY -- there is nothing to
# carry over.  The failure rate is flat (sp3 20.7 % / metal 22.4 %), so what is
# called for is GENERAL generation and not a special rule.  620 diastereomer
# sub-cases are covered by no mechanism today.
#
# THE GAP SITS IN A SINGLE BRANCH.  For k > KMAX the enumeration falls back from
# 2^k to k -- "single-centre flips only".  A single flip is by definition the case
# in which EXACTLY ONE centre differs.  A diastereomer that differs at TWO or more
# centres is no longer in the pot after that.  So the fallback does not throw away
# "a few rare" cases, but by construction EXACTLY the multi-centre diastereomers --
# and that is the class that was measured above as uncovered.
#
# THE ANSWER IS THINNING, NOT DISCARDING.  2^k is already 4096 at k = 12 and thus
# beyond any cap; the question is never "all or none", but WHICH subset.  This
# selection must have two properties:
#
#   1. DETERMINISTIC.  The project checks byte-determinism of whole pools.  A set
#      (`set`) in the OUTPUT ORDER makes the run fail, which is why `seen` here
#      appears only in the duplicate check, never in the ordering.
#
#   2. BALANCED -- and that is exactly what the obvious solution cannot do.  Take
#      `itertools.combinations` in lexicographic order and cut off at the cap, and
#      at k = 10 with cap 32 the pairs are (0,1)..(0,9), (1,2)..(1,9), (2,3)..(2,9),
#      (3,4)..(3,9), (4,5), (4,6).  Participation per centre, measured in the
#      self-test: [9, 9, 9, 9, 6, 5, 5, 4, 4, 4] -- the rear donors appear less
#      than half as often, and a triple never appears.  The cap would then be no
#      truncation but a silent PRE-SELECTION in favour of the low donor indices --
#      i.e. the same mistake that `dofs = dofs[:4]` and `_WELL_MAX_SIBLINGS = 12`
#      already made.
#
# CIRCULANT INSTEAD OF LEXICOGRAPHIC.  A combination is generated as a START INDEX i
# plus a STEP VECTOR (s_1..s_{r-1}) modulo k.  For a fixed step vector, i runs over
# all k centres, so every centre appears in exactly r combinations of this family
# -- the selection is exactly balanced after EVERY complete family and unbalanced
# by at most one combination within a family that has been started.  The cap thus
# truncates evenly instead of one-sidedly.  Measured in the self-test:
# k = 10, cap 32 -> participation per centre [8, 7, 7, 7, 8, 7, 7, 7, 7, 7]
# (the antipode contributes the one to every centre, the circulant pairs the rest).
#
# ORDER (the thinning always cuts from the BACK, so the best-justified targets
# stand at the front):
#   * first all k SINGLE FLIPS -- unchanged, so that every centre is still
#     guaranteed to be OFFERED both signs (the completeness law);
#   * then the ANTIPODE (all k centres flipped).  It is the only multi-combination
#     whose partner is CHEMICALLY GUARANTEED to exist: the antipodal sign pattern
#     is the enantiomeric configuration of the base, same energy, same sterics.
#     It stands before the pairs so that no cap can ever cut it away;
#   * then pairs, then triples, ... , each as circulant families by step sum
#     (tight neighbourhoods first) and lexicographic within the sum.
#     Small flip counts first, because a fold with fewer simultaneously flipped
#     centres survives the clash check `_flip_clash` more often.
#
# COMPLETENESS STAYS REACHABLE: the generation runs over ALL r-subsets (every
# r-subset of Z_k is a rotation of a step vector with sum <= k-1), so the cap is a
# TRUNCATION of a complete enumeration and not a restricted family.  The self-test
# checks exactly that: with a large cap, k = 5 yields exactly the
# 2^5 - 1 - 5 = 26 multi-combinations.
#
# ⚠ DEFAULT OFF.  With DELFIN_STEREOCENTER_MULTI_CENTRE = 0, `_fold_targets` is
# character-identical to the old code in BOTH branches, including the old warning.

def _step_vectors(m: int, span: int):
    """All step vectors (s_1..s_m) with s_j >= 1 and sum <= ``span``, generated
    lazily and stably sorted by (sum, lexicographic).

    The sum comes first because it corresponds to the span of the combination on
    the circle of centres: tight groups first.  Compositions are generated via cut
    points (`combinations` is lexicographic, so the output is too)."""
    if m < 1:
        return
    for total in range(m, span + 1):
        for cuts in itertools.combinations(range(1, total), m - 1):
            prev = 0
            parts = []
            for c in cuts:
                parts.append(c - prev)
                prev = c
            parts.append(total - prev)
            yield tuple(parts)


def _multi_flip_combos(k: int, limit: int) -> List[Tuple[int, ...]]:
    """Deterministic, balanced selection of MULTI-FLIP combinations over k centres:
    ascending-sorted index tuples of length >= 2, at most ``limit`` of them, in the
    order justified above (antipode, then circulant pairs, triples, ...).  Generated
    lazily -- for large k the full power set is never built up, generation stops at
    the cap."""
    out: List[Tuple[int, ...]] = []
    if k < 2 or limit <= 0:
        return out
    seen = set()                                   # duplicate check ONLY, never ordering

    def _emit_combo(idx) -> bool:
        """Takes in one combination; True means 'cap reached, stop'."""
        t = tuple(sorted(idx))
        if len(t) < 2 or t in seen:
            return False
        seen.add(t)
        out.append(t)
        return len(out) >= limit

    if _emit_combo(range(k)):                      # the antipode first -- cannot be cut
        return out
    for r in range(2, k):                          # r == k is the antipode, already in
        for steps in _step_vectors(r - 1, k - 1):
            for i in range(k):
                idx = [i]
                pos = i
                for s in steps:
                    pos = (pos + s) % k
                    idx.append(pos)
                if len(set(idx)) != r:
                    continue                       # cannot occur with sum <= k-1
                if _emit_combo(idx):
                    return out
    return out


def _fold_targets(base: List[str], kmax: int, multi_centre: int, multi_max: int):
    """The target sign patterns for ONE arrangement family.

    Returns ``(targets, multi_info)``.  ``multi_info`` is None when this call
    thinned nothing (k <= KMAX, or switch off -- then the output is
    character-identical to the behaviour before 18.08.2026); otherwise
    ``(k, built, possible)`` for the report line.  A cap that does not report
    reads afterwards like 'there were no more'."""
    k = len(base)
    if k <= kmax:
        return [["+" if (combo >> i) & 1 else "-" for i in range(k)]
                for combo in range(1 << k)], None
    # Rare: 2^k too large.  Still guarantee BOTH signs per centre via single-centre flips.
    targets = []
    for i in range(k):
        t = list(base)
        t[i] = "-" if base[i] == "+" else "+"
        targets.append(t)
    if not multi_centre:
        return targets, None
    combos = _multi_flip_combos(k, multi_max)
    for combo in combos:
        t = list(base)
        for i in combo:
            t[i] = "-" if base[i] == "+" else "+"
        targets.append(t)
    return targets, (k, len(combos), (1 << k) - 1 - k)


def expand_results(results):
    """ADDITIVE stereocentre-fold expansion of a finished ``results`` list of (xyz, label[, ...]).

    Groups the frames by coordination isomer (``_coord_iso_key``) so the fold axis is enumerated ONCE
    per isomer (both +/- at every X-H centre — the completeness law), NOT once per conformer.  Every
    fold already present among the base frames is seeded so it is never re-built; every MISSING
    buildable fold is appended (on the isomer's representative frame) under a ``<base>_stereo-<udud>``
    label.  Originals are preserved verbatim -> never-worse-safe.  Deterministic + fail-safe + bounded
    (``DELFIN_STEREOCENTER_MAX_ADDED`` backstop, logged when hit — no silent truncation)."""
    if not results:
        return results
    max_added = _env_int("DELFIN_STEREOCENTER_MAX_ADDED", 128)
    kmax = _env_int("DELFIN_STEREOCENTER_KMAX", 8)
    h_heavy_min = _env_float("DELFIN_STEREOCENTER_H_HEAVY_MIN", 1.45)
    h_h_min = _env_float("DELFIN_STEREOCENTER_H_H_MIN", 1.25)
    # Multi-centre folds beyond KMAX -- rationale in the block above
    # `_step_vectors`.  Default 0 -> byte-identical.
    multi_centre = _env_int("DELFIN_STEREOCENTER_MULTI_CENTRE", 0)
    multi_max = _env_int("DELFIN_STEREOCENTER_MULTI_MAX", 32)

    # ===== FAMILY PARTITION (DELFIN_STEREOCENTER_FAMILY_PARTITION, default OFF -> byte-identical) =====
    #
    # THE BUG, measured 2026-08-10 on YANLEG (archive_tpr6abl10k_champ, 276 frames):
    #     O-trans   60 frames   ALL anti,  0 syn
    #     N0-trans  60 frames   ALL anti,  0 syn
    #     all-cis   40 frames   ALL anti,  0 syn
    # 160 frames carrying exactly ONE of the two diastereomers, and not one `_stereo-` suffix
    # among them -- every fold this module added went to a family holding a SINGLE frame.
    # The crystal (YANLEG/YANLIK, same SMILES) is the anti fold: O-Mn-O 178.5 deg, the two
    # amine N cis at 84.6 deg, and each N-H pointing at a different phenolate O (H...O 2.65
    # and 2.57 A).  So the built member is right and its (R,R)/(S,S) partner is simply absent.
    #
    # WHY.  ``_coord_iso_key`` is a purely GEOMETRIC signature -- donor elements plus a
    # cis/trans bucket multiset plus a chirality sign.  Different arrangements can share it:
    # YANLEG collapses 38 label families into 17 groups, one of them holding 106 frames from
    # 22 families.  The module then sees the fold "already present" somewhere in the group and
    # adds nothing to the other 21 families.  It reports completeness on a partition COARSER
    # than the one it is completing -- the same failure shape as tier2 on 2026-08-08, where a
    # gate term compared two different frames.
    #
    # THE FIX is to make the group key what the manifold is actually partitioned by: the
    # coordination signature AND the arrangement family.  ``_arrangement_key`` is the codebase's
    # own definition of that family (it strips -confN, the Delta/Lambda hand and the duplicate
    # -N suffix), so this reuses it instead of inventing a second notion; it is imported lazily
    # because smiles_converter imports THIS module, and by call time it is fully loaded.
    #
    # It also fixes the second half of the same bug: ``reps[g]`` is the representative whose
    # LABEL the new fold inherits, so with the family in the key a fold built for O-trans is
    # finally labelled O-trans instead of borrowing some other family's name.
    #
    # ⚠ Default OFF, and it must stay off until measured: it strictly INCREASES the number of
    # folds built (more groups -> more missing targets), so it pushes against DELFIN_STEREOCENTER_
    # MAX_ADDED and against the hard-frame proportion.  It makes the manifold more COMPLETE,
    # which is not automatically more CLEAN.
    _fam_part = _env_int("DELFIN_STEREOCENTER_FAMILY_PARTITION", 0)
    _arrk = None
    if _fam_part:
        try:
            from delfin.smiles_converter import _arrangement_key as _arrk
        except Exception:
            _arrk = None                               # cannot partition -> behave exactly as before

    # Pass 1: analyse every base frame; group by coordination isomer; seed the folds already present.
    reps: dict = {}                                    # group_key -> (analysis, base_label, order)
    present: set = set()                               # (group_key, base_sign_tuple) already in manifold
    for order, entry in enumerate(results):
        try:
            xyz = entry[0]
            lbl = entry[1] if len(entry) > 1 else ""
        except Exception:
            continue
        try:
            A = _analyze_frame(xyz)
        except Exception as exc:
            _LOG.debug("stereocenter_enum: analyse failed for '%s': %s", lbl, exc)
            A = None
        if A is None:
            continue
        g = A["group"]
        if _arrk is not None:
            try:
                g = (g, _arrk(lbl))
            except Exception:
                pass                                   # unparseable label -> fall back to the geometric key
        present.add((g, tuple(A["base_signs"])))
        if g not in reps:                              # representative = FIRST (best-ranked) frame of the isomer
            reps[g] = (A, lbl, order)

    if not reps:
        return results                                 # no coordination-created X-H stereocentre -> bit-exact no-op

    # Pass 2: per isomer, attempt every fold not already present (both signs per centre).
    out = list(results)                                # additive: keep every original frame
    added = 0
    capped = False
    multi_built = 0                                    # tally line: multi-flips taken in
    multi_dropped = 0                                  # tally line: those cut away by the cap
    for g, (A, lbl, _order) in reps.items():
        base = A["base_signs"]
        k = len(base)
        targets, multi_info = _fold_targets(base, kmax, multi_centre, multi_max)
        if k > kmax and multi_info is None:
            _LOG.warning(
                "stereocenter_enum: %d centres > KMAX=%d -> single-centre flips only "
                "(full 2^k fold set not enumerated; raise DELFIN_STEREOCENTER_KMAX)", k, kmax)
        elif multi_info is not None:
            # ⚠ NO SILENT CAPPING.  A cap that does not report looks afterwards like
            # "there were no more" -- exactly the confusion that `dofs[:4]` and
            # `_WELL_MAX_SIBLINGS = 12` produced in this project.
            _k, _n_built, _n_possible = multi_info
            multi_built += _n_built
            multi_dropped += max(0, _n_possible - _n_built)
            if _n_built < _n_possible:
                _LOG.warning(
                    "stereocenter_enum: %d centres > KMAX=%d -> %d single flips + %d of %d "
                    "multi-centre flips (thinned, %d not enumerated; raise "
                    "DELFIN_STEREOCENTER_MULTI_MAX=%d or DELFIN_STEREOCENTER_KMAX=%d)",
                    _k, kmax, _k, _n_built, _n_possible, _n_possible - _n_built,
                    multi_max, kmax)
            else:
                _LOG.info(
                    "stereocenter_enum: %d centres > KMAX=%d -> %d single flips + all %d "
                    "multi-centre flips (nothing thinned)", _k, kmax, _k, _n_built)
        present_signs = {sv for (gg, sv) in present if gg == g}   # folds already in this isomer
        for target in targets:
            if tuple(target) in present_signs:
                continue                               # already built (base frame or a prior fold)
            if added >= max_added:
                capped = True
                break
            try:
                r = _build_fold(A, target, h_heavy_min, h_h_min)
            except Exception as exc:
                _LOG.debug("stereocenter_enum: build fold %s failed: %s", target, exc)
                r = None
            if not r:
                continue
            fxyz, tag = r
            new_lbl = f"{lbl}_stereo-{tag}" if lbl else f"stereo-{tag}"
            out.append((fxyz, new_lbl))
            present.add((g, tuple(target)))
            added += 1
        if capped:
            break
    if capped:
        _LOG.warning(
            "stereocenter_enum: added capped at %d folds (DELFIN_STEREOCENTER_MAX_ADDED); "
            "some feasible folds not emitted", max_added)
    if multi_built or multi_dropped:
        # ONE line for the whole run, so that a grep yields the overall tally and
        # not only the report of the individual family.
        _LOG.info("stereocenter_enum: Mehrzentren-Faltungen: %d angeboten, %d gedeckelt "
                  "(DELFIN_STEREOCENTER_MULTI_MAX=%d)", multi_built, multi_dropped, multi_max)
    return out


# ---------------------------------------------------------------------------
# Self-test:  python delfin/manta/_stereocenter_enum.py
# ---------------------------------------------------------------------------

def _synthetic_frame() -> str:
    """A synthetic octahedral fac-[Co(NHR2)3Cl3] with THREE coordination-created
    N-H stereocentres.  Every N carries the metal, one H and two C.

    ⚠ The azimuths of the three N substituents are DELIBERATELY asymmetric
    (0 / 100 / 215 degrees instead of 0 / 120 / 240).  With a symmetric arrangement,
    H, the centroid of the two C and the M-N axis lie in ONE plane; the triple
    product in `_center_sign` is then exactly 0, the centre is degenerate and is
    rightly rejected.  A symmetrically built test molecule would therefore have
    reported 'no centre found' and looked like a failure of the mechanism, although
    it was the test case itself -- the symmetry IS the absence of handedness.

    The C deliberately carry no H: what is tested is the branching, not the
    chemistry, and the thinner the scaffold, the more surely the folds are
    collision-free."""
    axes = [(np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
            (np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])),
            (np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))]
    tet = np.radians(70.53)                            # 109.47 degrees against the M direction
    rows = [("Co", np.zeros(3))]
    for u, e1, e2 in axes:
        d = 2.10 * u
        rows.append(("N", d))
        for sym, r, phi in (("H", 1.02, 0.0), ("C", 1.47, 100.0), ("C", 1.47, 215.0)):
            p = np.radians(phi)
            dirv = np.cos(tet) * u + np.sin(tet) * (np.cos(p) * e1 + np.sin(p) * e2)
            rows.append((sym, d + r * dirv))
    for u, _e1, _e2 in axes:
        rows.append(("Cl", -2.30 * u))
    out = [str(len(rows)), "synthetischer fac-[Co(NHR2)3Cl3] -- Selbsttest"]
    for s, v in rows:
        out.append(f"{s:4s} {float(v[0]):12.6f} {float(v[1]):12.6f} {float(v[2]):12.6f}")
    return "\n".join(out) + "\n"


def _targets_before_18_08(base, kmax):
    """VERBATIM replica of the branch before 18.08.2026 -- the yardstick for
    'byte-identical with switch 0'.  Needed only in the self-test."""
    k = len(base)
    if k <= kmax:
        return [["+" if (combo >> i) & 1 else "-" for i in range(k)] for combo in range(1 << k)]
    targets = []
    for i in range(k):
        t = list(base)
        t[i] = "-" if base[i] == "+" else "+"
        targets.append(t)
    return targets


def _self_test_folds() -> int:                         # pragma: no cover -- self-test
    class _Catcher(logging.Handler):
        """Catches the report lines -- a capping that does not report is the failure."""

        def __init__(self):
            super().__init__()
            self.records = []

        def emit(self, record):
            self.records.append(record)

    state = {"n": 0, "bad": 0}

    def _pruefe(name: str, ok: bool, extra: str = "") -> None:
        state["n"] += 1
        if not ok:
            state["bad"] += 1
        print(f"{state['n']:2d} {name}: {'OK' if ok else 'FEHLER'}"
              + (f"  [{extra}]" if extra else ""))

    # === 1. DEFAULT OFF = CHARACTER-IDENTICAL to the branch before 18.08. =========
    same = True
    detail = ""
    for k in (1, 3, 5, 8, 9, 10, 14):
        b = ["+" if (i % 3) else "-" for i in range(k)]
        for kmax in (2, 8):
            t_new, info = _fold_targets(b, kmax, 0, 32)
            t_old = _targets_before_18_08(b, kmax)
            if t_new != t_old or info is not None:
                same = False
                detail = f"k={k} kmax={kmax}"
    _pruefe("Schalter AUS ist zeichengleich mit dem alten Zweig (k=1..14, KMAX=2/8)",
            same, detail)

    # === 2. k <= KMAX stays the full 2^k enumeration, the switch changes nothing ==
    b3 = ["+", "-", "+"]
    t_off, i_off = _fold_targets(b3, 8, 0, 32)
    t_on, i_on = _fold_targets(b3, 8, 1, 32)
    _pruefe("k=3 unter KMAX=8: EIN und AUS identisch, 2^3 = 8 Ziele",
            t_off == t_on and len(t_off) == 8 and i_off is None and i_on is None,
            f"aus {len(t_off)} / ein {len(t_on)}")

    # === 3. k > KMAX: single flips PLUS thinned multi-flips =======================
    b10 = ["+" if i % 2 else "-" for i in range(10)]
    t10_off, _i = _fold_targets(b10, 8, 0, 32)
    t10_on, info10 = _fold_targets(b10, 8, 1, 32)
    uniq = len({tuple(t) for t in t10_on}) == len(t10_on)
    no_base = all(t != b10 for t in t10_on)
    _pruefe("k=10 ueber KMAX=8: 10 -> 42 Ziele, alle verschieden, keins ist die Basis",
            len(t10_off) == 10 and len(t10_on) == 42 and uniq and no_base
            and info10 == (10, 32, (1 << 10) - 1 - 10),
            f"aus {len(t10_off)} / ein {len(t10_on)} / Meldung {info10}")

    # === 4. DETERMINISM: two runs, the same order ================================
    det = (_multi_flip_combos(11, 40) == _multi_flip_combos(11, 40)
           and _multi_flip_combos(7, 200) == _multi_flip_combos(7, 200)
           and _fold_targets(b10, 8, 1, 32)[0] == t10_on)
    _pruefe("Determinismus: identische Reihenfolge bei Wiederholung", det)

    # === 5. TRUNCATION, NOT A RESTRICTED FAMILY ==================================
    full5 = _multi_flip_combos(5, 10 ** 6)
    want5 = {c for r in range(2, 6) for c in itertools.combinations(range(5), r)}
    _pruefe("mit grossem Deckel kommen bei k=5 ALLE 26 Mehrfachkombinationen",
            len(full5) == 26 and set(full5) == want5 and len(set(full5)) == len(full5),
            f"{len(full5)} von {len(want5)}")

    # === 6. BALANCE: the cap must not starve any centre ==========================
    combos10 = _multi_flip_combos(10, 32)
    load = [sum(1 for c in combos10 if i in c) for i in range(10)]
    lex = []
    for r in range(2, 11):
        for c in itertools.combinations(range(10), r):
            lex.append(c)
            if len(lex) >= 32:
                break
        if len(lex) >= 32:
            break
    load_lex = [sum(1 for c in lex if i in c) for i in range(10)]
    _pruefe("zirkulant: Beteiligung je Zentrum fast gleich -- lexikographisch nicht",
            max(load) - min(load) <= 1 and max(load_lex) - min(load_lex) > 1,
            f"zirkulant {load} / lexikographisch {load_lex}")

    # === 7. THE ANTIPODE STANDS AT THE FRONT and survives every cap ==============
    _pruefe("der Antipode (alle Zentren gekippt) ist die erste Mehrfachkombination",
            _multi_flip_combos(9, 1) == [tuple(range(9))]
            and _multi_flip_combos(12, 3)[0] == tuple(range(12)))

    # === 8. THE MECHANISM REACHES THE BUILD (end to end) =========================
    frame = _synthetic_frame()
    A = _analyze_frame(frame)
    _pruefe("synthetischer Komplex traegt 3 nicht-entartete N-H-Zentren",
            A is not None and len(A["centers"]) == 3,
            "keins gefunden" if A is None else f"{len(A['centers'])} Zentren, Vorzeichen "
            + "".join(A["base_signs"]))

    res = [(frame, "iso0")]
    prev = {kk: os.environ.get(kk) for kk in
            ("DELFIN_STEREOCENTER_KMAX", "DELFIN_STEREOCENTER_MULTI_CENTRE",
             "DELFIN_STEREOCENTER_MULTI_MAX")}
    os.environ["DELFIN_STEREOCENTER_KMAX"] = "2"       # k=3 is thereby the case k > KMAX
    os.environ["DELFIN_STEREOCENTER_MULTI_CENTRE"] = "0"
    out_off = expand_results(res)
    os.environ["DELFIN_STEREOCENTER_MULTI_CENTRE"] = "1"
    os.environ["DELFIN_STEREOCENTER_MULTI_MAX"] = "32"
    cat = _Catcher()
    _LOG.addHandler(cat)
    out_on = expand_results(res)
    _LOG.removeHandler(cat)
    os.environ["DELFIN_STEREOCENTER_KMAX"] = "8"       # k=3 <= KMAX -> full 2^k
    os.environ["DELFIN_STEREOCENTER_MULTI_CENTRE"] = "0"
    out_full = expand_results(res)

    _pruefe("k=3 > KMAX, Schalter AUS: nur Einzelflips, Original unberuehrt",
            out_off[0] == res[0] and len(out_off) - 1 == 3,
            f"{len(out_off) - 1} Faltungen")
    _pruefe("k=3 > KMAX, Schalter EIN: mehr Faltungen, Einzelflips als Praefix erhalten",
            len(out_on) > len(out_off) and out_on[:len(out_off)] == out_off,
            f"{len(out_on) - 1} statt {len(out_off) - 1} Faltungen")
    _pruefe("EIN bei KMAX=2 liefert GENAU den Satz der vollen 2^k-Aufzaehlung",
            sorted(l for _x, l in out_on) == sorted(l for _x, l in out_full),
            " ".join(sorted(l for _x, l in out_on)))

    # === 9. THE CAPPING REPORTS ITSELF ===========================================
    bilanz = [r for r in cat.records if "Mehrzentren-Faltungen" in r.getMessage()]
    _pruefe("ungekappt: genau eine Bilanzzeile im Log, keine Warnung",
            len(bilanz) == 1 and not [r for r in cat.records
                                      if r.levelno >= logging.WARNING],
            " | ".join(r.getMessage() for r in cat.records))

    os.environ["DELFIN_STEREOCENTER_KMAX"] = "2"
    os.environ["DELFIN_STEREOCENTER_MULTI_CENTRE"] = "1"
    os.environ["DELFIN_STEREOCENTER_MULTI_MAX"] = "2"  # 2 of 4 -> MUST warn
    cat2 = _Catcher()
    _LOG.addHandler(cat2)
    out_cap = expand_results(res)
    _LOG.removeHandler(cat2)
    warn_cap = [r for r in cat2.records
                if r.levelno >= logging.WARNING and "multi-centre flips" in r.getMessage()]
    _pruefe("Deckel MULTI_MAX=2 schlaegt zu und WARNT (keine stille Kappung)",
            len(warn_cap) == 1 and len(out_cap) < len(out_on),
            warn_cap[0].getMessage() if warn_cap else "keine Warnung")

    for kk, vv in prev.items():
        if vv is None:
            os.environ.pop(kk, None)
        else:
            os.environ[kk] = vv

    # === 10. THE NUMBER TABLE FOR THE REPORT =====================================
    print("\n   Zusatzziele je Anordnungsfamilie (KMAX=8, MULTI_MAX=32):")
    print("     k    AUS    EIN   Zusatz   moegliche Mehrfachflips")
    for k in (3, 5, 10):
        b = ["+" if i % 2 else "-" for i in range(k)]
        n_off = len(_fold_targets(b, 8, 0, 32)[0])
        n_on = len(_fold_targets(b, 8, 1, 32)[0])
        print(f"    {k:2d}   {n_off:4d}   {n_on:4d}   {n_on - n_off:+6d}   {(1 << k) - 1 - k:8d}")
    print("   (k <= KMAX ist per Bauart ein No-op: 2^k war dort schon vollstaendig.)")

    print(f"\n{state['n'] - state['bad']}/{state['n']} bestanden")
    return 1 if state["bad"] else 0


if __name__ == "__main__":                             # pragma: no cover -- self-test
    import sys as _sys
    logging.basicConfig(level=logging.INFO, format="   %(levelname)s %(message)s")
    _sys.exit(_self_test_folds())
