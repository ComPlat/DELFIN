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
    ref = (np.mean(pts[c["heavy"]], axis=0) if c["type"] == "XH" else pts[c["keyed"][0]])
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
    for d in range(len(syms)):
        if d in metal_set:
            continue
        if syms[d] not in _STEREO_DONOR_ELEMENTS:
            continue                                   # only pnictogen donors form a stable centre
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
        present.add((g, tuple(A["base_signs"])))
        if g not in reps:                              # representative = FIRST (best-ranked) frame of the isomer
            reps[g] = (A, lbl, order)

    if not reps:
        return results                                 # no coordination-created X-H stereocentre -> bit-exact no-op

    # Pass 2: per isomer, attempt every fold not already present (both signs per centre).
    out = list(results)                                # additive: keep every original frame
    added = 0
    capped = False
    for g, (A, lbl, _order) in reps.items():
        base = A["base_signs"]
        k = len(base)
        if k <= kmax:
            targets = [["+" if (combo >> i) & 1 else "-" for i in range(k)] for combo in range(1 << k)]
        else:
            # Rare: 2^k too large.  Still guarantee BOTH signs per centre via single-centre flips.
            targets = []
            for i in range(k):
                t = list(base)
                t[i] = "-" if base[i] == "+" else "+"
                targets.append(t)
            _LOG.warning(
                "stereocenter_enum: %d centres > KMAX=%d -> single-centre flips only "
                "(full 2^k fold set not enumerated; raise DELFIN_STEREOCENTER_KMAX)", k, kmax)
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
    return out
