"""Deterministic ring-pucker conformer CONSTRUCTION (Cremer-Pople + torsion-
constrained relaxation).

ETKDG samples only the global-minimum ring pucker — cyclohexane embeds
300/300 as the CHAIR, never the twist-boat/boat; a 5-ring never leaves its
lowest envelope.  The higher ring basins are genuine, populated, distinct
conformers that a COMPLETE manifold must contain, but no amount of seed
sampling reaches them.  They have to be CONSTRUCTED explicitly.

Mechanism (license-clean, no CCDC, deterministic):
  1. Displace the ring atoms onto a target Cremer-Pople puckering coordinate
     (Q, theta, phi) — chair / twist-boat / boat for a 6-ring, twist / envelope
     for a 5-ring — keeping their in-plane positions, moving only along the
     ring's mean-plane normal.
  2. Relax with UFF while HOLDING the pucker: every ring torsion is restrained
     to its just-set value +/- a window, so bond lengths and angles relax to
     physical values (C-C -> ~1.53 A, tetrahedral angles) but the pucker is
     preserved instead of collapsing back to the global-minimum chair.
  3. Keep the result iff its ring pucker is genuinely distinct from every
     conformer already in the pool (Cremer-Pople theta/phi separation) and the
     geometry is physical.

The primitive takes an optional ``frozen`` atom set (metal + coordinating
donors) so the SAME construction serves metal-containing chelate rings — the
coordination sphere stays put while the backbone puckers.  Metal-free rings
pass ``frozen=None``.
"""

import os as _os
from typing import List, Optional, Set, Tuple

try:
    import numpy as _np
except Exception:                                    # pragma: no cover
    _np = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    _RDKIT = True
except Exception:                                    # pragma: no cover
    _RDKIT = False


# --- GENERAL Cremer-Pople pucker candidates for ANY ring size ---------------
# A ring of size N has (N-3) puckering degrees of freedom.  Rather than hard-code
# named conformers per size (chair/boat/... only make sense for N=6), we sample
# the Cremer-Pople pucker sphere GENERICALLY: the polar caps (even N -> the
# chair-family "alternating" pucker) plus K phases around the m=2 pseudorotation
# circle (the low-energy boat/twist/envelope family).  A torsion-held relax then
# lets each candidate fall into the nearest genuine minimum, and a Cremer-Pople
# dedup distils the candidates to the DISTINCT populated conformers for THIS
# ring.  This works uniformly for N = 5, 6, 7, 8, 9, ... with no per-size table.
# Rings that get ONLY the flat state (conjugated metallacycles that
# `_is_puckerable` rejects).  Refilled on every `generate` call; a ring in
# here receives EXCLUSIVELY the Q=0 candidate in `_ring_pucker_states`,
# so that a conjugated ring is straightened and not folded.
_FLAT_ONLY: set = set()


def _pucker_candidates(n: int) -> List[Tuple[float, Optional[float], float]]:
    """(q_scale, theta, phi) per candidate.  `q_scale` multiplies `_amp(n)`.

    ⚠ WHY THE TUPLES NOW HAVE THREE VALUES (17.08.2026).  Up to here this function
    sampled the Cremer-Pople sphere at FIXED RADIUS: `_amp(n)` gives 0.40 A (5-ring)
    to 0.80 A (8-ring) and was passed through unchanged at both call sites.  The
    candidate `(0.0, 0.0)` below is **theta = 0**, i.e. the CHAIR (polar cap) -- **not**
    Q = 0.  **The centre of the sphere, the PLANE, was not a candidate.**  For odd
    rings it was even narrower: only the equatorial pseudorotation (`theta=None`), i.e.
    envelope and twist, never flat.

    MEASURED on 16./17.08. (`folds`, 965 systems): of 326 missing ring motifs,
    **218 are PLANAR** -- `5M:planar` 126, `6M:planar` 55, `4M:planar` 37 -- against `6M:boat` 30,
    `6:chair` 25, `5M:puckered` 19, all of which this generator can produce.  **Two thirds
    of the gap are exactly the one state it does not know by construction.**

    Chemically this is no edge case: a five-membered chelate ring with sp2 donors often lies
    FLAT; the generator treats it like cyclopentane.

    `_set_pucker` needs NO change for this: with Q = 0, q2 and q3 become zero, zj = 0,
    and every non-frozen ring atom is projected onto the mean plane -- metal and
    donors stay put because they are in `frozen`.  That is exactly the planar state.

    Default OFF -> the list is identical to before (all q_scale = 1.0), hence
    byte-identical.
    """
    cands: List[Tuple[float, Optional[float], float]] = []
    even = (n % 2 == 0)
    # equatorial pseudorotation ring — sample fine enough to hit both the boat
    # (phi = 0, 360/n, ...) and the twist (phi halfway between) positions.
    K = max(8, 2 * n)
    for k in range(K):
        phi = 360.0 * k / K
        cands.append((1.0, 90.0 if even else None, phi))
    if even:
        # polar caps: the chair / inverted-chair (alternating) puckers
        cands.append((1.0, 0.0, 0.0))
        cands.append((1.0, 180.0, 0.0))
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        # THE FLAT STATE.  Q = 0 -> all non-frozen ring atoms onto the
        # mean plane.  ONE candidate per ring, not K -- the plane has no phi.
        cands.append((0.0, 0.0, 0.0))
    return cands


def _amp(n: int) -> float:
    """Cremer-Pople puckering amplitude (Angstrom) as a LAW instead of a table.

    ⚠ THE OLD VERSION CONTRADICTED ITSELF.  It carried a calibrated table
    for 5..8 AND a linear fallback `0.45 + 0.06*n` for everything else -- and the
    fallback lies EVERYWHERE above the table:

        n     table     fallback
        5      0.40      0.75      <- almost double
        6      0.63      0.81
        7      0.72      0.87
        8      0.80      0.93

    This was never noticed because `_is_puckerable` rejected every ring outside 5..8:
    the fallback NEVER RAN.  A dark branch with a wrong number in it.

    The table values saturate (increments 0.23 / 0.09 / 0.08) -- a large ring
    does not fold arbitrarily deep, the amplitude approaches a bound.  A
    linear fallback is thus qualitatively wrong, not merely numerically off.

    Replacement: ONE saturating law for all N that reproduces the calibrated
    values (max. deviation 0.04 A):

        Q_max(N) = 1.15 * (N-4) / (N-4+1.6)
        N=5 0.44 · N=6 0.64 · N=7 0.75 · N=8 0.82 · N=12 0.98 · N=24 1.07

    ⚠ From N=9 on this is EXTRAPOLATION, not calibration -- stated honestly, not
      hidden.  And it is only the UPPER BOUND anyway: `_pucker_space_grid` samples
      the amplitude from 0 up to here, the relax decides what survives.

    ⛔ THE CALIBRATED VALUES STAY EXACTLY AS THEY ARE.  `_amp` is also read by the LEGACY
      path (`_set_pucker(..., _qs * _amp(n), ...)`).  If the law replaced them,
      the default would NOT be byte-identical -- at N=5 it would read 0.44 instead of
      0.40.  The law therefore applies only where the wrong fallback used to stand:
      outside 5..8.  Byte-identity is not a formality, it is the condition
      for an A/B measuring the mechanism and not the instrument.
    """
    _kal = {5: 0.40, 6: 0.63, 7: 0.72, 8: 0.80}
    if n in _kal:
        return _kal[n]
    if n <= 4:
        return 0.35
    return 1.15 * (n - 4.0) / (n - 4.0 + 1.6)


def _ring_normal_and_center(P, ring):
    n = len(ring)
    R = P[list(ring)]
    C = R.mean(0)
    Rc = R - C
    Rp = sum(Rc[j] * _np.sin(2 * _np.pi * j / n) for j in range(n))
    Rpp = sum(Rc[j] * _np.cos(2 * _np.pi * j / n) for j in range(n))
    nrm = _np.cross(Rp, Rpp)
    ln = float(_np.linalg.norm(nrm))
    if ln < 1e-9:
        return None, C
    return nrm / ln, C


def _cp_theta_phi(P, ring):
    """Cremer-Pople (Q, theta, phi) for the ring (atom coords in ring order)."""
    n = len(ring)
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return 0.0, 0.0, 0.0
    z = (P[list(ring)] - C) @ nrm
    q2c = _np.sqrt(2.0 / n) * sum(z[j] * _np.cos(2 * _np.pi * 2 * j / n) for j in range(n))
    q2s = -_np.sqrt(2.0 / n) * sum(z[j] * _np.sin(2 * _np.pi * 2 * j / n) for j in range(n))
    q2 = float(_np.hypot(q2c, q2s))
    if n % 2 == 0:
        q3 = _np.sqrt(1.0 / n) * sum(z[j] * ((-1) ** j) for j in range(n))
    else:
        q3 = 0.0
    Q = float(_np.sqrt(q2 ** 2 + q3 ** 2))
    theta = float(_np.degrees(_np.arctan2(q2, q3))) if n % 2 == 0 else 90.0
    phi = float(_np.degrees(_np.arctan2(q2s, q2c)) % 360.0)
    return Q, theta, phi


def _set_pucker(conf, ring, Q, theta, phi, frozen: Optional[Set[int]] = None):
    n = len(ring)
    P = conf.GetPositions()
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return
    ph = _np.radians(phi)
    even = (n % 2 == 0)
    th = _np.radians(theta) if (theta is not None) else None
    q2 = Q * (_np.sin(th) if th is not None else 1.0)
    q3 = Q * (_np.cos(th) if th is not None else 0.0)
    frozen = frozen or set()
    for j, idx in enumerate(ring):
        # a frozen ring atom (metal / coordinating donor of a chelate ring) keeps
        # its position -> only the backbone puckers, the coordination sphere is
        # preserved.
        if int(idx) in frozen:
            continue
        zj = _np.sqrt(2.0 / n) * q2 * _np.cos(ph + 2 * _np.pi * 2 * j / n)
        if even and th is not None:
            zj += _np.sqrt(1.0 / n) * q3 * ((-1) ** j)
        p = P[idx]
        inplane = p - ((p - C) @ nrm) * nrm
        newp = inplane + zj * nrm
        conf.SetAtomPosition(int(idx), (float(newp[0]), float(newp[1]), float(newp[2])))


def _prod_laenge(per_ring_states) -> int:
    """Size of the full cross product (incl. the ground state)."""
    n = 1
    for s in per_ring_states:
        n *= max(1, len(s))
    return n


def _ring_bahnen(mol, rings, max_aut: int = 20000):
    """Decompose the rings into ORBITS under the automorphism group of the molecule.

    Two rings lie in the same orbit if an automorphism maps one onto the other
    as an ATOM SET.  Only then are they indistinguishable, and only then may
    their state ordering be merged.

    ⛔ WHY NOT THE RANK MULTISET, which would be much cheaper.  Equal
    canonical ranks are NECESSARY for equivalence, but not SUFFICIENT.
    Recomputed by hand on 27.08.: of 27 ring pairs with equal rank set,
    26 were a genuine orbit and ONE was not (3.7 %).  Reducing via the rank
    set would have deleted a REAL state for that one pair.  One
    duplicate too many costs compute time; a missing state costs
    completeness -- and that is the north star.

    ⚠ CAP.  `GetSubstructMatches(mol, mol)` can explode for highly symmetric
    molecules (on 27.08. CONPUS ran into 200 000 matches).  If the cap is
    reached, the function returns None -> NO reduction, full
    product.  The fallback is always the LARGER set, never the smaller.

    Return: list of lists of ring indices, or None (do not reduce).
    """
    try:
        treffer = mol.GetSubstructMatches(mol, uniquify=False,
                                          useChirality=False,
                                          maxMatches=max_aut)
    except Exception:
        return None
    if not treffer or len(treffer) >= max_aut:
        return None                      # cap reached -> do not reduce
    if len(treffer) == 1:
        return None                      # only the identity -> nothing to gain

    ring_mengen = [frozenset(int(i) for i in r) for r in rings]
    index_von = {m: i for i, m in enumerate(ring_mengen)}
    if len(index_von) != len(ring_mengen):
        return None                      # duplicate ring sets -> hands off

    eltern = list(range(len(rings)))

    def _wurzel(x):
        while eltern[x] != x:
            eltern[x] = eltern[eltern[x]]
            x = eltern[x]
        return x

    for abb in treffer:
        for i, menge in enumerate(ring_mengen):
            try:
                bild = frozenset(int(abb[a]) for a in menge)
            except Exception:
                return None
            j = index_von.get(bild)
            if j is None or j == i:
                continue
            ri, rj = _wurzel(i), _wurzel(j)
            if ri != rj:
                eltern[ri] = rj

    gruppen = {}
    for i in range(len(rings)):
        gruppen.setdefault(_wurzel(i), []).append(i)
    # Deterministic order -- otherwise the enumeration depends on the
    # hash order and the run is not reproducible.
    return [sorted(v) for _k, v in sorted(gruppen.items())]


def _cp_abstand(a, b) -> float:
    """Angular distance of two pucker states ON the Cremer-Pople sphere (degrees).

    a, b are (Q, theta, phi).  The great-circle distance is used

        cos d = cos(th_a) cos(th_b) + sin(th_a) sin(th_b) cos(ph_a - ph_b)

    ⚠ WHY NOT SIMPLY |dtheta| + |dphi|.  At the POLE (theta = 0 or 180) phi is
      MEANINGLESS -- a ring in the chair has no phase.  A naive metric holds
      two chairs with phi = 136 and phi = 339 to be 200 degrees apart, although they
      are THE SAME state.  Exactly that is in the measurement of 26.08.:
          theta=180.0  phi=136.4
          theta=180.0  phi=338.6      <- identical, only the phase is noise
      The great-circle distance takes care of this by itself: at sin(theta) = 0 the
      phi term drops out.  Geometry solves the problem, not a special rule.
    """
    ta, pa = _np.radians(a[1]), _np.radians(a[2])
    tb, pb = _np.radians(b[1]), _np.radians(b[2])
    c = (_np.cos(ta) * _np.cos(tb)
         + _np.sin(ta) * _np.sin(tb) * _np.cos(pa - pb))
    return float(_np.degrees(_np.arccos(max(-1.0, min(1.0, float(c))))))


def _set_pucker_general(conf, ring, qs, phis, frozen: Optional[Set[int]] = None):
    """Cremer-Pople inversion in FULL generality -- for EVERY ring size.

    ``_set_pucker`` above covers only m = 2 plus the alternation term.  That is
    complete for N = 4, 5, 6 and INCOMPLETE from N = 7 on: a seven-ring has four
    puckering degrees of freedom (q2, phi2, q3, phi3), an eight-ring five.  The pairs with
    m >= 3 were missing there without replacement, and the alternation term was fixed as q3
    -- but for the eight-ring it is q4.

    THE SPACE, exactly.  A ring with N atoms has N-3 puckering degrees of freedom:
        N even:  pairs (q_m, phi_m) for m = 2 .. N/2-1,  plus ONE q_(N/2)
                   2*(N/2-2) + 1 = N-3
        N odd:   pairs (q_m, phi_m) for m = 2 .. (N-1)/2
                   2*((N-1)/2 - 1) = N-3
    The displacement of the j-th ring atom from the mean plane is

        z_j = sqrt(2/N) * SUM_m  q_m * cos(phi_m + 2*pi*m*j/N)
              + [N even]  sqrt(1/N) * q_(N/2) * (-1)^j

    The named forms are POINTS on it, not separate cases: chair at the poles
    (q2 = 0), boat and twist at the equator (q3 = 0), half-chair and envelope
    IN BETWEEN -- exactly the region the old candidate list never sampled.

    ``qs``/``phis``: mappings m -> value.  Reduces exactly to ``_set_pucker`` for
    N <= 6; the formula is the same, just no longer truncated to m = 2.
    ⚠ ``frozen`` remains untouched -- metal and donors stay put, only the backbone
    folds.  Pure generator.
    """
    n = len(ring)
    P = conf.GetPositions()
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return
    frozen = frozen or set()
    even = (n % 2 == 0)
    m_last = n // 2 if even else None
    for j, idx in enumerate(ring):
        if int(idx) in frozen:
            continue
        zj = 0.0
        for m, qm in qs.items():
            if not qm:
                continue
            if even and m == m_last:
                zj += _np.sqrt(1.0 / n) * qm * ((-1) ** j)
            else:
                ph = _np.radians(phis.get(m, 0.0))
                zj += _np.sqrt(2.0 / n) * qm * _np.cos(ph + 2.0 * _np.pi * m * j / n)
        p = P[idx]
        inplane = p - ((p - C) @ nrm) * nrm
        newp = inplane + zj * nrm
        conf.SetAtomPosition(int(idx), (float(newp[0]), float(newp[1]), float(newp[2])))


def _pucker_space_grid(n: int, n_amp: int, n_phase: int):
    """SYSTEMATIC grid over the ENTIRE pucker space of an N-ring.

    Yields candidates as ``(qs, phis)`` -- mappings m -> value -- for
    ``_set_pucker_general``.  Instead of named forms, the (N-3)-dimensional
    Cremer-Pople space is sampled; chair, boat, twist, half-chair and envelope
    arise by themselves as grid points.

    ⚠ COMPLETENESS IS A QUESTION OF RESOLUTION, not a yes/no question.  A
    continuous space cannot be sampled "entirely".  What stands here is the
    honest version: the space is covered COMPLETELY at the STATED resolution,
    and the resolution is in the trace.  No corner is left out,
    no direction preferred -- the difference from the old list, which knew only equator and
    poles and never varied the amplitude.

    ⚠ PRICE: the candidate count grows like (n_amp+1)^(#q) * n_phase^(#phi).
    Six-ring at n_amp=2, n_phase=8: 3 * 8 * 5 = 120 per ring.  Eight-ring: considerably
    more.  That is why both resolutions are env parameters and appear in the log.

    ===== THE RESOLUTION MUST FALL WITH THE DIMENSION (26.08.2026) =================

    COMPUTED, not estimated.  The number of candidates is

        prod over m_pairs of (1 + n_amp * n_phase)   times   (2*n_amp+1) for even n

    and thus at n_amp=2, n_phase=8:

        n= 8    845          n=10   24 565        n=12     417 605
        n= 9  4 913          n=11   83 521        n=16  120 687 845

    A 16-ring would thus have received 1.2e8 candidates per ring, each with relax and
    collision gate.  That is not a slow run, that is a run that dies and
    delivers ZERO folds.  Exactly so was `foldspace6k` queued (NAMP=2 NPHASE=8).
    ⇒ Infinite fineness is not completeness, it is infeasibility.

    WHAT DOES **NOT** HAPPEN HERE: no result is truncated.  The list
    remains the COMPLETE product of the chosen resolution -- what is reduced is the
    SAMPLING DENSITY, and first where it carries the least physically.
    Cremer-Pople amplitudes q_m fall with m: the high m are the fine ripple
    with small displacement, m=2 is the dominant fold.  Therefore the phase count
    is halved at the LARGEST m first and m=2 is touched last.

    ⚠ And it does NOT happen silently: `_grid_res` carries the per-m chosen resolution,
      the caller writes it into the log under DELFIN_FFFREE_PUCKER_TRACE.  Whether
      the chosen density suffices is not said by this code, but by
      `selbsttest_konvergenz` -- the number of DISTINGUISHABLE states, not the
      number of grid points, is the measure.
    """
    if n < 4:
        return []
    even = (n % 2 == 0)
    m_pairs = list(range(2, (n // 2) if even else ((n - 1) // 2) + 1))
    m_last = (n // 2) if even else None
    amp = _amp(n)
    _budget = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000") or 50000))
    # ===== WHICH MODES ARE EXCITED AT ALL?  (26.08.2026) ============================
    #
    # The self-test found the saving point itself: for the 21-ring (crown ether,
    # e.g. VEDCOA) the budget only sufficed if m=2 was reduced as well --
    # and m=2 is the DOMINANT fold.  Saved at the wrong end.
    #
    # The cause is not the phase count but the NUMBER OF PAIRS: it grows
    # like N/2, and the mere amplitude selection alone costs 3^(N/2-1).  At N=30
    # that is 8 million points BEFORE a single phase is sampled.
    #
    # Cremer-Pople amplitudes of real rings fall sharply with m: the low modes
    # carry the fold, the high ones are fine ripple near zero.  Large rings
    # are described in the literature by a few low modes for exactly this reason.
    # ⇒ Modes above M_MAX are set to amplitude 0 -- the statement is
    #   "this mode is NOT EXCITED", not "this mode was skipped".
    #   The degree of freedom stays in the parametrisation, it is merely at zero.
    #
    # ⚠ THIS IS A MODEL ASSUMPTION, NOT A MEASUREMENT.  It is testable and MUST be
    #   tested: raise M_MAX and ask `selbsttest_konvergenz` whether the number of
    #   DISTINGUISHABLE states changes.  If it changes, M_MAX is too small.
    #   Until then it is in the trace and carries its name.
    _mmax = max(2, int(_os.environ.get("DELFIN_FFFREE_PUCKER_MMAX", "4") or 4))
    _aktiv = [m for m in m_pairs if m <= _mmax]
    _ruhend = [m for m in m_pairs if m > _mmax]
    # Phase count per m; start equal everywhere, then halve from the top down.
    _ph_m = {m: max(1, int(n_phase)) for m in _aktiv}
    m_pairs = _aktiv

    def _zahl():
        t = (2 * n_amp + 1) if m_last is not None else 1
        for _m in m_pairs:
            t *= (1 + n_amp * _ph_m[_m])
        return t

    while _zahl() > _budget:
        _kand = [m for m in m_pairs if _ph_m[m] > 1]
        if not _kand:
            break                                  # already at phase 1 -- nothing more to gain
        _hoch = max(_kand)                         # largest m first: smallest amplitude
        _ph_m[_hoch] = max(1, _ph_m[_hoch] // 2)
    # ⚠ `m_last` (the alternation term q_{N/2}) ALWAYS stays active and is never
    #   set dormant: for the six-ring it IS the chair.  It costs nothing either --
    #   a single signed amplitude factor (2*n_amp+1), not exponential.
    _pucker_space_grid._grid_res = {"n": n, "n_amp": n_amp, "n_phase_je_m": dict(_ph_m),
                                    "kandidaten": _zahl(), "budget": _budget,
                                    "m_max": _mmax, "ruhende_moden": list(_ruhend),
                                    "m_last": m_last,
                                    "reduziert": (any(v < n_phase for v in _ph_m.values())
                                                  or bool(_ruhend))}
    # Amplitude steps per pair: 0 (axis flat) up to n_amp * amp.  The zero MUST be
    # included -- it is the planar state, and exactly that was missing (218 of 326 motifs).
    lv_pair = [amp * k / max(1, n_amp) for k in range(0, n_amp + 1)]
    # The alternation term runs SIGNED: +q is the chair, -q the inverted one.
    lv_last = [amp * k / max(1, n_amp) for k in range(-n_amp, n_amp + 1)]
    # Phases PER m -- same formula, just with the density chosen for this m.
    _phasen = {m: [360.0 * k / _ph_m[m] for k in range(_ph_m[m])] for m in m_pairs}

    out = []

    def _rek(i, qs, phis):
        if i < len(m_pairs):
            m = m_pairs[i]
            for q in lv_pair:
                if q == 0.0:                      # amplitude 0 -> phase meaningless
                    _rek(i + 1, {**qs, m: 0.0}, {**phis, m: 0.0})
                else:
                    for ph in _phasen[m]:
                        _rek(i + 1, {**qs, m: q}, {**phis, m: ph})
            return
        if m_last is not None:
            for q in lv_last:
                out.append(({**qs, m_last: q}, dict(phis)))
        else:
            out.append((dict(qs), dict(phis)))

    _rek(0, {}, {})
    # keep the zero point (everything flat) exactly ONCE -- it is the planar state
    _seen = set()
    uniq = []
    for qs, phis in out:
        key = tuple(sorted((m, round(q, 6), round(phis.get(m, 0.0), 3) if q else 0.0)
                           for m, q in qs.items()))
        if key in _seen:
            continue
        _seen.add(key)
        uniq.append((qs, phis))
    return uniq


def _relax_hold_pucker(mol, ring, frozen: Set[int], window: float = 18.0, iters: int = 1200) -> bool:
    """UFF-relax that frees bonds+angles but HOLDS the pucker: each ring torsion
    restrained to its current value +/- ``window``; any ``frozen`` atom fixed."""
    conf = mol.GetConformer()
    n = len(ring)
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for i in range(n):
        a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
        try:
            t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
            ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
        except Exception:
            pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


def _relax_hold_pucker_multi(mol, rings, frozen: Set[int], window: float = 18.0,
                             iters: int = 1500) -> bool:
    """UFF relax holding the pucker of EVERY ring simultaneously (each ring
    torsion restrained to its current value +/- window); frozen atoms fixed."""
    conf = mol.GetConformer()
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for ring in rings:
        n = len(ring)
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            try:
                t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
                ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
            except Exception:
                pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


_VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
        "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10}


def _has_bad_angles(mol, tol: float = 25.0, skip: Optional[Set[int]] = None) -> bool:
    """True if any heavy centre's VSEPR angle is off its hybridisation ideal by
    more than ``tol`` — i.e. the (multi-)ring pucker left a LOCAL geometry broken
    even though nothing clashes.  Essential for FUSED / BRIDGED ring systems
    (ACEQAC-type lactams): puckering one ring independently strains the shared
    fusion atoms, distorting their angles; the clash gate is blind to it.  A
    frame is realistic only if EVERY VSEPR body is correct, so any such pucker is
    rejected.  2-coordinate centres are hybridisation-ambiguous (sp/sp2/sp3) ->
    skipped; >=5 is non-molecular -> skipped.

    ``skip``: centres exempt from the VSEPR ideal, because they do not HAVE one.
    A COORDINATION centre is the case this exists for: its angles are set by the
    polyhedron, not by hybridisation.  A CN4 metal has nh == 4, so this function
    would demand 109.5 deg of it -- and a square-planar d8 has two 180 deg trans
    angles, i.e. a 70.5 deg "error" that no pucker caused and no pucker can fix.
    Without the exemption EVERY combination is rejected for every SP-4 and T-3
    complex, which reads as "the lever has no reach" for entirely the wrong reason.
    Default None -> empty -> byte-identical for every existing caller."""
    _skip = skip or frozenset()
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        for c in range(mol.GetNumAtoms()):
            if c in _skip:
                continue
            a = mol.GetAtomWithIdx(c)
            if a.GetSymbol() == "H":
                continue
            hv = [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetSymbol() != "H"]
            nh = len(hv)
            if nh < 3 or nh >= 5:
                continue
            exp = 120.0 if nh == 3 else 109.5
            nbset = {x.GetIdx() for x in a.GetNeighbors()}
            for i in range(len(hv)):
                for j in range(i + 1, len(hv)):
                    # skip 3-membered rings (real ~60deg geometry)
                    if hv[j] in {x.GetIdx() for x in mol.GetAtomWithIdx(hv[i]).GetNeighbors()}:
                        continue
                    v1 = P[hv[i]] - P[c]
                    v2 = P[hv[j]] - P[c]
                    d = float(_np.linalg.norm(v1) * _np.linalg.norm(v2))
                    if d < 1e-9:
                        continue
                    ang = _np.degrees(_np.arccos(max(-1.0, min(1.0, float(_np.dot(v1, v2) / d)))))
                    if abs(ang - exp) > tol:
                        return True
    except Exception:
        return False
    return False


def _has_clash(mol, frac: float = 0.60) -> bool:
    """True if any non-bonded heavy-atom pair (topological distance > 3 bonds)
    overlaps below ``frac`` x sum-of-vdW-radii — i.e. the combined ring puckers
    left the whole molecule sterically unrealistic despite the relax."""
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        dm = Chem.GetDistanceMatrix(mol)
        heavy = [i for i, s in enumerate(syms) if s != "H"]
        for a in range(len(heavy)):
            i = heavy[a]
            ri = _VDW.get(syms[i], 1.7)
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                if dm[i][j] <= 3:            # bonded / 1-3 / 1-4 -> expected close
                    continue
                d = float(_np.linalg.norm(P[i] - P[j]))
                if d < frac * (ri + _VDW.get(syms[j], 1.7)):
                    return True
    except Exception:
        return False
    return False


def _bindungs_ausreisser(mol, tol_lang: float = 1.30) -> frozenset:
    """The bonds of the GRAPH whose LENGTH no longer describes a bond.

    ===== WHY THE THIRD GATE IS MISSING AT ALL (26.08.2026) ==========================

    `generate` has two gates: `_has_clash` sees NON-bonded pairs that stand too
    close, `_has_bad_angles` sees angles.  The BOND LENGTH itself is checked by
    neither.  A fold that pulls a bond apart thus gets through.
    Exactly that happened on 18.08. and is on record: "in one frame a bond
    broke".

    ⚠ AND THE SELF-GATE CANNOT CATCH IT -- by construction, not out of
      carelessness.  `assemble_complex._collapsed_heavy_bonds_strict` runs over
      all heavy-atom PAIRS and decides from the DISTANCE whether they are bonded:

          if d > 1.30 * ideal:   continue        # "so not bonded at all"

      A bond stretched to 1.4 times its length thus DROPS OUT OF THE CHECK.
      It is not reported as broken but as absent.  The collapse
      (too short) is seen, the rupture (too long) is a blind spot.

    HERE the bond graph is available.  Thus "bonded" is no longer a distance question,
    and the same number 1.30 turns from an EXCLUSION criterion into the RUPTURE threshold.  Nothing
    is invented: floor (`_bd.COLLAPSE_FLOOR`, 0.82) and ceiling (1.30) are
    exactly the two numbers the self-gate already computes with anyway, and
    `_ideal_bond` is the same source.

    ⚠ SET INSTEAD OF BOOLEAN, and that is the whole difference between filter and
      verdict.  A nitrile sits at 1.20 A against a single-bond ideal of 1.52 --
      ratio 0.79, below the floor.  An absolute yes/no would reject EVERY fold
      of EVERY nitrile-containing molecule, and the finding would read "the gate has
      no reach" for the wrong reason.  The caller therefore subtracts the set
      of the GROUND STATE: only what the fold NEWLY INTRODUCES is rejected.  That is
      the same never-worse form that `converter_backend` imposes on its siblings.
    """
    try:
        from delfin.manta import _bond_decollapse as _bd
        P = mol.GetConformer().GetPositions()
    except Exception:
        return frozenset()
    try:
        floor = float(_bd.COLLAPSE_FLOOR)
    except Exception:
        floor = 0.82
    aus = set()
    try:
        for b in mol.GetBonds():
            i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
            si = mol.GetAtomWithIdx(i).GetSymbol()
            sj = mol.GetAtomWithIdx(j).GetSymbol()
            if si == "H" or sj == "H":
                continue                       # H as in the self-gate: not judged
            try:
                if _bd._is_metal(si) or _bd._is_metal(sj):
                    continue                   # M-D ideal is made up, see `_ideal_bond`
            except Exception:
                pass
            try:
                ideal = float(_bd._ideal_bond(si, sj, bool(b.GetIsAromatic())))
            except Exception:
                continue
            if ideal <= 1e-6:
                continue
            r = float(_np.linalg.norm(P[i] - P[j])) / ideal
            if r > tol_lang or r < floor:
                aus.add((min(i, j), max(i, j)))
    except Exception:
        return frozenset()
    return frozenset(aus)


def _is_puckerable(mol, ring) -> bool:
    """A ring is puckerable iff it is saturated enough to have out-of-plane
    minima: non-aromatic, size 5-8, and >= 3 sp3 ring atoms (an aromatic /
    fully-conjugated ring is planar and rigid -> no pucker conformers).

    Robust for BOTH a sanitised RDKit mol (uses hybridisation) AND a distance-
    perceived metal-complex mol (no hybridisation, all bonds single -> cannot
    tell benzene from cyclohexane from the graph, so read sp3 from the 3D SHAPE:
    a saturated centre is 4-coordinate tetrahedral or 3-coordinate pyramidal,
    an aromatic/sp2 centre is 3-coordinate planar)."""
    n = len(ring)
    # ===== THE SIZE WINDOW 5..8 IS A SHACKLE, NOT A LAW (26.08.2026) ================
    #
    # Cremer-Pople holds for EVERY ring from N = 4: the number of puckering degrees of freedom
    # is N-3, and `_set_pucker_general` by now carries all of them.  The window here
    # nevertheless cut off at 8 -- a four-ring (1 DOF, genuine butterfly fold)
    # and EVERY macrocycle from 9 on were thus unfoldable by construction.
    # Porphyrins, calixarenes, crown ethers, large chelate rings: zero folds, not
    # because the mathematics is missing, but because a number stood in the way.
    # ⚠ Side finding: `_amp` carries a table for 5..8 AND a fallback for the
    #   rest -- and the fallback lies EVERYWHERE above the table (n=5: 0.75 against
    #   0.40, almost double).  Because this window rejected every other ring, the
    #   fallback NEVER ran.  It is corrected together with the window.
    # ⛔ Default OFF -> old window -> byte-identical.
    if _os.environ.get("DELFIN_FFFREE_PUCKER_SPACE", "0") == "1":
        if n < 4:
            return False
    elif n < 5 or n > 8:
        return False
    try:
        P = mol.GetConformer().GetPositions()
    except Exception:
        P = None
    # ===== THE METALLACYCLE FAILS ON A CRITERION THAT DOES NOT APPLY TO IT ==========
    #
    # MEASURED 26.08. on 400 systems (`find_conformer_coverage`):
    #     rings total 1707 · METAL skipped 352 = 20.6 %
    #     six-rings 176, reaching the chair 57  ->  67.6 % NEVER
    #
    # The two conditions below are right for an ORGANIC ring and wrong for
    # a METALLACYCLE, and for the same reason: they look for the
    # softness at the RING ATOMS.  In a chelate ring it sits in the M-D BONDS
    # -- 2.0 to 2.4 A long, soft, with a low torsion barrier.  A salen ring
    # M-N=C-C(ar)-C(ar)-O really folds at exactly these two bonds (the "step"
    # or umbrella fold), although its organic part is rigid and aromatic.
    #
    #   * `GetIsAromatic() -> return False` tips the WHOLE ring as soon as ONE atom
    #     is aromatic.  In the fused salen metallacycle those are the
    #     phenolate carbons -> immediate exclusion.
    #   * `n_sat >= 3` demands three sp3 ring atoms.  A conjugated chelate ring does
    #     not have them and does not need them either.
    #
    # ⇒ For a ring WITH a metal: aromatic is an exclusion only if
    #   ALL non-metal ring atoms are aromatic (then the ring really is
    #   planar-rigid, e.g. a metallabenzene).  And the sp3 threshold drops to 1,
    #   because the metal itself provides the hinge, not an sp3 centre.
    # ⚠ The coordination sphere remains untouched: `frozen` holds metal AND donors
    #   fixed, only the ring atoms in between move.  Pure generator.
    # ⛔ Default OFF -> ring set unchanged -> byte-identical.
    _mc = False
    if _os.environ.get("DELFIN_FFFREE_PUCKER_MC", "0") == "1":
        try:
            from delfin.manta import _elements as _EL
            _mc = any(_EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol()) for i in ring)
        except Exception:
            _mc = False
        if _mc:
            _nonmetal = [int(i) for i in ring
                         if not _EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol())]
            if _nonmetal and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in _nonmetal):
                return False          # fully aromatic metallacycle: rigid
    n_sat = 0
    for idx in ring:
        a = mol.GetAtomWithIdx(int(idx))
        if a.GetIsAromatic():
            if not _mc:
                return False
            continue
        hyb = a.GetHybridization()
        if hyb == Chem.HybridizationType.SP3:
            n_sat += 1
            continue
        if hyb in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):
            continue
        # unspecified (perceived mol) -> geometric sp3 test
        if P is None:
            continue
        nbrs = [nb.GetIdx() for nb in a.GetNeighbors()]
        if len(nbrs) >= 4:
            n_sat += 1
        elif len(nbrs) == 3:
            c = P[int(idx)]
            q0, q1, q2 = P[nbrs[0]], P[nbrs[1]], P[nbrs[2]]
            nrm = _np.cross(q1 - q0, q2 - q0)
            ln = float(_np.linalg.norm(nrm))
            if ln > 1e-9 and abs(float(_np.dot(c - q0, nrm / ln))) > 0.25:
                n_sat += 1   # pyramidal -> sp3-like
    # ⚠ In the metallacycle the METAL provides the hinge, not an sp3 centre --
    #   there the threshold 3 is an organic criterion applied to the wrong object.
    return n_sat >= (1 if _mc else 3)


def _ring_order(mol, ring_set):
    """Return the ring atoms in connectivity (traversal) order."""
    ring = list(ring_set)
    adj = {i: [] for i in ring}
    rs = set(ring)
    for i in ring:
        for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            j = nb.GetIdx()
            if j in rs:
                adj[i].append(j)
    order = [ring[0]]
    prev = None
    cur = ring[0]
    for _ in range(len(ring) - 1):
        nxts = [x for x in adj[cur] if x != prev]
        if not nxts:
            return ring  # fall back to arbitrary order
        nxt = nxts[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def _conf_to_xyz(mol) -> str:
    conf = mol.GetConformer()
    out = []
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i)
        p = conf.GetAtomPosition(i)
        out.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(out) + "\n"


def _tfd(acc_mol, id_a: int, id_b: int) -> float:
    """Torsion-Fingerprint-Deviation between two conformers of ``acc_mol``.
    TFD is the field-standard conformer discriminator — it compares ALL ring +
    rotatable-bond torsions with the molecule's topological symmetry folded in,
    so pucker families (chair vs twist-boat), rotamers and axial/equatorial
    substituents separate cleanly where heavy-atom RMSD conflates them."""
    try:
        from rdkit.Chem import TorsionFingerprints as _TF
        return float(_TF.GetTFDBetweenConformers(acc_mol, [id_a], [id_b])[0])
    except Exception:
        return 1.0     # no torsions / failure -> treat as distinct (keep)


# ===== THE RING-LOCAL TFD (26.08.2026) =============================================
#
# MEASURED (`selbsttest_trennschaerfe`, dilution series): the same cyclohexane ring on
# a growing rigid acene, states PER RING --
#     cyclohexylbenzene      ring share 50.0 %   TFD 0.05 -> 11   TFD 0.005 -> 62
#     cyclohexylnaphthalene             37.5 %                6                49
#     cyclohexylanthracene              30.0 %                1                23
#     cyclohexyltetracene               25.0 %                1                10
# The folds are THERE -- at 0.005 they come back.  TFD merges them.
#
# THE CAUSE IS IN RDKIT'S OWN FORMULA, and it is SHARPER than "averages".
# `CalculateTFD` forms sum(d_i * w_i) / sum(w_i) over ALL torsions of the molecule.
# If only the one ring moves, d_i = 0 for every other torsion, and what remains is
#     TFD_global = d_ring * w_ring / sum(w)
# `CalculateTorsionWeights` sets w = exp(-beta * d^2) with d = topological distance to the
# MOST CENTRAL bond of the molecule.  An appended cyclohexyl ring slides further to the
# edge with every additional acene ring -- its weight share falls EXPONENTIALLY, not
# like 1/N.  Measured on the same four MMFF-optimised samples (w_ring / sum(w)):
#     cyclohexylbenzene      0.1875   ->   5.3-fold dilution
#     cyclohexylnaphthalene  0.0548   ->  18.3-fold
#     cyclohexylanthracene   0.0166   ->  60.2-fold
#     cyclohexyltetracene    0.0069   -> 146.0-fold
# 0.05 / 146 = 0.00034 is thus the threshold the ring in the tetracene EFFECTIVELY sees.
# That is exactly why the states only come back at 0.005, and exactly why the
# number falls monotonically with the scaffold size.
#
# ⇒ Word for word the RMSD error one level up -- and the effect GROWS with the
#   ligand size, so it hits hardest the systems that matter.
#
# THE REPAIR TAKES RDKIT'S OWN ROUTE, no re-implementation: `CalculateTorsionLists` returns
# non-ring and ring torsions SEPARATELY, `CalculateTorsionAngles` and
# `CalculateTFD` accept exactly such lists.  So it is only FILTERED: the
# entry of the ring under consideration stays, everything else drops out.
#
# ⚠ THE SYMMETRY FOLDING SURVIVES -- the point at which the attempt of 26.08. died,
#   to REPLACE TFD by a pure CP distance (n=5 went from 3,3,3 to 9,13,14,
#   because phi depends on the atom numbering).  Nothing is re-implemented here: RDKit's
#   ring entry is the MEAN of |torsion| over the whole ring, i.e. a number
#   that is invariant under rotation AND reflection of the ring numbering.  This
#   invariance IS the folding.  Measured: unsubstituted five-ring yields ring-locally
#   3 states -- exactly as globally.
#
# ⚠ THE THRESHOLD STAYS 0.05, and that is not a setting but an IDENTITY.  An
#   unsubstituted single-ring molecule has NO non-ring torsion and EXACTLY ONE ring entry;
#   sum(w) is then w_ring, the fraction cancels, and global TFD = ring-local TFD on
#   every conformer pair.  MEASURED, denominator 760 conformer pairs (n=5,6,7,8 with 190 each):
#   largest difference 5.6e-17 -- that is floating-point noise, not a small
#   difference.  On the calibration sample the two measures are thus not similarly
#   calibrated but THE SAME; the meaning of 0.05 changes there by exactly zero.
#   The state counts confirm it line by line (3/9/11/16 at 0.05 and
#   6/22/64/103 at 0.005, global as ring-local).  `selbsttest_tfd_lokal` measures both.
#
# ⚠ WHAT THIS VERSION PAYS FOR IT, and it is stated here so nobody finds it later as a
#   surprise: RDKit's ring entry is ONE number per ring.  The same
#   invariance that folds the symmetry makes the measure one-dimensional -- two
#   genuinely different folds with the same mean of |torsion| are
#   merged.  The fold axis itself (theta, phi) is NOT seen by the measure.
#   ⇒ Ring-local TFD is a sharper dedup, NOT a complete
#     fold descriptor.  Whoever needs the axis needs `_cp_theta_phi` in addition -- and
#     that does NOT fold the symmetry (measurement 26.08.), so it is fit only as a SECOND
#     instrument alongside this one, never as a replacement.
# ⚠ WITHOUT WEIGHTS when several rings are selected.  The weights ARE the
#   dilution mechanism (distance to the most central bond); bringing them back in
#   ring-locally would bring back the effect this version removes.  For ONE ring
#   it is irrelevant anyway -- a weight cancels against itself.
# ⛔ Default OFF -> `_tfd_distinct` runs the old line -> byte-identical.


def _tfd_lokal_listen(mol, ringe):
    """RDKit's ring torsion list, FILTERED to the rings passed in.

    ⚠ ASSIGNMENT IS BY ATOM SET, not by list index.
      `tors_list_rings` comes from `Chem.GetSymmSSSR`, the caller's rings from
      `RingInfo.AtomRings()`.  Both deliver the same rings -- but RELYING on their
      index equality would be an assumption, and a wrongly
      assigned ring would not be recognisable here as an error, only as
      "the other ring simply did not move".  The caller therefore checks the
      LENGTH of the return value against the number of requested rings.
    """
    from rdkit.Chem import TorsionFingerprints as _TF
    _tl, _tlr = _TF.CalculateTorsionLists(mol)
    ziel = {frozenset(int(a) for a in r) for r in ringe}
    # Ring entry k consists of the N consecutive quadruples of the ring; the
    # FIRST atoms of these groups are exactly the N ring atoms (RDKit builds them so).
    return [(q, d) for q, d in _tlr if frozenset(int(t[0]) for t in q) in ziel]


def _tfd_lokal(acc_mol, listen, id_a: int, id_b: int) -> float:
    """TFD over ONLY the torsion entries passed in -- no scaffold in the denominator."""
    from rdkit.Chem import TorsionFingerprints as _TF
    t_a = _TF.CalculateTorsionAngles(acc_mol, [], listen, confId=id_a)
    t_b = _TF.CalculateTorsionAngles(acc_mol, [], listen, confId=id_b)
    return float(_TF.CalculateTFD(t_a, t_b, weights=None))


# ⚠ TWO SWITCHES FOR ONE MEASURE, and that is not knob proliferation.  The same measure
#   moves the number at the two call sites in OPPOSITE directions:
#     `_ring_pucker_states` (per ring)         cyclohexyltetracene  1 -> 11 states
#     `generate` (per combination)             decalin             22 ->  8 frames
#   If both hang on ONE switch, an A/B measures their SUM and nobody can say
#   which share came from where -- exactly the design on which verdicts in this
#   project have already failed.  Switched separately they are two measurements.
# ⚠ THE FINDING HANGS ON THE FIRST.  The dilution was measured on the states PER
#   RING; the combination level is an EXTRAPOLATION from that and therefore stands under
#   its own switch, likewise switched off.
def _tfd_distinct(acc_mol, cid: int, kept_ids, thr: float, ringe=None,
                  schalter: str = "DELFIN_FFFREE_PUCKER_TFD_LOCAL") -> bool:
    if ringe and _os.environ.get(schalter, "0") == "1":
        try:
            _listen = _tfd_lokal_listen(acc_mol, ringe)
            # ⚠ ALL OR NONE.  If the assignment recovers only PART of the rings,
            #   the ring-local comparison silently measures fewer rings
            #   than the caller meant -- and the missing ones would take place
            #   nowhere.  Half a measurement looks like a whole one from outside; that is
            #   exactly the design that has already passed as a finding several times
            #   in this project.  Better to fall back entirely to the global measure.
            if len(_listen) == len({frozenset(int(a) for a in r) for r in ringe}):
                # ⚠ OWN THRESHOLD ONLY IF SOMEONE SETS IT.  The measurement says: on
                #   the unsubstituted ring both measures are identical, so 0.05 keeps
                #   its meaning.  The knob is there for RE-MEASURING, not for
                #   re-tuning -- empty means "unchanged".
                _s = _os.environ.get("DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR", "")
                _thr = float(_s) if _s.strip() else thr
                # The list ONCE per candidate, not per pair: `GetTFDBetweenConformers`
                # rebuilds it in the global path on EVERY call -- ring-local is thus
                # also cheaper, not only sharper.
                return all(_tfd_lokal(acc_mol, _listen, k, cid) >= _thr
                           for k in kept_ids)
        except Exception:
            pass          # fall back to the global measure -- never silently empty
    return all(_tfd(acc_mol, k, cid) >= thr for k in kept_ids)


# ===== DISTINGUISHABILITY IS A MAXIMUM, NOT A MEAN (26.08.2026) =====================
#
# The error of RMSD is NOT that it measures geometry.  It is that it AVERAGES --
# and a ring fold is a LOCAL event in a large molecule.  Worked through with our
# own numbers on a folded six-ring:
#
#     ring atoms move 0.203 A in the median (largest single displacement 0.33)
#     ring share of the heavy atoms 13.2 %
#     ⇒ total RMSD = sqrt(0.132) * 0.203 = 0.086 A
#
# 0.086 lies BELOW every dedup threshold this project carries (`_DEDUP_RMSD`
# 0.30 · `rmsd_dedup` 0.5).  So the fold does not vanish because it is small,
# but because it is divided by 87 % unmoved atoms.  The largest displacement
# after Kabsch alignment stays at 0.33 -- FACTOR 4 between the two numbers, measured on
# the same geometry.
#
# WHAT A CRYSTALLOGRAPHER READS INSTEAD.  In the difference Fourier map the
# LARGEST unmodelled deviation appears as the residual-density maximum; the mean over all
# atoms does not occur in it.  Two models whose largest atomic displacement lies below the
# resolution would be INDISTINGUISHABLE on the same data -- they are ONE
# entry in the manifold, not two.
#
# THE THRESHOLD, justified instead of set:
#   * Coordinate esd of a routine structure lies at 0.002 to 0.01 A.  That is the
#     LOWER bound -- below it every difference is refinement noise.
#   * Disorder is only modelled as TWO sites from about 0.3 to 0.5 A on.
#     That is the UPPER bound -- above it the crystallographer sees two conformers.
#   0.15 A lies in between: factor 15 to 75 above the esd, factor 2 to 3 below the
#   disorder limit.  Large enough not to count noise; small enough not to
#   merge anything a crystallographer would still model separately.
# ⚠ IT IS AN ENV PARAMETER because it is a CONVENTION and not a natural constant -- the
#   resolution depends on the dataset, and whoever shifts it should be able to measure that.
#
# ⚠ WHAT THIS METRIC DOES **NOT** DO: it does not rank.  It says "indistinguishable"
#   or "distinguishable", never "better".  A rank would need an energy; that is
#   stage (3) and is OFF for good reason.


def _kabsch_max_rmsd(A, B) -> Tuple[float, float]:
    """(LARGEST displacement, RMSD) of two point sets after Kabsch alignment.

    Both numbers from THE SAME alignment -- otherwise the comparison of the two
    measures would be none.  Kabsch minimises the RMSD; the maximum is thus measured against
    the superposition MOST FAVOURABLE for RMSD and is therefore rather too small than too large.

    ⚠ HEAVY ATOMS ONLY, and that is not economising.  X-ray diffraction sees ELECTRON
      DENSITY; a hydrogen carries one electron and in a routine structure is
      CALCULATED, not found.  A metric that pretends to distinguish H positions
      judges something that is not in the data at all.  The caller therefore passes
      already filtered coordinates.

    ⚠ NO SYMMETRY FOLDING, deliberately.  Both point sets come from THE SAME
      molecule instance in THE SAME atom order; the degeneracy "two numberings
      of the same conformer" is handled in `generate` by the TFD before it, and TFD can do that
      because it folds in the topological symmetry (the CP distance could not -- the
      attempt of 26.08. failed on exactly that).  Whoever additionally minimised over
      automorphisms here would pay N! and measure the same thing.
    """
    v = _kabsch_abweichungen(A, B)
    return float(v.max()), float(_np.sqrt(float((v ** 2).mean())))


def _kabsch_abweichungen(A, B):
    """The deviation PER ATOM after Kabsch alignment -- the common raw material.

    ⚠ ONE alignment, then both measures from it.  If maximum and RMSD were each
      aligned individually, they would compare two different superpositions and
      the factor between them would be partly instrument, partly alignment.  Kabsch
      minimises the RMSD -- the maximum is thus measured against the superposition
      MOST FAVOURABLE for the opponent and is rather too small than too large.
    """
    Am = A.mean(0)
    Bm = B.mean(0)
    Ac = A - Am
    Bc = B - Bm
    try:
        U, _S, Vt = _np.linalg.svd(Bc.T @ Ac)
        d = 1.0 if float(_np.linalg.det(Vt.T @ U.T)) > 0.0 else -1.0
        R = Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T
        Bd = Bc @ R.T
    except Exception:
        Bd = Bc
    return _np.linalg.norm(Ac - Bd, axis=1)


def _add_conf(acc_mol, src_mol) -> int:
    return acc_mol.AddConformer(Chem.Conformer(src_mol.GetConformer()), assignId=True)


def _ring_pucker_states(mol_with_conf, ring, frozen: Set[int],
                        tfd_thr: float) -> List[Optional[Tuple[Optional[float], float]]]:
    """DISTINCT pucker settings for ONE ring, TFD-deduped.  ``None`` = the base
    pucker; each other entry is a ``(theta, phi)`` SETTING that, after a torsion-
    held relax, gives a conformer whose torsion fingerprint differs from every
    kept one (cyclohexane -> {base chair, the twist-boat(s)}, not 7 relabelled
    pseudorotation copies)."""
    states: List[Optional[Tuple[float, Optional[float], float]]] = [None]
    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    n = len(ring)
    # ===== DEDUP IN CP INSTEAD OF IN TFD (26.08.2026) ================================
    #
    # MEASURED on our own convergence test.  The TFD threshold 0.05 over-splits:
    #     n=6, NPHASE=16 -> 10 states, smallest pairwise CP distance  5.0 degrees
    #     n=7, NPHASE=16 -> 12 states, smallest pairwise CP distance  1.8 degrees
    # Two folds that lie 1.8 degrees apart are THE SAME.  And two
    # entries stood at theta = 180 with phi = 136 and phi = 339 -- at the pole phi is
    # meaningless, so provably the same state.
    #
    # ⚠️ THIS IS NOT COSMETICS.  The state count PER RING is the base of the
    #    cross product over all rings (measured: 4.3 rings per system):
    #         3 states, 4 rings ->      81 combinations   computable
    #        10 states, 4 rings ->  10 000                not computable
    #    Over-splitting makes the COMPLETE combinatorics unaffordable.  Whoever wants the
    #    whole pucker space must first stop counting noise as a basin
    #    -- otherwise the cap comes back through the back door.
    #
    # Chemical yardstick: cyclohexane has chair + twist-boat family, after
    # symmetry folding 2-3 classes.  The five-ring converges to 3 by itself.
    #
    # Dedup therefore happens on the SPHERE, with the great-circle distance -- the
    # pole degeneracy resolves itself there (see `_cp_abstand`).
    # ⛔ Default OFF -> TFD as before -> byte-identical.
    _cpd = _os.environ.get("DELFIN_FFFREE_PUCKER_CPDEDUP", "0") == "1"
    _cp_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_CPTOL", "15") or 15.0)
    _cp_qtol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_CPQTOL", "0.15") or 0.15)
    _cp_kept: List[Tuple[float, float, float]] = []
    # ===== (2) THE SAME INDISTINGUISHABILITY, BUT PER RING (26.08.2026) ==============
    #
    # ⚠ THE LEVER IS HERE, NOT IN THE CROSS PRODUCT.  The state count PER RING enters
    #   the combinatorics as a POWER (measured 4.3 rings per system): one state
    #   fewer per ring saves more than any rule further down, because down there the
    #   relax has already happened.  A gate behind the relax kills the result, not the
    #   cost -- that has been in the docstring of `selbsttest_kombinatorik` since 26.08.
    #   and applies to the defect filter just as to the indistinguishability.
    # ⚠ WHETHER IT REALLY REDUCES IS A MEASUREMENT AND NOT A HOPE.  The maximum is
    #   insensitive to dilution (a maximum knows no denominator), so nothing
    #   suggests it merges genuine ring basins -- those lie far
    #   apart.  `selbsttest_trennschaerfe` step 3 counts it.
    # ⛔ Default OFF -> byte-identical.
    _xrd_r = _os.environ.get("DELFIN_FFFREE_PUCKER_XRD", "0") == "1"
    _xrd_r_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_XRDTOL", "0.15") or 0.15)
    _schwer_r: List[int] = []
    _xrd_r_kept: List = []
    if _xrd_r:
        try:
            _schwer_r = [i for i in range(mol_with_conf.GetNumAtoms())
                         if mol_with_conf.GetAtomWithIdx(i).GetSymbol() != "H"]
            _xrd_r_kept = [mol_with_conf.GetConformer().GetPositions()[_schwer_r]]
        except Exception:
            _xrd_r = False
    # ===== THE WHOLE PUCKER SPACE INSTEAD OF THREE SPOTS ON IT (26.08.2026) =========
    #
    # The old candidate list samples the Cremer-Pople sphere at exactly three places:
    # the EQUATOR (theta = 90, K phases), and for even rings the two POLES.
    # `q_scale` is constant 1.0 throughout.
    #   ⇒ theta between 0 and 90 is NEVER sampled -- that is where half-chair
    #     (theta ~50) and envelope (theta ~55) lie.
    #   ⇒ the amplitude is NEVER varied -- only ONE spherical shell.
    #   ⇒ odd rings get `theta=None`, i.e. pure pseudorotation.
    # The missing forms are thus not "not implemented" but NOT
    # SAMPLED -- a difference that makes the repair cheap.
    #
    # With `DELFIN_FFFREE_PUCKER_SPACE=1` the (N-3)-dimensional space is instead
    # covered systematically (`_pucker_space_grid`), for EVERY ring size and
    # via `_set_pucker_general`, which also knows the pairs m >= 3 -- without those
    # every ring from N = 7 on was incompletely parametrised.
    # ⚠ Completeness here is a question of RESOLUTION: the space is continuous.
    #   It is covered completely at the stated resolution, and that is
    #   in the trace -- no corner left out, no direction preferred.
    # ⛔ Default OFF -> old list -> byte-identical.
    _raum = _os.environ.get("DELFIN_FFFREE_PUCKER_SPACE", "0") == "1"
    if _raum:
        _namp = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_NAMP", "2") or 2))
        _nph = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_NPHASE", "8") or 8))
        _cands = _pucker_space_grid(n, _namp, _nph)
        if _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
            # ⚠ THE SAMPLING DENSITY GOES INTO THE LOG TOO.  Without it a
            # reduced resolution could not later be told from a full one --
            # and exactly that would be a silent cap.
            _res = getattr(_pucker_space_grid, "_grid_res", None)
            print("[pucker] RAUM n=%d: %d Kandidaten (%d-dim, Amplitudenstufen %d, "
                  "Phasen %d)" % (n, len(_cands), max(0, n - 3), _namp, _nph))
            if isinstance(_res, dict) and _res.get("n") == n and _res.get("reduziert"):
                print("[pucker] RAUM n=%d: Dichte REDUZIERT auf Budget %d -- Phasen je m %s"
                      % (n, _res.get("budget"), _res.get("n_phase_je_m")))
    else:
        _cands = _pucker_candidates(n)
    if frozenset(ring) in _FLAT_ONLY:
        # ONLY straighten, do not fold -- see the block in `generate`.
        _cands = [({}, {})] if _raum else [(0.0, 0.0, 0.0)]
    for _cand in _cands:
        try:
            m2 = Chem.Mol(mol_with_conf)
            if _raum:
                _qs, _phis = _cand
                _set_pucker_general(m2.GetConformer(), ring, _qs, _phis, frozen)
                theta = phi = None
            else:
                _qs, theta, phi = _cand
                _set_pucker(m2.GetConformer(), ring, _qs * _amp(n), theta, phi, frozen)
            if not _relax_hold_pucker(m2, ring, frozen):
                continue
            _Pr = None
            if _xrd_r:
                # Indistinguishable from the ground state OR from an already kept
                # state -> the same entry, not a second ring state.
                _Pr = m2.GetConformer().GetPositions()[_schwer_r]
                if any(_kabsch_max_rmsd(_Pa, _Pr)[0] < _xrd_r_tol
                       for _Pa in _xrd_r_kept):
                    continue
            if _cpd:
                # Measure AFTER the relax, do not compare the TARGET values: the relax
                # pulls the starting point into the nearest genuine basin, and exactly its
                # location decides whether it is a new one.
                _cp = _cp_theta_phi(m2.GetConformer().GetPositions(), ring)
                if any(_cp_abstand(_cp, _k) < _cp_tol and abs(_cp[0] - _k[0]) < _cp_qtol
                       for _k in _cp_kept):
                    continue
                _cp_kept.append(_cp)
                if _xrd_r and _Pr is not None:
                    _xrd_r_kept.append(_Pr)
                states.append(_cand if _raum else (_qs, theta, phi))
                continue
            cid = _add_conf(acc, m2)
            # ⚠ HERE IS THE MEASURING POINT OF THE FINDING.  `ringe` names the ONE ring
            #   that is folded here; with the switch ON only its torsion
            #   counts, the scaffold is no longer in the denominator.  Switch OFF -> the
            #   argument is not even looked at in `_tfd_distinct`.
            if _tfd_distinct(acc, cid, kept_ids, tfd_thr, ringe=(ring,)):
                kept_ids.append(cid)
                if _xrd_r and _Pr is not None:
                    _xrd_r_kept.append(_Pr)
                # In space mode the state is the coordinate pair itself; the
                # legacy form stays a 3-tuple.  `generate` only indexes, it does not
                # read the content -- both forms are equivalent there.
                states.append(_cand if _raum else (_qs, theta, phi))
            else:
                acc.RemoveConformer(cid)
        except Exception:
            continue
    return states


def _neuer_zaehler() -> dict:
    """Fresh counter set for ``generate(..., _zaehler=...)``.

    ⚠ WHY THE MEASURING POINT SITS IN `generate` AND NOT IN A COPY.  The question
      of how large the cross product is AFTER the physics can only be answered on the
      code that actually builds the frames.  A re-implemented loop measures
      the re-implementation -- in this project exactly that has already passed as a finding
      several times and was none.  The price is an `if _zaehler is not None`
      at six places; the default path (`_zaehler is None`) runs unchanged.
    """
    return {"ringgroessen": [], "zustaende_je_ring": [], "kreuzprodukt": 0,
            "gemeinsame_atome": 0, "gem_max": 0, "aufzaehlung": 0, "gebaut": 0,
            "relax_fehler": 0, "kollision": 0, "winkel": 0, "tor_ueberlebt": 0,
            "tfd_doppelt": 0, "ausnahme": 0, "energien": [],
            # (1) defect filter and (2) crystallographic indistinguishability get
            # counters OF THEIR OWN.  Lumping them into `kollision`/`tfd_doppelt` would be
            # exactly the mistake this project has already made three times: a
            # detector name that covers two mechanisms is not a measurement.
            "bindung": 0, "xrd_doppelt": 0}


def generate(mol_with_conf, frozen: Optional[Set[int]] = None,
             budget: int = 64, tfd_thr: float = 0.05,
             angle_skip: Optional[Set[int]] = None,
             _zaehler: Optional[dict] = None) -> List[Tuple[str, str]]:
    """Construct the COMBINATORIAL ring-pucker conformers from a base conformer.

    ``_zaehler``: optional counter set (`_neuer_zaehler()`).  If set, this run
    records how many combinations were enumerated, built, rejected at the collision or
    angle gate and merged by TFD -- the measurement that
    `selbsttest_kombinatorik` evaluates.  ``None`` (default) = not a single counter
    is touched, the build is byte-identical to before.

    ``mol_with_conf`` carries ONE embedded conformer (a chain/rotamer pose whose
    rings sit at their base pucker).  Every puckerable ring's distinct pucker
    states are enumerated, then the CARTESIAN PRODUCT across all rings is
    constructed (Cy3P: 3 rings x {chair, twist-boat, ...} -> 3xchair, 2xchair+
    twist, ...).  Each combination sets all rings' puckers simultaneously and is
    relaxed with EVERY ring's pucker HELD but all inter-ring bonds/torsions FREE,
    so the free degrees of freedom relieve any inter-ring steric clash while the
    puckers themselves survive.  A whole-molecule clash gate then drops any
    combination that stayed sterically unrealistic, and TFD dedup keeps one
    representative per distinct torsion fingerprint (so Cy3P's three equivalent
    rings collapse correctly).  ``frozen`` fixes metal + donor atoms so metal
    chelate rings pucker without disturbing the coordination sphere.  Returns the
    NEW distinct, clash-free conformers ``[(xyz, label), ...]``; never raises.
    """
    if not (_RDKIT and _np is not None):
        return []
    try:
        if mol_with_conf.GetNumConformers() == 0:
            return []
        ri = mol_with_conf.GetRingInfo()
        rings_raw = [set(r) for r in ri.AtomRings()]
    except Exception:
        return []
    frozen = frozen or set()
    rings = [_ring_order(mol_with_conf, r) for r in rings_raw
             if _is_puckerable(mol_with_conf, r)]
    # ===== THE FLAT STATE FOR CONJUGATED RINGS (18.08.2026) ===========================
    # MEASURED (16./17.08., `folds`, 965 systems): of 326 missing ring motifs,
    # **218 are PLANAR** -- 5M 126 - 6M 55 - 4M 37 -- and ALL three classes are
    # METALLACYCLES (`find_conformer_completeness:254` builds the name as
    # f"{sz}{'M' if is_metallacycle else ''}:{basin}").
    #
    # `_is_puckerable` does not let exactly these in: it demands "no aromatic
    # ring atom" and ">= 3 sp3 ring atoms", and justifies that by claiming a conjugated ring
    # IS planar and rigid anyway.  The eye measures the opposite: the crystal
    # planar state is MISSING in the build.  Both together mean -- the ring is folded by
    # something else, and the only module that could deliberately set it flat
    # is not allowed to touch it.  `planar138` and `pktrace` confirmed this
    # (affected 0, even with the reach barrier deliberately bypassed).
    #
    # ⚠ ONLY THE Q=0 STATE, no pucker.  These rings are not to be folded but
    # STRAIGHTENED; `_ring_pucker_states` therefore offers them exclusively the
    # projection onto the mean plane.  Metal and donors are in `frozen` and
    # do not move -- the coordination sphere remains untouched.
    #
    # ⚠ CLASS: ENUMERATOR, not a repairer (module census 18.08.).  It PROJECTS with
    # a frozen core instead of generating anew -- the same class as the
    # mirror closure (+1.0 pp), not that of BACKBONE_REEMBED (+11.9 pp).
    #
    # Default OFF -> ring set unchanged -> byte-identical.
    # ⚠️ HERE, NOT IN _set_pucker -- a mistake of mine, measured on 18.08.
    # `_FLAT_ONLY.clear()` stood in `_set_pucker`, i.e. in the function that runs
    # once per RING.  Consequence: ring 0 was correctly held flat, after that the
    # set was empty, and EVERY further ring got the full fold ladder at
    # Q = 0.63 Angstrom.  So the mechanism FOLDED the conjugated chelate ring
    # instead of STRAIGHTENING it -- the exact opposite of its purpose.
    #
    # The proof was in the labels: 246 of 256 pucker labels begin with
    # `r0:base` (only ring 0 stayed flat), and the state indices on rings from 1 on
    # run up to 15 -- a flat-only ring can have at most index 1, and 15 is
    # exactly the full candidate count of a six-ring.
    # Chemically measured on four salicylaldiminato chelates: Walsh angle at the
    # azomethine carbon 13.1 to 17.5 degrees, while the crystal holds the same
    # centres at no more than 1.4 degrees; in one frame a bond broke.
    #
    # ⇒ The gate term `smiles_ccdc_regressed`, on which `planarA` was blocked, was
    # RIGHT.  Loosening it would have cemented the build error in.
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        _FLAT_ONLY.clear()
        _have = {frozenset(r) for r in rings}
        for _r in rings_raw:
            if frozenset(_r) in _have or not (5 <= len(_r) <= 8):
                continue
            try:
                _ro = _ring_order(mol_with_conf, _r)
            except Exception:
                continue
            if _ro:
                _FLAT_ONLY.add(frozenset(_ro))
                rings.append(_ro)
    if not rings:
        return []
    # per-ring distinct pucker states (index 0 == base pucker for every ring),
    # TFD-deduped so an unsubstituted ring yields only its genuine minima.
    if _zaehler is not None:
        import time as _t
        _zaehler["t0"] = _t.perf_counter()
    per_ring_states = [_ring_pucker_states(mol_with_conf, ring, frozen, tfd_thr)
                       for ring in rings]
    if _zaehler is not None:
        # ⚠ THE TIME IS THE ACTUAL ANSWER to "affordable?".  It splits into
        #   two items that scale completely differently: the states per ring
        #   cost LINEARLY in the ring count, the cross product costs EXPONENTIALLY.  Whoever
        #   measures only the total time does not see the difference.
        _zaehler["t_zustaende"] = _t.perf_counter() - _zaehler["t0"]
        # (a) and (b) of the measurement -- and the COUPLING NUMBER with them.  Shared atoms
        # between two rings are the independent variable of the whole test:
        # fused = 2, spiro = 1, independent = 0.  It is counted here from THE SAME
        # ring list the build uses -- not looked up from the SMILES.
        _zaehler["ringgroessen"] = [len(r) for r in rings]
        _zaehler["zustaende_je_ring"] = [len(s) for s in per_ring_states]
        _p = 1
        for _s in per_ring_states:
            _p *= len(_s)
        _zaehler["kreuzprodukt"] = _p
        # ⚠ SUM AND MAXIMUM ARE TWO DIFFERENT STATEMENTS, and only the MAXIMUM
        #   names the KIND of coupling.  A pair shares 0 atoms (separate), 1 (spiro),
        #   2 (fused, one shared bond) or >= 3 (bridged).  The sum
        #   over all pairs, by contrast, simply grows with the ring count and confuses
        #   three loose rings with a cage.
        _gem, _gmax = 0, 0
        for _i in range(len(rings)):
            for _j in range(_i + 1, len(rings)):
                _n_ij = len(set(rings[_i]) & set(rings[_j]))
                _gem += _n_ij
                _gmax = max(_gmax, _n_ij)
        _zaehler["gemeinsame_atome"] = _gem
        _zaehler["gem_max"] = _gmax

    # cartesian product of state indices, deterministic order, budget-capped;
    # skip the all-base (identity) combination; fewest-changed rings first.
    import itertools as _it
    # ===== SYMMETRY REDUCTION (27.08.2026, DELFIN_FFFREE_PUCKER_SYMM) ===============
    #
    # WHY HERE.  Two rings that lie in ONE orbit under the automorphism group of the
    # molecule are indistinguishable.  Then (s1,s2) is THE SAME STRUCTURE
    # as (s2,s1) -- the cross product generates duplicates there, which the
    # RMSD dedup later throws away again, after they have paid for relax and clash gate.
    #
    # MEASURED 27.08. with TWO independent instruments:
    #   build side     534 systems with >=2 foldable rings -> 371 (69.5 %) with orbit
    #   crystal side  (CCDC clean_v2, DELFIN-independent)
    #                1757 with >=2 folded rings -> 1393 (79.3 %) with orbit
    #
    # THE ACTUAL GAIN IS NOT COMPUTE TIME BUT COVERAGE.  The cap
    # thirty lines further down (`combos[:budget]`) cuts off by fold depth
    # -- with >2 flexible rings the state "all folded simultaneously" is NEVER
    # built.  If the product shrinks below the budget, the cap stops biting,
    # and exactly the deep states arise again.  So the reduction takes
    # duplicates away and gives back GENUINE states in return.
    #
    # ⛔ NEVER-WORSE, AND STRICTLY SO.  Only what a GENUINE automorphism maps
    # onto each other is merged -- not what merely carries the same rank multiset.
    # The difference is measured: of 27 hand-checked ring pairs,
    # 26 were a genuine orbit and ONE was merely rank-equal (3.7 %).  Had I taken the
    # rank set, this one pair would have been merged and a real state
    # WOULD HAVE BEEN MISSING.  With completeness as the north star that is the worse
    # mistake of the two.
    # ⛔ AND IF THE AUTOMORPHISM SEARCH DOES NOT HOLD, there is NO reduction (full
    # product branch).  Fallback is always the larger set, never the smaller.
    # ⛔ Default OFF -> byte-identical.
    _symm = _os.environ.get("DELFIN_FFFREE_PUCKER_SYMM", "0") == "1"
    _bahnen = None
    if _symm and len(rings) > 1:
        _bahnen = _ring_bahnen(mol_with_conf, rings)
    if _bahnen and any(len(b) > 1 for b in _bahnen):
        # Per orbit MULTISETS instead of ordered tuples: n^k becomes C(n+k-1, k).
        # The state lists of an orbit have equal length (same ring size,
        # same environment) -- checked, otherwise the orbit falls back to the product.
        _pro_bahn = []
        for _b in _bahnen:
            _n = len(per_ring_states[_b[0]])
            if any(len(per_ring_states[i]) != _n for i in _b):
                _pro_bahn.append(list(_it.product(*[range(len(per_ring_states[i]))
                                                    for i in _b])))
            else:
                _pro_bahn.append(list(_it.combinations_with_replacement(
                    range(_n), len(_b))))
        combos = []
        for _wahl in _it.product(*_pro_bahn):
            _c = [0] * len(rings)
            for _b, _w in zip(_bahnen, _wahl):
                for _ri, _st in zip(_b, _w):
                    _c[_ri] = _st
            if any(_c):
                combos.append(tuple(_c))
        if _zaehler is not None:
            _zaehler["symm_bahnen"] = [len(b) for b in _bahnen]
            _zaehler["symm_voll"] = _prod_laenge(per_ring_states)
            _zaehler["symm_reduziert"] = len(combos)
        if _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
            print("[pucker] SYMMETRIE: %d Ringe in %d Bahnen %s -- "
                  "Kombinationen %d statt %d"
                  % (len(rings), len(_bahnen), [len(b) for b in _bahnen],
                     len(combos), _prod_laenge(per_ring_states) - 1))
    else:
        combos = [c for c in _it.product(*[range(len(s)) for s in per_ring_states])
                  if any(c)]
    combos.sort(key=lambda c: (sum(1 for x in c if x), c))
    # ===== THE TRUNCATION WAS THE DAMAGE, NOT THE CAP (26.08.2026) ==================
    #
    # `combos[:budget]` takes the FIRST 48 of a list sorted by "how many rings
    # deviate from the ground state".  That means: first ALL
    # one-ring changes, then all two-ring combinations -- and then it stops.
    # The state in which ALL rings are folded simultaneously is NEVER built
    # with more than two flexible rings.  Four rings with three basins each are 80
    # combinations; 48 of them cover depth 1 and 2, depth 3 and 4 drop out.
    #
    # MEASURED 26.08. on 400 systems -- the distribution is exactly that of a
    # truncated product:
    #     realised conformers per system: 1:25 · 2:57 · 3:26 · 4:7 · 5:4 · 6:2 · 7:1 · 9:1
    #     82 of 123 systems have ONE OR TWO conformers, at 2.18 basins per ring.
    #
    # ⇒ With `DELFIN_FFFREE_PUCKER_FULL=1` there is NO truncation: the complete
    #   product is built, every fold depth arises.  The cap stays as the
    #   default (byte-identical), but from now on it is NEVER SILENT -- what it
    #   throws away is in the trace.  A silently capped coverage reads
    #   like completeness and is none.
    # ⚠ PRICE, honestly: the product grows exponentially with the ring count (8 rings with
    #   3 basins each = 6561 combinations, each with relax and clash gate).  That is why
    #   the full version is behind a switch and not in the default.
    _n_voll = len(combos)
    # ══ ADD, NEVER REPLACE -- AT THE CAP, NOT AT THE REDUCTION ═════════════════
    # (04.09.2026, user decision.)
    #
    # FINDING (register #305/#315).  `symmfold6k` blocks at RONSIW with
    # `broken_regressed`: 12 frames before, 12 after -- and yet ONE new.
    #     off   SP-4-chelate-1-pucker r0:base+r1:1+r2:3+r3:base
    #     on    SP-4-chelate-1-pucker r0:3+r1:base+r2:4+r3:base
    #
    # ⛔ MY FIRST EXPLANATION WAS WRONG.  I took this for a
    #    replacement decision of the symmetry reduction.  It is a CAP
    #    consequence: the reduction frees up places under `combos[:budget]`, so
    #    OTHER combinations slide in.  The reduction replaces nothing --
    #    the cap does.  That is why the guard sits HERE and not there.
    #
    # THE CONSEQUENCE.  The champion set is not a subset of the reduced set:
    # a state the champion built drops out, and if the one that slid in
    # embeds worse, `broken_frac` rises (RONSIW
    # 0.250 -> 0.333).
    #
    # THE GUARD.  First the set the champion would have built, then filled up
    # with the reduced one.  Thus the champion frame set is contained BY
    # CONSTRUCTION and the axis is never-worse.
    # ⚠️ PRICE, honestly: the union becomes up to twice as large as the
    #    cap.  That IS the point -- the reduction is meant to make DEEPER states
    #    reachable, not to displace shallower ones.  Whoever wants to hold the cap
    #    leaves the guard off.
    # ⚠️ Default OFF -> byte-identical.
    if (_symm and _bahnen and any(len(b) > 1 for b in _bahnen)
            and _os.environ.get("DELFIN_FFFREE_PUCKER_SYMM_ADD", "0") == "1"
            and _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") != "1"):
        _champ = [c for c in _it.product(*[range(len(s)) for s in per_ring_states])
                  if any(c)]
        _champ.sort(key=lambda c: (sum(1 for x in c if x), c))
        _champ = _champ[:max(0, int(budget))]
        _gesehen = set(_champ)
        combos = _champ + [c for c in combos[:max(0, int(budget))]
                           if c not in _gesehen]
        if _zaehler is not None:
            _zaehler["symm_add_champion"] = len(_champ)
            _zaehler["symm_add_gesamt"] = len(combos)
    elif _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") == "1":
        pass                                    # all folds, no cap
    else:
        combos = combos[:max(0, int(budget))]
    if _zaehler is not None:
        # ⚠ TWO DIFFERENT NUMBERS, and confusing them would be the whole error.
        #   `aufzaehlung` is the cross product WITHOUT the ground state -- what WOULD HAVE
        #   to be enumerated.  `gebaut` is what really goes through relax and
        #   gate after the cap.  Only the second number costs compute time, only the first is the
        #   completeness question.
        _zaehler["aufzaehlung"] = _n_voll
        _zaehler["gebaut"] = len(combos)
    if len(combos) < _n_voll and _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
        print("[pucker] KOMBINATIONEN GEKAPPT: %d von %d gebaut, %d verworfen "
              "(%d Ringe, Mulden %s) -- DELFIN_FFFREE_PUCKER_FULL=1 baut alle"
              % (len(combos), _n_voll, _n_voll - len(combos), len(rings),
                 "x".join(str(len(s)) for s in per_ring_states)))

    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    out: List[Tuple[str, str]] = []
    # ===== (1) DEFECT FILTER: THE THIRD GATE, AS A SET DIFFERENCE (26.08.2026) =======
    #
    # Collision and angle are already below.  What is missing is the BOND LENGTH --
    # see `_bindungs_ausreisser`: the self-gate takes a stretched bond for an
    # absent one and never reports it.
    # ⚠ A RANK WOULD BE WRONG HERE.  A defect is not a "worse", it is a "not
    #   real".  Hence filter, not score -- and hence measured against the
    #   GROUND STATE: only what the fold NEWLY introduces is rejected.  A molecule
    #   that already carries an unusual bond before the fold (nitrile, carbene,
    #   a badly embedded core) thus does not lose all its folds.
    # ⛔ Default OFF -> the gate is never consulted -> byte-identical.
    _defekt = _os.environ.get("DELFIN_FFFREE_PUCKER_DEFEKT", "0") == "1"
    _basis_bind = _bindungs_ausreisser(mol_with_conf) if _defekt else frozenset()
    # ===== (2) CRYSTALLOGRAPHIC INDISTINGUISHABILITY (26.08.2026) ====================
    #
    # Two folds whose LARGEST atomic displacement lies below the resolution are
    # ONE entry in the manifold -- see the block at `_kabsch_max_rmsd`.
    # ⚠ IN ADDITION TO TFD, not instead of it, and that is not caution but a
    #   division of labour: TFD folds in the molecular symmetry and thereby kills the
    #   numbering duplicates; the maximum cannot do that (it has no topology) and
    #   instead kills the sub-threshold duplicates that TFD does not see.  The attempt
    #   of 26.08. to REPLACE TFD by a pure geometric distance failed on exactly
    #   that (n=5: 3,3,3 -> 9,13,14).
    # ⚠ ORDER: the maximum comes BEFORE the TFD because it is cheaper -- Kabsch on
    #   the heavy atoms versus a torsion fingerprint against all kept ones.  The
    #   order changes nothing in the result: what passes BOTH checks is kept,
    #   and the set of kept ones grows identically in both orders.
    # ⛔ Default OFF -> byte-identical.
    _xrd = _os.environ.get("DELFIN_FFFREE_PUCKER_XRD", "0") == "1"
    _xrd_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_XRDTOL", "0.15") or 0.15)
    _schwer: List[int] = []
    _xrd_kept: List = []
    if _xrd:
        try:
            _schwer = [i for i in range(mol_with_conf.GetNumAtoms())
                       if mol_with_conf.GetAtomWithIdx(i).GetSymbol() != "H"]
            _xrd_kept = [mol_with_conf.GetConformer().GetPositions()[_schwer]]
        except Exception:
            _xrd = False
    if _zaehler is not None:
        _zaehler["t0"] = _t.perf_counter()
    for combo in combos:
        try:
            m2 = Chem.Mol(mol_with_conf)
            conf = m2.GetConformer()
            active = False
            for ri_i, st_i in enumerate(combo):
                if st_i == 0:
                    continue
                st = per_ring_states[ri_i][st_i]
                if st is None:
                    continue
                # ⚠ TWO STATE FORMS, and the blind unpacking was a trap.
                #   Legacy: (q_scale, theta, phi) -- three values.
                #   Space:  (qs, phis) -- two mappings m -> value.
                #   A `_qs2, theta, phi = st` on the space form throws ValueError, and
                #   the surrounding `except Exception: continue` would have swallowed that
                #   SILENTLY: every multi-ring combination would have failed without a sound
                #   and the run would have reported "no effect".  Exactly the design
                #   that has already produced a null measurement three times today.
                if len(st) == 2 and isinstance(st[0], dict):
                    _set_pucker_general(conf, rings[ri_i], st[0], st[1], frozen)
                else:
                    _qs2, theta, phi = st
                    _set_pucker(conf, rings[ri_i],
                                _qs2 * _amp(len(rings[ri_i])), theta, phi, frozen)
                active = True
            if not active:
                continue
            # relax holding EVERY ring's pucker; inter-ring bonds/torsions free so
            # a clash between two puckered rings is relieved without collapsing
            # the puckers.
            if not _relax_hold_pucker_multi(m2, rings, frozen):
                if _zaehler is not None:
                    _zaehler["relax_fehler"] = _zaehler.get("relax_fehler", 0) + 1
                continue
            # realism gate: a combination that stayed clashed OR left any VSEPR
            # body distorted (fused/bridged rings strain their shared atoms) is
            # not a physical ensemble member -> drop it.  Everything must be
            # right, or the frame is unrealistic.
            if _zaehler is None:
                if _has_clash(m2) or _has_bad_angles(m2, skip=angle_skip):
                    continue
                if _defekt and (_bindungs_ausreisser(m2) - _basis_bind):
                    continue                    # NEWLY introduced bond defect
            else:
                # ⚠ IN MEASUREMENT MODE BOTH GATES ARE ASKED, in the default path not.
                #   `or` short-circuits: if the collision fires, the angle gate is
                #   NEVER consulted -- the two causes could then not be separated.
                #   Exactly their separation is the whole finding for fused rings:
                #   there two rings share atoms, the angle gate sees the strain at
                #   the fusion centres, and the collision gate is blind to it by
                #   construction (see `_has_bad_angles`).  A second gate call costs
                #   time -- hence only when someone is measuring.
                _kl = _has_clash(m2)
                _wk = _has_bad_angles(m2, skip=angle_skip)
                _bi = bool(_bindungs_ausreisser(m2) - _basis_bind) if _defekt else False
                if _kl:
                    _zaehler["kollision"] = _zaehler.get("kollision", 0) + 1
                if _wk:
                    _zaehler["winkel"] = _zaehler.get("winkel", 0) + 1
                if _bi:
                    _zaehler["bindung"] = _zaehler.get("bindung", 0) + 1
                if _kl or _wk or _bi:
                    continue
                _zaehler["tor_ueberlebt"] = _zaehler.get("tor_ueberlebt", 0) + 1
            # (2) INDISTINGUISHABILITY before the dedup -- see the block above.
            _Pk = None
            if _xrd:
                try:
                    _Pk = m2.GetConformer().GetPositions()[_schwer]
                    if any(_kabsch_max_rmsd(_Pa, _Pk)[0] < _xrd_tol for _Pa in _xrd_kept):
                        if _zaehler is not None:
                            _zaehler["xrd_doppelt"] = _zaehler.get("xrd_doppelt", 0) + 1
                        continue
                except Exception:
                    _Pk = None
            cid = _add_conf(acc, m2)
            # At the combination level SEVERAL rings fold simultaneously -- ring-local
            # here means "all folded rings, but only those".  ⚠ That is a
            # DIFFERENT statement than one level up: here the exocyclic torsion
            # also drops out of the comparison, so two combinations that differ ONLY
            # in it are merged.  That is intended (this
            # module folds rings) and is stated here so nobody discovers it later as
            # a side effect.
            # ⚠ OWN SWITCH, because the effect here has the opposite sign from
            #   one level up -- see the block at `_tfd_distinct`.
            if not _tfd_distinct(acc, cid, kept_ids, tfd_thr, ringe=rings,
                                 schalter="DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"):
                acc.RemoveConformer(cid)
                if _zaehler is not None:
                    _zaehler["tfd_doppelt"] = _zaehler.get("tfd_doppelt", 0) + 1
                continue
            kept_ids.append(cid)
            if _xrd and _Pk is not None:
                _xrd_kept.append(_Pk)
            if _zaehler is not None:
                # ⚠ THE ENERGY IS THE SECOND ANSWER to the same question.  If (d)
                #   stays large, the selection must run on ENERGY and not on RMSD
                #   (user rule).  So that this sentence is not merely an intention, the
                #   number is here: UFF energy of the finished frame, in kcal/mol, in the
                #   order of the kept conformers.
                #   ⚠ UFF is an ORDERING here, not thermochemistry -- the same force
                #     that also relaxed, so at least self-consistent.  Whoever
                #     turns it into populations overstretches it.
                try:
                    _zaehler.setdefault("energien", []).append(
                        float(AllChem.UFFGetMoleculeForceField(m2).CalcEnergy()))
                except Exception:
                    pass
            label = "pucker " + "+".join(
                f"r{ri_i}:{'base' if combo[ri_i] == 0 else combo[ri_i]}"
                for ri_i in range(len(rings)))
            out.append((_conf_to_xyz(m2), label))
        except Exception:
            if _zaehler is not None:
                # The silent failure gets a number.  Without it,
                # "the gate rejected" could not be told from "something blew up"
                # -- two completely different findings.
                _zaehler["ausnahme"] = _zaehler.get("ausnahme", 0) + 1
            continue
    if _zaehler is not None:
        _zaehler["t_kombis"] = _t.perf_counter() - _zaehler["t0"]
    return out


def selbsttest_raum() -> int:
    """Proves that the CP grid really covers the pucker space.

    Invocation:  python -m delfin.manta._ring_pucker
    Without this test `_pucker_space_grid` would be an assertion -- and the error
    falls in `generate` into an `except Exception: continue`, i.e. SILENTLY.
    """
    import itertools as _itt
    fehler = 0
    print("=== Selbsttest: der Faltungsraum ===")

    # 1 DIMENSION COUNT.  An N-ring has exactly N-3 puckering degrees of freedom.
    for n in range(4, 9):
        even = (n % 2 == 0)
        n_paare = len(range(2, (n // 2) if even else ((n - 1) // 2) + 1))
        dof = 2 * n_paare + (1 if even else 0)
        if dof != n - 3:
            print("  ✗ 1 DIMENSION n=%d: %d statt %d" % (n, dof, n - 3)); fehler += 1
    if not fehler:
        print("  ✓ 1 DIMENSION: N-3 Freiheitsgrade fuer N=4..8 (1,2,3,4,5)")

    # 2 THE PLANAR STATE is in the grid -- the old list lacked it (218 motifs).
    for n in (5, 6, 7):
        g = _pucker_space_grid(n, 2, 8)
        if not any(all(v == 0.0 for v in qs.values()) for qs, _ in g):
            print("  ✗ 2 PLANAR fehlt bei n=%d" % n); fehler += 1
    if fehler == 0:
        print("  ✓ 2 PLANAR: Q=0 ist Gitterpunkt fuer n=5,6,7")

    # 3 CHAIR AND INVERTED CHAIR.  For the six-ring the two poles: q2=0, q3=+/-.
    g6 = _pucker_space_grid(6, 2, 8)
    pole = [qs for qs, _ in g6 if qs.get(2, 0.0) == 0.0 and qs.get(3, 0.0) != 0.0]
    if len([1 for qs in pole if qs[3] > 0]) < 1 or len([1 for qs in pole if qs[3] < 0]) < 1:
        print("  ✗ 3 POLE: Sessel/inv. Sessel nicht beide im Gitter"); fehler += 1
    else:
        print("  ✓ 3 POLE: Sessel UND inverser Sessel (q3 mit beiden Vorzeichen)")

    # 4 THE INTERMEDIATE REGION -- exactly what the old list NEVER sampled.
    #   Half-chair/envelope lie between pole and equator: q2>0 AND q3!=0.
    zwischen = [qs for qs, _ in g6 if qs.get(2, 0.0) > 0 and qs.get(3, 0.0) != 0]
    if not zwischen:
        print("  ✗ 4 ZWISCHENBEREICH leer -- Half-Chair/Envelope unerreichbar"); fehler += 1
    else:
        print("  ✓ 4 ZWISCHENBEREICH: %d Punkte mit q2>0 UND q3!=0" % len(zwischen))

    # 5 HIGHER PAIRS from n=7 on -- without them the seven-ring is incomplete.
    g7 = _pucker_space_grid(7, 2, 8)
    if not any(qs.get(3, 0.0) != 0.0 for qs, _ in g7):
        print("  ✗ 5 m=3 fehlt beim Siebenring"); fehler += 1
    else:
        print("  ✓ 5 HOEHERE PAARE: m=3 wird beim Siebenring belegt")

    # 6 NO DUPLICATE grid points (otherwise the product bloats without gain).
    for n in (5, 6, 7, 8):
        g = _pucker_space_grid(n, 2, 6)
        keys = [tuple(sorted((m, round(q, 6)) for m, q in qs.items())) for qs, _ in g]
        print("     n=%d: %4d Kandidaten (%d-dim)" % (n, len(g), n - 3))

    # 7 THE MACROCYCLE MUST NOT BRING THE BRANCH TO A HALT.
    #   Without a budget there would be about 1.2e8 candidates PER RING at n=16, each with relax
    #   and collision gate -- the run dies and delivers ZERO folds.  Infeasibility
    #   is the opposite of completeness.  Two things are checked:
    #   the number stays under the budget, AND the ring is folded nevertheless
    #   (m=2 keeps full phase resolution, only the fine ripple drops out).
    print("     -- Makrozyklen (Budget %s) --"
          % _os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000"))
    _budget_soll = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000") or 50000))
    for n in (12, 16, 21):
        g = _pucker_space_grid(n, 2, 8)
        res = getattr(_pucker_space_grid, "_grid_res", {}) or {}
        if len(g) > _budget_soll:
            print("  ✗ 7 n=%d: %d Kandidaten UEBER Budget %d" % (n, len(g), _budget_soll))
            fehler += 1
            continue
        # m=2 is the dominant fold and must keep its full phase count
        if res.get("n_phase_je_m", {}).get(2) != 8:
            print("  ✗ 7 n=%d: m=2 wurde reduziert (%s) -- die dominante Falte"
                  % (n, res.get("n_phase_je_m", {}).get(2)))
            fehler += 1
            continue
        if not any(qs.get(2, 0.0) > 0 for qs, _ in g):
            print("  ✗ 7 n=%d: keine einzige gefaltete Konfiguration" % n)
            fehler += 1
            continue
        # the alternation term is the chair for the even ring -- it must never be missing
        if n % 2 == 0 and not any(qs.get(n // 2, 0.0) != 0 for qs, _ in g):
            print("  ✗ 7 n=%d: Alternierungsterm q_%d fehlt -- kein Sessel" % (n, n // 2))
            fehler += 1
            continue
        print("     n=%2d: %6d Kandidaten (%2d-dim), Phasen je m %s, ruhend %s"
              % (n, len(g), n - 3, res.get("n_phase_je_m"),
                 res.get("ruhende_moden") or "keine"))
    if fehler == 0:
        print("  ✓ 7 MAKROZYKLEN: unter Budget, m=2 voll aufgeloest, Faltung erreichbar")

    print("=== Faltungsraum: %s ===" % ("BESTANDEN" if fehler == 0 else "%d FEHLER" % fehler))
    return 1 if fehler else 0


def selbsttest_tfd_sweep(sizes=(5, 6, 7, 8)) -> int:
    """THE THRESHOLD, NOT THE INSTRUMENT.  Does TFD over-split at 0.05?

    BACKGROUND.  `selbsttest_konvergenz` reports for n=6 and n=7 that the number
    of states keeps growing with the resolution (8->9->10 and 11->11->12 respectively), at
    smallest CP distances of 5.0 and 1.8 degrees.  Two causes are possible and
    have OPPOSITE repairs:
        (a) genuine basins, grid too coarse   -> sample finer
        (b) over-splitting by TFD             -> raise threshold
    The first attempt to settle this via CP dedup FAILED and the
    measurement stands: n=5 went from 3,3,3 to 9,13,14.  TFD folds in the MOLECULAR
    SYMMETRY, the CP distance does not -- phi depends on the ring numbering.  CP is no
    replacement.  So exactly this route remains: the same instrument, different threshold.

    WHY THIS DECIDES THE COMBINATORICS.  Measured 4.3 rings per system.
    The cross product over the rings grows like (states per ring)^(rings):
        3 states, 4 rings  ->      81 combinations   computable
       10 states, 4 rings  ->  10 000                not computable
    If 0.05 is too fine, the complete combinatorics becomes affordable through this -- without
    truncating anywhere.  If it is right, the price is real and the
    cap would otherwise come back through the back door.

    ⚠ THIS TEST DOES NOT JUDGE CHEMISTRY.  It measures on UNSUBSTITUTED rings,
      whose symmetry is high; a substituted ring legitimately has more states.
      What it shows is the UPPER BOUND of the over-splitting, not the production number.
    """
    if not (_RDKIT and _np is not None):
        print("=== TFD-Sweep: RDKit fehlt, uebersprungen ==="); return 0
    print("=== Selbsttest: TFD-Schwellensweep (NPHASE=8, unsubstituierte Ringe) ===")
    schwellen = (0.02, 0.05, 0.10, 0.15, 0.20, 0.30)
    print("    Ring  " + "".join("%7.2f" % t for t in schwellen))
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NPHASE",
             "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_TRACE",
             "DELFIN_FFFREE_PUCKER_CPDEDUP")}
    tabelle = {}
    try:
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        _os.environ.pop("DELFIN_FFFREE_PUCKER_CPDEDUP", None)   # measure pure TFD
        for n in sizes:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles("C1" + "C" * (n - 1) + "1"))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    n=%d  Einbettung fehlgeschlagen" % n); continue
                AllChem.MMFFOptimizeMolecule(m)
                ri = m.GetRingInfo().AtomRings()
                if not ri:
                    print("    n=%d  kein Ring gefunden" % n); continue
                ring = _ring_order(m, set(ri[0]))
            except Exception as e:
                print("    n=%d  Aufbau fehlgeschlagen: %s" % (n, type(e).__name__)); continue
            zeile = []
            for thr in schwellen:
                try:
                    zeile.append(len(_ring_pucker_states(m, ring, set(), thr)))
                except Exception:
                    zeile.append(-1)
            tabelle[n] = zeile
            print("    n=%-3d " % n + "".join("%7d" % z for z in zeile))
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    if not tabelle:
        print("=== TFD-Sweep: nichts gemessen ==="); return 1

    # ---- WHAT THIS COSTS.  4.3 rings per system, measured on 400 systems.
    print()
    print("    Kreuzprodukt bei 4 Ringen je System (Zustaende^4):")
    i05 = schwellen.index(0.05)
    for n, zeile in sorted(tabelle.items()):
        s05, s10 = zeile[i05], zeile[schwellen.index(0.10)]
        print("      n=%-3d  Schwelle 0,05 -> %8d      Schwelle 0,10 -> %8d"
              % (n, max(0, s05) ** 4, max(0, s10) ** 4))

    # ---- VERDICT.  Over-splitting means: the number falls sharply and then STAYS flat.
    #      If it keeps falling uniformly, the higher threshold merges genuine
    #      basins -- then 0.05 is not too fine, rather the threshold is the wrong
    #      tool.  Exactly this distinction is the point of the sweep.
    print()
    print("    URTEIL je Ringgroesse:")
    verdacht = 0
    for n, zeile in sorted(tabelle.items()):
        if min(zeile) < 0 or zeile[i05] <= 0:
            print("      n=%-3d  nicht messbar" % n); continue
        _sturz = 1.0 - (zeile[i05 + 1] / float(zeile[i05]))       # 0.05 -> 0.10
        _rest = 1.0 - (zeile[-1] / float(max(1, zeile[i05 + 1])))  # 0.10 -> 0.30
        if _sturz >= 0.34 and _rest <= _sturz:
            print("      n=%-3d  UEBERSPLITTUNG: %d -> %d bei 0,05 -> 0,10 (%.0f %%), "
                  "danach nur noch %.0f %% -- der Sturz sitzt AN der Schwelle"
                  % (n, zeile[i05], zeile[i05 + 1], 100 * _sturz, 100 * _rest))
            verdacht += 1
        elif _sturz < 0.15:
            print("      n=%-3d  STABIL: %d -> %d (%.0f %%) -- 0,05 splittet NICHT ueber"
                  % (n, zeile[i05], zeile[i05 + 1], 100 * _sturz))
        else:
            print("      n=%-3d  GLEITEND: %.0f %% dann %.0f %% -- die Schwelle verschmilzt "
                  "fortlaufend, also auch ECHTE Mulden.  Kein sauberer Schnittpunkt."
                  % (n, 100 * _sturz, 100 * _rest))
    print("=== TFD-Sweep: %d von %d Ringgroesse(n) mit Uebersplittungsverdacht ==="
          % (verdacht, len(tabelle)))
    return 0


def selbsttest_konvergenz(sizes=(5, 6, 7)) -> int:
    """CONVERGENCE instead of assertion: does the number of minima still grow with the resolution?

    THE QUESTION this answers.  The parameter space (q_m, phi_m) is CONTINUOUS
    -- every real combination is a valid geometry.  The CONFORMER space is
    not: a ring has finitely many energy minima.  The grid is therefore not a
    result but a distribution of STARTING POINTS; `_relax_hold_pucker` pulls every
    point into the nearest genuine minimum, TFD dedups.

    ⇒ Completeness is REACHABLE, not merely approachable: the grid must be fine
      enough that every basin of attraction is hit at least once.  Whether that
      is the case is told by exactly one measurement -- the number of distinct states
      against the resolution.  If it no longer grows, the space is covered.

    ⚠ WHY THIS IS HERE.  `conformer_enum.py:7-9` asserts the same ("finer grid
    stops adding distinct minima") and has NO measuring point for it.  An assertion
    without evidence is exactly the design that has already carried a
    wrong number several times in this project.
    """
    if not (_RDKIT and _np is not None):
        print("=== Konvergenz: RDKit fehlt, uebersprungen ==="); return 0
    print("=== Selbsttest: Konvergenz des Faltungsgitters ===")
    print("    Ring   NPHASE=4   8   16     konvergiert?")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NPHASE",
             "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_TRACE")}
    fehler = 0
    try:
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        # Cross-check with CP dedup if the caller has set it --
        # otherwise the test re-measures the old over-splitting.
        # ⚠ READ ONCE, NAME ONCE.  The first version read the switch here
        #   and called it `_cpd_an` below -- a name that never existed.  The
        #   NameError fell into the `except Exception` of the CP scatter and was
        #   printed as "not measurable": the mode-corrected verdict thus ran NOT
        #   a single time, and the test still reported something.  A swallowed
        #   error is a null measurement that looks like a finding.
        _cpd_an = _os.environ.get("DELFIN_FFFREE_PUCKER_CPDEDUP") == "1"
        if _cpd_an:
            print("    (CP-Entdopplung AN, Toleranz %s Grad / %s A)"
                  % (_os.environ.get("DELFIN_FFFREE_PUCKER_CPTOL", "15"),
                     _os.environ.get("DELFIN_FFFREE_PUCKER_CPQTOL", "0.15")))
        for n in sizes:
            smi = "C1" + "C" * (n - 1) + "1"
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    n=%d  Einbettung fehlgeschlagen" % n); continue
                AllChem.MMFFOptimizeMolecule(m)
                ri = m.GetRingInfo().AtomRings()
                if not ri:
                    print("    n=%d  kein Ring gefunden" % n); continue
                ring = _ring_order(m, set(ri[0]))
            except Exception as e:
                print("    n=%d  Aufbau fehlgeschlagen: %s" % (n, type(e).__name__)); continue
            zahlen = []
            for nph in (4, 8, 16):
                _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = str(nph)
                try:
                    st = _ring_pucker_states(m, ring, set(), 0.05)
                    zahlen.append(len(st))
                except Exception as e:
                    zahlen.append(-1)
            ok = (len(zahlen) == 3 and zahlen[1] > 0 and zahlen[2] <= zahlen[1])
            print("    n=%-3d  %8d %3d %4d      %s"
                  % (n, zahlen[0], zahlen[1], zahlen[2],
                     "JA" if ok else "NEIN -- siehe CP-Streuung"))
            if not ok:
                fehler += 1
                # ---- WHY does the number grow?  Two causes, OPPOSITE fixes --
                # (1) space not yet covered -> the states lie FAR apart in CP
                #     coordinates -> sample finer.
                # (2) TFD threshold separates chemically IDENTICAL states -> they lie
                #     CLOSE together -> the magnifier is too fine, not the grid coarse.
                #
                # ⚠️ THIS IS NOT AN ACADEMIC QUESTION.  The number of states PER RING
                #    enters the cross product over all rings as the base:
                #        3 states, 4 rings ->      81 combinations  (computable)
                #       10 states, 4 rings ->  10 000                (not computable)
                #    Measured are 4.3 rings per system.  Over-splitting makes the
                #    COMPLETE combinatorics unaffordable -- convergence anomaly and
                #    computability are the same problem.
                # Chemical yardstick: cyclohexane has chair + twist-boat family,
                # after symmetry folding 2-3 classes.  The five-ring converges to 3.
                _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "16"
                try:
                    cps = []
                    for _s in _ring_pucker_states(m, ring, set(), 0.05):
                        if _s is None or not (isinstance(_s, tuple) and len(_s) == 2
                                              and isinstance(_s[0], dict)):
                            continue
                        m3 = Chem.Mol(m)
                        _set_pucker_general(m3.GetConformer(), ring, _s[0], _s[1], set())
                        _relax_hold_pucker(m3, ring, set())
                        cps.append(_cp_theta_phi(m3.GetConformer().GetPositions(), ring))
                    cps.sort(key=lambda t: (round(t[1], 0), round(t[2], 0)))
                    print("        CP der ueberlebenden Zustaende (Q, theta, phi):")
                    for _Q, _th, _ph in cps:
                        print("          Q=%.3f  theta=%6.1f  phi=%6.1f" % (_Q, _th, _ph))
                    dmin = None
                    for _i in range(len(cps)):
                        for _j in range(_i + 1, len(cps)):
                            _, ti, pi_ = cps[_i]
                            _, tj, pj = cps[_j]
                            dph = min(abs(pi_ - pj), 360.0 - abs(pi_ - pj))
                            d = ((ti - tj) ** 2 + dph ** 2) ** 0.5
                            dmin = d if dmin is None else min(dmin, d)
                    if dmin is not None:
                        print("        kleinste paarweise CP-Distanz: %.1f Grad" % dmin)
                        # ⚠ THE VERDICT MUST KNOW WHICH MODE RAN.  The first version
                        #   was hard-coupled to `dmin` and claimed "TFD too fine"
                        #   even when the CP dedup was running -- i.e. a verdict
                        #   about an instrument that was not in use at all.
                        if _cpd_an:
                            print("        URTEIL: CP-Entdopplung laeuft und liefert MEHR"
                                  " Zustaende als TFD.  Grund: TFD faltet die MOLEKUEL-"
                                  "SYMMETRIE mit, die CP-Distanz nicht.  phi haengt an der"
                                  " Ringnummerierung -- bei einem unsubstituierten Ring"
                                  " sind alle phi bei gleichem (Q, theta) DERSELBE"
                                  " Konformer.  CP allein ist KEIN Ersatz fuer TFD.")
                        else:
                            print("        URTEIL: %s" % (
                                "ECHTE Mulden -- Gitter zu grob, feiner abtasten"
                                if dmin > 20.0 else
                                "Zustaende liegen dicht -- Uebersplittung moeglich; "
                                "PRUEFEN durch TFD-Schwellensweep, NICHT durch Ersetzen "
                                "von TFD (siehe CP-Modus)"))
                except Exception as _e:
                    print("        CP-Streuung nicht messbar: %s" % type(_e).__name__)
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v
    print("=== Konvergenz: %s ===" %
          ("BESTANDEN -- der Raum ist bei NPHASE=8 ueberdeckt" if fehler == 0
           else "%d Ringgroesse(n) NICHT konvergiert" % fehler))
    return 1 if fehler else 0


# ===== THE SAMPLES: A COUPLING GRADIENT, NOT A COLLECTION =========================
#
# The question is not "how many folds does molecule X have", but WHAT it depends on
# how much of the cross product remains.  Only a VARIABLE can answer that, one that
# stands systematically different from sample to sample -- here the number of SHARED ATOMS
# between two rings.  It runs down the list from 2 to 0:
#
#   2 shared atoms, three bridges       bridged (bicyclo[2.2.2]octane, norbornane)
#                                       -- the stiffest case there is
#   2 shared atoms, one bond            fused (decalin, perhydroanthracene)
#   1 shared atom                       spiro (spiro[5.5]undecane)
#   0, direct ring-ring bond            only STERICALLY coupled (bicyclohexyl)
#   0, two CH2 in between               practically independent (1,2-dicyclohexylethane)
#   0, three rings on one P             tricyclohexylphosphine -- the case the
#                                       `generate` docstring itself cites as an
#                                       example, and a real ligand
#
# Cyclohexane stands as the ZERO POINT: ONE ring, so no cross product at all.  Without it
# one could not separate what the COUPLING costs from what the single ring already costs.
# Perhydroanthracene and the phosphine are the only THREE-ring samples -- only at
# three rings does it show whether the curve runs exponentially or capped.
_KOMBI_PROBEN = (
    ("Cyclohexan",            "C1CCCCC1",                      "1 Ring -- Nullpunkt"),
    ("Bicyclo[2.2.2]octan",   "C1CC2CCC1CC2",                  "verbrueckt, 3 Bruecken"),
    ("Norbornan",             "C1CC2CCC1C2",                   "verbrueckt, 1 Bruecke"),
    ("Decalin",               "C1CCC2CCCCC2C1",                "kondensiert, 1 Bindung"),
    ("Perhydroanthracen",     "C1CCC2CC3CCCCC3CC2C1",          "3 Ringe, kondensiert"),
    ("Spiro[5.5]undecan",     "C1CCC2(CC1)CCCCC2",             "spiro, 1 Atom"),
    ("Bicyclohexyl",          "C1CCCCC1C1CCCCC1",              "0 Atome, 1 Bindung"),
    ("1,2-Dicyclohexylethan", "C1CCCCC1CCC1CCCCC1",            "0 Atome, 2 CH2"),
    ("Tricyclohexylphosphin", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1", "3 Ringe, unabhaengig"),
)


def selbsttest_kombinatorik(proben=None, deckel_s: float = 900.0,
                            max_kombis: int = 4000) -> int:
    """HOW LARGE IS THE CROSS PRODUCT **AFTER** THE PHYSICS?

    THE QUESTION that decides the complete ring folding.  Measured are 4.3
    rings per system and -- with the CP space grid -- 9 to 16 pucker states per ring
    (`selbsttest_tfd_sweep`, threshold 0.05).  The naive cross product is thus 6500 to
    65000 combinations per system.  But rings of a molecule are NOT independent:
    fused and bridged rings share atoms; fold one, and the
    other is fixed.  The cross product is an upper bound of the ENUMERATION -- the
    question is what of it survives the realism gate.

    ⚠ THE DECISIVE DIFFERENCE this test makes visible: the gate kills the
      RESULT, but not the COST.  Every combination is first set, then relaxed with
      held folds and ONLY THEN rejected.  Whoever reads "(c) is small, hence
      cheap" has confused the order.  That is why TWO numbers stand side by side
      here: (b) is the price, (d) is the yield.

    ⚠ TWO GATES, NOT ONE.  `_has_clash` sees only overlapping vdW spheres.  For
      FUSED rings, however, the contradiction arises at the shared atoms, and
      there the VSEPR angle is no longer right without anything colliding --
      exactly for that `_has_bad_angles` exists.  Which of the two gates fires is
      therefore itself a finding and is counted separately.

    ⚠ WHAT THIS TEST DOES NOT MEASURE.  It runs on METAL-FREE hydrocarbons with
      high symmetry.  A substituted ring legitimately has more states, a real
      DELFIN system is larger and thus more expensive per combination.  The times here are
      a LOWER BOUND of the cost, not the production number.

    ``max_kombis``: test limit.  (a) and (b) are ALWAYS determined -- they cost only the
    linear item.  If (b) lies above it, (c) and (d) are NOT measured and exactly
    that is printed, instead of outputting a capped number.  ⚠ A cap would be
    especially insidious here: `combos` is sorted by FOLD DEPTH, a prefix of it
    contains only shallow combinations and would have systematically too high survival rates.
    ``0`` = no limit (then a single multi-ring molecule can run for hours).
    """
    if not (_RDKIT and _np is not None):
        print("=== Kombinatorik: RDKit fehlt, uebersprungen ==="); return 0
    import time as _time
    proben = proben or _KOMBI_PROBEN
    print("=== Selbsttest: das Kreuzprodukt NACH dem Realismustor ===")
    print("    Raumgitter AN (NAMP=2, NPHASE=8), Kappe AUS -- die VOLLE Kombinatorik.")
    print("    (a) Zustaende je Ring · (b) Kreuzprodukt · (c) ueberlebt das Tor"
          " · (d) davon TFD-distinkt")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE")}
    zeilen = []
    fehler = 0
    try:
        # ===== 0 THE DEFAULT PATH MUST REMAIN UNCHANGED -- MEASURED, NOT ASSERTED
        #
        # `generate` now carries an optional counter set.  The assertion "with
        # `_zaehler=None` nothing changes" is exactly the kind of assertion that has
        # already been wrong several times in this project.  So it is measured: the same
        # molecule, the same default path (all switches OFF, cap ON), once without and
        # once with counters -- the returned frames must be CHARACTER-IDENTICAL.
        # ⚠ In measurement mode both gates are asked instead of one short-circuited; if
        #   `_has_bad_angles` changed anything, it would show up exactly here.
        for _k in _alt:
            _os.environ[_k] = "0"
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))    # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) == 0:
            AllChem.MMFFOptimizeMolecule(_mv)
            _ohne = generate(_mv, budget=48)
            _mit = generate(_mv, budget=48, _zaehler=_neuer_zaehler())
            if _ohne == _mit:
                print("    ✓ 0 VORGABE UNVERAENDERT: Decalin, Schalter AUS, %d Frames "
                      "mit und ohne Zaehler identisch" % len(_ohne))
            else:
                print("    ✗ 0 VORGABE VERAENDERT: %d Frames ohne Zaehler, %d mit -- "
                      "die Messstelle ist nicht folgenlos" % (len(_ohne), len(_mit)))
                fehler += 1
        else:
            print("    ? 0 VORGABE: Decalin nicht einbettbar, Identitaet NICHT gemessen")
            fehler += 1

        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"     # no cap -- whole product
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        print()
        print("    %-22s %2s %4s  %-14s %8s %7s %7s %7s %8s"
              % ("Molekuel", "R", "gmax", "(a) je Ring", "(b)", "(c)", "(d)",
                 "c/b", "Zeit/s"))
        print("    (gmax = Atome, die sich das ENGSTE Ringpaar teilt: 0 getrennt · "
              "1 spiro · 2 kondensiert · >=3 verbrueckt)")
        for name, smi, klasse in proben:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    %-22s Einbettung fehlgeschlagen" % name); continue
                try:
                    AllChem.MMFFOptimizeMolecule(m)
                except Exception:
                    AllChem.UFFOptimizeMolecule(m)
            except Exception as e:
                print("    %-22s Aufbau fehlgeschlagen: %s" % (name, type(e).__name__))
                continue
            # ---- THE TEST LIMIT: (a) and (b) FIRST, separate from the build.  The build costs
            #      one relax plus two gates PER combination; whether it is affordable
            #      is decided by (b) -- so (b) must be known BEFORE building.
            #      ⚠ The pre-pass determines the states a second time (`generate` does
            #        it again right after).  That is the LINEAR item, i.e. the cheap one --
            #        but it is paid nonetheless, which is why it runs ONLY if the limit
            #        is set at all.  With `max_kombis=0` there is no pre-pass
            #        and thus no duplicated work in the timing either.
            if max_kombis:
                _t_vor = _time.perf_counter()
                try:
                    _rings = [_ring_order(m, set(r)) for r in m.GetRingInfo().AtomRings()
                              if _is_puckerable(m, r)]
                    _stv = [_ring_pucker_states(m, r, set(), 0.05) for r in _rings]
                except Exception as e:
                    print("    %-22s Zustaende nicht bestimmbar: %s"
                          % (name, type(e).__name__))
                    continue
                _prod = 1
                for _s in _stv:
                    _prod *= len(_s)
                _b_vor = max(0, _prod - 1)
                _gem_vor = max([len(set(_rings[i]) & set(_rings[j]))
                                for i in range(len(_rings))
                                for j in range(i + 1, len(_rings))] or [0])
                if _b_vor > max_kombis:
                    print("    %-22s %2d %4d  %-14s %8d %7s %7s %7s %8.1f"
                          % (name, len(_rings), _gem_vor,
                             "x".join(str(len(s)) for s in _stv), _b_vor,
                             "-", "-", "-", _time.perf_counter() - _t_vor))
                    print("        %-28s (c) und (d) NICHT GEMESSEN: (b) = %d ueber der "
                          "Testgrenze %d.  Kein gekappter Ersatzwert -- `combos` ist nach "
                          "Faltungstiefe sortiert, ein Praefix waere systematisch zu flach."
                          % (klasse, _b_vor, max_kombis))
                    zeilen.append({"name": name, "klasse": klasse, "ringe": len(_rings),
                                   "gem": _gem_vor, "gmax": _gem_vor,
                                   "b": 0, "c": 0, "d": 0, "dt": 0.0,
                                   "zustaende": [len(s) for s in _stv],
                                   "t_kombis": 0.0, "t_zust": 0.0,
                                   "kollision": 0, "winkel": 0, "energien": [],
                                   "ungemessen": _b_vor})
                    continue
            z = _neuer_zaehler()
            _t0 = _time.perf_counter()
            try:
                out = generate(m, budget=10 ** 9, _zaehler=z)
            except Exception as e:
                print("    %-22s generate() geplatzt: %s" % (name, type(e).__name__))
                continue
            _dt = _time.perf_counter() - _t0
            _b = int(z.get("aufzaehlung", 0))
            _c = int(z.get("tor_ueberlebt", 0))
            _d = len(out)
            # ⚠ NEVER A PERCENTAGE WITHOUT A DENOMINATOR.  The denominator here is (b), the number of
            #   enumerated combinations without the ground state -- not the
            #   cross product itself, because the ground state is never built.
            _cb = ("%6.1f%%" % (100.0 * _c / _b)) if _b else "   n/a"
            print("    %-22s %2d %4d  %-14s %8d %7d %7d %7s %8.1f"
                  % (name, len(z.get("ringgroessen") or []),
                     int(z.get("gem_max", 0)),
                     "x".join(str(v) for v in (z.get("zustaende_je_ring") or [])) or "-",
                     _b, _c, _d, _cb, _dt))
            print("        %-28s verworfen: Kollision %d · Winkel %d · TFD %d · "
                  "Relax %d · Ausnahme %d"
                  % (klasse, int(z.get("kollision", 0)), int(z.get("winkel", 0)),
                     int(z.get("tfd_doppelt", 0)), int(z.get("relax_fehler", 0)),
                     int(z.get("ausnahme", 0))))
            _en = sorted(z.get("energien") or [])
            if len(_en) >= 2:
                _e0 = _en[0]
                _in10 = sum(1 for e in _en if e - _e0 <= 10.0)
                print("        UFF-Energie der (d): Spanne %.1f kcal/mol · "
                      "innerhalb 10 kcal/mol %d von %d"
                      % (_en[-1] - _e0, _in10, len(_en)))
            zeilen.append({"name": name, "klasse": klasse,
                           "ringe": len(z.get("ringgroessen") or []),
                           "gem": int(z.get("gemeinsame_atome", 0)),
                           "gmax": int(z.get("gem_max", 0)),
                           "b": _b, "c": _c, "d": _d, "dt": _dt,
                           "zustaende": list(z.get("zustaende_je_ring") or []),
                           "t_kombis": float(z.get("t_kombis", 0.0)),
                           "t_zust": float(z.get("t_zustaende", 0.0)),
                           "kollision": int(z.get("kollision", 0)),
                           "winkel": int(z.get("winkel", 0)),
                           "energien": _en})
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    if not zeilen:
        print("=== Kombinatorik: nichts gemessen ==="); return 1

    # ---- 1 COUPLING VERSUS SURVIVAL.
    # ⚠ PARTITION, NOT MEAN.  "coupled versus independent" would be the wrong
    #   dichotomy: it throws a fused ring pair (one shared BOND) into one pot with
    #   a cage (four shared atoms), and their survival rates
    #   lie two orders of magnitude apart.  A mean cannot see a dead class
    #   -- only a partition can.  Therefore the split is by the number of
    #   atoms the TIGHTEST ring pair shares; that is at the same time the chemical name
    #   of the coupling.
    mehr = [r for r in zeilen if r["ringe"] >= 2 and r["b"] > 0]

    def _klasse(r):
        g = r.get("gmax", 0)
        return 0 if g == 0 else (1 if g == 1 else (2 if g == 2 else 3))

    _NAMEN = {0: "0 Atome  getrennt", 1: "1 Atom   spiro",
              2: "2 Atome  kondensiert", 3: ">=3      verbrueckt"}

    def _quote(gruppe):
        _b = sum(r["b"] for r in gruppe)
        _c = sum(r["c"] for r in gruppe)
        return _b, _c, (100.0 * _c / _b if _b else 0.0)

    print()
    print("    ===== 1 KOPPLUNG GEGEN UEBERLEBEN (nur Mehrringer) =====")
    print("      engstes Ringpaar teilt ...")
    _quoten = {}
    for _kl in (0, 1, 2, 3):
        _g = [r for r in mehr if _klasse(r) == _kl]
        if not _g:
            print("      %-22s keine Probe" % _NAMEN[_kl]); continue
        _b, _c, _q = _quote(_g)
        _quoten[_kl] = _q
        # ⚠ (a) MUST BE SHOWN ALONGSIDE, otherwise c/b is not interpretable.  A bridged
        #   ring may already have FEWER states -- then (b) is small because the
        #   coupling acted earlier, and not because the gate kills more.  Two
        #   different routes to the same small product, and only both together say
        #   which one it was.
        _zust = [v for r in _g for v in r["zustaende"]]
        print("      %-22s %d Probe(n) · %6d von %6d ueberleben = %5.1f %% · (a) im "
              "Mittel %.1f je Ring (%d Ringe) · %s"
              % (_NAMEN[_kl], len(_g), _c, _b, _q,
                 (sum(_zust) / float(len(_zust))) if _zust else 0.0, len(_zust),
                 ", ".join(r["name"] for r in _g)))
    if len(_quoten) >= 2:
        _hi = max(_quoten.values())
        _lo = min(_quoten.values())
        # ⚠ The finding is the RANGE across the partition, not a group mean.
        print("      ⇒ Spanne ueber die Kopplungsklassen: %.1f %% bis %.1f %% -- %s"
              % (_lo, _hi,
                 "die Kopplungsart entscheidet, nicht die Kopplung an sich"
                 if _hi - _lo > 20.0 else
                 "die Kopplungsart macht kaum einen Unterschied"))

    # ---- 2 WHICH GATE FIRES.  Collision and angle separately, otherwise "the gate"
    #      is a name for two different mechanisms (detector name != measurement).
    print()
    print("    ===== 2 WELCHES TOR TOETET =====")
    for _kl in (0, 1, 2, 3):
        _g = [r for r in mehr if _klasse(r) == _kl]
        if not _g:
            continue
        _b = sum(r["b"] for r in _g)
        _k = sum(r["kollision"] for r in _g)
        _w = sum(r["winkel"] for r in _g)
        print("      %-22s von %6d Kombinationen: Kollision %6d (%5.1f %%) · "
              "Winkel %6d (%5.1f %%)"
              % (_NAMEN[_kl], _b, _k, 100.0 * _k / _b if _b else 0.0,
                 _w, 100.0 * _w / _b if _b else 0.0))
    # ⚠ A GATE NAME IS NOT A MEASUREMENT.  If "the collision gate" in truth never fires
    #   and the whole selection comes from the angle gate, then every statement about "the
    #   sterics cut the product" rests on the wrong mechanism -- and a repair
    #   at the collision gate would be ineffective before it is written.
    _kges = sum(r["kollision"] for r in mehr)
    _wges = sum(r["winkel"] for r in mehr)
    _bges = sum(r["b"] for r in mehr)
    if _bges:
        if _kges == 0 and _wges > 0:
            print("      ⇒ DAS KOLLISIONSTOR HAT NULL REICHWEITE: 0 von %d Kombinationen."
                  % _bges)
            print("        Der Filter ist AUSSCHLIESSLICH das WINKELTOR (%d von %d = "
                  "%.1f %%).  Wer die Kombinatorik am Kollisionstor beschneiden will, "
                  "greift den Mechanismus an, der gar nicht feuert."
                  % (_wges, _bges, 100.0 * _wges / _bges))
        else:
            print("      ⇒ Kollision %d von %d (%.1f %%) · Winkel %d von %d (%.1f %%) "
                  "-- beide Tore tragen."
                  % (_kges, _bges, 100.0 * _kges / _bges,
                     _wges, _bges, 100.0 * _wges / _bges))

    # ---- 3 AFFORDABILITY.  The cost hangs on (b), not on (d).
    _sum_b = sum(r["b"] for r in zeilen)
    _sum_t = sum(r["t_kombis"] for r in zeilen)
    print()
    print("    ===== 3 BEZAHLBARKEIT =====")
    if _sum_b <= 0 or _sum_t <= 0.0:
        print("      Kosten je Kombination NICHT MESSBAR (b=%d, t=%.3f s)"
              % (_sum_b, _sum_t))
        print("=== Kombinatorik: unvollstaendig ==="); return 1
    _ms = 1000.0 * _sum_t / _sum_b
    print("      Jede der (b) Kombinationen wird GEBAUT und RELAXIERT, bevor das Tor")
    print("      sie verwirft -- das Tor spart nichts, es waehlt nur aus.")
    print("      Gemessen: %.1f ms je Kombination (Nenner: %d Kombinationen ueber %d "
          "Proben, %.1f s gesamt)" % (_ms, _sum_b, len(zeilen), _sum_t))
    print("      ⚠ UNTERGRENZE: metallfreie Kohlenwasserstoffe, 7 bis 21 Schweratome. "
          "Ein echtes System ist groesser und je Kombination teurer.")
    # ⚠ THE COST PER COMBINATION IS NOT A CONSTANT, and the scatter must be
    #   stated with it.  It rises with the SURVIVAL RATE: what passes the gate is
    #   checked by TFD against EVERY already kept conformer, i.e. quadratically.
    #   A molecule whose combinations all survive is thus doubly expensive --
    #   more candidates AND a more expensive check per candidate.
    _je = sorted((1000.0 * r["t_kombis"] / r["b"], r["name"])
                 for r in zeilen if r["b"] > 0 and r["t_kombis"] > 0.0)
    if len(_je) >= 2:
        print("      Streuung je Kombination: %.1f ms (%s) bis %.1f ms (%s) -- sie "
              "steigt mit der Ueberlebensquote, weil TFD gegen alle Behaltenen prueft."
              % (_je[0][0], _je[0][1], _je[-1][0], _je[-1][1]))
    # ⚠ TWO COST ITEMS WITH DIFFERENT GROWTH.  The states per ring cost
    #   LINEARLY in the ring count (each ring once), the cross product EXPONENTIALLY.  If
    #   the enumeration item still looks small today, that means nothing for 6 rings --
    #   the other one is the one that explodes.
    _t_zust = sum(r["t_zust"] for r in zeilen)
    print("      Aufteilung: Zustaende je Ring %.1f s (linear in der Ringzahl) · "
          "Kreuzprodukt %.1f s (exponentiell) -- Summe %.1f s ueber %d Proben"
          % (_t_zust, _sum_t, _t_zust + _sum_t, len(zeilen)))
    # ⚠ THE EXTRAPOLATION MUST NOT RUN WITH AN INVENTED STATE COUNT.  What was
    #   measured here stands beside it -- and the measured range extends beyond
    #   what the sweep table finds on UNSUBSTITUTED rings: a ring in
    #   a cage is symmetry-poor and splits further.
    _az = [v for r in zeilen for v in r["zustaende"]]
    if _az:
        print("      Gemessene (a): %d bis %d Zustaende je Ring, Mittel %.1f (Nenner: "
              "%d Ringe ueber %d Proben)"
              % (min(_az), max(_az), sum(_az) / float(len(_az)), len(_az), len(zeilen)))
    print()
    print("      Hochrechnung auf 4,3 Ringe je System (gemessen) -- EIN Kern, EIN System:")
    _stufen = sorted({3, 9, 16} | ({max(_az)} if _az else set()))
    for _z in _stufen:
        _n = _z ** 4.3
        _s = _n * _ms / 1000.0
        # ⚠ THE CAP DOES NOT BELONG TO THIS MECHANISM ALONE.  `deckel_s` is the
        #   deadline for the WHOLE build of a system; the ring folding is one of
        #   many steps in it.  "Under the cap" is therefore not yet a "works" --
        #   only the SHARE says whether anything else still has room beside it.
        print("        %2d Zustaende je Ring -> %10.0f Kombinationen -> %10.0f s "
              "= %6.1f h  = %6.1f %% des Arm-Deckels (%.0f s)%s"
              % (_z, _n, _s, _s / 3600.0, 100.0 * _s / deckel_s, deckel_s,
                 "" if _s <= deckel_s else "   UEBER dem Deckel"))

    # ---- 4 VERDICT.  Two sides, and they come out differently.
    # ⚠ THE VERDICT COMPUTES WITH THE MEASURED STATE COUNT, not with the assumed one.
    #   The assumption this measurement arose from was "9 to 16 states per ring"
    #   -- a number from the SWEEP on UNSUBSTITUTED rings.  In real multi-ring molecules
    #   this test measures 4 to 26 with a mean around 14: the environment breaks the ring symmetry,
    #   and TFD then separates more.  Judging with the assumed number while our
    #   own stands beside it would be cooking the books in its purest form.
    _zmit = (sum(_az) / float(len(_az))) if _az else 9.0
    _smess = (_zmit ** 4.3) * _ms / 1000.0
    _s3 = (3 ** 4.3) * _ms / 1000.0
    _dmax = max(r["d"] for r in zeilen)
    print()
    print("    ===== 4 URTEIL =====")
    if _smess <= deckel_s:
        print("      AUFZAEHLUNG: BEZAHLBAR -- bei der GEMESSENEN Zustandszahl %.1f je "
              "Ring und 4,3 Ringen %.0f s je System = %.0f %% des Arm-Deckels (%.0f s), "
              "den sich die Ringfaltung mit jedem anderen Bauschritt teilt."
              % (_zmit, _smess, 100.0 * _smess / deckel_s, deckel_s))
    elif _s3 <= deckel_s:
        print("      AUFZAEHLUNG: NICHT BEZAHLBAR bei der GEMESSENEN Aufloesung -- %.1f "
              "Zustaende je Ring, 4,3 Ringe: %.0f s je System = %.0f %% des Arm-Deckels "
              "(%.0f s)." % (_zmit, _smess, 100.0 * _smess / deckel_s, deckel_s))
        print("        Bezahlbar wird es erst weit darunter: bei 3 Zustaenden je Ring "
              "%.0f s = %.0f %% des Deckels.  Der Weg dahin ist WENIGER ZUSTAENDE JE "
              "RING (groebere Entdopplung), nicht eine Kappe auf dem Produkt -- eine "
              "Kappe schneidet nach Faltungstiefe und laesst die tiefen Faltungen weg."
              % (_s3, 100.0 * _s3 / deckel_s))
    else:
        print("      AUFZAEHLUNG: NICHT BEZAHLBAR -- selbst bei 3 Zustaenden je Ring "
              "%.0f s je System gegen einen Deckel von %.0f s." % (_s3, deckel_s))
    print("      ERGEBNIS: groesste gemessene Ausbeute (d) einer einzelnen Probe: %d "
          "Konformere." % _dmax)
    # ⚠ WHAT WAS NOT MEASURED MUST APPEAR IN THE VERDICT.  A sample that was skipped
    #   because of its size is the strongest case against affordability --
    #   leaving it silently out of the balance would be exactly the book-cooking
    #   this test is built against.
    _uv = [r for r in zeilen if r.get("ungemessen")]
    if _uv:
        print("      ⚠ %d von %d Proben UNGEMESSEN, weil (b) ueber der Testgrenze %d "
              "lag: %s" % (len(_uv), len(zeilen), max_kombis,
                           " · ".join("%s (b=%d)" % (r["name"], r["ungemessen"])
                                      for r in _uv)))
        print("        Das ist selbst ein Befund: bei diesen Systemen ist das volle "
              "Produkt schon zu gross, um es ueberhaupt einmal zu bauen.")
    # ---- THE THREE FILTERS IN SEQUENCE, each with its own denominator.
    # ⚠ This is the core statement of the whole test, and it is readable only as a CHAIN:
    #   which of the three filters actually makes the product small is a measurement
    #   and not a guess -- and the guess was that it is the first one.
    _sum_c = sum(r["c"] for r in zeilen)
    _sum_d = sum(r["d"] for r in zeilen)
    _alle_en = [r for r in zeilen if len(r["energien"]) >= 2]
    _ges = sum(len(r["energien"]) for r in _alle_en)
    _in10 = sum(sum(1 for e in r["energien"] if e - r["energien"][0] <= 10.0)
                for r in _alle_en)
    print()
    print("      DREI FILTER HINTEREINANDER (alle Proben zusammen, jeder mit Nenner):")
    print("        1 PHYSIK  (Kollision + Winkel)  (b)->(c)  %6d von %6d = %5.1f %%"
          % (_sum_c, _sum_b, 100.0 * _sum_c / _sum_b if _sum_b else 0.0))
    print("        2 TFD     (Entdopplung)         (c)->(d)  %6d von %6d = %5.1f %%"
          % (_sum_d, _sum_c, 100.0 * _sum_d / _sum_c if _sum_c else 0.0))
    if _ges:
        print("        3 ENERGIE (<= 10 kcal/mol)      (d)->(e)  %6d von %6d = %5.1f %%"
              % (_in10, _ges, 100.0 * _in10 / _ges))
        # The sharpest filter is the one with the SMALLEST pass rate.
        _kette = (("die PHYSIK", 100.0 * _sum_c / max(1, _sum_b)),
                  ("die TFD-Entdopplung", 100.0 * _sum_d / max(1, _sum_c)),
                  ("die ENERGIE", 100.0 * _in10 / _ges))
        _eng = min(_kette, key=lambda t: t[1])
        print("      ⇒ Der schaerfste Filter ist %s (%.1f %% Durchlass)."
              % (_eng[0], _eng[1]))
    print("=== Kombinatorik: %s ==="
          % ("gemessen" if fehler == 0 else "gemessen, aber %d Pruefung(en) FEHLGESCHLAGEN"
             % fehler))
    return 1 if fehler else 0


# The reference value against which step 0 checks byte-identity.  It comes from
# `selbsttest_kombinatorik` step 0 (decalin, all switches OFF, budget=48) and is
# thus the number BEFORE the changes of 26.08.  Carrying it here as a constant is
# the difference between "computed the same thing twice" and "computed against the
# previous state": two identical runs of the NEW code prove nothing at all.
_REF_DECALIN_FRAMES = 22
# The same role for the candidate grid -- `_pucker_space_grid(n, 2, 6)`.
_REF_GITTER = {5: 13, 6: 65, 7: 169, 8: 845}

# ===== THE SAMPLES FOR THE DISCRIMINATION TEST ====================================
# Only MULTI-RING molecules, and for a reason that carries the whole measurement: the question
# is whether two frames of the same molecule that DIFFER IN FOLD are judged the same
# by the two measures.  A single-ring molecule delivers too few frames to carry a
# pair statistic, and above all its ring share of the heavy atoms is
# near 1 -- exactly the case in which RMSD and maximum do NOT diverge.
# The effect at stake is a DILUTION effect; it needs atoms that do
# not move.  Measuring it on cyclohexane would mean defining it away.
_TRENN_PROBEN = tuple(p for p in _KOMBI_PROBEN if p[0] != "Cyclohexan")

# ===== THE DILUTION SERIES: THE INDEPENDENT VARIABLE OF THE WHOLE DESIGN ===========
#
# ⚠ THE FIRST RUN REFUTED ITS OWN SAMPLE CHOICE.  `_KOMBI_PROBEN` are pure
#   ring hydrocarbons -- measured ring share of the heavy atoms: 100 % for
#   six of eight samples.  But the effect at stake is a DILUTION
#   effect: RMSD divides the ring displacement by ALL atoms, the maximum by none.
#   At ring share 1 there is nothing to dilute, and the measurement consequently sees
#   only a factor 1.6 to 2.1 instead of the expected 4.  It did not refute the effect,
#   it DEFINED IT AWAY -- on samples in which it does not occur by construction.
#
# THE REPAIR is a series in which exactly ONE quantity varies: the same folded
# cyclohexane ring, attached to an ever larger RIGID scaffold, which moreover
# is FROZEN.  The ring share falls from 50 % to 25 %, the fold stays
# the same.  What then opens up between maximum and RMSD is the effect.
#   cyclohexyl + acene:  6 / (6 + C_acene) heavy atoms
#   benzene 50.0 % · naphthalene 37.5 % · anthracene 30.0 % · tetracene 25.0 %
#
# ⚠ WHY ACENES AND NOT OLIGOPHENYLS.  The first attempt was cyclohexyl-oligophenyl up
#   to sexiphenyl (down to 14.3 % ring share).  FAILED, and measurably so: the
#   maximum grew over the series from 1.19 to 2.41 A, although all members
#   contain THE SAME fold.  A maximum that grows with the scaffold measures the
#   scaffold -- the biaryl torsions are free and the chain flips over during the relax.
#   A fused acene does not have these degrees of freedom.
# ⚠ THE SERIES DOES NOT REACH DOWN TO 13.2 %, and that is not forced with an even larger
#   molecule (heptacene would be geometrically usable and chemically nonsense).
#   Instead the LAW is checked at these four points -- RMSD falls like
#   sqrt(ring share), the maximum stays put -- and then applied to the 1227 pairs of the
#   coupling samples.  Applying a law confirmed at four points to measured
#   pairs is something other than extending a curve.
_VERD_PROBEN = (
    ("Cyclohexylbenzol",     "C1CCCCC1c1ccccc1"),
    ("Cyclohexylnaphthalin", "C1CCCCC1c1ccc2ccccc2c1"),
    ("Cyclohexylanthracen",  "C1CCCCC1c1ccc2cc3ccccc3cc2c1"),
    ("Cyclohexyltetracen",   "C1CCCCC1c1ccc2cc3cc4ccccc4cc3cc2c1"),
)


def _xyz_schwer(txt: str):
    """Heavy-atom coordinates from a frame as `generate` returns it.

    ⚠ FROM THE OUTPUT TEXT, not from a conformer kept in parallel.  What the
      module delivers is this text; any metric that computes on something else
      measures an intermediate stage that never reaches the caller in this form.
    """
    P = []
    for ln in txt.splitlines():
        t = ln.split()
        if len(t) < 4 or t[0] == "H":
            continue
        try:
            P.append([float(t[1]), float(t[2]), float(t[3])])
        except Exception:
            continue
    return _np.array(P, dtype=float)


def _greedy_eintraege(frames, index: int, tol: float) -> int:
    """How many MANIFOLD ENTRIES remain when deduping with ``tol``.

    ``index`` 0 = largest displacement, 1 = RMSD.  Greedy and in EMISSION ORDER
    -- exactly how `generate` dedups, and exactly how the project's RMSD filters
    dedup.  An optimal covering would be a different number and a different question.
    """
    kept = []
    for P in frames:
        if any(_kabsch_max_rmsd(K, P)[index] < tol for K in kept):
            continue
        kept.append(P)
    return len(kept)


def selbsttest_trennschaerfe(proben=None, tol: float = 0.15,
                             rmsd_projekt: float = 0.30,
                             max_kombis: int = 800) -> int:
    """DOES THE MAXIMUM SEPARATE WHAT THE MEAN MERGES?  With denominator.

    THE QUESTION.  The design claims: RMSD averages a ring fold away, the largest
    displacement after Kabsch alignment does not.  The calculation for it is at
    `_kabsch_max_rmsd` (0.203 A ring displacement, 13.2 % ring share, 0.086 A RMSD versus
    0.33 A maximum -- factor 4).  But a calculation is not a measurement: it assumes
    a ring share and a displacement that can come out differently on real multi-ring
    molecules.  This test recomputes them on built frames.

    MEASURED ON PAIRS, not on frames.  "Discrimination" is a statement about
    two states, not about one; the denominator is therefore the number of PAIRS
    of frames differing in fold, and it is stated everywhere.

    ⚠ TWO COMPARISONS, and only the first isolates the INSTRUMENT:
        (A) same threshold, both 0.15 -- measures solely the difference between
            maximum and mean.
        (B) maximum 0.15 against the RMSD threshold 0.30 carried in the project -- measures
            what really happens today, but mixes instrument and threshold.
      Whoever shows only (B) can produce any desired effect through the choice of
      threshold.  Whoever shows only (A) talks past practice.

    ⚠ THE REVERSE DIRECTION IS MEASURED TOO, although it MUST be zero: the maximum is
      never smaller than the root mean square of the same deviations.  A pair that
      the maximum merges is therefore necessarily also merged by the RMSD at the same
      threshold.  If this number does NOT come out zero, a calculation error is in play
      and not a finding -- that is why it is there.
    """
    if not (_RDKIT and _np is not None):
        print("=== Trennschaerfe: RDKit fehlt, uebersprungen ==="); return 0
    import time as _time
    proben = proben or _TRENN_PROBEN
    fehler = 0
    print("=== Selbsttest: Trennschaerfe -- groesste Auslenkung gegen RMSD ===")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE", "DELFIN_FFFREE_PUCKER_DEFEKT",
             "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_XRDTOL")}
    zeilen = []
    try:
        # ===== 0 DEFAULT OFF -> BYTE-IDENTICAL.  Against the state BEFORE 26.08. ======
        for _k in _alt:
            _os.environ[_k] = "0"
        _os.environ["DELFIN_FFFREE_PUCKER_XRDTOL"] = "0.15"
        for _n, _soll in sorted(_REF_GITTER.items()):
            _ist = len(_pucker_space_grid(_n, 2, 6))
            if _ist != _soll:
                print("    ✗ 0 GITTER n=%d: %d Kandidaten statt %d" % (_n, _ist, _soll))
                fehler += 1
        if not fehler:
            print("    ✓ 0 GITTER unveraendert: n=5,6,7,8 -> %s"
                  % ", ".join(str(_REF_GITTER[k]) for k in (5, 6, 7, 8)))
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))        # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) != 0:
            print("    ? 0 VORGABE: Decalin nicht einbettbar, NICHT gemessen"); fehler += 1
        else:
            AllChem.MMFFOptimizeMolecule(_mv)
            _aus = generate(_mv, budget=48)
            _aus2 = generate(_mv, budget=48, _zaehler=_neuer_zaehler())
            if len(_aus) != _REF_DECALIN_FRAMES:
                print("    ✗ 0 VORGABE VERAENDERT: Decalin liefert %d Frames statt %d "
                      "-- der Vorgabepfad ist NICHT mehr byte-identisch"
                      % (len(_aus), _REF_DECALIN_FRAMES))
                fehler += 1
            elif _aus != _aus2:
                print("    ✗ 0 MESSMODUS VERAENDERT: %d gegen %d Frames"
                      % (len(_aus), len(_aus2)))
                fehler += 1
            else:
                print("    ✓ 0 VORGABE UNVERAENDERT: Decalin %d Frames (Referenz %d), "
                      "mit und ohne Zaehler zeichengleich"
                      % (len(_aus), _REF_DECALIN_FRAMES))
            # ... and the proof that the new switches have any REACH at all.
            # A switch that changes nothing cannot be told from one that is not
            # wired up -- has happened five times in one day in this project.
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "1"
            _os.environ["DELFIN_FFFREE_PUCKER_XRDTOL"] = "0.15"
            _mit_xrd = generate(_mv, budget=48)
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
            _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"
            _mit_def = generate(_mv, budget=48)
            _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "0"
            print("    ✓ 0 REICHWEITE der Schalter (Decalin, budget=48): Vorgabe %d "
                  "-> XRD %d -> DEFEKT %d Frames"
                  % (len(_aus), len(_mit_xrd), len(_mit_def)))
            if len(_mit_xrd) == len(_aus) and len(_mit_def) == len(_aus):
                print("      ⚠ BEIDE Schalter ohne Wirkung auf DIESER Probe -- das ist "
                      "kein Fehler, aber es beweist an Decalin nichts.  Die Reichweite "
                      "muss dann aus den Proben unten kommen.")

        # ===== 1 THE PAIR STATISTICS ==================================================
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"     # no cap: all depths
        _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"   # (1) defect filter ON
        _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"      # (2) NOT yet here
        print()
        print("    Kandidatenliste wie in der Vorgabe (SPACE=0), Kappe AUS, "
              "(1) Defektfilter AN, (2) noch AUS.")
        print("    %-22s %5s %7s %8s %8s %8s %8s"
              % ("Molekuel", "Frames", "Paare", "max~", "rmsd~", "max/rmsd", "s"))
        for name, smi, klasse in proben:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    %-22s Einbettung fehlgeschlagen" % name); continue
                try:
                    AllChem.MMFFOptimizeMolecule(m)
                except Exception:
                    AllChem.UFFOptimizeMolecule(m)
            except Exception as e:
                print("    %-22s Aufbau fehlgeschlagen: %s" % (name, type(e).__name__))
                continue
            # Pre-pass: know (b) BEFORE the build, otherwise a sample runs for hours.
            try:
                _rings = [_ring_order(m, set(r)) for r in m.GetRingInfo().AtomRings()
                          if _is_puckerable(m, r)]
                _stv = [_ring_pucker_states(m, r, set(), 0.05) for r in _rings]
            except Exception as e:
                print("    %-22s Zustaende nicht bestimmbar: %s"
                      % (name, type(e).__name__)); continue
            _prod = 1
            for _s in _stv:
                _prod *= len(_s)
            if max_kombis and (_prod - 1) > max_kombis:
                print("    %-22s UEBERSPRUNGEN: (b) = %d ueber der Testgrenze %d"
                      % (name, _prod - 1, max_kombis))
                continue
            _z = _neuer_zaehler()
            _t0 = _time.perf_counter()
            try:
                _out = generate(m, budget=10 ** 9, _zaehler=_z)
            except Exception as e:
                print("    %-22s generate() geplatzt: %s" % (name, type(e).__name__))
                continue
            _dt = _time.perf_counter() - _t0
            # The GROUND STATE is itself a manifold entry and belongs in the
            # pair set: a fold that cannot be told from the starting frame
            # is just as much a duplicate as two indistinguishable folds.
            _frames = [_xyz_schwer(_conf_to_xyz(m))] + [_xyz_schwer(x) for x, _l in _out]
            _frames = [P for P in _frames if P.size and P.shape == _frames[0].shape]
            # Ring share of the heavy atoms -- the dilution at stake.
            try:
                _ring_at = set()
                for _r in m.GetRingInfo().AtomRings():
                    _ring_at |= {int(x) for x in _r}
                _hv = [i for i in range(m.GetNumAtoms())
                       if m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _ring_pos = [k for k, i in enumerate(_hv) if i in _ring_at]
                _anteil = len(_ring_pos) / float(len(_hv)) if _hv else 0.0
            except Exception:
                _ring_pos, _anteil = [], 0.0
            _paare = []
            for _i in range(len(_frames)):
                for _j in range(_i + 1, len(_frames)):
                    _mx, _rm = _kabsch_max_rmsd(_frames[_i], _frames[_j])
                    _paare.append((_mx, _rm))
            if not _paare:
                print("    %-22s %5d %7d   -- kein Paar, keine Aussage"
                      % (name, len(_frames), 0))
                continue
            _mxs = sorted(p[0] for p in _paare)
            _rms = sorted(p[1] for p in _paare)
            _med = lambda v: v[len(v) // 2]
            print("    %-22s %5d %7d %8.3f %8.3f %8.1f %8.1f"
                  % (name, len(_frames), len(_paare), _med(_mxs), _med(_rms),
                     (_med(_mxs) / _med(_rms)) if _med(_rms) > 1e-9 else float("inf"),
                     _dt))
            print("        %-26s Ringanteil an den schweren Atomen %4.1f %% (%d von %d)"
                  " · verworfen: Kollision %d · Winkel %d · Bindung %d · TFD %d"
                  % (klasse, 100.0 * _anteil, len(_ring_pos),
                     len(_frames[0]) if _frames else 0,
                     int(_z.get("kollision", 0)), int(_z.get("winkel", 0)),
                     int(_z.get("bindung", 0)), int(_z.get("tfd_doppelt", 0))))
            # ---- RUN B: THE SAME SAMPLE WITH (2) ON.  That is the number for (c), and
            #      it is BUILT and not extrapolated from run A: with (2) on,
            #      `_ring_pucker_states` already dedups per ring, so the cross product (b)
            #      is a different one.  Whoever greedily re-creates that from the frames
            #      of run A measures the re-creation -- the mistake this project
            #      has already let pass as a finding several times.
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "1"
            _zb = _neuer_zaehler()
            _t1 = _time.perf_counter()
            try:
                _outb = generate(m, budget=10 ** 9, _zaehler=_zb)
            except Exception:
                _outb, _zb = [], _neuer_zaehler()
            _dtb = _time.perf_counter() - _t1
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
            zeilen.append({"name": name, "klasse": klasse, "paare": _paare,
                           "frames": _frames, "anteil": _anteil, "dt": _dt,
                           "b": int(_z.get("aufzaehlung", 0)),
                           "c": int(_z.get("tor_ueberlebt", 0)),
                           "d": len(_out),
                           "zustaende": list(_z.get("zustaende_je_ring") or []),
                           "bindung": int(_z.get("bindung", 0)),
                           "kollision": int(_z.get("kollision", 0)),
                           "winkel": int(_z.get("winkel", 0)),
                           "b2": int(_zb.get("aufzaehlung", 0)),
                           "c2": int(_zb.get("tor_ueberlebt", 0)),
                           "d2": len(_outb), "dt2": _dtb,
                           "xrd2": int(_zb.get("xrd_doppelt", 0)),
                           "zustaende2": list(_zb.get("zustaende_je_ring") or []),
                           "mol": m, "smi": smi})
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    # ===== 1b THE DILUTION SERIES ==================================================
    # Run separately, reported separately -- it answers a DIFFERENT question than
    # the coupling samples above (there: which pairs are separated by whom; here: WHAT the
    # difference depends on at all).  Thrown together, both would be unreadable.
    verd, _verd_tfd = [], []
    _altv = {k: _os.environ.get(k) for k in
             ("DELFIN_FFFREE_PUCKER_FULL", "DELFIN_FFFREE_PUCKER_DEFEKT",
              "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_SPACE",
              "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_NPHASE")}
    try:
        # ⚠ SPACE GRID ON so that the ONE ring delivers enough states -- a
        #   pair statistic from three pairs would be none.  It is the same ring in all
        #   six molecules, so this changes nothing about the independent variable.
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
        print()
        print("    ===== 1b VERDUENNUNGSREIHE: derselbe Ring, wachsendes Geruest =====")
        print("      Geruest EINGEFROREN (`frozen`) -- s. den Kommentar bei "
              "_VERD_PROBEN: ohne das misst die Reihe die Geruestbewegung mit.")
        print("      %-24s %6s %7s %7s %7s %6s %8s %8s %9s"
              % ("Molekuel", "schwer", "Anteil", "(a)0,05", "(a)0,005", "Paare",
                 "max~", "rmsd~", "max/rmsd"))
        for name, smi in _VERD_PROBEN:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("      %-24s Einbettung fehlgeschlagen" % name); continue
                AllChem.MMFFOptimizeMolecule(m)
            except Exception as e:
                print("      %-24s Aufbau ausgefallen: %s"
                      % (name, type(e).__name__)); continue
            try:
                _ring_at, _ring_ord = set(), []
                for _r in m.GetRingInfo().AtomRings():
                    if _is_puckerable(m, _r):
                        _ring_at |= {int(x) for x in _r}
                        if not _ring_ord:
                            _ring_ord = _ring_order(m, set(_r))
                _hv = [i for i in range(m.GetNumAtoms())
                       if m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _rp = [k for k, i in enumerate(_hv) if i in _ring_at]
                # free = the ring and ITS hydrogens; everything else is fixed.
                _frei = set(_ring_at)
                for _i in list(_ring_at):
                    for _nb in m.GetAtomWithIdx(int(_i)).GetNeighbors():
                        if _nb.GetSymbol() == "H":
                            _frei.add(int(_nb.GetIdx()))
                _fr = set(range(m.GetNumAtoms())) - _frei
            except Exception as e:
                print("      %-24s Ringmenge unbestimmbar: %s"
                      % (name, type(e).__name__)); continue
            if not _ring_at or not _hv:
                print("      %-24s kein faltbarer Ring" % name); continue
            # ===== TFD AVERAGES TOO -- AND IT IS THE INSTRUMENT IN USE ================
            #
            # This series uncovered it itself during the build: anthracene and tetracene
            # delivered ZERO folds, and not at the realism gate but already
            # at the states PER RING -- (b) was 0, so the cross product was 1x nothing.
            # The ring is the same as in cyclohexylbenzene, which delivers 55 pairs.
            #
            # The cause is the same disease one level up: TFD compares ALL
            # torsions of the molecule and AVERAGES over them.  A large rigid
            # scaffold brings along many torsions that do not change -- the contribution
            # of the six ring torsions is divided by them and falls below the
            # threshold 0.05.  The fold vanishes in the mean, exactly as with RMSD.
            #
            # ⚠ THIS IS NOT A SIDE FINDING.  TFD is the dedup measure that runs in the
            #   build TODAY.  If it gets blurrier with ligand size, then
            #   the manifold loses folds exactly on the systems that matter
            #   -- large ligands, small ring share.
            # It is measured with the only means that separates the two causes:
            # the same candidates, two thresholds.  If (a) rises sharply at 0.005, it was
            # the threshold (i.e. the dilution); if it stays the same, the
            # states are really not there.
            try:
                _a05 = len(_ring_pucker_states(m, _ring_ord, _fr, 0.05))
                _a005 = len(_ring_pucker_states(m, _ring_ord, _fr, 0.005))
            except Exception:
                _a05 = _a005 = -1
            _zv = _neuer_zaehler()
            try:
                _out = generate(m, frozen=_fr, budget=10 ** 9, _zaehler=_zv)
            except Exception as e:
                print("      %-24s generate() geplatzt: %s"
                      % (name, type(e).__name__)); continue
            _frames = [_xyz_schwer(_conf_to_xyz(m))] + [_xyz_schwer(x) for x, _l in _out]
            _frames = [P for P in _frames if P.size and P.shape == _frames[0].shape]
            if len(_frames) < 2:
                # ⚠ A ZERO GETS ITS REASON.  "generated and rejected" looks from
                #   outside exactly like "never built" -- the fallacy that on
                #   14.08. rendered the firing census worthless.
                print("      %-24s %6d %6.1f%% %7d %7d      -- nur %d Frame(e), (b)=%d"
                      % (name, len(_hv),
                         100.0 * len(_rp) / float(len(_hv)) if _hv else 0.0,
                         _a05, _a005, len(_frames), int(_zv.get("aufzaehlung", 0))))
                _verd_tfd.append({"name": name, "n": len(_hv),
                                  "anteil": len(_rp) / float(len(_hv)) if _hv else 0.0,
                                  "a05": _a05, "a005": _a005})
                continue
            _anteil = len(_rp) / float(len(_hv))
            _mx, _rm, _pa = [], [], []
            _geruest = 0.0
            _nicht_ring = [k for k in range(len(_hv)) if k not in set(_rp)]
            for _i in range(len(_frames)):
                for _j in range(_i + 1, len(_frames)):
                    _v = _kabsch_abweichungen(_frames[_i], _frames[_j])
                    _a = float(_v.max())
                    _b = float(_np.sqrt(float((_v ** 2).mean())))
                    _mx.append(_a); _rm.append(_b); _pa.append((_a, _b))
                    # ⚠ FREEZE CHECK WITHOUT KABSCH, and that is the point.  Both frames
                    #   stand in the SAME reference frame -- nothing was re-embedded,
                    #   only relaxed.  "Does the scaffold stand?" is thus a question to the
                    #   RAW coordinates.  After Kabsch it would be unanswerable: the
                    #   alignment minimises the total RMSD and distributes the error over
                    #   ALL atoms, including held ones -- a first version reported
                    #   exactly from that 1.4 A of scaffold motion that did not exist.
                    if _nicht_ring:
                        _roh = _np.linalg.norm(_frames[_i][_nicht_ring]
                                               - _frames[_j][_nicht_ring], axis=1)
                        _geruest = max(_geruest, float(_roh.max()))
            _mx.sort(); _rm.sort()
            _m1, _r1 = _mx[len(_mx) // 2], _rm[len(_rm) // 2]
            print("      %-24s %6d %6.1f%% %7d %7d %6d %8.3f %8.3f %9.1f"
                  % (name, len(_hv), 100.0 * _anteil, _a05, _a005, len(_pa), _m1, _r1,
                     (_m1 / _r1) if _r1 > 1e-9 else float("inf")))
            _verd_tfd.append({"name": name, "n": len(_hv), "anteil": _anteil,
                              "a05": _a05, "a005": _a005})
            verd.append({"name": name, "n": len(_hv), "anteil": _anteil,
                         "max": _m1, "rmsd": _r1, "paare": _pa, "geruest": _geruest})
    finally:
        for k, v in _altv.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v
    if verd:
        # ---- CHECK: DID THE FREEZING HOLD?  Raw, without alignment.
        _gmax = max(r["geruest"] for r in verd)
        if _gmax < 1e-6:
            print("      ✓ EINFROSTPROBE (rohe Koordinaten): groesste Geruestauslenkung "
                  "%.2e A -- das Geruest steht, die Reihe isoliert die Verduennung."
                  % _gmax)
        else:
            print("      ✗ EINFROSTPROBE: Geruest bewegt sich um bis zu %.3f A -- die "
                  "Reihe misst NICHT nur Verduennung." % _gmax)
            fehler += 1
    if _verd_tfd:
        # ---- THE ACTUAL RESULT OF THIS SERIES, and it was not the one sought.
        print()
        print("      ===== TFD MITTELT GENAUSO -- und TFD laeuft heute im Bau =====")
        for r in _verd_tfd:
            print("        %-24s %2d schwere Atome, Ringanteil %5.1f %% -> (a) %d bei "
                  "Schwelle 0,05 · %d bei 0,005"
                  % (r["name"], r["n"], 100.0 * r["anteil"], r["a05"], r["a005"]))
        _v0, _vn = _verd_tfd[0], _verd_tfd[-1]
        if _vn["a05"] < _v0["a05"] and _vn["a005"] > _vn["a05"]:
            print("      ⇒ BEFUND: mit wachsendem Geruest faellt (a) bei Schwelle 0,05 "
                  "von %d auf %d -- DERSELBE Ring, dieselben Kandidaten.  Bei 0,005 "
                  "kommen die Zustaende zurueck (%d).  Die Faltungen sind also DA und "
                  "werden von TFD verschmolzen, nicht vom Generator ausgelassen."
                  % (_v0["a05"], _vn["a05"], _vn["a005"]))
            print("        URSACHE: TFD vergleicht ALLE Torsionen und mittelt ueber sie. "
                  "Ein grosses starres Geruest bringt unbewegte Torsionen mit; der "
                  "Beitrag der Ringtorsionen wird durch sie geteilt.  Das ist Wort fuer "
                  "Wort der RMSD-Fehler, eine Ebene hoeher -- und TFD ist das Mass, das "
                  "im Bau ENTSCHEIDET.")
            print("        ⚠ TRAGWEITE: der Effekt waechst mit der LIGANDGROESSE.  Er "
                  "trifft also am haertesten die Systeme, um die es geht -- grosse "
                  "Liganden, kleiner Ringanteil.  Die 59,7 % Ringidentitaeten mit nur "
                  "EINER Faltung (16.08.) haben hier eine kandidatenfaehige Ursache, "
                  "die nichts mit dem Generator zu tun hat.")
        else:
            print("      ⇒ Kein Verduennungsmuster in (a) -- TFD verschmilzt hier nicht.")
        _vp = [p for r in verd for p in r["paare"]]
        _vN = len(_vp)
        if _vN:
            _vA = sum(1 for a, b in _vp if b < tol <= a)
            _vB = sum(1 for a, b in _vp if b < rmsd_projekt <= a)
            print("      ⇒ Trennschaerfe auf dieser Reihe (Nenner %d Paare, Ringanteil "
                  "%.0f bis %.0f %%): (A) %d = %.1f %% · (B) %d = %.1f %%"
                  % (_vN, 100.0 * min(r["anteil"] for r in verd),
                     100.0 * max(r["anteil"] for r in verd),
                     _vA, 100.0 * _vA / _vN, _vB, 100.0 * _vB / _vN))

    if not zeilen:
        print("=== Trennschaerfe: nichts gemessen ==="); return 1

    _alle = [p for r in zeilen for p in r["paare"]]
    _N = len(_alle)
    _A = sum(1 for mx, rm in _alle if rm < tol <= mx)
    _B = sum(1 for mx, rm in _alle if rm < rmsd_projekt <= mx)
    _G = sum(1 for mx, rm in _alle if mx < tol <= rm)
    print()
    print("    ===== 1 TRENNSCHAERFE (Nenner: %d Paare faltungsverschiedener Frames "
          "ueber %d Molekuele) =====" % (_N, len(zeilen)))
    print("      (A) gleiche Schwelle %.2f A -- nur das INSTRUMENT:" % tol)
    print("          Maximum trennt, RMSD verschmilzt:  %6d von %6d = %5.1f %%"
          % (_A, _N, 100.0 * _A / _N if _N else 0.0))
    print("      (B) Maximum %.2f gegen die Projekt-RMSD-Schwelle %.2f -- die PRAXIS:"
          % (tol, rmsd_projekt))
    print("          Maximum trennt, RMSD verschmilzt:  %6d von %6d = %5.1f %%"
          % (_B, _N, 100.0 * _B / _N if _N else 0.0))
    print("      GEGENRICHTUNG (muss 0 sein, das Maximum ist nie kleiner als das "
          "quadratische Mittel): %d" % _G)
    # ---- THE NUMBER THAT APPLIES TO REAL SYSTEMS.
    # ⚠ THE SAMPLES ABOVE HAVE 86 TO 100 %% RING SHARE -- the dilution effect does
    #   not occur in them by construction.  What they deliver is the LOWER BOUND of the
    #   discriminating power.  A real DELFIN system has a metal core, aromatic
    #   backbones and substituents; the premise names 13.2 % ring share.
    #   What is applied is the law CONFIRMED above at four points: the maximum stays,
    #   the RMSD falls like sqrt(share).  No new molecule, no curve extended
    #   -- the same measured pairs, with the denominator of a real system.
    # ---- THE RING SHARE IS THE INDEPENDENT VARIABLE, and both points are MEASURED.
    # ⚠ AN EXTRAPOLATION TO 13.2 % ONCE STOOD HERE -- twice, and wrong both times.
    #   (i) "RMSD falls like sqrt(ring share)" holds only under IDENTICAL
    #   alignment; Kabsch, however, aligns and redistributes the error.  (ii) Padding
    #   the pair with rigid copies and re-superimposing only shifted the error:
    #   the scaffold lay on the molecule and bound the alignment so hard
    #   that the maximum rose from 0.94 to 3.95 A -- what was measured was the change
    #   of the ALIGNMENT REGIME, not the dilution.  Both attempts are removed.
    #   What remains are two MEASURED points on real molecules; series 1b
    #   does not reach 13.2 %, and the reason for that is itself the finding (TFD
    #   no longer delivers any frames there).
    _vp2 = [p for r in verd for p in r["paare"]]
    if _vp2:
        _vA2 = sum(1 for a, b in _vp2 if b < tol <= a)
        _vB2 = sum(1 for a, b in _vp2 if b < rmsd_projekt <= a)
        print("      ⇒ GEGEN DEN RINGANTEIL, beide Punkte gemessen:")
        print("        ~100 %% Ringanteil (Kopplungsproben, %4d Paare): (A) %5.1f %% · "
              "(B) %5.1f %%" % (_N, 100.0 * _A / _N, 100.0 * _B / _N))
        print("        38-50 %% Ringanteil (Verduennungsreihe, %4d Paare): (A) %5.1f %% "
              "· (B) %5.1f %%" % (len(_vp2), 100.0 * _vA2 / len(_vp2),
                                  100.0 * _vB2 / len(_vp2)))
        print("        ⇒ die Trennschaerfe WAECHST mit der Verduennung (A: %.1f -> "
              "%.1f %%).  Die Kopplungsproben sind damit die UNTERGRENZE, nicht die "
              "Antwort." % (100.0 * _A / _N, 100.0 * _vA2 / len(_vp2)))
    if _G:
        print("      ✗ GEGENRICHTUNG NICHT NULL -- Rechenfehler, kein Befund."); fehler += 1

    # ---- 2 WHAT THAT MEANS FOR THE NUMBER OF MANIFOLD ENTRIES.
    # ⚠ A percentage over pairs does not yet say how many ENTRIES arise:
    #   dedup is greedy and not cleanly transitive, three pairwise-close frames
    #   can become one or two.  So count instead of extrapolating.
    print()
    print("    ===== 2 EINTRAEGE JE MOLEKUEL, je Kriterium =====")
    print("      %-22s %7s %8s %8s %8s %8s"
          % ("Molekuel", "roh", "max%.2f" % tol, "rmsd%.2f" % tol,
             "rmsd%.2f" % rmsd_projekt, "rmsd0.50"))
    _sum = {"roh": 0, "max": 0, "r_gleich": 0, "r_proj": 0, "r_50": 0}
    for r in zeilen:
        _f = r["frames"]
        _e = (len(_f), _greedy_eintraege(_f, 0, tol), _greedy_eintraege(_f, 1, tol),
              _greedy_eintraege(_f, 1, rmsd_projekt), _greedy_eintraege(_f, 1, 0.50))
        print("      %-22s %7d %8d %8d %8d %8d" % ((r["name"],) + _e))
        for _k, _v in zip(("roh", "max", "r_gleich", "r_proj", "r_50"), _e):
            _sum[_k] += _v
    print("      %-22s %7d %8d %8d %8d %8d"
          % ("SUMME", _sum["roh"], _sum["max"], _sum["r_gleich"], _sum["r_proj"],
             _sum["r_50"]))
    if _sum["max"] > 0:
        print("      ⇒ RMSD bei %.2f A behaelt %d von %d Eintraegen, die das Maximum "
              "bei derselben Schwelle als UNTERSCHEIDBAR fuehrt (%.0f %%).  Die "
              "Differenz %d sind Faltungen, die eine RMSD-Entdopplung LOESCHT."
              % (tol, _sum["r_gleich"], _sum["max"],
                 100.0 * _sum["r_gleich"] / _sum["max"],
                 _sum["max"] - _sum["r_gleich"]))
        print("      ⇒ bei der Projektschwelle %.2f A: %d von %d (%.0f %%), "
              "geloescht %d." % (rmsd_projekt, _sum["r_proj"], _sum["max"],
                                 100.0 * _sum["r_proj"] / _sum["max"],
                                 _sum["max"] - _sum["r_proj"]))
    # ---- 3 WHAT (1)+(2) COST THE COMBINATORICS -- OR SAVE.
    # ⚠ THE DECISIVE DISTINCTION, and it is easy to miss: (1) and (2)
    #   can act at TWO places, and only one of them saves compute time.
    #     PER RING  (`_ring_pucker_states`)  -> lowers (a), hence (b) RAISED TO A POWER:
    #                                           the only place where anything gets cheaper.
    #     PER COMBINATION (`generate`)       -> lowers only the number of ENTRIES.  The
    #                                           relax is already paid for there; the gate
    #                                           selects, it saves nothing.
    #   A filter that acts only at the bottom does NOT make the complete fold
    #   affordable, no matter how sharp it is.  That is why (a) and (b) stand side by side here.
    print()
    print("    ===== 3 KOMBINATORIK MIT (1)+(2) -- (a) je Ring und (b) das Produkt =====")
    print("      %-22s %-12s %-12s %8s %8s %7s %7s"
          % ("Molekuel", "(a) ohne (2)", "(a) mit (2)", "(b) ohne", "(b) mit",
             "(d) ohne", "(d) mit"))
    _sb = _sb2 = _sd = _sd2 = 0
    for r in zeilen:
        print("      %-22s %-12s %-12s %8d %8d %7d %7d"
              % (r["name"],
                 "x".join(str(v) for v in r["zustaende"]) or "-",
                 "x".join(str(v) for v in r["zustaende2"]) or "-",
                 r["b"], r["b2"], r["d"], r["d2"]))
        _sb += r["b"]; _sb2 += r["b2"]; _sd += r["d"]; _sd2 += r["d2"]
        # ⚠ A ZERO IN (d) MUST NAME ITS REASON, otherwise it reads like a
        #   defect.  For norbornane it is the INTENDED result: its largest
        #   fold displacement lies at 0.139 A, i.e. BELOW the resolution threshold.
        #   A rigidly bridged bicycle HAS no second fold -- (2) says exactly
        #   that, and the manifold keeps the ground state (which never passes through this
        #   gate).  From the crystallographer's point of view one entry is right, not three.
        if r["d"] > 0 and r["d2"] == 0:
            print("        %-20s (d) faellt auf 0: alle %d Faltungen liegen unter %.2f A "
                  "Maximalauslenkung -- ununterscheidbar vom Grundzustand, EIN Eintrag."
                  % (r["name"], r["d"], tol))
    print("      %-22s %-12s %-12s %8d %8d %7d %7d"
          % ("SUMME", "", "", _sb, _sb2, _sd, _sd2))
    if _sb:
        print("      ⇒ (b), die AUFZAEHLUNG und damit der PREIS: %d -> %d = %+.1f %%."
              % (_sb, _sb2, 100.0 * (_sb2 - _sb) / _sb))
        print("      ⇒ (d), die EINTRAEGE und damit das ERGEBNIS: %d -> %d = %+.1f %%."
              % (_sd, _sd2, 100.0 * (_sd2 - _sd) / max(1, _sd)))
        # ⚠ PARTITION, NOT AN AVERAGE.  An overall percentage over eight samples cannot
        #   say whether (2) saves a little everywhere or a lot on two samples
        #   and nothing at all on six -- and those are completely different mechanisms.
        #   The second case would be NO general cost lever, but a finding
        #   about a class.  The split is by the number of samples with an effect.
        _wirkt = [r for r in zeilen if r["b2"] < r["b"]]
        _still = [r for r in zeilen if r["b2"] >= r["b"]]
        print("      ⇒ PARTITION: (2) senkt (b) bei %d von %d Proben (%s); bei den "
              "anderen %d aendert sie (a) um keinen einzigen Zustand (%s)."
              % (len(_wirkt), len(zeilen),
                 ", ".join(r["name"] for r in _wirkt) or "keiner", len(_still),
                 ", ".join(r["name"] for r in _still) or "keine"))
        if len(_wirkt) <= len(zeilen) // 2:
            print("      ⇒ URTEIL: (2) ist KEIN allgemeiner Kostenhebel.  Sie greift "
                  "dort, wo Ringzustaende ohnehin fast entartet sind (die VERBRUECKTEN "
                  "Proben -- ein verbrueckter Ring KANN kaum falten), und nirgends "
                  "sonst.  Das Maximum kennt keinen Nenner: genau die Eigenschaft, die "
                  "es gegen Verduennung unempfindlich macht, hindert es daran, echte "
                  "Ringmulden zusammenzuziehen.  Es macht die Auswahl RICHTIG, nicht "
                  "BILLIG.")
        else:
            print("      ⇒ URTEIL: (2) senkt (a) und damit (b) bei der MEHRHEIT der "
                  "Proben -- ein echter Kostenhebel, nicht nur eine Korrektur.")
        _sbi = sum(r["bindung"] for r in zeilen)
        _swi = sum(r["winkel"] for r in zeilen)
        _skl = sum(r["kollision"] for r in zeilen)
        print("      ⇒ (1) DEFEKTFILTER, aufgeschluesselt (Nenner %d Kombinationen): "
              "Kollision %d (%.1f %%) · Winkel %d (%.1f %%) · Bindung %d (%.1f %%)."
              % (_sb, _skl, 100.0 * _skl / _sb, _swi, 100.0 * _swi / _sb,
                 _sbi, 100.0 * _sbi / _sb))
        if _skl == 0:
            print("        ⚠ Das Kollisionstor feuert NULL mal -- dieselbe Nullreichweite "
                  "wie am 26.08. (0 von 17 754).  Ein zweites Mal gemessen, ein zweites "
                  "Mal null: der Name `Kollisionstor` beschreibt keinen wirksamen Filter.")
        if _sbi:
            print("        ✓ Das NEUE Bindungstor feuert %d mal -- es ist verdrahtet und "
                  "hat Reichweite; es sieht genau den Bruch, den das Selbstgate per "
                  "Konstruktion fuer 'nicht gebunden' haelt." % _sbi)
    # ---- 4 IS STAGE (3), THE ENERGY, NEEDED?
    # ⚠ THE QUESTION IS NOT "would energy be nice", but "does it solve the problem that
    #   (1) and (2) leave open".  And the problem is the PRICE (b), not the number
    #   of entries (d).  An energy is -- like every other gate here -- evaluated AFTER the
    #   relax; so it cannot lower (b) at all.  A mechanism that
    #   by construction does not reach the bottleneck is not built, but named.
    print()
    print("    ===== 4 BRAUCHT ES DIE ENERGIE? =====")
    print("      Der Engpass ist (b) = %d Kombinationen, jede mit einem Relax BEVOR "
          "irgendein Tor sie sieht." % _sb)
    print("      Eine Energieauswahl wird an derselben Stelle ausgewertet wie (1) und "
          "(2) -- nach dem Relax.  Sie kann (b) also per Konstruktion nicht senken.")
    print("      ⇒ ENERGIE NICHT GEBAUT.  Sie wuerde das Ergebnis weiter ausduennen "
          "(%d Eintraege) und den Preis unveraendert lassen.  Der einzige Ort, an dem "
          "etwas zu sparen ist, ist (a) -- die Zustaende JE RING, vor dem Kreuzprodukt."
          % _sd2)
    print("=== Trennschaerfe: %s ==="
          % ("gemessen" if fehler == 0 else "%d Pruefung(en) FEHLGESCHLAGEN" % fehler))
    return 1 if fehler else 0


# The unsubstituted calibration rings.  They are the ONLY place where the
# meaning of the threshold can be checked: there the molecule has no scaffold that
# could dilute, and both measures must therefore say THE SAME.  If that fails,
# the ring-local version is not "calibrated differently", but a different instrument.
_LOKAL_KALIBER = (("Cyclopentan", "C1CCCC1"), ("Cyclohexan", "C1CCCCC1"),
                  ("Cycloheptan", "C1CCCCCC1"), ("Cyclooctan", "C1CCCCCCC1"))


def _lokal_ringlage(mol):
    """(first puckerable ring in ring order, ring atoms, frozen scaffold).

    The same preparation as in the dilution series of `selbsttest_trennschaerfe`
    -- ⚠ and that is the purpose: the repair must be checked on THE SAME measurement
    that showed the error.  A second, slightly different preparation
    would compare two measurements instead of two measures.
    """
    ring_at, ring_ord = set(), []
    for r in mol.GetRingInfo().AtomRings():
        if _is_puckerable(mol, r):
            ring_at |= {int(x) for x in r}
            if not ring_ord:
                ring_ord = _ring_order(mol, set(r))
    if not ring_ord:
        return None, None, None
    frei = set(ring_at)
    for i in list(ring_at):
        for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            if nb.GetSymbol() == "H":
                frei.add(int(nb.GetIdx()))
    return ring_ord, ring_at, set(range(mol.GetNumAtoms())) - frei


def selbsttest_tfd_lokal() -> int:
    """RE-RUNS THE DILUTION TABLE -- global against ring-local, the same samples.

    HISTORY.  `selbsttest_trennschaerfe` threw off a finding while being built
    that it had not been looking for at all: the same cyclohexane ring yields fewer
    and fewer states on a growing rigid acene (11 -> 6 -> 1 -> 1 at TFD 0.05), and at
    0.005 they come back (62 -> 49 -> 23 -> 10).  So the folds ARE THERE and
    are merged by the dedup measure.

    THIS TEST ANSWERS FOUR QUESTIONS, and in this order, because each
    next one would be pointless if the one before fails:
        0  Is the default path unchanged?             (grid, decalin)
        1  WHAT exactly causes the dilution?          (RDKit's weights, computed)
        2  Does the symmetry folding survive?         (unsubstituted ring)
        3  Does the gradient vanish?                  (the four acenes)

    ⚠ STEP 2 IS THE ABORT CRITERION, not step 3.  On 26.08. a replacement for TFD
      already died because it did not fold in the molecular symmetry
      -- n=5 went from 3,3,3 to 9,13,14.  A version that removes the gradient
      and in doing so splits the five-ring is NOT a repair, but the same
      fallacy with a different sign.

    Invocation:  python -m delfin.manta._ring_pucker tfdlokal
    """
    if not (_RDKIT and _np is not None):
        print("=== Ringlokale TFD: RDKit fehlt, uebersprungen ==="); return 0
    from rdkit.Chem import TorsionFingerprints as _TF
    fehler = 0
    # ⚠ BEFORE the `try`, not inside it.  If step 0 blows up, the verdict below would
    #   otherwise run into a NameError -- a crash that looks like "no finding".
    _reihe = []
    print("=== Selbsttest: die ringlokale TFD ===")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE", "DELFIN_FFFREE_PUCKER_DEFEKT",
             "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_CPDEDUP",
             "DELFIN_FFFREE_PUCKER_TFD_LOCAL", "DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI",
             "DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR")}
    try:
        # ===== 0 DEFAULT OFF -> BYTE-IDENTICAL =======================================
        for _k in _alt:
            _os.environ[_k] = "0"
        _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR"] = ""
        for _n, _soll in sorted(_REF_GITTER.items()):
            _ist = len(_pucker_space_grid(_n, 2, 6))
            if _ist != _soll:
                print("  ✗ 0 GITTER n=%d: %d Kandidaten statt %d" % (_n, _ist, _soll))
                fehler += 1
        if not fehler:
            print("  ✓ 0 GITTER unveraendert: n=5,6,7,8 -> %s"
                  % ", ".join(str(_REF_GITTER[k]) for k in (5, 6, 7, 8)))
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))        # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) != 0:
            print("  ? 0 VORGABE: Decalin nicht einbettbar, NICHT gemessen"); fehler += 1
        else:
            AllChem.MMFFOptimizeMolecule(_mv)
            _aus = generate(_mv, budget=48)
            if len(_aus) != _REF_DECALIN_FRAMES:
                print("  ✗ 0 VORGABE VERAENDERT: Decalin liefert %d Frames statt %d -- "
                      "der Vorgabepfad ist NICHT mehr byte-identisch"
                      % (len(_aus), _REF_DECALIN_FRAMES))
                fehler += 1
            else:
                print("  ✓ 0 VORGABE UNVERAENDERT: Decalin %d Frames (Referenz %d)"
                      % (len(_aus), _REF_DECALIN_FRAMES))
            # ⚠ A SWITCH WITHOUT REACH CANNOT BE DISTINGUISHED FROM AN UNWIRED
            #   ONE.  Happened five times in one day in this project --
            #   that is why the counter-check stands right next to the identity.
            # ⚠ THE TWO SWITCHES INDIVIDUALLY, never together.  Their effects have
            #   opposite signs; measured jointly, the sum would give a
            #   number from which no share can be computed back.
            _za = _neuer_zaehler()
            generate(_mv, budget=48, _zaehler=_za)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
            _zm = _neuer_zaehler()
            _mit = generate(_mv, budget=48, _zaehler=_zm)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"] = "1"
            _zk = _neuer_zaehler()
            _kom = generate(_mv, budget=48, _zaehler=_zk)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"] = "0"
            _zur = generate(_mv, budget=48)
            print("  %s 0 REICHWEITE (Decalin, budget=48), die Schalter EINZELN:"
                  % ("✓" if (len(_mit) != len(_aus) or len(_kom) != len(_aus)) else "⚠"))
            print("      AUS                %2d Frames, Zustaende je Ring %s"
                  % (len(_aus), _za["zustaende_je_ring"]))
            print("      nur je RING        %2d Frames, Zustaende je Ring %s"
                  % (len(_mit), _zm["zustaende_je_ring"]))
            print("      nur je KOMBINATION %2d Frames, Zustaende je Ring %s"
                  % (len(_kom), _zk["zustaende_je_ring"]))
            print("      zurueck auf AUS    %2d Frames" % len(_zur))
            if len(_zur) != len(_aus):
                print("      ✗ NICHT ZURUECKSCHALTBAR -- ein Schalter hinterlaesst Zustand")
                fehler += 1
            if len(_kom) == len(_aus):
                # ⚠ A SWITCH WITHOUT A MEASURED EFFECT IS NAMED AS SUCH.  On
                #   decalin the zero is even PREDICTABLE -- two equivalent rings,
                #   no acyclic torsion, weights 1:1: global and ring-local compute
                #   literally the same number there.  That explains the zero, but it does
                #   not substantiate the switch.  Whoever uses it measures it first.
                print("      ⚠ DER KOMBINATIONSSCHALTER ist auf Decalin wirkungslos, und "
                      "das ist vorhersagbar: zwei gleichwertige Ringe, keine acyclische "
                      "Torsion, Gewichte 1:1 -- beide Masse rechnen dieselbe Zahl.  Er "
                      "ist damit in dieser Datei NICHT BELEGT; kein Test zeigt bisher "
                      "eine Wirkung von ihm.")
            if len(_mit) == len(_aus) and len(_kom) == len(_aus):
                print("      ⚠ BEIDE ohne Wirkung -- die Reichweite muss dann aus "
                      "Schritt 3 kommen.")
            else:
                # ⚠ HERE THE NUMBER FALLS, WHILE IN STEP 3 IT RISES, and that is
                #   no contradiction, but THE SAME statement from two sides.
                #   Ring-local means "only the torsion of THIS ring" -- and that removes
                #   TWO contaminations at once:
                #     (a) the rigid scaffold in the DENOMINATOR -> acenes, states RISE
                #     (b) the motion of the NEIGHBOURING RING -> decalin, states FALL
                #   In (b) the global comparison counted states of ring 1 as
                #   states of ring 0; the cross product counts them AFTERWARDS once
                #   more.  Decalin has no scaffold to dilute (ring share 100 %,
                #   no acyclic torsion, weights 1:1), so only (b) remains here.
                # ⚠ NOT PROVEN by this is that the dropped frames WERE
                #   duplicates -- shown is only WHERE the number changes.  Whoever wants
                #   the verdict needs the eye, not this test.
                print("      ⚠ HIER FAELLT die Zahl, in Schritt 3 STEIGT sie.  Dieselbe "
                      "Aussage von zwei Seiten: ringlokal entfernt das Geruest aus dem "
                      "Nenner (Acene: mehr Zustaende) UND die Bewegung des Nachbarrings "
                      "aus der Zustandszahl eines Rings (Decalin: weniger).  Decalin hat "
                      "kein Geruest -- Ringanteil 100 %, keine acyclische Torsion, "
                      "Gewichte 1:1 -- also bleibt hier nur der zweite Anteil.")

        # ===== 1 WHAT CAUSES THE DILUTION -- RDKIT'S OWN WEIGHTS =====================
        # ⚠ COMPUTED, NOT ESTIMATED.  `CalculateTFD` forms sum(d_i*w_i)/sum(w_i).
        #   If only ONE ring moves, d_ring * w_ring / sum(w) remains -- the quotient
        #   w_ring/sum(w) therefore IS the dilution factor, without any model assumption.
        print()
        print("  ===== 1 DER VERDUENNUNGSFAKTOR STEHT IN RDKITS GEWICHTEN =====")
        print("    %-24s %7s %6s %6s %10s %12s"
              % ("Molekuel", "Anteil", "nring", "ring", "w_R/sum(w)", "Verduennung"))
        _wfak = {}
        for _name, _smi in _LOKAL_KALIBER + _VERD_PROBEN:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    print("    %-24s Einbettung fehlgeschlagen" % _name); continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                if not _ro:
                    print("    %-24s kein faltbarer Ring" % _name); continue
                _tl, _tlr = _TF.CalculateTorsionLists(_m)
                _w = _TF.CalculateTorsionWeights(_m)
                _ziel = frozenset(int(a) for a in _ro)
                _k = next((i for i, (_q, _d) in enumerate(_tlr)
                           if frozenset(int(t[0]) for t in _q) == _ziel), None)
                if _k is None:
                    print("    %-24s Ring NICHT in RDKits Ringliste -- "
                          "die Zuordnung ueber die Atommenge greift nicht" % _name)
                    fehler += 1
                    continue
                _hv = [i for i in range(_m.GetNumAtoms())
                       if _m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _wr = _w[len(_tl) + _k] / sum(_w)
                _wfak[_name] = _wr
                print("    %-24s %6.1f%% %6d %6d %10.4f %11.1fx"
                      % (_name, 100.0 * len(_rat) / len(_hv), len(_tl), len(_tlr),
                         _wr, 1.0 / _wr))
            except Exception as _e:
                print("    %-24s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        for _name, _ in _LOKAL_KALIBER:
            if _name in _wfak and abs(_wfak[_name] - 1.0) > 1e-9:
                print("    ✗ 1 KALIBER %s hat Gewichtsanteil %.6f statt 1 -- die "
                      "Schwellenherleitung in Schritt 2 traegt dann nicht"
                      % (_name, _wfak[_name]))
                fehler += 1
        if all(abs(_wfak.get(n, 1.0) - 1.0) <= 1e-9 for n, _ in _LOKAL_KALIBER):
            print("    ✓ 1 KALIBER: unsubstituierter Einringer hat GENAU EINEN "
                  "Torsionseintrag, Gewichtsanteil 1,0000 -- dort gibt es per "
                  "Konstruktion nichts zu verduennen.")

        # ===== 2 SYMMETRY FOLDING AND THRESHOLD ON THE CALIBRATION RING ===============
        # (2a) THE NUMBERS THEMSELVES: are global and ring-local THE SAME on the
        #      unsubstituted ring?  Not "similar" -- the derivation claims equality,
        #      so equality is measured, with a denominator.
        print()
        print("  ===== 2 SYMMETRIEFALTUNG UND SCHWELLE (unsubstituierte Ringe) =====")
        _paare_ges, _dmax_ges = 0, 0.0
        for _name, _smi in _LOKAL_KALIBER:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                _ids = list(AllChem.EmbedMultipleConfs(_m, numConfs=20, randomSeed=42))
                if len(_ids) < 2:
                    print("    %-14s nur %d Konformer -- nicht messbar"
                          % (_name, len(_ids))); continue
                AllChem.MMFFOptimizeMoleculeConfs(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                _li = _tfd_lokal_listen(_m, (_ro,))
                _d, _np_ = 0.0, 0
                for _i in range(len(_ids)):
                    for _j in range(_i + 1, len(_ids)):
                        _g = _tfd(_m, _ids[_i], _ids[_j])
                        _l = _tfd_lokal(_m, _li, _ids[_i], _ids[_j])
                        _d = max(_d, abs(_g - _l)); _np_ += 1
                _paare_ges += _np_; _dmax_ges = max(_dmax_ges, _d)
                print("    %-14s %4d Konformerpaare, groesste Differenz "
                      "|global - ringlokal| = %.3e" % (_name, _np_, _d))
            except Exception as _e:
                print("    %-14s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        if _paare_ges and _dmax_ges <= 1e-9:
            print("    ✓ 2a IDENTISCH auf %d Konformerpaaren (groesste Differenz %.1e). "
                  "⇒ DIE SCHWELLE BLEIBT 0,05: auf dem Kaliberring sind die beiden "
                  "Masse nicht aehnlich geeicht, sondern DASSELBE." % (_paare_ges, _dmax_ges))
        elif _paare_ges:
            print("    ✗ 2a NICHT identisch: groesste Differenz %.3e ueber %d Paare -- "
                  "die Schwelle muesste dann neu geeicht werden, und die Herleitung "
                  "in Schritt 1 ist falsch." % (_dmax_ges, _paare_ges))
            fehler += 1

        # (2b) THE STATE COUNT -- the test on which the CP replacement died.
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        print("    Zustaende je Ring, Raumgitter NAMP=2 NPHASE=8:")
        print("    %-14s %8s %8s %9s %9s" % ("Ring", "gl 0,05", "lo 0,05",
                                             "gl 0,005", "lo 0,005"))
        _n5 = None
        for _name, _smi in _LOKAL_KALIBER:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                _g05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _g005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
                _l05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _l005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                print("    %-14s %8d %8d %9d %9d" % (_name, _g05, _l05, _g005, _l005))
                if (_g05, _g005) != (_l05, _l005):
                    print("      ✗ 2b %s: ringlokal weicht auf dem KALIBERRING ab -- "
                          "0,05 bedeutet dort dann nicht mehr dasselbe" % _name)
                    fehler += 1
                if len(_ro) == 5:
                    _n5 = _l05
            except Exception as _e:
                print("    %-14s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        if _n5 is None:
            print("    ✗ 2b FUENFRING nicht gemessen -- die Symmetrieprobe fehlt")
            fehler += 1
        elif _n5 == 3:
            print("    ✓ 2b SYMMETRIEFALTUNG: unsubstituierter Fuenfring gibt ringlokal "
                  "3 Zustaende.  Der CP-Ersatz vom 26.08. gab hier 9/13/14 -- die "
                  "Faltung ueberlebt, weil RDKits Ringeintrag der MITTELWERT von "
                  "|Torsion| ueber den Ring ist und damit nummerierungsinvariant.")
        else:
            print("    ✗ 2b SYMMETRIEFALTUNG ZERSTOERT: Fuenfring gibt %d statt 3 "
                  "Zustaende -- derselbe Fehlschluss wie beim CP-Ersatz." % _n5)
            fehler += 1

        # ===== 3 THE DILUTION SERIES RE-RUN ==========================================
        print()
        print("  ===== 3 DIESELBE TABELLE, GLOBAL GEGEN RINGLOKAL =====")
        print("    Geruest EINGEFROREN, Raumgitter NAMP=2 NPHASE=8 -- exakt die "
              "Vorbereitung der Verduennungsreihe in `selbsttest_trennschaerfe`.")
        print("    %-24s %7s %10s %8s %9s %10s"
              % ("Molekuel", "Anteil", "w_R/sum(w)", "gl 0,05", "gl 0,005", "lo 0,05"))
        _reihe = []
        for _name, _smi in _VERD_PROBEN:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    print("    %-24s Einbettung fehlgeschlagen" % _name); continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                if not _ro:
                    print("    %-24s kein faltbarer Ring" % _name); continue
                _hv = [i for i in range(_m.GetNumAtoms())
                       if _m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                _g05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _g005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
                _l05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                print("    %-24s %6.1f%% %10.4f %8d %9d %10d"
                      % (_name, 100.0 * len(_rat) / len(_hv), _wfak.get(_name, 0.0),
                         _g05, _g005, _l05))
                _reihe.append({"name": _name, "g05": _g05, "g005": _g005, "l05": _l05,
                               "anteil": len(_rat) / float(len(_hv))})
            except Exception as _e:
                print("    %-24s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
    finally:
        for _k, _v in _alt.items():
            if _v is None:
                _os.environ.pop(_k, None)
            else:
                _os.environ[_k] = _v

    # ---- THE VERDICT.  ⚠ IT MAY ALSO GO AGAINST THE REPAIR -- a ring-local
    #      version that does NOT remove the gradient is a finding and not an error.
    if len(_reihe) >= 2:
        _g = [r["g05"] for r in _reihe]
        _l = [r["l05"] for r in _reihe]
        _gf = _g[-1] < _g[0]                       # global falls over the series
        _lf = _l[-1] < _l[0]                       # ring-local too?
        print()
        if not _gf:
            print("  ⇒ KEIN GRADIENT IN DER GLOBALEN SPALTE (%s) -- der Befund, den "
                  "dieser Test pruefen soll, tritt auf dieser Reihe gar nicht auf.  "
                  "Das Urteil ueber die Reparatur haengt in der Luft."
                  % " -> ".join(str(x) for x in _g))
            fehler += 1
        elif _lf:
            print("  ⇒ RINGLOKAL HILFT NICHT.  global %s, ringlokal %s -- der Abfall "
                  "bleibt.  Die Verduennung war dann nicht (oder nicht allein) die "
                  "Ursache; der naechste Verdaechtige ist der Relax, nicht das Mass."
                  % (" -> ".join(str(x) for x in _g), " -> ".join(str(x) for x in _l)))
        else:
            print("  ⇒ DER GRADIENT IST WEG.  global %s (Ringanteil %.0f -> %.0f %%), "
                  "ringlokal %s -- DERSELBE Ring, DIESELBEN Kandidaten, DIESELBE "
                  "Schwelle 0,05."
                  % (" -> ".join(str(x) for x in _g), 100.0 * _reihe[0]["anteil"],
                     100.0 * _reihe[-1]["anteil"], " -> ".join(str(x) for x in _l)))
            print("    Das kleinste Glied der Reihe gewinnt %d -> %d Zustaende je Ring. "
                  "⚠ ZUSTAENDE JE RING GEHEN POTENZIERT ins Kreuzprodukt ein "
                  "(gemessen 4,3 Ringe je System) -- das ist der Hebel, nicht die "
                  "Kombinatorik dahinter." % (_g[-1], _l[-1]))
            print("    ⚠ WAS DAMIT NICHT BEWIESEN IST: dass diese Zustaende das "
                  "Realismustor ueberleben.  Diese Reihe misst (a), die Zustaende JE "
                  "RING -- Kollision, Winkel und Bindungstor sitzen dahinter.")
            print("    Getragen wird dieser Befund von DELFIN_FFFREE_PUCKER_TFD_LOCAL "
                  "allein.  DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI steht in dieser Reihe "
                  "NICHT im Spiel und ist mit ihr auch nicht belegt -- seine einzige "
                  "Messung ist die Decalinzeile in Schritt 0.")
    print("=== Ringlokale TFD: %s ==="
          % ("gemessen" if fehler == 0 else "%d Pruefung(en) FEHLGESCHLAGEN" % fehler))
    return 1 if fehler else 0


if __name__ == "__main__":
    # ⚠ AT THE END OF THE FILE, and that is not cosmetic.  At MODULE LEVEL the
    #   order counts: if this block stands before one of the test functions, its name
    #   is still unbound at execution time -> NameError.  (Inside a
    #   function that does not apply -- exactly the confusion that produced the
    #   [Z4] death certificate on 09.08., only the other way round.)
    import sys as _sys
    # Without an argument everything runs -- a name runs exactly one test.  That is
    # not a convenience: `selbsttest_kombinatorik` BUILDS, and whoever wants to re-measure
    # it during a change should not have to pay for the TFD sweep three times for that.
    _TESTS = (("raum", selbsttest_raum),
              ("konvergenz", selbsttest_konvergenz),
              # The sweep stands AFTER the convergence, because it answers its open
              # question: it reports "states lie close together", the sweep measures whether
              # that is due to the threshold.
              ("sweep", selbsttest_tfd_sweep),
              # The discriminating power stands before the combinatorics: it decides WHICH
              # measure the combinatorics is supposed to dedup with at all.  A cost verdict with
              # the wrong dedup measure would be a verdict about the instrument.
              ("trennschaerfe", selbsttest_trennschaerfe),
              # Directly BEHIND IT, because the discriminating power throws off the finding that
              # this test repairs: it measures that the same ring on a growing scaffold
              # gets fewer and fewer states, this one measures the same series once more
              # with ring-local TFD.  Run separately they would be two measurements; this way
              # it is one measurement and its counter-check.
              ("tfdlokal", selbsttest_tfd_lokal),
              # Last the combinatorics: it builds relax + gate PER combination and is
              # thereby the most expensive of the four.  It answers what the three before
              # raise -- the state count per ring is only interesting because it
              # is raised to a power.
              ("kombinatorik", selbsttest_kombinatorik))
    _wahl = [a for a in _sys.argv[1:] if not a.startswith("-")]
    _unbekannt = [a for a in _wahl if a not in dict(_TESTS)]
    if _unbekannt:
        print("unbekannter Test: %s -- bekannt: %s"
              % (", ".join(_unbekannt), ", ".join(n for n, _ in _TESTS)))
        _sys.exit(2)
    # `--ohne-grenze` lifts the test limit of the combinatorics.  Then EVERY sample is
    # built completely -- including those whose cross product is five digits.  That is the
    # measurement, not the default: it runs for hours and does not belong in a regular run.
    _ohne_grenze = "--ohne-grenze" in _sys.argv[1:]
    _rc = 0
    for _name, _fn in _TESTS:
        if _wahl and _name not in _wahl:
            continue
        _rc = (_fn(max_kombis=0) if (_name == "kombinatorik" and _ohne_grenze)
               else _fn()) or _rc
    _sys.exit(_rc)
