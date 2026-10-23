"""delfin.manta._frame_assertions — WHAT THE CONSTRUCTION HAS ASSERTED.

THE PROBLEM, measured three times independently (17.08.2026)
------------------------------------------------------------
The builder seats atoms WITH INTENT: this ring is conjugated, hence flat; this donor
is sp2, hence trigonal; this M=O bond is 1.63 A.  Whatever comes afterwards -- UFF, MMFF,
reembeds, correctors -- knows NOTHING of that and destroys it:

  * `pyramidal_sp2` 15.7 % versus 1.0 % (31.07.) and `dhyb1` **463 planar_pyramidalised**
    (16.08.): the builder seats planar, the force field pyramidalises.  It knows no
    pi conjugation.
  * **218 planar metallacycles** are missing (17.08.), and `_is_puckerable` even believes
    explicitly that conjugated rings ARE already flat -- a builder conviction that the
    optimiser never received.
  * `me42d` (17.08.): a SUBSEQUENT single-bond correction pulled a terminal
    Re=O from 2.0 to 1.63 A and cost **2 CCDC isomers and 2 polyhedra** -- because the
    SEATING had already built the coordination sphere around the old length.
  * `bbrefix` (17.08.): appended reembed frames cost **+11.9 pp hard frames**,
    because the generator folds per ligand WITHOUT knowing what the seating asserted.

Four findings, ONE root: **the assertion does not travel with the frame.**

WHAT ALREADY EXISTS -- and why that is not enough
--------------------------------------------------
Carriers exist, but SCATTERED and reinvented per mechanism:
  `torsion_relax._md_pairs` / `._md_ok`  -- M-D lengths, only there
  `_ring_pucker(frozen=...)`             -- frozen indices, only there
  `converter_backend(donors=...)`        -- donor indices, only there
No module knows the assertions of the others, and none can CHECK whether a
foreign assertion was violated.  That is exactly what this module provides.

WHAT IT DOES -- AND WHAT IT DELIBERATELY DOES NOT
--------------------------------------------------
It carries an **assertion** (what the construction claims) and **measures violations**.
It enforces NOTHING.  That is intentional and follows the rule "measure first, then build":
before any relaxer is shackled, the number must be on the table of HOW OFTEN and
HOW STRONGLY the construction's own assertions are broken downstream.  That number does
not exist to this day.

Deterministic, geometry-only, no force field, no `mol` -- the assertion is derived from the
FRAME, so that the atom-order trap (on which the ring-pucker emitter and the
sp2 planarisers have already failed) never arises in the first place.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _is_metal_sym,
    _parse_xyz,
)

# Tolerances.  Deliberately GENEROUS: this measures whether an assertion was BROKEN,
# not whether it is kept perfectly.  Too tight a value produces noise and makes
# the number useless -- the same mistake as with eye thresholds that are too sharp.
_TRANS_MIN_DEG = 150.0    # from here on a donor pair counts as opposite each other
_TRANS_DROP_DEG = 15.0    # permissible drop of the largest D-M-D angle
_MD_TOL_A = 0.10          # Angstrom, M-D length
_PLANAR_TOL_A = 0.20      # Angstrom, distance of the centre from the plane of its neighbours
_FROZEN_TOL_A = 0.05      # Angstrom, movement of a frozen atom


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def derive(xyz, frozen=None) -> Optional[Dict]:
    """The assertion that a fully constructed frame CARRIES.

    Derived from the frame itself -- no `mol`, no atom indices from outside:
      * ``frozen``  : metal + all of its donors (the coordination sphere that EVERY
                      downstream pass leaves standing, according to its own docstring)
      * ``md``      : (metal, donor, length) for every M-D bond
      * ``planar``  : every three-coordinate C/N/O/B that stands flat HERE -- the
                      builder seated it that way, so it is a claim
    """
    if not xyz:
        return None
    try:
        # ``(syms, P)`` directly as well: the builder holds its frames as arrays and would
        # otherwise have to build an XYZ string and parse it again for every check --
        # formatting twice per frame, just so that the signature fits.
        if isinstance(xyz, tuple):
            syms, P = xyz[0], np.asarray(xyz[1], float)
        else:
            syms, P, _lines = _parse_xyz(xyz)
        nbrs, _bl = _build_geometric_adjacency(syms, P)
    except Exception:
        return None
    if not syms or P is None:
        return None
    P = np.asarray(P, float)

    metals = [i for i, s in enumerate(syms) if _is_metal_sym(s)]
    md: List[Tuple[int, int, float]] = []
    donors: set = set()
    for m in metals:
        for d in nbrs[m] if m < len(nbrs) else []:
            if syms[d] == "H" or _is_metal_sym(syms[d]):
                continue
            donors.add(d)
            md.append((m, d, float(np.linalg.norm(P[d] - P[m]))))

    planar: List[Tuple[int, Tuple[int, int, int]]] = []
    for i, s in enumerate(syms):
        if s not in ("C", "N", "O", "B"):
            continue
        heavy = [k for k in (nbrs[i] if i < len(nbrs) else []) if syms[k] != "H"]
        if len(heavy) != 3:
            continue
        oop = _oop(P, i, heavy[0], heavy[1], heavy[2])
        if oop is not None and oop <= _env_float("DELFIN_ASSERT_PLANAR_TOL", _PLANAR_TOL_A):
            planar.append((i, (heavy[0], heavy[1], heavy[2])))

    # ===== THE ARRANGEMENT OF THE DONORS IS ALSO A CLAIM ==========================
    # Measured 2026-08-18 on 30921 systems: a net 988 systems flow from the octahedron
    # into the trigonal prism (McNemar X2 = 860.8).  The builder produces 2.99 times too many
    # prisms, while every other shape stays between 0.86 and 1.29 -- and the
    # inflow into the prism is exactly as large as the outflow from the octahedron.
    #
    # But they are NOT prisms: the CShM mass against OC-6 is unimodal at 8 to 12
    # (median 11.03), whereas an ideal TPR-6 would lie at 16.7.  They are octahedra
    # halfway along -- an incomplete Bailar twist.  And it is not a selection: in
    # 1061 of 1061 cases poly_match is false, although the eye reads the BEST frame over
    # the whole manifold.  In the whole manifold there is no octahedron.
    #
    # The signal is the INTERLOCKING, monotonic: 2.49 % at zero chelate rings,
    # 16.12 % at five.  Metal and d count are flat.  The chelate pull turns the polyhedron
    # after the seating has put it right.
    #
    # THE INVARIANT that pins this down is cheap and sharp in direction: the number of
    # TRANS PAIRS per metal.  An octahedron has three, a trigonal prism zero.  A
    # half twist lowers it.  No CShM, no reference shape, no threshold with
    # fine-tuning -- only "what the builder seated as opposite donors,
    # the relaxation must not lose".
    trans: Dict[int, Tuple[int, float]] = {}
    _tmin = _env_float("DELFIN_ASSERT_TRANS_MIN_DEG", _TRANS_MIN_DEG)
    for m in metals:
        _dn = [d for d in (nbrs[m] if m < len(nbrs) else [])
               if syms[d] != "H" and not _is_metal_sym(syms[d])]
        trans[m] = _trans_stats(P, m, _dn, _tmin)

    # ⚠️ WHICH SET IS FROZEN -- THE DERIVED ONE OR THE PROMISED ONE?
    # Measured 18.08. on the octahedron cases: the first hit of the assertion was
    # `frozen=2` at a movement of 0.194 Angstrom.  That looks like a defect --
    # but can just as well be a DEFINITION GAP: derive() derives metal plus
    # geometric neighbours, while refine() receives the set `fixed` from the
    # builder.  If those diverge, the assertion reports a break for an
    # atom that was never promised -- and the rollback switch would fire on a
    # measurement error.
    # That is why derive() now accepts the builder's set.  Then the
    # invariant is exactly the builder's own contract ("these atoms you hold fixed") and not
    # my reconstruction of it.  Without the argument, everything stays as before.
    _frozen = (sorted({int(i) for i in frozen if 0 <= int(i) < len(syms)})
               if frozen is not None else sorted({*metals, *donors}))

    return {"n": len(syms), "frozen": _frozen,
            "frozen_source": "builder" if frozen is not None else "derived",
            "md": sorted(md), "planar": planar, "trans": trans,
            "P": P.copy()}


def _trans_stats(P, m: int, donors, tmin_deg: float) -> Tuple[int, float]:
    """``(number of pairs above the threshold, largest D-M-D angle)`` at metal m.

    ⚠ WHY BOTH AND NOT ONLY THE COUNT.  The first draft counted only pairs
    above 150 degrees -- and its own self-test refuted it: a HALF
    Bailar twist (30 degrees) moves a trans angle from 180 to 155.6 degrees, so it
    stays above the threshold, whereas a full twist (60 degrees) drops to 131.8.
    But what is measured is precisely the half one -- the CShM mass of the 1061 cases
    lies at 8 to 12 instead of at 16.7.  A threshold count would have missed the measured
    defect completely and reported "no violation" while doing so.
    The largest angle drops by 24.4 degrees under the half twist and is therefore the
    right measure; the count stays alongside as a second, coarse signal.
    """
    n, mx = 0, 0.0
    for a in range(len(donors)):
        for b in range(a + 1, len(donors)):
            va = P[donors[a]] - P[m]
            vb = P[donors[b]] - P[m]
            na, nb = float(np.linalg.norm(va)), float(np.linalg.norm(vb))
            if na < 1e-9 or nb < 1e-9:
                continue
            c = float(np.dot(va, vb)) / (na * nb)
            c = max(-1.0, min(1.0, c))
            ang = math.degrees(math.acos(c))
            if ang >= tmin_deg:
                n += 1
            mx = max(mx, ang)
    return n, mx


def _oop(P, c: int, a: int, b: int, d: int) -> Optional[float]:
    """Distance of the centre c from the plane through a,b,d (Angstrom)."""
    try:
        nrm = np.cross(P[b] - P[a], P[d] - P[a])
        ln = float(np.linalg.norm(nrm))
        if ln < 1e-9:
            return None
        return abs(float(np.dot(P[c] - P[a], nrm / ln)))
    except Exception:
        return None


def violations(assertion: Optional[Dict], xyz_after: str) -> Optional[Dict]:
    """Which assertions has a later frame BROKEN, and how strongly?

    Returns counters and the largest deviations.  ``None`` if not comparable
    (different atom count -- then it is not a "later frame of the same build").
    """
    if not assertion or xyz_after is None or (
            not isinstance(xyz_after, tuple) and not xyz_after):
        return None
    try:
        if isinstance(xyz_after, tuple):
            syms2, P2 = xyz_after[0], np.asarray(xyz_after[1], float)
        else:
            syms2, P2, _l = _parse_xyz(xyz_after)
    except Exception:
        return None
    if not syms2 or P2 is None or len(syms2) != assertion.get("n"):
        return None
    P2 = np.asarray(P2, float)
    P1 = assertion["P"]

    md_tol = _env_float("DELFIN_ASSERT_MD_TOL", _MD_TOL_A)
    pl_tol = _env_float("DELFIN_ASSERT_PLANAR_TOL", _PLANAR_TOL_A)
    fr_tol = _env_float("DELFIN_ASSERT_FROZEN_TOL", _FROZEN_TOL_A)

    n_md, worst_md = 0, 0.0
    for (m, d, l0) in assertion["md"]:
        l1 = float(np.linalg.norm(P2[d] - P2[m]))
        dev = abs(l1 - l0)
        if dev > md_tol:
            n_md += 1
            worst_md = max(worst_md, dev)

    n_pl, worst_pl = 0, 0.0
    for (c, (a, b, d)) in assertion["planar"]:
        oop = _oop(P2, c, a, b, d)
        if oop is not None and oop > pl_tol:
            n_pl += 1
            worst_pl = max(worst_pl, oop)

    n_fr, worst_fr = 0, 0.0
    for i in assertion["frozen"]:
        mv = float(np.linalg.norm(P2[i] - P1[i]))
        if mv > fr_tol:
            n_fr += 1
            worst_fr = max(worst_fr, mv)

    # TRANS PAIRS: only the LOSS counts.  A pass that turns a distorted frame into a
    # more regular one gains pairs -- that is not a break, that is the
    # purpose of the relaxation.  Only the opposite direction is the measured defect.
    _tmin = _env_float("DELFIN_ASSERT_TRANS_MIN_DEG", _TRANS_MIN_DEG)
    _tdrop = _env_float("DELFIN_ASSERT_TRANS_DROP_DEG", _TRANS_DROP_DEG)
    n_tr, worst_tr = 0, 0.0
    try:
        _nb2, _ = _build_geometric_adjacency(syms2, P2)
    except Exception:
        _nb2 = None
    for m, (c0, mx0) in (assertion.get("trans") or {}).items():
        if _nb2 is None or m >= len(_nb2):
            continue
        _dn = [d for d in _nb2[m]
               if syms2[d] != "H" and not _is_metal_sym(syms2[d])]
        c1, mx1 = _trans_stats(P2, m, _dn, _tmin)
        _drop = mx0 - mx1
        if c1 < c0 or _drop > _tdrop:
            n_tr += 1
            worst_tr = max(worst_tr, _drop)

    return {"md_broken": n_md, "md_total": len(assertion["md"]), "md_worst": round(worst_md, 3),
            "planar_broken": n_pl, "planar_total": len(assertion["planar"]),
            "planar_worst": round(worst_pl, 3),
            "frozen_moved": n_fr, "frozen_total": len(assertion["frozen"]),
            "frozen_worst": round(worst_fr, 3),
            "trans_lost_metals": n_tr, "trans_total_metals": len(assertion.get("trans") or {}),
            "trans_worst": worst_tr,
            "any_broken": bool(n_md or n_pl or n_fr or n_tr)}


def holds(assertion: Optional[Dict], xyz_after: str) -> bool:
    """True if the later frame breaks NO assertion.

    For a reshaping pass that wants to check itself BEFORE it returns its
    result -- instead of every mechanism rebuilding its own partial invariant.
    Not comparable -> True (a checker that cannot check must not block anything).
    """
    v = violations(assertion, xyz_after)
    return True if v is None else not v["any_broken"]


# ---------------------------------------------------------------------------
# Self-test:  python delfin/manta/_frame_assertions.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "t"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    fails = 0
    # A square-planar Pt with two Cl and two N, plus a planar sp2 C.
    base = _xyz([("Pt", 0.0, 0.0, 0.0),
                 ("Cl", 2.30, 0.0, 0.0), ("Cl", -2.30, 0.0, 0.0),
                 ("N", 0.0, 2.05, 0.0), ("N", 0.0, -2.05, 0.0),
                 ("C", 0.0, 3.45, 0.0),
                 ("C", 1.20, 4.15, 0.0), ("C", -1.20, 4.15, 0.0)])
    a = derive(base)
    ok = a is not None and len(a["md"]) == 4 and 0 in a["frozen"]
    print(f"1 Zusicherung abgeleitet: {len(a['md']) if a else '-'} M-D, "
          f"{len(a['frozen']) if a else '-'} eingefroren, {len(a['planar']) if a else '-'} planar"
          f"  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 2) The same frame breaks nothing.
    v = violations(a, base)
    ok = v is not None and not v["any_broken"]
    print(f"2 identischer Frame bricht nichts: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 3) ONE M-D bond shortened by 0.4 A -> md_broken AND frozen_moved.
    #    (Exactly the me42d case: one atom pulled, the rest stays put.)
    pulled = _xyz([("Pt", 0.0, 0.0, 0.0),
                   ("Cl", 1.90, 0.0, 0.0), ("Cl", -2.30, 0.0, 0.0),
                   ("N", 0.0, 2.05, 0.0), ("N", 0.0, -2.05, 0.0),
                   ("C", 0.0, 3.45, 0.0),
                   ("C", 1.20, 4.15, 0.0), ("C", -1.20, 4.15, 0.0)])
    v = violations(a, pulled)
    ok = v["md_broken"] == 1 and abs(v["md_worst"] - 0.4) < 1e-6 and v["frozen_moved"] == 1
    print(f"3 me42d-Fall (0,4 A gezogen): md_broken={v['md_broken']} worst={v['md_worst']} "
          f"frozen_moved={v['frozen_moved']}  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 4) A planar sp2 C pyramidalises -> planar_broken.
    #    (Exactly the dhyb1 case: 463 planar_pyramidalised.)
    pyr = _xyz([("Pt", 0.0, 0.0, 0.0),
                ("Cl", 2.30, 0.0, 0.0), ("Cl", -2.30, 0.0, 0.0),
                ("N", 0.0, 2.05, 0.0), ("N", 0.0, -2.05, 0.0),
                ("C", 0.0, 3.45, 0.55),
                ("C", 1.20, 4.15, 0.0), ("C", -1.20, 4.15, 0.0)])
    v = violations(a, pyr)
    ok = v["planar_broken"] >= 1
    print(f"4 dhyb1-Fall (sp2 pyramidalisiert): planar_broken={v['planar_broken']} "
          f"worst={v['planar_worst']}  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 5) A RIGID ROTATION of the whole frame breaks NOTHING -- the assertion must
    #    react not to coordinates but only to geometry.  ⚠ Exception:
    #    `frozen` compares positions, so the rotation has to trigger there; that is
    #    correct and intended (a frozen atom is SUPPOSED not to move).
    import math
    th = math.radians(37.0)
    rot = []
    for line in base.splitlines()[2:]:
        p = line.split()
        x, y, z = float(p[1]), float(p[2]), float(p[3])
        rot.append((p[0], x * math.cos(th) - y * math.sin(th),
                    x * math.sin(th) + y * math.cos(th), z))
    v = violations(a, _xyz(rot))
    ok = v["md_broken"] == 0 and v["planar_broken"] == 0
    print(f"5 starre Drehung: md_broken={v['md_broken']} planar_broken={v['planar_broken']} "
          f"(frozen_moved={v['frozen_moved']}, erwartet >0)  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 6) holds() as a shield for a reshaping pass.
    ok = holds(a, base) and not holds(a, pulled)
    print(f"6 holds() trennt sauber: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 7) Different atom count -> not comparable -> None, and holds() does NOT block.
    ok = violations(a, _xyz(rot[:5])) is None and holds(a, _xyz(rot[:5]))
    print(f"7 nicht vergleichbar -> None, holds() blockiert nicht: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 8) THE MEASURED DEFECT: octahedron -> half Bailar twist.  A pure OC-6 has three
    #    trans pairs; turning the lower triangle by 30 degrees towards the prism loses
    #    all three.  That is exactly the movement that carries a net 988 systems into the
    #    prism (McNemar X2 = 860.8 on 30921 systems), and its CShM mass lies at 8 to 12
    #    instead of 16.7 -- so a HALF twist, not a full one.
    def _oct(twist_deg):
        rows = [("Fe", 0.0, 0.0, 0.0)]
        r, zh = 2.00, 1.1547           # r*cos(54.7 degrees); yields exactly 90/180 degrees
        rho = 1.63299                  # r*sin(54.7 degrees)
        for k, (ph, zs) in enumerate([(0.0, +1), (120.0, +1), (240.0, +1),
                                      (60.0, -1), (180.0, -1), (300.0, -1)]):
            a = math.radians(ph + (twist_deg if zs < 0 else 0.0))
            rows.append(("N", rho * math.cos(a), rho * math.sin(a), zs * zh))
        return _xyz(rows)
    a8 = derive(_oct(0.0))
    v8 = violations(a8, _oct(30.0))            # 30 degrees = halfway to the prism
    v8b = violations(a8, _oct(0.0))
    ok = (a8 is not None and a8["trans"].get(0, (0, 0.0))[0] == 3
          and abs(a8["trans"][0][1] - 180.0) < 0.5
          and v8 is not None and v8["trans_lost_metals"] == 1 and v8["trans_worst"] > 20.0
          and v8b["trans_lost_metals"] == 0 and not v8b["any_broken"])
    print(f"8 halber Bailar-Twist: max D-M-D vorher="
          f"{a8['trans'][0][1]:.1f} Grad, Abfall={v8['trans_worst']:.1f} Grad, "
          f"gebrochen={v8['trans_lost_metals']} (unveraendert: "
          f"{v8b['trans_lost_metals']})  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    print(f"\n{8 - fails}/8 bestanden")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
