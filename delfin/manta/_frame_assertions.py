"""delfin.manta._frame_assertions — WAS DIE KONSTRUKTION ZUGESICHERT HAT.

DAS PROBLEM, dreimal unabhaengig gemessen (17.08.2026)
------------------------------------------------------
Der Bauer setzt Atome MIT ABSICHT: dieser Ring ist konjugiert, also flach; dieser Donor
ist sp2, also trigonal; diese M=O-Bindung ist 1,63 A.  Was danach kommt -- UFF, MMFF,
Reembeds, Korrektoren -- weiss davon NICHTS und zerstoert es:

  * `pyramidal_sp2` 15,7 % gegen 1,0 % (31.07.) und `dhyb1` **463 planar_pyramidalised**
    (16.08.): der Bauer setzt planar, das Kraftfeld pyramidalisiert.  Es kennt keine
    pi-Konjugation.
  * **218 planare Metallacyclen** fehlen (17.08.), und `_is_puckerable` glaubt sogar
    ausdruecklich, konjugierte Ringe SEIEN schon flach -- eine Bauer-Ueberzeugung, die
    der Optimierer nie erhalten hat.
  * `me42d` (17.08.): eine NACHTRAEGLICHE Einzelbindungs-Korrektur zog ein terminales
    Re=O von 2,0 auf 1,63 A und kostete **2 CCDC-Isomere und 2 Polyeder** -- weil die
    SETZUNG die Koordinationssphaere bereits um die alte Laenge herum gebaut hatte.
  * `bbrefix` (17.08.): angehaengte Reembed-Frames kosten **+11,9 pp harte Frames**,
    weil der Erzeuger je Ligand faltet OHNE zu wissen, was die Setzung zugesichert hat.

Vier Befunde, EINE Wurzel: **die Zusicherung reist nicht mit dem Frame.**

WAS ES SCHON GIBT -- und warum das nicht reicht
------------------------------------------------
Traeger existieren, aber VERSTREUT und je Mechanismus neu erfunden:
  `torsion_relax._md_pairs` / `._md_ok`  -- M-D-Laengen, nur dort
  `_ring_pucker(frozen=...)`             -- eingefrorene Indizes, nur dort
  `converter_backend(donors=...)`        -- Donorindizes, nur dort
Kein Modul kennt die Zusicherungen der anderen, und keines kann PRUEFEN, ob eine
fremde Zusicherung verletzt wurde.  Genau das leistet dieses Modul.

WAS ES TUT -- UND WAS BEWUSST NICHT
------------------------------------
Es traegt eine **Zusicherung** (was die Konstruktion behauptet) und **misst Verletzungen**.
Es erzwingt NICHTS.  Das ist Absicht und folgt der Regel "erst messen, dann bauen":
bevor irgendein Relaxer gefesselt wird, muss die Zahl auf dem Tisch liegen, WIE OFT und
WIE STARK die eigenen Zusicherungen stromabwaerts gebrochen werden.  Diese Zahl gibt es
bis heute nicht.

Deterministisch, geometrieonly, kein Kraftfeld, kein `mol` -- die Zusicherung wird aus dem
FRAME abgeleitet, damit die Atomreihenfolge-Falle (an der der Ring-Pucker-Emitter und die
sp2-Planarisierer schon gescheitert sind) gar nicht erst entsteht.
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

# Toleranzen.  Bewusst GROSSZUEGIG: dies misst, ob eine Zusicherung GEBROCHEN wurde,
# nicht ob sie perfekt eingehalten ist.  Ein zu enger Wert erzeugt Rauschen und macht
# die Zahl unbrauchbar -- derselbe Fehler wie bei zu scharfen Augen-Schwellen.
_TRANS_MIN_DEG = 150.0    # ab hier gilt ein Donorpaar als gegenueberstehend
_TRANS_DROP_DEG = 15.0    # zulaessiger Abfall des groessten D-M-D-Winkels
_MD_TOL_A = 0.10          # Angstroem, M-D-Laenge
_PLANAR_TOL_A = 0.20      # Angstroem, Abstand des Zentrums von der Ebene seiner Nachbarn
_FROZEN_TOL_A = 0.05      # Angstroem, Bewegung eines eingefrorenen Atoms


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def derive(xyz: str) -> Optional[Dict]:
    """Die Zusicherung, die ein fertig konstruierter Frame TRAEGT.

    Abgeleitet aus dem Frame selbst -- kein `mol`, keine Atomindizes von aussen:
      * ``frozen``  : Metall + alle seine Donoren (die Koordinationssphaere, die JEDER
                      nachgelagerte Pass laut seinem eigenen Docstring stehen laesst)
      * ``md``      : (metall, donor, laenge) fuer jede M-D-Bindung
      * ``planar``  : jedes dreifach koordinierte C/N/O/B, das HIER flach steht -- der
                      Bauer hat es so gesetzt, also ist es eine Behauptung
    """
    if not xyz:
        return None
    try:
        # Auch ``(syms, P)`` direkt: der Bauer haelt seine Frames als Arrays und muesste
        # sonst fuer jede Pruefung eine XYZ-Zeichenkette bauen und wieder zerlegen --
        # zweimal Formatierung pro Frame, nur damit die Signatur stimmt.
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

    # ===== DIE ANORDNUNG DER DONOREN IST AUCH EINE BEHAUPTUNG =====================
    # Gemessen 2026-08-18 auf 30921 Systemen: netto 988 Systeme fliessen vom Oktaeder
    # ins trigonale Prisma (McNemar X2 = 860,8).  Der Bauer erzeugt 2,99 mal zu viele
    # Prismen, waehrend jede andere Form zwischen 0,86 und 1,29 bleibt -- und der
    # Zufluss ins Prisma ist exakt so gross wie der Abfluss aus dem Oktaeder.
    #
    # Es sind aber KEINE Prismen: die CShM-Masse gegen OC-6 liegt unimodal bei 8 bis 12
    # (Median 11,03), waehrend ein ideales TPR-6 bei 16,7 laege.  Es sind Oktaeder auf
    # halbem Weg -- ein unvollstaendiger Bailar-Twist.  Und es ist keine Auswahl: in
    # 1061 von 1061 Faellen ist poly_match false, obwohl das Auge den BESTEN Frame ueber
    # den ganzen Manifold liest.  Im ganzen Manifold gibt es kein Oktaeder.
    #
    # Das Signal ist die VERZAHNUNG, monoton: 2,49 % bei null Chelatringen, 16,12 % bei
    # fuenf.  Metall und d-Zahl sind flach.  Der Chelatzug dreht das Polyeder, nachdem
    # die Setzung es richtig gestellt hat.
    #
    # DIE INVARIANTE, die das festhaelt, ist billig und richtungsscharf: die Zahl der
    # TRANS-PAARE je Metall.  Ein Oktaeder hat drei, ein trigonales Prisma null.  Ein
    # halber Twist senkt sie.  Kein CShM, keine Referenzform, keine Schwelle mit
    # Feinabstimmung -- nur "was der Bauer an gegenueberliegenden Donoren gesetzt hat,
    # darf die Relaxation nicht verlieren".
    trans: Dict[int, Tuple[int, float]] = {}
    _tmin = _env_float("DELFIN_ASSERT_TRANS_MIN_DEG", _TRANS_MIN_DEG)
    for m in metals:
        _dn = [d for d in (nbrs[m] if m < len(nbrs) else [])
               if syms[d] != "H" and not _is_metal_sym(syms[d])]
        trans[m] = _trans_stats(P, m, _dn, _tmin)

    return {"n": len(syms), "frozen": sorted({*metals, *donors}),
            "md": sorted(md), "planar": planar, "trans": trans,
            "P": P.copy()}


def _trans_stats(P, m: int, donors, tmin_deg: float) -> Tuple[int, float]:
    """``(Zahl der Paare ueber der Schwelle, groesster D-M-D-Winkel)`` am Metall m.

    ⚠ WARUM BEIDES UND NICHT NUR DIE ZAHL.  Der erste Entwurf zaehlte nur Paare
    ueber 150 Grad -- und der eigene Selbsttest hat ihn widerlegt: ein HALBER
    Bailar-Twist (30 Grad) bewegt einen trans-Winkel von 180 auf 155,6 Grad, bleibt
    also ueber der Schwelle, waehrend ein ganzer Twist (60 Grad) auf 131,8 faellt.
    Gemessen ist aber genau der halbe -- die CShM-Masse der 1061 Faelle liegt bei
    8 bis 12 statt bei 16,7.  Eine Schwellenzaehlung haette den gemessenen Defekt
    komplett uebersehen und dabei "keine Verletzung" gemeldet.
    Der groesste Winkel faellt beim halben Twist um 24,4 Grad und ist damit das
    richtige Mass; die Zaehlung bleibt als zweites, grobes Signal daneben.
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
    """Abstand des Zentrums c von der Ebene durch a,b,d (Angstroem)."""
    try:
        nrm = np.cross(P[b] - P[a], P[d] - P[a])
        ln = float(np.linalg.norm(nrm))
        if ln < 1e-9:
            return None
        return abs(float(np.dot(P[c] - P[a], nrm / ln)))
    except Exception:
        return None


def violations(assertion: Optional[Dict], xyz_after: str) -> Optional[Dict]:
    """Welche Zusicherungen hat ein spaeterer Frame GEBROCHEN, und wie stark?

    Gibt Zaehler und die groessten Abweichungen.  ``None``, wenn nicht vergleichbar
    (andere Atomzahl -- dann ist es kein "spaeterer Frame desselben Baus").
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

    # TRANS-PAARE: nur der VERLUST zaehlt.  Ein Pass, der aus einem verzerrten Frame ein
    # regelmaessigeres macht, gewinnt Paare hinzu -- das ist kein Bruch, das ist der
    # Zweck der Relaxation.  Nur die Gegenrichtung ist der gemessene Defekt.
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
    """True, wenn der spaetere Frame KEINE Zusicherung bricht.

    Fuer einen umformenden Pass, der sich selbst pruefen will, BEVOR er sein Ergebnis
    zurueckgibt -- statt dass jeder Mechanismus seine eigene Teil-Invariante nachbaut.
    Nicht vergleichbar -> True (ein Pruefer, der nicht pruefen kann, darf nichts blockieren).
    """
    v = violations(assertion, xyz_after)
    return True if v is None else not v["any_broken"]


# ---------------------------------------------------------------------------
# Selbsttest:  python delfin/manta/_frame_assertions.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "t"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    fails = 0
    # Ein quadratisch-planarer Pt mit zwei Cl und zwei N, dazu ein planares sp2-C.
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

    # 2) Derselbe Frame bricht nichts.
    v = violations(a, base)
    ok = v is not None and not v["any_broken"]
    print(f"2 identischer Frame bricht nichts: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 3) EINE M-D-Bindung um 0,4 A verkuerzt -> md_broken UND frozen_moved.
    #    (Genau der me42d-Fall: ein Atom gezogen, Rest steht.)
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

    # 4) Ein planares sp2-C pyramidalisiert -> planar_broken.
    #    (Genau der dhyb1-Fall: 463 planar_pyramidalised.)
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

    # 5) Eine STARRE DREHUNG des ganzen Frames bricht NICHTS -- die Zusicherung darf
    #    nicht auf Koordinaten, sondern nur auf Geometrie reagieren.  ⚠ Ausnahme:
    #    `frozen` vergleicht Positionen, also muss die Drehung dort anschlagen; das ist
    #    korrekt und gewollt (ein eingefrorenes Atom SOLL sich nicht bewegen).
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

    # 6) holds() als Schutzschild fuer einen umformenden Pass.
    ok = holds(a, base) and not holds(a, pulled)
    print(f"6 holds() trennt sauber: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 7) Andere Atomzahl -> nicht vergleichbar -> None, und holds() blockiert NICHT.
    ok = violations(a, _xyz(rot[:5])) is None and holds(a, _xyz(rot[:5]))
    print(f"7 nicht vergleichbar -> None, holds() blockiert nicht: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 8) DER GEMESSENE DEFEKT: Oktaeder -> halber Bailar-Twist.  Ein reines OC-6 hat drei
    #    Trans-Paare; dreht man das untere Dreieck um 30 Grad Richtung Prisma, verliert es
    #    alle drei.  Genau das ist die Bewegung, die 988 Systeme netto ins Prisma traegt
    #    (McNemar X2 = 860,8 auf 30921 Systemen), und ihre CShM-Masse liegt bei 8 bis 12
    #    statt 16,7 -- also HALBER Twist, nicht ganzer.
    def _oct(twist_deg):
        rows = [("Fe", 0.0, 0.0, 0.0)]
        r, zh = 2.00, 1.1547           # r*cos(54,7 Grad); ergibt exakt 90/180 Grad
        rho = 1.63299                  # r*sin(54,7 Grad)
        for k, (ph, zs) in enumerate([(0.0, +1), (120.0, +1), (240.0, +1),
                                      (60.0, -1), (180.0, -1), (300.0, -1)]):
            a = math.radians(ph + (twist_deg if zs < 0 else 0.0))
            rows.append(("N", rho * math.cos(a), rho * math.sin(a), zs * zh))
        return _xyz(rows)
    a8 = derive(_oct(0.0))
    v8 = violations(a8, _oct(30.0))            # 30 Grad = halbe Strecke zum Prisma
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
