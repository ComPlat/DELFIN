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

    return {"n": len(syms), "frozen": sorted({*metals, *donors}),
            "md": sorted(md), "planar": planar,
            "P": P.copy()}


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
    if not assertion or not xyz_after:
        return None
    try:
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

    n_md = worst_md = 0.0, 0.0
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

    return {"md_broken": n_md, "md_total": len(assertion["md"]), "md_worst": round(worst_md, 3),
            "planar_broken": n_pl, "planar_total": len(assertion["planar"]),
            "planar_worst": round(worst_pl, 3),
            "frozen_moved": n_fr, "frozen_total": len(assertion["frozen"]),
            "frozen_worst": round(worst_fr, 3),
            "any_broken": bool(n_md or n_pl or n_fr)}


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

    print(f"\n{7 - fails}/7 bestanden")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
