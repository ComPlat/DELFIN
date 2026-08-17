"""delfin.manta._me_bond_snap — terminale M=E Mehrfachbindungen auf dem FF-FREIEN Pfad.

WOZU (gemessen 17.08.2026).  Der Art-Parameter ``kind="me"`` wurde am 16.08. an drei Stellen
"verdrahtet".  Nachgezaehlt sind die drei nicht drei:

  * ``_clamp_metalloid_md_xyz``   -- METALLOID-only (Sb/As/Bi/Te/Se/Ge/Sn/Pb).  Sieht kein
    einziges Oxo-Metall; Re/W/V/Mo/Os stehen nicht in ``_METALLOID_MD_DONORS``.
  * ``_md_distance_in_tolerance`` -- per Signatur ``-> bool``.  Ein PRAEDIKAT, es bewegt
    kein Atom.
  * ``_manual_metal_embed``       -- der CN4-Pfad.  Terminale Oxo-Komplexe sind CN5-7.

Uebrig bleibt ``_snap_md_distances_to_ideal`` -- mit GENAU EINER Aufrufstelle
(``smiles_converter.py:33066``), und die liegt HINTER dem FF-freien Return (``:32399``,
``return _ff, None``); beide im selben Rumpf ``_smiles_to_xyz_isomers_impl:32156``.
=> Der einzige echte Setzer ist LEGACY-ONLY, und FF-frei ist 98,8 % der Bauten.
Gemessene Folge: ``me42b`` erreichte 2 von 42, ``byte_identical`` 40.

Vierte Instanz des 10.08.-Musters ("drei additive Module nur auf LEGACY -> AUFRUFSTELLEN
zaehlen, nie Zeilen").  Der Fehler war, FUNDSTELLEN des Art-Parameters zu zaehlen statt
SETZSTELLEN auf dem lebenden Pfad.

WARUM DIESES MODUL KEIN ``mol`` NIMMT.  ``_ffree_shared_tail`` warnt woertlich vor der
Atomreihenfolge: ein aus SMILES geparstes mol traegt RDKits Reihenfolge, FF-freie Frames
tragen Metall-auf-0 plus AddHs-Bloecke.  Genau daran ist der Ring-Pucker-Emitter an dieser
Stelle schon einmal zum Nulltest geworden (185 von 187 byte-identisch).  Das Kriterium fuer
ein terminales M=E braucht den Graphen aber gar nicht -- es ist STRUKTURELL und aus dem
Frame selbst ablesbar.  Damit entfaellt die Falle, statt umgangen zu werden.

DAS KRITERIUM ist identisch zu ``smiles_converter._ml_bond_kind``, nur aus der Geometrie
statt aus dem Graphen gelesen: Donor aus {O, N, C}, KEIN gebundener Wasserstoff, und sein
EINZIGER schwerer Nachbar ist ein Metall.  Ein OH/NH2 ist damit korrekt kein Oxo/Imido; ein
verbrueckendes mu-Oxo (zwei Metalle) faellt ueber "genau ein schwerer Nachbar" heraus --
ebenfalls wie im Original, wo eine einzelne Translation zwei M-D-Ideale ohnehin nicht
erfuellen koennte.

WARUM DIE KORREKTUR HIER BESONDERS SICHER IST.  Ein terminaler Donor hat per Definition
keinen weiteren schweren Nachbarn und keinen Wasserstoff.  Verschoben wird deshalb genau EIN
Atom -- kein BFS-Fragment, kein Ligandrumpf.  Es kann keine Bindung zerreissen, weil das Atom
ausser der M-D-Bindung keine hat.

Vorgabe AUS -> nicht gerufen -> byte-identisch.  Liegt fuer ein Paar kein kalibriertes Band
vor, gibt der Aufrufer ``None`` und dieses Modul ruehrt nichts an.  LIZENZ: die Werte sind
CCDC-abgeleitet, stehen NICHT in diesem Repo und werden hier auch nicht gelesen -- der
Aufrufer reicht sie herein.  Deterministisch (sortierte Reihenfolge, kein RNG), nie eine
nicht-endliche Koordinate.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple

import numpy as np

from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

# Elemente, fuer die es ueberhaupt terminale M=E-Chemie gibt (Oxo / Nitrido / Imido / Carbin).
_ME_ELEMENTS = frozenset({"O", "N", "C"})

# Harte Untergrenze fuer einen Schwer-Schwer-Kontakt.  Kommt der Donor durch die Verkuerzung
# einem DRITTEN Atom naeher als das, wird er zurueckgerollt -- die Verkuerzung zieht ihn ja
# in die Koordinationssphaere hinein.
_CLASH_FLOOR_A = 1.30

# Relative Mindestaenderung, damit ueberhaupt gesetzt wird.  Verhindert Rauschen auf Paaren,
# deren Band praktisch 1.0 ist (62 % der gemessenen Paare liegen zwischen 0.95 und 1.05).
_MIN_REL_DELTA = 0.02


def terminal_me_pairs(syms: List[str], nbrs: List[List[int]]) -> List[Tuple[int, int]]:
    """(metall_idx, donor_idx) fuer jeden strukturell terminalen M=E-Donor, sortiert."""
    out: List[Tuple[int, int]] = []
    for d, s in enumerate(syms):
        if s not in _ME_ELEMENTS:
            continue
        nb = nbrs[d] if d < len(nbrs) else []
        if any(syms[x] == "H" for x in nb):
            continue                                  # OH / NH2 / CH ist kein Oxo/Imido
        heavy = [x for x in nb if syms[x] != "H"]
        if len(heavy) != 1:
            continue                                  # terminal: genau EIN schwerer Nachbar
        m = heavy[0]
        if not _is_metal_sym(syms[m]):
            continue                                  # und der muss das Metall sein
        out.append((m, d))
    out.sort()
    return out


def snap_me_bonds(
    xyz_str: str,
    target_for: Callable[[str, str], Optional[float]],
) -> str:
    """Setze jede terminale M=E-Bindung auf ihre kalibrierte Laenge.

    ``target_for(metall_symbol, donor_symbol)`` gibt die Ziellaenge in Angstroem, oder
    ``None``, wenn fuer das Paar KEIN kalibriertes Band vorliegt.  ``None`` heisst
    ausdruecklich "nichts tun" -- ohne Tabelle ist der Durchlauf byte-identisch.

    Gibt bei jedem Fehlschlag den unveraenderten Eingang zurueck.
    """
    if not xyz_str:
        return xyz_str
    try:
        syms, pts, lines = _parse_xyz(xyz_str)
    except Exception:
        return xyz_str
    if not syms or pts is None or len(syms) < 2:
        return xyz_str
    try:
        nbrs, _blen = _build_geometric_adjacency(syms, pts)
    except Exception:
        return xyz_str

    pairs = terminal_me_pairs(syms, nbrs)
    if not pairs:
        return xyz_str

    new_pts = np.array(pts, dtype=float, copy=True)
    moved = False
    for m, d in pairs:
        try:
            target = target_for(syms[m], syms[d])
        except Exception:
            target = None
        if target is None:
            continue
        try:
            target = float(target)
        except Exception:
            continue
        if not np.isfinite(target) or target <= 0.0:
            continue
        v = new_pts[d] - new_pts[m]
        cur = float(np.linalg.norm(v))
        if cur < 1e-8:
            continue
        if abs(cur - target) / target < _MIN_REL_DELTA:
            continue
        cand = new_pts[m] + v * (target / cur)
        if not np.all(np.isfinite(cand)):
            continue
        # ROLLBACK: der Donor darf keinem DRITTEN Atom naeher kommen als der harten
        # Untergrenze.  Das Metall ist ausgenommen -- zu ihm ist die neue Distanz das Ziel.
        others = [k for k in range(len(syms)) if k != d and k != m]
        if others:
            before = float(np.min(np.linalg.norm(new_pts[others] - new_pts[d], axis=1)))
            after = float(np.min(np.linalg.norm(new_pts[others] - cand, axis=1)))
            if after < before and after < _CLASH_FLOOR_A:
                continue                              # verworfen, dieser Donor bleibt
        new_pts[d] = cand
        moved = True

    if not moved:
        return xyz_str
    try:
        return _format_xyz(lines, syms, new_pts)
    except Exception:
        return xyz_str


# ---------------------------------------------------------------------------
# Selbsttest:  python delfin/manta/_me_bond_snap.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    def _dist(xyz, i, j):
        s, p, _l = _parse_xyz(xyz)
        return float(np.linalg.norm(p[i] - p[j]))

    fails = 0

    # 1) TERMINALES OXO wird auf die Ziellaenge gesetzt.
    base = _xyz([("W", 0.0, 0.0, 0.0),
                 ("O", 2.10, 0.0, 0.0),      # terminal -> soll wandern
                 ("Cl", 0.0, 2.30, 0.0),
                 ("Cl", 0.0, -2.30, 0.0),
                 ("Cl", 0.0, 0.0, 2.30)])
    got = snap_me_bonds(base, lambda m, d: 1.905 if (m, d) == ("W", "O") else None)
    d1 = _dist(got, 0, 1)
    ok = abs(d1 - 1.905) < 1e-4
    print(f"1 terminales W=O 2.100 -> {d1:.4f} (Ziel 1.905)  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 2) OHNE BAND (target None) passiert NICHTS -- byte-identisch.
    same = snap_me_bonds(base, lambda m, d: None)
    ok = (same == base)
    print(f"2 ohne Band byte-identisch: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 3) HYDROXO wird NICHT angefasst (O traegt ein H).
    oh = _xyz([("W", 0.0, 0.0, 0.0),
               ("O", 2.10, 0.0, 0.0),
               ("H", 2.70, 0.90, 0.0),
               ("Cl", 0.0, 2.30, 0.0),
               ("Cl", 0.0, -2.30, 0.0)])
    ok = (snap_me_bonds(oh, lambda m, d: 1.905) == oh)
    print(f"3 Hydroxo unangetastet: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 4) VERBRUECKENDES mu-Oxo (zwei Metalle) wird NICHT angefasst.
    mu = _xyz([("W", 0.0, 0.0, 0.0),
               ("O", 1.95, 0.0, 0.0),
               ("W", 3.90, 0.0, 0.0),
               ("Cl", 0.0, 2.30, 0.0),
               ("Cl", 3.90, 2.30, 0.0)])
    ok = (snap_me_bonds(mu, lambda m, d: 1.60) == mu)
    print(f"4 mu-Oxo unangetastet: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 5) Der ZWEITE Aufruf aendert nichts mehr (idempotent).
    twice = snap_me_bonds(got, lambda m, d: 1.905 if (m, d) == ("W", "O") else None)
    ok = (twice == got)
    print(f"5 idempotent: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # 6) Nur der Donor bewegt sich -- alle uebrigen Atome stehen exakt still.
    s0, p0, _ = _parse_xyz(base)
    s1, p1, _ = _parse_xyz(got)
    moved = [i for i in range(len(s0)) if float(np.linalg.norm(p0[i] - p1[i])) > 1e-9]
    ok = (moved == [1])
    print(f"6 bewegte Atome = {moved} (erwartet [1])  {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    print(f"\n{6 - fails}/6 bestanden")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
