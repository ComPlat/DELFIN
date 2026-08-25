"""_refine_gate.py — DAS RUECKNAHMETOR.  Ein FILTER, kein Korrektor.

WOZU (24.08.2026, aus der Frage des Users).  Ein Verfeinerer, der Atome bewegt, gewinnt auf
vielen Systemen und verliert auf wenigen.  Das Landetor kennt 33 Verbotsterme, alle auf der
VERLUSTSEITE, und `capability_gained` wird berechnet, gedruckt -- und fliesst in keine
einzige Blockierentscheidung ein.  Ein netto positiver Verfeinerer kann damit nie landen.

DIE ANTWORT IST NICHT, DAS TOR ZU LOCKERN.  Ein Mittelwert kann keine tote Klasse sehen, und
das Auge selbst aendert sich -- never-worse ist das einzige Kriterium, das einen Augenwechsel
ueberlebt.  Die Antwort ist, dem Verfeinerer den Verlust WEGZUNEHMEN:

    nach der Kette ist KEIN Frame schlechter als davor.

Dann ist "Atome so bewegen, dass sie besser stehen" genau das, was passiert -- ohne
Kollateralschaden, und never-worse gilt PER KONSTRUKTION statt im Schnitt.

⛔ DREI REGELN, DIE DIE BAUART BESTIMMEN (User, 24.08.):

 1. **KEIN KORREKTOR EINES KORREKTORS.**  Dieses Modul schreibt NIE Koordinaten.  Es waehlt
    zwischen zwei bereits vorhandenen Frames -- Ergebnis oder Original.  Ein Korrektor, der
    einen Korrektor gutmacht, ist ein Zeichen, dass einer von beiden falsch ist; dann gehoert
    der eine richtig gemacht, nicht ein dritter daneben gestellt.

 2. **KEINE AUFBLAEHUNG.**  EINE oeffentliche Funktion, EIN Modul, und ZWEI Zeilen je
    Verfeinerungskette -- nicht je Verfeinerer.  Es gibt ueber zwanzig `_apply_*`-Stellen;
    jede einzeln zu umwickeln waere genau die Aufblaehung, die hier ausgeschlossen ist.
    Und es werden KEINE Radien, Schwellen oder Graphen neu definiert: alles kommt aus
    `_h_placement`, wo es kalibriert steht und gegen die Detektoren begruendet ist.

 3. **GESCHWINDIGKEIT.**  Der Punkt, an dem das billig wird: ein Verfeinerer laesst die
    meisten Frames BYTE-IDENTISCH.  Ein Stringvergleich davor kostet nichts und schliesst
    sie aus.  Bewertet wird nur, was sich wirklich geaendert hat -- typisch eine Handvoll
    Frames je System statt aller.  Ohne diesen Vorfilter waere das Tor teurer als die
    Verfeinerer, die es bewacht.

DAS KRITERIUM -- referenzfrei, lokal, ZWEI GANZE ZAHLEN, beide "kleiner ist besser":

    n_clash      Atompaare, die ihren Boden unterschreiten (H...H 1,50 A; H...schwer
                 0,85 x vdW-Summe; schwer...schwer 0,70 x vdW-Summe), nur fuer wirklich
                 nichtgebundene Paare ab Graphabstand `_NONBONDED_MIN_HOPS`
    n_bond_out   kovalente Bindungen ausserhalb [0,85 .. 1,15] x Kovalenzradiensumme

Ein Frame gilt als NICHT SCHLECHTER, wenn KEINE der beiden Zahlen steigt.  Zahl gegen Zahl,
also dimensionsrein -- kein Vergleich einer Schwere mit einer Framezahl, wie er dieses Repo
schon dreimal in die Irre gefuehrt hat.

⚠ DAS KRITERIUM IST NICHT DAS AUGE, und das ist Absicht.  Wer im Bauer das Auge nachbaut,
uebergibt Goodhart die Schluessel: der Bauer optimiert dann die Messung statt die Geometrie.
Hier stehen nur physikalische Invarianten, die der Bauer aus sich selbst pruefen kann.

⚠ NEUE ETIKETTEN GEHEN UNBERUEHRT DURCH.  Enumeratoren (Spiegel, Faltung, Atropisomer)
HAENGEN AN; sie koennen per Konstruktion nichts verlieren, und ein angehaengtes Geschwister
DARF schlechter sein als das Original -- das ist Vollstaendigkeit, keine Regression.  Die
Zuordnung laeuft ueber das ETIKETT, damit dieses Modul die beiden Mechanismusklassen ohne
Kenntnis der einzelnen Paesse auseinanderhaelt.

⚠ NAMENSGEBUNG `_rg_`: der Wachhund `check_exists_first` hat `_enabled` und `_score` als
Kollision gemeldet (`joint_declash`, `sphere_flex`, `smiles_converter`).  Dieselbe Lehre wie
bei `_hp_` und `_atrop_` -- ein eindeutiger Name macht die Kopie sichtbar, statt sie zu tarnen.

Vorgabe AUS (``DELFIN_FFFREE_REFINE_GATE``) -> byte-identisch.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

# KEINE eigenen Radien, Boeden oder Graphen: alles aus `_h_placement`, wo es kalibriert steht.
from delfin.manta._h_placement import (
    _HH_FLOOR,
    _H_HEAVY_FRAC,
    _NONBONDED_MIN_HOPS,
    _XH_FIRE_HI,
    _XH_FIRE_LO,
    _hp_cov,
    _hp_graph,
    _hp_metal,
    _hp_read,
    _hp_vdw,
    _hp_within,
)

_LOG = logging.getLogger(__name__)

FLAG = "DELFIN_FFFREE_REFINE_GATE"

# Schwer-gegen-schwer ist der EINZIGE Boden, den `_h_placement` nicht braucht (es bewegt nur
# H).  0,70 x vdW-Summe ist bewusst LOCKER: eine echte Bindung liegt weit darunter und wird
# ohnehin ueber den Graphabstand ausgeschlossen; gemeint sind nur grobe Durchdringungen.
_RG_HEAVY_HEAVY_FRAC = 0.70


def _rg_on() -> bool:
    """DIE eine Lesestelle.  Ein Schalter, der zweimal gelesen wird, driftet."""
    return os.environ.get(FLAG, "0") == "1"


# ── DIE DRITTE KENNZAHL: sp2-Zentren, die pyramidal gebaut wurden ────────────────
#
# WOZU (25.08.2026).  Die Defektrangliste ueber 79 495 Frames (aromrad6k_off, 1618
# Systeme) sagt, dass die zwei GROESSTEN Klassen dieselbe Chemie sind:
#     smiles_hyb-angle   22,18 %   Kristall 0,00 %   Masse 3911
#     pyramidal_sp2      20,02 %   Kristall 0,98 %   Masse 3028
# Zusammen 6939 -- mehr als die naechsten drei Klassen zusammen, und der Kristall
# liegt bei null: reiner Baufehler, kein Detektorrauschen.
#
# `_rg_score` war fuer beide BLIND: es kannte nur Kollisionen und Bindungslaengen.
# Damit konnte das Ruecknahmetor jede Reparatur schuetzen AUSSER der an der
# groessten Klasse.  Ein Verflacher, der 300 Zentren richtet und 30 verbiegt, wurde
# vom Tor durchgewinkt und vom Landungstor erschlagen.
#
# ⚠ DIE GROESSE IST DIE DES AUGES, NICHT MEINE EIGENE.  Erst wollte ich
# `_frame_assertions._oop` nehmen -- den Abstand des Zentrums von der Ebene seiner
# drei Nachbarn, in ANGSTROEM.  Das Auge misst aber den WALSH-WINKEL in GRAD
# (`find_pyramidalization.py:515`):
#       walsh_deg = atan2(d_oop, r_mean)
# also genau `_oop`, NORMIERT auf die mittlere Zentrum-Nachbar-Bindungslaenge.  Ein
# Angstroem-Mass mit einer Grad-Schwelle waere wieder ein Etikett statt einer
# Messung gewesen -- der Fehler, an dem der Atropisomer-Enumerator bis zum 23.08.
# lahmlag.  Deshalb steht hier die Formel des Auges, Zeichen fuer Zeichen.
#
# ⚠ WAS HIER BEWUSST GROEBER IST ALS DAS AUGE, und warum das reicht: das Auge waehlt
# seine Decke UMGEBUNGSABHAENGIG (`_PYR_ENV`, p99.9 aus 307k Kristallen).  Dieses
# Tor VERGLEICHT dagegen dasselbe Frame vor und nach der Kette -- dieselben Atome,
# dieselbe Umgebung.  Fuer ein "ist es schlechter geworden" genuegt ein KONSISTENTES
# Mass; die absolute Eichung entscheidet nur, welche Zentren ueberhaupt mitzaehlen.
# Der Boden 12,0 Grad ist trotzdem der des Auges (`find_graph_geometry._PYR_BUILD_MIN`),
# damit die Zahlen anschlussfaehig bleiben.
#
# ⚠ mu-BRUECKEN SIND NICHT PYRAMIDAL, SONDERN VERBRUECKEND -- die R3-Regel des Auges
# (`full_verdict.py:514-525`, auf CCDC belegt).  Ein Leichtatom an ZWEI oder mehr
# Metallen wird uebersprungen; sonst zaehlte jede mu2-Oxo-Ecke als Defekt.
#
# ⚠ EIGENER SCHALTER, und das ist keine Vorsicht, sondern Pflicht: `picopgate6k` und
# `hplacegate6k` laufen SEIT HEUTE MITTAG mit dem Zwei-Tupel, und der Bauer liest
# seine Datei bei JEDEM System neu.  Ohne eigenen Schalter haette diese Zeile zwei
# laufende A/B mitten im Lauf veraendert.  Vorgabe AUS -> Zwei-Tupel, byte-identisch.
FLAG_PYR = "DELFIN_FFFREE_REFINE_GATE_PYR"
_RG_PYR_FLOOR = 12.0                      # Grad -- find_graph_geometry._PYR_BUILD_MIN
_RG_PYR_ELEMS = frozenset(("C", "N", "O", "S", "P"))


def _rg_pyr_on() -> bool:
    """DIE eine Lesestelle fuer die dritte Kennzahl."""
    return os.environ.get(FLAG_PYR, "0") == "1"


def _rg_walsh(P, c: int, nb: Sequence[int]) -> Optional[float]:
    """Walsh-Winkel in Grad, Formel aus `find_pyramidalization.py:501-515`.
    Planares sp2 -> ~0; je pyramidaler, desto groesser."""
    a, b, d = int(nb[0]), int(nb[1]), int(nb[2])
    nrm = np.cross(P[b] - P[a], P[d] - P[a])
    ln = float(np.linalg.norm(nrm))
    if ln < 1e-9:
        return None
    d_oop = abs(float(np.dot(P[c] - P[a], nrm / ln)))
    r_mean = (float(np.linalg.norm(P[c] - P[a]))
              + float(np.linalg.norm(P[c] - P[b]))
              + float(np.linalg.norm(P[c] - P[d]))) / 3.0
    if r_mean < 1e-6:
        return None
    return math.degrees(math.atan2(d_oop, r_mean))


def _rg_pyr_count(syms: Sequence[str], P, adj) -> int:
    """Wieviele DREIBINDIGE Leichtatom-Zentren sind ueber den Boden pyramidalisiert?

    Dreibindig, weil ein sp2-Zentrum genau drei sigma-Nachbarn hat -- H zaehlt mit,
    wie beim Auge ("including every sp2 C-H", `full_verdict.py:507`).  Ohne Bindungs-
    ordnungen ist "drei Nachbarn an C/N/O/S/P" die beste verfuegbare sp2-Naeherung;
    ein echtes sp3-Zentrum hat vier und faellt damit von selbst heraus."""
    n_pyr = 0
    for i in range(len(syms)):
        if syms[i] not in _RG_PYR_ELEMS:
            continue
        nb = list(adj[i])
        if len(nb) != 3:
            continue
        if sum(1 for j in nb if _hp_metal(syms[j])) >= 2:
            continue                      # mu-Bruecke, keine Fehlplanarisierung (R3)
        w = _rg_walsh(P, i, nb)
        if w is not None and w > _RG_PYR_FLOOR:
            n_pyr += 1
    return n_pyr


def _rg_score(xyz: str) -> Tuple[int, ...]:
    """(n_clash, n_bond_out[, n_pyr]) -- alle "kleiner ist besser".

    Die dritte Stelle kommt nur mit `DELFIN_FFFREE_REFINE_GATE_PYR=1` dazu; ohne sie
    ist der Rueckgabewert das alte Zwei-Tupel.  Der Vergleich in `keep_better` ist
    ein Tupelvergleich und traegt beide Laengen -- solange BEIDE Seiten aus demselben
    Lauf stammen, und das tun sie per Konstruktion (derselbe Prozess, dieselbe
    Umgebung).  `(-1, -1)` bzw. `(-1, -1, -1)` heisst unlesbar."""
    _pyr = _rg_pyr_on()
    try:
        syms, P, _lines = _hp_read(xyz)
    except Exception:
        return (-1, -1, -1) if _pyr else (-1, -1)
    n = len(syms)
    if n < 2:
        return (0, 0, 0) if _pyr else (0, 0)
    adj = _hp_graph(syms, P)
    n_bond_out = 0
    for i in range(n):
        for j in adj[i]:
            if j <= i:
                continue
            tgt = _hp_cov(syms[i]) + _hp_cov(syms[j])
            if tgt <= 1e-6:
                continue
            r = float(np.linalg.norm(P[i] - P[j])) / tgt
            if r < _XH_FIRE_LO or r > _XH_FIRE_HI:
                n_bond_out += 1
    # Kollisionen: nur wirklich nichtgebundene Paare, sonst meldet jeder normale
    # 1,3-Kontakt eine Verletzung.
    nah: List[Set[int]] = [_hp_within(adj, i, _NONBONDED_MIN_HOPS - 1) for i in range(n)]
    n_clash = 0
    for i in range(n):
        si = syms[i]
        if _hp_metal(si):
            continue                          # die M-D-Sphaere ist keine Kollisionsfrage
        for j in range(i + 1, n):
            if j in nah[i]:
                continue
            sj = syms[j]
            if _hp_metal(sj):
                continue
            if si == "H" and sj == "H":
                thr = _HH_FLOOR
            elif si == "H" or sj == "H":
                thr = _H_HEAVY_FRAC * (_hp_vdw("H") + _hp_vdw(sj if si == "H" else si))
            else:
                thr = _RG_HEAVY_HEAVY_FRAC * (_hp_vdw(si) + _hp_vdw(sj))
            if float(np.linalg.norm(P[i] - P[j])) < thr:
                n_clash += 1
    if _pyr:
        return (n_clash, n_bond_out, _rg_pyr_count(syms, P, adj))
    return (n_clash, n_bond_out)


def keep_better(before: Sequence, after: Sequence):
    """Nach der Verfeinerungskette ist KEIN Frame schlechter als davor.

    `before` ist die Liste VOR der Kette, `after` die danach; beide (xyz, etikett).
    Zurueck kommt eine Liste derselben Laenge wie `after`:
      * Etikett war vorher da UND der Text hat sich geaendert UND die Bewertung ist
        schlechter -> das ORIGINAL,
      * sonst unveraendert das Ergebnis der Kette.

    Byte-identische Frames werden ueber einen Stringvergleich ausgeschlossen, BEVOR
    irgendetwas gerechnet wird -- der Grund, warum dieses Tor im Bau fast nichts kostet.
    """
    if not _rg_on() or not before or not after:
        return after

    def _rg_lbl(e):
        try:
            return e[1] if len(e) > 1 else ""
        except Exception:
            return ""

    def _rg_xyz(e):
        try:
            return e[0]
        except Exception:
            return None

    vorher: Dict[str, str] = {}
    for e in before:
        l, x = _rg_lbl(e), _rg_xyz(e)
        if l and x is not None and l not in vorher:
            vorher[l] = x

    out = list(after)
    n_rueck = n_geprueft = 0
    for k, e in enumerate(out):
        l, x = _rg_lbl(e), _rg_xyz(e)
        alt = vorher.get(l)
        if alt is None or x is None or alt == x:
            continue                          # neu (Enumerator) oder unveraendert -> frei
        n_geprueft += 1
        s_alt, s_neu = _rg_score(alt), _rg_score(x)
        # ⚠ BEIDE ZEILEN WAREN AUF DAS ZWEI-TUPEL FESTGENAGELT (bemerkt 25.08., bevor
        # die dritte Kennzahl je lief).  `== (-1, -1)` ist fuer `(-1, -1, -1)` FALSCH,
        # also waere ein unlesbares Frame nicht mehr uebersprungen worden; und der
        # Vergleich las nur Index 0 und 1, also waere `n_pyr` berechnet und dann
        # WEGGEWORFEN worden -- ein Zaehler ohne Wirkzeile, genau die Bauart, die hier
        # schon mehrfach als "dunkler Schalter" geendet ist.  Jetzt laengenunabhaengig.
        if -1 in s_alt or -1 in s_neu:
            continue                          # unlesbar -> nicht urteilen, durchlassen
        if len(s_alt) != len(s_neu):
            continue                          # kann nur bei Schalterwechsel MITTEN im
            # Lauf passieren; dann ist kein Vergleich moeglich und Durchlassen richtig.
        if any(b > a for a, b in zip(s_alt, s_neu)):
            try:
                out[k] = ((alt,) + tuple(e[1:])) if isinstance(e, tuple) else ([alt] + list(e[1:]))
            except Exception:
                continue
            n_rueck += 1
    if n_rueck:
        _LOG.info("refine-gate: %d von %d veraenderten Frames zurueckgenommen "
                  "(Kollisionen oder Bindungslaengen wurden schlechter)", n_rueck, n_geprueft)
    return out


def _rg_selbsttest() -> int:  # pragma: no cover - Werkzeug, kein Produktivpfad
    """Behauptungen dieses Moduls gegen Zahlen, nicht gegen Zuversicht.

        python -m delfin.manta._refine_gate

    Geprueft wird, was schiefgehen KANN, nicht was bequem ist:
      1. Vorgabe AUS  -> `_rg_score` liefert ein ZWEI-Tupel.  Ohne das waeren die
         laufenden A/B `picopgate6k` und `hplacegate6k` mitten im Lauf veraendert.
      2. Schalter AN  -> DREI-Tupel, und die dritte Zahl UNTERSCHEIDET flach von
         pyramidal.  Ein Zaehler, der auf beiden Formen dasselbe sagt, ist eine
         Tabellenkonstante und kein Detektor.
      3. Der Walsh-Winkel trifft die Zahlen des AUGES: planares Formaldehyd ~0 Grad,
         und die Formel ist `atan2(d_oop, r_mean)`, nicht `d_oop` allein.
      4. `keep_better` NIMMT eine Verschlechterung der dritten Zahl ZURUECK -- sonst
         waere sie berechnet und weggeworfen.
    """
    import numpy as _np

    def _xyz(rows):
        return "%d\ntest\n" % len(rows) + "\n".join(
            f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, x, y, z in rows) + "\n"

    # Formaldehyd-artig: C mit drei Nachbarn, exakt planar (z = 0 fuer alle)
    flach = _xyz([("C", 0.0, 0.0, 0.0), ("O", 0.0, 1.21, 0.0),
                  ("H", 0.94, -0.54, 0.0), ("H", -0.94, -0.54, 0.0)])
    # dasselbe C, aber 0.35 A aus der Ebene gezogen -> deutlich pyramidal
    pyr = _xyz([("C", 0.0, 0.0, 0.35), ("O", 0.0, 1.21, 0.0),
                ("H", 0.94, -0.54, 0.0), ("H", -0.94, -0.54, 0.0)])

    fehler = 0

    def _urteil(name, ok, zusatz=""):
        nonlocal fehler
        print(f"  {'OK  ' if ok else 'FEHL'} {name}{('  ' + zusatz) if zusatz else ''}")
        if not ok:
            fehler += 1

    _alt = os.environ.get(FLAG_PYR)
    _alt_haupt = os.environ.get(FLAG)
    try:
        os.environ[FLAG_PYR] = "0"
        s_aus = _rg_score(flach)
        _urteil("Vorgabe AUS liefert ein Zwei-Tupel", len(s_aus) == 2, f"-> {s_aus}")

        os.environ[FLAG_PYR] = "1"
        s_flach = _rg_score(flach)
        s_pyr = _rg_score(pyr)
        _urteil("Schalter AN liefert ein Drei-Tupel", len(s_flach) == 3, f"-> {s_flach}")
        _urteil("flaches sp2 zaehlt NICHT als pyramidal",
                len(s_flach) == 3 and s_flach[2] == 0, f"n_pyr={s_flach[2] if len(s_flach)>2 else '?'}")
        _urteil("gezogenes sp2 zaehlt SEHR WOHL",
                len(s_pyr) == 3 and s_pyr[2] >= 1, f"n_pyr={s_pyr[2] if len(s_pyr)>2 else '?'}")
        _urteil("die beiden Formen sind UNTERSCHEIDBAR (kein konstanter Zaehler)",
                len(s_flach) == 3 and len(s_pyr) == 3 and s_flach[2] != s_pyr[2])

        # Walsh-Winkel gegen die Formel des Auges nachgerechnet
        P = _np.array([[0.0, 0.0, 0.35], [0.0, 1.21, 0.0],
                       [0.94, -0.54, 0.0], [-0.94, -0.54, 0.0]])
        w = _rg_walsh(P, 0, [1, 2, 3])
        _urteil("Walsh-Winkel liegt ueber dem Boden", w is not None and w > _RG_PYR_FLOOR,
                f"{w:.1f} Grad, Boden {_RG_PYR_FLOOR}")
        Pf = _np.array([[0.0, 0.0, 0.0], [0.0, 1.21, 0.0],
                        [0.94, -0.54, 0.0], [-0.94, -0.54, 0.0]])
        wf = _rg_walsh(Pf, 0, [1, 2, 3])
        _urteil("planar ergibt ~0 Grad", wf is not None and wf < 1.0, f"{wf:.2f} Grad")

        # keep_better muss die Verschlechterung der DRITTEN Zahl zurueckholen.
        # ⚠ ZWEI SCHALTER, und der erste Testentwurf setzte nur EINEN -- die Probe
        # schlug fehl und sah wie ein Codefehler aus.  `keep_better` steigt in Zeile
        # 258 als ERSTES an `_rg_on()` aus, also am HAUPTSCHALTER; `FLAG_PYR` allein
        # bewirkt nichts.  Das ist keine Testkosmetik: ein Lauf, der nur
        # DELFIN_FFFREE_REFINE_GATE_PYR setzt, misst NICHTS und meldete es als
        # "keine Wirkung".  Die Warteschlangenzeile muss BEIDE tragen.
        os.environ[FLAG] = "1"
        vor = [(flach, "f0")]
        nach = [(pyr, "f0")]
        zurueck = keep_better(vor, nach)
        _urteil("keep_better nimmt die Pyramidalisierung ZURUECK",
                bool(zurueck) and zurueck[0][0] == flach)

        # und bei AUS darf es genau das NICHT tun (der Beweis, dass es an der
        # dritten Zahl lag und nicht an Kollisionen oder Bindungslaengen)
        os.environ[FLAG_PYR] = "0"
        zurueck_aus = keep_better(vor, nach)
        _urteil("bei Vorgabe AUS bleibt dieselbe Aenderung STEHEN",
                bool(zurueck_aus) and zurueck_aus[0][0] == pyr,
                "sonst kaeme die Ruecknahme woanders her")
    finally:
        for _f, _v in ((FLAG_PYR, _alt), (FLAG, _alt_haupt)):
            if _v is None:
                os.environ.pop(_f, None)
            else:
                os.environ[_f] = _v

    print(f"\n  {'ALLE PROBEN BESTANDEN' if fehler == 0 else str(fehler) + ' PROBE(N) GESCHEITERT'}")
    return 1 if fehler else 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys
    _sys.exit(_rg_selbsttest())
