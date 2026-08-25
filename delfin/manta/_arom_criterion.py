"""EINE QUELLE fuer das Aromatizitaetstor (19.08.2026).

DAS PROBLEM WAR NICHT DIE HOEHE DER SCHWELLE, SONDERN IHRE EINHEIT.
Drei Module entschieden bis heute unabhaengig voneinander, ob ein 5-/6-Ring aus
C/N/O/S aromatisch ist, und alle drei taten es mit derselben Zeile:

    mittlere Ringbindung < 1.46 A

Gemessen an 910 Kristallen (3237 Aromaten gegen 505 Nicht-Aromaten,
``agent_workspace/FORENSIK_2026_08_19_hapto42/``, Schritte E/G/I) trifft dieses
Tor die Aromaten zu 99,35 % -- und laesst gleichzeitig **54,46 %** aller
Nicht-Aromaten durch.  Der ausgewogene Fehler betraegt 27,55 %.  Je Klasse:

    Oxazolin 92,6 %   Imidazolin 96,7 %   THF 28,8 %
    Piperidin 30,0 %  Cyclohexan 18,3 %   Cyclopentan 0 %

DER GRUND IST DIE ABSOLUTE ANGSTROEM-SKALA.  1.46 A ist eine C-C-Zahl.  Ein Ring
mit O oder N hat von sich aus kuerzere Bindungen (C-O 1.43, C-N 1.47 im
GESAETTIGTEN Zustand), also faellt ein gesaettigtes Oxazolin unter eine Schwelle,
die fuer Kohlenstoff gebaut wurde -- ohne dass daran etwas aromatisch waere.

UND DIESELBE SKALA HAT EINEN ZWEITEN, BISHER UNBENANNTEN DEFEKT: THIOPHEN wird
von JEDER absoluten Schwelle zu **0 %** erkannt.  C-S ist 1.72 A lang, das
Ringmittel eines echten Thiophens betraegt 1.518 A -- es faellt schon heute durch
1.46 und wurde nie verflacht.

DIE ELEMENTNORMIERUNG REPARIERT BEIDES MIT EINER GROESSE.  Statt der rohen Laenge
wird das Mittel von ``d_ij / (r_i + r_j)`` bewertet -- die Bindung wird an der
Kovalenzradiensumme ihrer EIGENEN Elemente gemessen.  Vier Groessen wurden auf
denselben Kristallringen gegeneinander geprueft (Schritt G):

    Groesse                    Schwelle   Aromaten   Nicht-Arom.   ausgewogen
    roher Mittelwert (IST)     < 1.460     99,35 %     54,46 %       27,55 %
    roher Mittelwert           < 1.410     97,62 %      4,95 %        3,66 %
    Streuung                   < 0.035     94,19 %     34,06 %       19,93 %
    Verhaeltnis d/(r_i+r_j)    < 0.939     99,26 %      5,94 %        3,34 %

⚠ DIE STREUUNG IST NICHT DIE ANTWORT (fuenfmal schlechter): sie trennt Oxazolin
tadellos, ist aber BLIND fuer gesaettigte Carbocyclen -- Cyclohexan streut 0,014,
so gleichmaessig wie Benzol.

Preis der Normierung, ehrlich: Thiophen 0 -> 95,2 %, aber Pyrrol 100 -> 84,8 %
und Furan 100 -> 66,7 % (n=6).  Ein nicht erkannter Aromat wird nicht verflacht --
das ist ENTGANGENE VERFLACHUNG, kein neuer Defekt.

⚠ WARUM DIESE DATEI UEBERHAUPT EXISTIERT.  Am 18.08. wurde genau an diesem Tor
schon einmal repariert (das sp3-Verbot), und der Selbsttest deckte nur EINE der
drei Kopien ab.  Eine Zahl, die dreimal im Baum steht, laeuft frueher oder
spaeter auseinander.  Die drei Aufrufer sind:

    delfin/manta/_arom_planarize.py            ``_detect_aromatic_rings``
    delfin/manta/_aromatic_ring_flattener.py   ``_detect_aromatic_rings``
    delfin/manta/_bond_decollapse.py           ``_aromatic_ring_bonds``

⚠ IMPORTKOSTEN NULL, UND ZWAR ABSICHTLICH.  Dieses Modul zieht auf Modulebene NUR
``os`` und ``typing``.  Die Radientabelle wird erst beim ERSTEN Ring geholt, den
das normierte Tor bewertet -- also nie, solange der Schalter aus ist.  Damit kann
der Umbau auch die IMPORTREIHENFOLGE nicht verschieben: ``polyhedra`` liest beim
Import ``DELFIN_FFFREE_COV_COMPLETE``, und ein frueherer Import haette diesen
Lesezeitpunkt vor ``_apply_construction_env`` ziehen koennen.  Genau das passiert
hier nicht.  ``_bond_decollapse`` darf deshalb auch auf Modulebene importieren,
ohne seine Zirkelfreiheit zu verlieren.

Schalter ``DELFIN_FFFREE_AROM_CRITERION_RADII``, Vorgabe 0 -> byte-identisch.
"""
from __future__ import annotations

import os
from typing import Dict, Optional, Sequence, Tuple

# Der historische Ist-Stand: roher Mittelwert der Ringbindungen in Angstroem.
# Diese Zahl stand bis 19.08.2026 dreimal im Baum; sie steht jetzt EINMAL hier.
_AROMATIC_BOND_MAX: float = 1.46

# Elementnormiert: Mittel von d_ij / (r_i + r_j).  Optimum des AUSGEWOGENEN
# Fehlers auf 3237 Aromaten + 505 Nicht-Aromaten aus 910 Kristallen (Schritt G).
#
# WAS 0.939 IN ANGSTROEM HEISST -- die Schwelle je Bindungstyp (Cordero-Summen):
#     C-C 1.427   C-N 1.380   C-O 1.333   C-S 1.700   N-N 1.333   O-O 1.239
# ⚠ FUER REINE CARBOCYCLEN IST DAS SCHAERFER ALS HEUTE (1.427 statt 1.460), und
# auf BAU-Frames -- deren Bindungen laenger sind als im Kristall -- faellt dadurch
# ein spuerbarer Anteil echter Benzolringe heraus.  Ein nicht erkannter Aromat wird
# nicht verflacht: entgangene Verflachung, kein neuer Defekt.  Fuer S-haltige Ringe
# ist es umgekehrt viel WEITER (1.700 statt 1.460) -- genau die Thiophen-Reparatur.
# ⚠ 25.08.2026: 0.939 -> 0.963.  AUF BAU-FRAMES NEU GEEICHT, nicht auf Kristallen.
#
# WARUM DIE ALTE ZAHL FALSCH ANGEWANDT WAR.  0.939 ist das Optimum auf 910
# KRISTALLEN.  `aromrad6k` hat sie am 24.08. auf BAU-Frames gemessen, und das
# Verdikt zeigt beide Seiten der Medaille:
#     pyramidal_sp2         15911 -> 15360 Frames   (Masse 3028 -> 2837)
#     smiles_ccdc_regressed    35 Systeme           <- der Preis
# Also im Mittel besser, im Einzelfall zu scharf -- genau die Unsicherheit, die in
# der Vorregistrierung als DIE Frage des Laufs benannt war.
#
# DIE NEUE EICHUNG (`harness/arom_eichung_baurahmen.py`, 1500 Systeme):
#     Schwelle   Aromat verpasst   Nicht-Ar. durch   ausgewogen
#     0.939 alt      13,34 %            5,41 %         9,37 %
#     0.963 neu       2,14 %           10,81 %         6,48 %
# In Angstroem fuer C-C: 1.427 -> 1.464.  Damit liegt die Schwelle fuer reine
# Carbocyclen wieder dort, wo die historische rohe Zahl 1.460 stand -- sie war fuer
# Carbocyclen NIE das Problem; die Normierung hat sie versehentlich um 0.033 A
# verschaerft.  Fuer S-haltige Ringe bleibt der Thiophen-Gewinn erhalten (C-S jetzt
# 1.743 statt 1.700, weiterhin weit ueber 1.46).
#
# ⚠ NICHT ZIRKULAER GEMESSEN, und das war die eigentliche Schwierigkeit: die
# WAHRHEIT ("dieser Ring ist aromatisch") darf nicht aus der Bindungslaenge kommen,
# denn das ist die zu eichende Groesse.  Sie kommt aus der PLANARITAET DES RINGES
# IM KRISTALL -- eine Groesse, die mit der Bindungslaenge nichts zu tun hat und die
# der Bau nicht beeinflusst.  Gemessen wird dann dasselbe Ringmotiv im BAU-Frame.
#
# ⚠ EMPFINDLICHKEITSPROBE, weil die Flachheitsschwelle 0.08 A frei gewaehlt ist:
#     flach< 0.04 -> Optimum 0.954     flach< 0.06 -> 0.963     flach< 0.08 -> 0.963
#     flach< 0.12 -> 0.963             flach< 0.18 -> 0.963
# Vier von fuenf ergeben dieselbe Zahl; nur die schaerfste weicht ab, und dort ist
# die Nicht-Aromaten-Klasse mit n=46 am duennsten.  Die Zahl haengt also nicht an
# der gewaehlten Flachheit.
#
# ⚠ WAS AN DER MESSUNG DUENN IST, ehrlich benannt: 1402 Aromaten gegen 37
# Nicht-Aromaten.  Die Seite "Aromat verpasst" -- die den SCHADEN bestimmt -- steht
# auf ueber 1400 Ringen und traegt.  Die Seite "Nicht-Aromat faelschlich durch"
# steht auf 37 und ist grob; sie kostet aber nur eine unnoetige Verflachung, keinen
# Defekt.  Deshalb ist der Fehler hier bewusst asymmetrisch gewichtet worden.
#
# 🔴 DER PREIS HAT EINEN NAMEN, und der Selbsttest zeigt ihn: OXAZOLIN liegt bei
# 0.965.  Mit 0.939 wurde es mit einem Abstand von 0.026 verworfen; mit 0.963 ist
# der Abstand nur noch **0.002**.  Die Probe 11-13 haelt (0.965 >= 0.963), aber sie
# haelt knapp -- und sie rechnet mit einer IDEALGEOMETRIE.  Auf einem Bau-Frame
# streut derselbe Ring, und ein etwas kuerzer geratenes Oxazolin wird jetzt
# faelschlich verflacht.  Genau das ist die gemessene Verdopplung der falsch
# durchgelassenen Nicht-Aromaten (5,41 % -> 10,81 %), hier an einem benannten
# Stoff statt an einer Prozentzahl.
# ⇒ WORAUF IM VERDIKT ZU SCHAUEN IST: `smiles_ccdc_regressed` und
#   `pyramidal_sp2` auf Systemen mit Oxazolin/Imidazolin.  Steigt dort etwas,
#   ist 0.963 fuer diese Ringfamilie zu weit und die Schwelle gehoert
#   ELEMENTABHAENGIG gestaffelt (C-C anders als C-N-O), nicht global gesenkt.
#   Das waere die naechste Verfeinerung -- eine Zahl fuer alle Ringe ist selbst
#   eine Naeherung, und diese Messung zeigt ihre Grenze.
#
# Vorgabe des SCHALTERS unveraendert AUS -> byte-identisch.
_AROM_RATIO_MAX: float = 0.963

# Nur falls ein Ringatom kein Kovalenzradius hat.  Kann bei den Aufrufern nicht
# feuern (alle drei filtern die Ringe vorher auf C/N/O/S), steht aber auf
# demselben Wert wie in der Messung, damit die Zahl anschlussfaehig bleibt.
_COV_FALLBACK: float = 0.90

_COV_TABLE: Optional[Dict[str, float]] = None


def radii_criterion_enabled() -> bool:
    """Vorgabe AUS -> byte-identisch zum Auslieferungsstand.

    Beim AUFRUF gelesen, nicht beim Import: ``_bond_decollapse`` liest seine
    eigenen Schalter zwar beim Import (und ``_apply_construction_env`` laeuft in
    ``cli_manta`` bei :428 vor dem Import bei :435), aber ein Tor, das seine
    Wirkung nur bei richtiger Importreihenfolge entfaltet, ist genau die Bauart,
    die in diesem Projekt schon mehrfach als "dunkler Schalter" geendet hat.
    """
    return os.environ.get("DELFIN_FFFREE_AROM_CRITERION_RADII", "0") == "1"


def _cov_table() -> Dict[str, float]:
    """Die VORHANDENE Cordero-Tabelle aus ``polyhedra``, keine zweite Kopie.

    Erst hier geholt (siehe Modul-Docstring: Importreihenfolge).  ``polyhedra.COV``
    traegt C 0.76, N 0.71, O 0.66, S 1.05 -- also alle vier Elemente, die die
    Aufrufer als ringfaehig zulassen -- und stimmt fuer diese vier ziffernweise
    mit ``_bond_decollapse._COV`` ueberein, der Tabelle, auf der die Schwelle
    0.939 gemessen wurde.  Der Schalter ``DELFIN_FFFREE_COV_COMPLETE`` fuellt
    dort nur LUECKEN per ``setdefault`` und kann C/N/O/S nicht verschieben.
    """
    global _COV_TABLE
    if _COV_TABLE is None:
        from delfin.manta.polyhedra import COV
        _COV_TABLE = COV
    return _COV_TABLE


def ring_rejected_by_length(
    syms: Sequence[str],
    edges: Sequence[Tuple[int, int]],
    lens: Sequence[float],
) -> bool:
    """Faellt der Ring durch das Laengentor, ist also KEIN Aromat?

    ``edges`` und ``lens`` sind parallel: ``lens[k]`` ist die Laenge der Bindung
    ``edges[k]``; ``syms`` ist die volle Elementliste des Frames.

    ⚠ DER VERGLEICH IST BEWUSST ``>=``, NICHT ``not <``.  Die drei Aufrufer
    schrieben bisher ``if mittel >= 1.46: continue``; mit einem NaN (kaputter
    Frame) verhaelt sich ``>=`` anders als ``<``, und ein Umbau, der das
    stillschweigend dreht, waere nicht byte-identisch.
    """
    if not lens or len(edges) != len(lens):
        return True
    if radii_criterion_enabled():
        cov = _cov_table()
        acc = 0.0
        for (i, j), d in zip(edges, lens):
            acc += d / (cov.get(syms[i], _COV_FALLBACK)
                        + cov.get(syms[j], _COV_FALLBACK))
        return (acc / len(lens)) >= _AROM_RATIO_MAX
    return (sum(lens) / len(lens)) >= _AROMATIC_BOND_MAX


# ---------------------------------------------------------------------------
# Selbsttest -- prueft ALLE DREI Kopien, nicht nur eine:
#   PYTHONPATH=<worktree> python -m delfin.manta._arom_criterion
# (ueber -m, nicht ueber den Dateipfad: die editable-Installation zoege sonst
#  eine ANDERE Modulkopie und der Test masse den falschen Baum.)
# ---------------------------------------------------------------------------
def _selbsttest() -> int:  # pragma: no cover - Werkzeug, kein Produktivpfad
    import math

    from delfin.manta._pi_h_projector import (
        _parse_xyz, _build_geometric_adjacency,
    )
    from delfin.manta import _arom_planarize as AP
    from delfin.manta import _aromatic_ring_flattener as AF
    from delfin.manta import _bond_decollapse as BD

    stand = {"n": 0, "fehl": 0}

    def _urteil(name, ok):
        stand["n"] += 1
        if not ok:
            stand["fehl"] += 1
        print(f"{stand['n']:2d} {name}: {'OK' if ok else 'FEHLER'}")

    def _setze(v):
        os.environ["DELFIN_FFFREE_AROM_CRITERION_RADII"] = v

    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}")
        return "\n".join(out) + "\n"

    def _ring(el, laengen, z_of=None):
        """Ebener Ring mit VORGEGEBENEN Bindungslaengen, H nach aussen.

        Die Laengen werden ueber die Innenwinkel gesetzt: bei n gleichen Sehnen
        waere der Radius r = d/(2 sin(pi/n)); fuer ungleiche Sehnen wird der
        Radius so gewaehlt, dass die Summe der Zentriwinkel 2*pi ergibt.
        """
        n = len(el)

        def _summe(r):
            return sum(2.0 * math.asin(min(1.0, d / (2.0 * r))) for d in laengen)

        lo, hi = max(laengen) / 2.0 + 1e-6, 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _summe(mid) > 2 * math.pi:
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)
        rows, ang = [], 0.0
        winkel = []
        for k in range(n):
            winkel.append(ang)
            ang += 2.0 * math.asin(min(1.0, laengen[k] / (2.0 * r)))
        for k in range(n):
            zz = 0.0 if z_of is None else float(z_of.get(k, 0.0))
            rows.append((el[k], r * math.cos(winkel[k]), r * math.sin(winkel[k]), zz))
        for k in range(n):
            if el[k] == "S":
                continue                     # Ring-S traegt kein H
            dh = {"C": 1.08, "N": 1.01, "O": 0.96}.get(el[k], 1.08)
            rows.append(("H", (r + dh) * math.cos(winkel[k]),
                         (r + dh) * math.sin(winkel[k]), 0.0))
        return rows

    # --- Testringe.  Laengen sind Kristallmediane (Schritt G/I). --------------
    BENZOL = _ring(["C"] * 6, [1.39] * 6, {2: 0.12})
    # Thiophen: C-S 1.71, S-C 1.71, C-C 1.37, C=C 1.42, C-C 1.37
    #   -> roher Mittelwert 1.516 > 1.46, also vom Ist-Stand VERWORFEN.
    THIOPHEN = _ring(["S", "C", "C", "C", "C"], [1.71, 1.37, 1.42, 1.37, 1.71],
                     {2: 0.12})
    # Oxazolin (gesaettigt an C4/C5): O-C 1.35, C=N 1.28, N-C 1.47, C-C 1.54,
    #   C-O 1.44  -> roher Mittelwert 1.416 < 1.46, also vom Ist-Stand DURCHGELASSEN.
    OXAZOLIN = _ring(["O", "C", "N", "C", "C"], [1.35, 1.28, 1.47, 1.54, 1.44],
                     {3: 0.20})
    CYCLOHEXAN = _ring(["C"] * 6, [1.52] * 6, {1: 0.25})

    def _kopie1(x):
        s, p, _ = _parse_xyz(x)
        return len(AP._detect_aromatic_rings(s, p, _build_geometric_adjacency(s, p)))

    def _kopie2(x):
        s, p, _ = _parse_xyz(x)
        return len(AF._detect_aromatic_rings(s, p, _build_geometric_adjacency(s, p)))

    def _kopie3(x):
        # ⚠ Kopie 3 liefert BINDUNGEN, nicht Ringe -- ein erkannter 5-Ring sind
        # dort fuenf Eintraege.  Die Erwartungswerte unten tragen das mit.
        s, p, _ = BD._parse(x)
        return len(BD._aromatic_ring_bonds(s, p, BD._geometric_bonds(s, p)))

    KOPIEN = [("_arom_planarize", _kopie1),
              ("_aromatic_ring_flattener", _kopie2),
              ("_bond_decollapse", _kopie3)]

    # --- 1) das PRAEDIKAT selbst, ohne Umweg ueber die Ringsucher -------------
    #     Thiophen als RING: die fuenf Bindungen mitteln sich roh auf 1.516 A
    #     (> 1.46, also verworfen), normiert auf 0.925 (< 0.939, also gehalten).
    #     ⚠ Eine EINZELNE C-S-Bindung (1.71/1.81 = 0.945) faellt fuer sich
    #     genommen durch -- entschieden wird der RING, nicht die Bindung.
    t_syms = ["S", "C", "C", "C", "C"]
    t_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    t_lens = [1.71, 1.37, 1.42, 1.37, 1.71]
    _setze("0")
    a_off = ring_rejected_by_length(t_syms, t_edges, t_lens)
    _setze("1")
    a_on = ring_rejected_by_length(t_syms, t_edges, t_lens)
    _urteil("Thiophenring: absolut verworfen, normiert gehalten "
            f"(aus={a_off}, an={a_on})", a_off is True and a_on is False)

    _setze("0")
    c6 = ["C"] * 6
    e6 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    b_off = ring_rejected_by_length(c6, e6, [1.52] * 6)
    _setze("1")
    b_on = ring_rejected_by_length(c6, e6, [1.52] * 6)
    _urteil(f"Cyclohexanring bleibt in beiden verworfen (aus={b_off}, an={b_on})",
            b_off is True and b_on is True)

    # --- 2) die Radien kommen aus der VORHANDENEN Quelle ----------------------
    tab = _cov_table()
    _urteil("polyhedra.COV traegt C/N/O/S mit den Messwerten",
            all(abs(tab.get(e, -1) - v) < 1e-12 for e, v in
                (("C", 0.76), ("N", 0.71), ("O", 0.66), ("S", 1.05))))
    _urteil("polyhedra.COV stimmt fuer C/N/O/S mit _bond_decollapse._COV ueberein",
            all(abs(tab[e] - BD._COV[e]) < 1e-12 for e in ("C", "N", "O", "S")))

    # --- 3) bis 6) die vier Sonden durch ALLE DREI Kopien, in BEIDEN
    #        Schalterstellungen.  Erwartung je Kopie als (aus, an); Kopie 3
    #        zaehlt BINDUNGEN, deshalb 6 bzw. 5 statt 1.
    #
    #   BENZOL     roh 1.390 / normiert 0.914  -> beide Tore: Aromat
    #   THIOPHEN   roh 1.516 / normiert 0.925  -> nur das normierte Tor: Aromat
    #                                             (der ZWEITE, bisher unbenannte
    #                                              Defekt: 0 % bei jeder Angstroem-Schwelle)
    #   OXAZOLIN   roh 1.416 / normiert 0.965  -> nur das absolute Tor laesst durch
    #   CYCLOHEXAN roh 1.520 / normiert 1.000  -> beide Tore: kein Aromat
    SONDEN = [
        ("Benzol", BENZOL, [(1, 1), (1, 1), (6, 6)],
         "roh 1.390 / normiert 0.914 -- beide Tore halten ihn"),
        ("Thiophen", THIOPHEN, [(0, 1), (0, 1), (0, 5)],
         "roh 1.516 -- absolut BLIND, normiert 0.925 erkannt"),
        ("Oxazolin", OXAZOLIN, [(1, 0), (1, 0), (5, 0)],
         "roh 1.416 laesst durch, normiert 0.965 haelt"),
        ("Cyclohexan", CYCLOHEXAN, [(0, 0), (0, 0), (0, 0)],
         "in beiden Stellungen verworfen"),
    ]
    for probe_name, rows, erwartet, warum in SONDEN:
        x = _xyz(rows)
        for (name, fn), (e0, e1) in zip(KOPIEN, erwartet):
            _setze("0")
            n0 = fn(x)
            _setze("1")
            n1 = fn(x)
            _urteil(f"{probe_name} in {name}: {warum} "
                    f"(aus={n0}/{e0}, an={n1}/{e1})", n0 == e0 and n1 == e1)

    # --- 7) SCHALTER AUS ist byte-identisch: die Korrektoren, nicht nur die
    #        Detektoren.  Ein Detektor kann uebereinstimmen und der Pass
    #        trotzdem eine andere Datei schreiben.
    for probe, name in ((BENZOL, "Benzol"), (OXAZOLIN, "Oxazolin"),
                        (THIOPHEN, "Thiophen"), (CYCLOHEXAN, "Cyclohexan")):
        x = _xyz(probe)
        _setze("0")
        p1, p2, p3 = (AP.correct_xyz(x),
                      AF.correct_results(None, [(x, "l")]),
                      BD.correct_xyz(None, x))
        _setze("0")
        q1, q2, q3 = (AP.correct_xyz(x),
                      AF.correct_results(None, [(x, "l")]),
                      BD.correct_xyz(None, x))
        _urteil(f"Schalter AUS reproduzierbar auf {name} (alle drei Korrektoren)",
                p1 == q1 and p2 == q2 and p3 == q3)

    # --- 8) und die drei Kopien sehen NICHT dieselbe Ringmenge (Befund, kein
    #        Fehler): Kopie 3 perzipiert Bindungen mit 1.30*ideal statt
    #        Sigma_r_cov + 0.25 und wirft jeden Ring aus der Koordinationssphaere.
    fe = list(THIOPHEN) + [("Fe", 0.0, 0.0, 1.70)]
    x_fe = _xyz(fe)
    _setze("1")
    k1, k3 = _kopie1(x_fe), _kopie3(x_fe)
    _urteil(f"eta5-Thiophen am Fe: Kopie 1 sieht {k1}, Kopie 3 sieht {k3} "
            "(Koordinationsausschluss NUR in Kopie 3)", k1 == 1 and k3 == 0)

    _setze("0")
    print(f"\n{stand['n'] - stand['fehl']}/{stand['n']} bestanden")
    return 1 if stand["fehl"] else 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys
    _sys.exit(_selbsttest())
