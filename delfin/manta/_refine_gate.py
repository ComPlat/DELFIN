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
import os
from typing import Dict, List, Sequence, Set, Tuple

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


def _rg_score(xyz: str) -> Tuple[int, int]:
    """(n_clash, n_bond_out) -- beide "kleiner ist besser".  (-1, -1) heisst unlesbar."""
    try:
        syms, P, _lines = _hp_read(xyz)
    except Exception:
        return (-1, -1)
    n = len(syms)
    if n < 2:
        return (0, 0)
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
        if s_alt == (-1, -1) or s_neu == (-1, -1):
            continue                          # unlesbar -> nicht urteilen, durchlassen
        if s_neu[0] > s_alt[0] or s_neu[1] > s_alt[1]:
            try:
                out[k] = ((alt,) + tuple(e[1:])) if isinstance(e, tuple) else ([alt] + list(e[1:]))
            except Exception:
                continue
            n_rueck += 1
    if n_rueck:
        _LOG.info("refine-gate: %d von %d veraenderten Frames zurueckgenommen "
                  "(Kollisionen oder Bindungslaengen wurden schlechter)", n_rueck, n_geprueft)
    return out
