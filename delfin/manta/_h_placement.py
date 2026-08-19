"""_h_placement -- EINE Wasserstoff-Reparatur, gebaut fuer den FF-freien Pfad.

DIE MESSUNG, DIE DIESES MODUL BEGRUENDET (19.08.2026, 9 Pools, 248 099 Frames).
Vier Detektorfamilien tragen zusammen rund 45 000 harte Frames = 28 % der
Haertemasse, und keine von ihnen hatte einen aktiven Mechanismus:

    xh_hh_clash + core_hhclash   16 881   10,7 %   H...H unter dem Kristallboden
    xh_stretch  + xh_orphan      14 165    9,0 %   X-H-Laenge / kein Elternteil
    methyl_broken                10 049    6,4 %   H-C-H mehr als 15 Grad daneben
    h_axis_H_proximal_via_donor   7 205    4,6 %   Donor-H zeigt AUF das Metall

Es ist EIN physikalischer Fehler unter sechs Namen: das Schweratomgeruest wird
gesetzt, gedreht, radial korrigiert -- und der Wasserstoff bleibt liegen, wo die
Einbettung ihn hingelegt hatte.  Beleg ist der garantierte Frame-Ueberlapp
(Schubfachschluss, Untergrenze): methyl_broken & smiles_hyb-angle 68,4 %,
smiles_topology & xh_hh_clash 68,1 %, smiles_hyb-angle & xh_hh_clash 67,8 %.

WARUM DIESE BEWEGUNGSKLASSE UND KEINE ANDERE.  Das Kostengesetz ist an vier
Punkten gemessen: Ordnung ~0 · Isometrie +0,98 pp · starre Drehung mit neuer
Konformation +6,57 pp · Neueinbettung +11,9 pp harte Frames.  Alle drei Stufen
hier ruehren AUSSCHLIESSLICH Wasserstoffatome an; jedes Schweratom und jedes
Metall bleibt byte-genau stehen (am Ende geprueft, nicht behauptet).  Das ist die
billigste Klasse ueberhaupt -- eine Neueinbettung waere die teuerste und muesste
nach demselben Gesetz scheitern.

DIE DREI STUFEN, in dieser Reihenfolge:

  A  DACH (methyl_broken).  Delegiert an ``_vsepr_repair.repair_terminal_groups``
     -- den vorhandenen, XYZ-only, nie widerlegten EX3-Reparateur -- und
     UEBERNIMMT danach ausschliesslich die Wasserstoffe.  Jede Schweratom-
     Verschiebung (CF3, SO3) wird verworfen.  Stufe A ist damit KEIN Neubau,
     sondern das Anschliessen eines vorhandenen Mechanismus.

  B  LAENGE (xh_stretch, xh_orphan).  Setzt die X-H-Distanz auf die
     Kovalenzradiensumme, RADIAL entlang der bestehenden Richtung.  Es aendert
     sich kein einziger Winkel, also auch kein Vorzeichenvolumen: Stufe B kann
     Stereochemie mathematisch nicht anfassen.

  C  ROTOR (xh_hh_clash, h_axis_H_proximal_via_donor).  Dreht die H-Gruppe eines
     Rotors (Schweratom mit GENAU EINEM schweren Nachbarn: CH3, NH2, OH, SH)
     starr um die X-Y-Achse auf den Rasterwinkel, der den engsten Kontakt
     maximiert.  Bindungslaengen und Bindungswinkel bleiben exakt erhalten.

STEREOCHEMIE IST HIER BEWEISBAR UNVERSEHRT, nicht hoffentlich.  Am 19.08. hat
eine Ebenenprojektion 61 sp3-Zentren plattgedrueckt und zwei Vorzeichen gekippt;
genau das kann hier nicht passieren:
  * Stufe A feuert nur auf Zentren mit DREI terminalen Nachbarn desselben
    Elements.  Drei konstitutionell gleiche Substituenten -> das Zentrum ist per
    Definition kein Stereozentrum.
  * Stufe B aendert nur Radien, keine Richtungen -> das Vorzeichen der
    Determinante ist invariant.
  * Stufe C feuert nur auf Zentren mit GENAU EINEM schweren Nachbarn; alle
    uebrigen Substituenten sind Wasserstoffe, also gleich -> kein Stereozentrum.
Gemessen wird es trotzdem: ``stereo_signature`` / ``stereo_delta`` vergleichen
das Vorzeichenvolumen vorher und nachher, und der Zensus druckt es mit.

SCHWELLEN.  Die AUSLOESEschwellen liegen absichtlich INNERHALB der Detektorbaender
des Auges (0,85/1,15 gegen die Detektorgrenzen 0,75/1,25), damit nie auf der
Kante repariert wird.  Die ZIELwerte sind physikalisch -- die Cordero-Kovalenz-
radiensumme aus ``_elements.COV_R``, dieselbe Konvention wie im Auge
(``find_xh_integrity``), und die referenzfreie Messung ueber 296 696 X-H-Bindungen
des Champion-Archivs bestaetigt sie als MEDIAN (aromatisch C-H 1,080 gegen Summe
1,07; N-H 1,034 gegen 1,02; O-H 0,990 gegen 0,97).  Falsch ist nicht die Regel,
falsch ist der SCHWANZ -- und genau den korrigiert Stufe B.  Es wird kein
Detektorwert nachgebaut und keine Schwelle angepeilt.

NAMENSGEBUNG.  Alle modulinternen Helfer tragen das Praefix ``_hp_``.  Das ist
kein Stil, sondern eine Lehre: die generischen Namen (`_adjacency`, `_is_metal`,
`_parse_xyz`) existieren im Doppelbaum ein Dutzend Mal mit LEICHT
verschiedener Semantik, und genau daran ist am 16.08. der Radienzensus haengen
geblieben.  Ein eindeutiger Name macht die Kopie sichtbar, statt sie zu tarnen.

VERDRAHTUNG.  Vorgabe AUS (``DELFIN_FFFREE_H_PLACEMENT``), und dieses Modul hat
heute NULL Aufrufstellen -- der Bau ist damit byte-identisch, weil im Baupfad
keine einzige Zeile anders ist.  Die noetige Aufrufstelle steht in
``smiles_converter.py`` unmittelbar vor dem FF-freien ``return _ff, None``, neben
``_apply_mirror_enum_if_enabled``; sie ist im Bericht benannt und wegen des
Dateireviers bewusst NICHT gebaut.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from delfin.manta import _elements as _el

# --- Schalter -------------------------------------------------------------
# DIE EINE LESESTELLE.  Ein zweiter os.environ.get auf denselben Namen ist genau
# die Bauart, an der der Feuerzensus am 14.08. gescheitert ist.
FLAG = "DELFIN_FFFREE_H_PLACEMENT"

# --- Schwellen ------------------------------------------------------------
# Ausloesung INNERHALB der Detektorbaender (Auge: collision < 0,75 · stretch > 1,25).
_XH_FIRE_LO = 0.85
_XH_FIRE_HI = 1.15
# Kristallboden fuer H...H; identisch zu find_xh_integrity (1,50 A).
_HH_FLOOR = 1.50
# Echte Ueberlappung H gegen schweres Nichtmetall: 0,85 x vdW-Summe, aber NUR
# fuer wirklich nichtgebundene Paare (Graphabstand >= 4).  Ohne diese Bedingung
# meldet jeder normale 1,3-Kontakt eine Verletzung.
_H_HEAVY_FRAC = 0.85
_NONBONDED_MIN_HOPS = 4
# Elternerkennung: identisch zu metric_h_axis._HEAVY_PARENT_BOND.
_PARENT_FACTOR = 1.45
# Bindungsgraph (nur Nichtmetalle), Faktor wie _bond_criterion / _isolated_reseat.
_BOND_FACTOR = 1.30
# Ein H weiter weg als das ist nicht "liegengeblieben", sondern verloren;
# es 1,5 A weit zu ziehen waere keine Isometrie mehr.
_MAX_PULL = 2.20
# Kriterium des Auges fuer "Donor-H zeigt auf das Metall" (metric_h_axis).
_PROX_H_M_MAX = 2.00
_PROX_DELTA = 0.05
# Hydriderkennung: NICHT "irgendein Metall in Reichweite", sondern "das NAECHSTE
# Schweratom IST ein Metall" -- exakt die Regel, mit der find_xh_integrity ein
# M-H als Hydrid fuehrt und aus der X-H-Achse herausnimmt.
#
# WARUM DIE ERSTE FASSUNG FALSCH WAR.  Sie schloss jedes H aus, dem IRGENDEIN
# Metall naeher als 2,00 A kam -- und der Befund ``h_axis_H_proximal_via_donor``
# IST per Definition ein H mit d(M,H) < 2,00.  Die Regel schloss also genau die
# Klasse aus, fuer die Stufe C gebaut ist; der modulinterne Zensus konnte sie
# nie anders als mit 0 melden, und das sah wie "kommt nicht vor" aus.
#
# ⚠ EIN REST BLEIBT, und er gehoert hierher und nicht in eine Fussnote.  Auf den
# 2144 gemessenen Frames faellt ``xh_n_m_hydride`` 21 -> 18, und die Ursache ist
# NICHT diese Regel, sondern eine METALLDEFINITION, die auseinandergeht:
# ``_elements.METALS`` fuehrt Te und Sb als METALLOIDE, ``find_xh_integrity``
# fuehrt sie als METALLE.  In ALEKOS sitzt ein H bei 1,249 A von Te und 1,281 A
# von C -- 0,03 A entscheiden, wer Elternteil ist.  Das Modul liest Te-H, setzt
# es auf die Kovalenzsumme 1,690 (Lehrbuch Te-H 1,65-1,70), und danach ist C mit
# 1,150 A das naechste Schweratom, also kein Hydrid mehr.  ``n_m_hydride_bad``
# bleibt in BEIDEN Armen 13, Stretch/Kollision/Waise sinken -- zerstoert wurde
# nichts.  Aber "naechstes Schweratom" ist auf gequetschten Frames FRAGIL, und
# das ist eine offene Schwaeche, keine erledigte Frage.
# Rotorraster: 24 Schritte a 15 Grad.  Fest, deterministisch, keine Optimierung.
_ROTOR_STEPS = 24
# Eine Verbesserung unter diesem Wert ist Rauschen und wird verworfen.
_EPS_GAIN = 1e-3
# Detektorschwellen fuer den modulinternen Zensus (siehe _hp_frame_census).
_DET_ORPHAN_MAX = 1.6
_DET_COLL_FRAC = 0.75
_DET_COLL_ABS = 0.90
_DET_STRETCH_FRAC = 1.25
_DET_STRETCH_ABS = 1.30
_DET_METHYL_DEG = 15.0
_TETRA_DEG = 109.471

_VDW: Dict[str, float] = {}
try:                                     # pragma: no cover - Importweg
    from delfin.manta._vdw_radii import VDW_RADII as _VDW  # type: ignore
except Exception:                        # pragma: no cover
    _VDW = {}
_VDW_DEFAULT = 1.80


# ---------------------------------------------------------------------------
# XYZ-Ein/Ausgabe -- Kopfzeilen und Atomreihenfolge bleiben unangetastet.
# ---------------------------------------------------------------------------
def _hp_atom_line(line: str) -> bool:
    parts = line.split()
    if len(parts) < 4 or not parts[0][:1].isalpha():
        return False
    try:
        float(parts[1]); float(parts[2]); float(parts[3])
    except (ValueError, IndexError):
        return False
    return True


def _hp_read(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Symbole, Koordinaten und die ORIGINALZEILEN.

    Es wird nichts umsortiert und nichts weggeworfen: ``_hp_write`` schreibt
    ausschliesslich die Zahlen der Atomzeilen neu, jede andere Zeile (Anzahl,
    Kommentar, Leerzeile) geht unveraendert durch.
    """
    syms: List[str] = []
    pts: List[List[float]] = []
    lines = xyz_str.splitlines()
    for line in lines:
        if not _hp_atom_line(line):
            continue
        parts = line.split()
        syms.append(parts[0])
        pts.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if not pts:
        return syms, np.zeros((0, 3), dtype=float), lines
    return syms, np.asarray(pts, dtype=float), lines


def _hp_write(orig_lines: Sequence[str], syms: Sequence[str],
              P: np.ndarray) -> str:
    out: List[str] = []
    k = 0
    for line in orig_lines:
        if _hp_atom_line(line) and k < len(syms):
            x, y, z = P[k]
            out.append(f"{syms[k]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            k += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if orig_lines else "")


# ---------------------------------------------------------------------------
# Elementwissen -- KEINE weitere Radientabelle.  _elements ist der Kanon.
# ---------------------------------------------------------------------------
def _hp_cov(sym: str) -> float:
    return _el.covalent_radius(sym)


def _hp_metal(sym: str) -> bool:
    """Delegiert an den Kanon ``_elements.is_metal`` -- kein eigenes Praedikat."""
    return _el.is_metal(sym)


def _hp_vdw(sym: str) -> float:
    return float(_VDW.get(_el.normalise(sym), _VDW_DEFAULT))


def _hp_xh_target(parent_sym: str) -> float:
    """Sollaenge einer X-H-Bindung = Kovalenzradiensumme (Cordero)."""
    return _hp_cov("H") + _hp_cov(parent_sym)


# ---------------------------------------------------------------------------
# Graph.
# ---------------------------------------------------------------------------
def _hp_hydrogens(syms: Sequence[str]) -> List[int]:
    return [i for i, s in enumerate(syms) if _el.normalise(s) == "H"]


def _hp_metals(syms: Sequence[str]) -> List[int]:
    return [i for i, s in enumerate(syms) if _hp_metal(s)]


def _hp_graph(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Bindungsgraph ueber Kovalenzradien; Metalle bleiben ausserhalb.

    Metalle werden bewusst NICHT verknuepft: die M-D-Bindung ist keine kovalente
    Bindung im Sinne dieser Abstandsregel, und sie mitzunehmen wuerde die halbe
    Koordinationssphaere zum 1,3-Nachbarn machen.  Anders als
    ``_isolated_reseat._adjacency`` (Metalle drin) und
    ``_vsepr_repair._adjacency`` (15 Elemente, Rueckfall 0,90) liest dieser Graph
    den Kanon ``_elements.COV_R``.
    """
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        si = _el.normalise(syms[i])
        if _hp_metal(si):
            continue
        for j in range(i + 1, n):
            sj = _el.normalise(syms[j])
            if _hp_metal(sj):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * (_hp_cov(si) + _hp_cov(sj)):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _hp_link_parents(adj: List[List[int]], parents: Dict[int, int]) -> None:
    """Die Elternbindung IN den Graphen eintragen -- und warum das noetig ist.

    GEMESSEN im Selbsttest, und es ist kein Sonderfall, sondern der Regelfall
    dieses Moduls: ein gestrecktes C-H bei 1,60 A liegt UEBER der Bindungsgrenze
    des Graphen (1,30 x 1,07 = 1,39), also gilt das H dort als ungebunden.  Dann
    ist sein 1,3-Nachbar nicht ausgeschlossen, der Nachbar zaehlt als echter
    nichtgebundener Kontakt -- und die Rueckrollung verwirft genau die Korrektur,
    die den Fehler behoben haette.  Der Reparateur haette sich an seinem eigenen
    Defekt blockiert.  ``parents_of_h`` kennt die Bindung (es zieht die weitere
    Grenze), also traegt sie sie hier nach.
    """
    for h, p in parents.items():
        if p not in adj[h]:
            adj[h].append(p)
        if h not in adj[p]:
            adj[p].append(h)


# ---------------------------------------------------------------------------
# DAS LOCH IM EIGENEN BEWEIS -- und der Weg, der NICHT genommen wurde.
#
# Gemessen 19.08. auf 2144 Frames: EIN echter Stereokandidat unter 1076 kippte
# sein Vorzeichen, bei unveraenderter Nachbarschaft.  Der Beweis "Stufe B
# aendert nur Radien, also kann das Vorzeichen nicht kippen" gilt nur, wenn das
# H genau EIN Schweratom als Nachbarn hat.  Ein gequetschtes H zwischen zwei
# Schweratomen ist im Abstandsgraphen an BEIDE gebunden; radial zu seinem
# Elternteil ist dann NICHT radial zum anderen Zentrum, und dessen Determinante
# darf kippen.
#
# ⛔ DER NAHELIEGENDE WEG WURDE GEMESSEN UND VERWORFEN.  "Solche H einfach
# auslassen" (ein Filter auf genau einen schweren Graphnachbarn) macht das Modul
# sicher und gleichzeitig fast wirkungslos -- denn die mehrdeutigen H SIND die
# Kollisionen.  Auf denselben 2144 Frames, Auge ``find_xh_integrity``:
#     ohne Filter   xh_n_collision 80 -> 48,  hh_n_hard 225 -> 142
#     mit  Filter   xh_n_collision 80 -> 77,  hh_n_hard 225 -> 212
# Der Filter kostet also rund 90 % der Wirkung auf der haertesten Stufe, um
# EINEN Frame von 2144 zu schuetzen.  Das ist der falsche Handel.
#
# ✅ STATTDESSEN EIN TOR STATT EINES VERBOTS: das Vorzeichenvolumen wird vor und
# nach der Reparatur gemessen, und ein Frame, in dem ein ECHTER Stereokandidat
# kippt oder plattgedrueckt wird, wird als GANZES verworfen.  Reichweite bleibt,
# der Schaden wird unmoeglich statt unwahrscheinlich -- dieselbe Bauart wie die
# Rueckrollungen in _isolated_reseat und _fix_sp3_h_tetrahedrality.
# ---------------------------------------------------------------------------


def _hp_within(adj: Sequence[Sequence[int]], start: int, hops: int) -> Set[int]:
    """Alle Atome bis ``hops`` Bindungen entfernt (inklusive ``start``)."""
    seen = {int(start)}
    front = [int(start)]
    for _ in range(hops):
        nxt: List[int] = []
        for a in front:
            for b in adj[a]:
                if b not in seen:
                    seen.add(b)
                    nxt.append(b)
        front = nxt
        if not front:
            break
    return seen


def parents_of_h(syms: Sequence[str], P: np.ndarray) -> Dict[int, int]:
    """H-Index -> Index des schweren Elternatoms (Nichtmetall).

    Drei Regeln, alle drei aus dem Auge uebernommen statt neu erfunden:
      * Ist das NAECHSTE Schweratom ein METALL, ist das H ein Hydrid und bekommt
        gar keinen Elternteil -- exakt die Klassifikation von
        ``find_xh_integrity``.  Hydride sind nicht unser Fehler.
      * Sonst ist der Elternteil das naechste schwere Nichtmetall.
      * Es zaehlt, wenn d < 1,45 x Sollaenge (normale Bindung) ODER
        d <= 2,20 A (liegengebliebenes H, das noch erreichbar ist).
    """
    out: Dict[int, int] = {}
    n = len(syms)
    if n == 0 or P.shape[0] != n:
        return out
    for h in _hp_hydrogens(syms):
        best_j = -1
        best_d = 1e9
        near_j = -1
        near_d = 1e9
        for j in range(n):
            if j == h:
                continue
            sj = _el.normalise(syms[j])
            if sj == "H":
                continue
            d = float(np.linalg.norm(P[h] - P[j]))
            if d < near_d:                 # naechstes Schweratom, Metall inklusive
                near_d = d
                near_j = j
            if _hp_metal(sj):
                continue
            if d < best_d:                 # naechstes schweres NICHTmetall
                best_d = d
                best_j = j
        if near_j >= 0 and _hp_metal(syms[near_j]):
            continue                       # Hydrid -> unberuehrt
        if best_j < 0:
            continue
        if (best_d < _PARENT_FACTOR * _hp_xh_target(syms[best_j])
                or best_d <= _MAX_PULL):
            out[h] = best_j
    return out


# ---------------------------------------------------------------------------
# Freiraum -- die EINE Zielgroesse aller drei Stufen.
# ---------------------------------------------------------------------------
def _hp_clearance(syms: Sequence[str], P: np.ndarray, h: int, parent: int,
                  skip: Set[int]) -> float:
    """Kleinster NORMIERTER Abstand des H zu allem, was ihn angeht.

    Normiert heisst d / Boden: >= 1 ist sauber, < 1 ist eine Verletzung.
    Drei Boeden, jeder aus dem Detektor, den er bedient:
      H...H            1,50 A                  (find_xh_integrity)
      H...Schweratom   0,85 x vdW-Summe        (nur bei Graphabstand >= 4)
      H...Metall       min(2,00, d(M,X)-0,05)  (metric_h_axis: das H darf dem
                       Metall nicht naeher sein als sein eigener Elternteil)
    Ein weit entferntes Metall bindet die Zielgroesse von selbst nicht -- der
    Quotient wird dann gross, ohne dass es dafuer eine Sonderregel braucht.
    """
    worst = 1e9
    for j in range(len(syms)):
        if j == h or j in skip:
            continue
        sj = _el.normalise(syms[j])
        d = float(np.linalg.norm(P[h] - P[j]))
        if _hp_metal(sj):
            d_mx = (float(np.linalg.norm(P[parent] - P[j]))
                    if parent >= 0 else 1e9)
            thr = min(_PROX_H_M_MAX, d_mx - _PROX_DELTA)
            if thr <= 1e-6:
                continue
        elif sj == "H":
            thr = _HH_FLOOR
        else:
            thr = _H_HEAVY_FRAC * (_hp_vdw("H") + _hp_vdw(sj))
        worst = min(worst, d / thr)
    return worst if worst < 1e9 else 1e9


def _hp_metal_clear(syms: Sequence[str], P: np.ndarray, h: int,
                    parent: int) -> float:
    """NUR die Metallachse, als EIGENE Zahl -- und warum sie eine sein muss.

    GEMESSEN 19.08. auf 2144 Frames: mit der Metallachse INNERHALB des Minimums
    von ``_hp_clearance`` stieg ``h_prox_donor`` von 48 auf 51.  Der Grund ist
    kein Chemiefehler, sondern Arithmetik: ein Minimum verdeckt seine
    Unterachsen.  Stand der engste Kontakt eines H bei 0,60 (H...H), durfte eine
    Korrektur die Metallachse von 1,20 auf 0,80 druecken -- das Minimum stieg
    trotzdem von 0,60 auf 0,80, die Rueckrollung sah eine Verbesserung, und der
    Donor-H-am-Metall-Befund war neu erzeugt.  Dieselbe Form wie
    "ein Mittelwert kann keine tote Klasse sehen": die Achse braucht ihr eigenes
    Tor, nicht einen Platz in einer Summe.
    """
    worst = 1e9
    for j in range(len(syms)):
        if j == h or not _hp_metal(syms[j]):
            continue
        d = float(np.linalg.norm(P[h] - P[j]))
        d_mx = (float(np.linalg.norm(P[parent] - P[j]))
                if parent >= 0 else 1e9)
        thr = min(_PROX_H_M_MAX, d_mx - _PROX_DELTA)
        if thr <= 1e-6:
            continue
        worst = min(worst, d / thr)
    return worst if worst < 1e9 else 1e9


def _hp_skip(syms: Sequence[str], adj: Sequence[Sequence[int]], h: int,
             parent: int, extra: Sequence[int] = ()) -> Set[int]:
    """Was fuer dieses H nicht als Kontakt zaehlt.

    Der eigene Elternteil, die mitgegebenen Gruppenmitglieder -- und alle
    SCHWERATOME innerhalb von drei Bindungen (nichtgebunden zaehlt ab vier).
    Wasserstoffe bleiben immer im Test: 1,3- und 1,4-H...H sind genau die
    Kontakte, die der Rotor entdrehen soll, und ihr Boden (1,50 A) liegt weit
    unter jedem gesunden geminalen Abstand (~1,78 A).
    """
    out: Set[int] = set(int(x) for x in extra)
    if parent >= 0:
        out.add(int(parent))
    for j in _hp_within(adj, h, _NONBONDED_MIN_HOPS - 1):
        if _el.normalise(syms[j]) != "H":
            out.add(int(j))
    return out


# ---------------------------------------------------------------------------
# STUFE A -- das Dach (methyl_broken).  Vorhandener Mechanismus, H-only gemacht.
# ---------------------------------------------------------------------------
def _hp_stage_umbrella(syms: List[str], P: np.ndarray,
                       tol_deg: float = _DET_METHYL_DEG) -> int:
    """``_vsepr_repair.repair_terminal_groups`` anwenden -- NUR die H behalten.

    Der vorhandene Reparateur setzt eine verzerrte terminale EX3-Gruppe auf ihre
    ideale VSEPR-Lage und laesst Zentrum und Anker stehen.  Bewegt werden also
    ausschliesslich Gruppenmitglieder, und die sind entweder alle H (CH3, NH3)
    oder alle schwer (CF3, SO3).  Indem hier jedes Nicht-H verworfen wird,
    bleibt exakt der H-Anteil uebrig; das Schweratomgeruest ist unberuehrt.

    ``tol_deg`` = 15 ist die Schwelle des Detektors ``methyl_broken`` selbst;
    die Modulvorgabe des Reparateurs (20) wuerde einen Teil der gemeldeten
    Faelle nicht anfassen.
    """
    try:
        from delfin.manta import _vsepr_repair as _vr
    except Exception:
        return 0
    try:
        block = "".join(
            f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n"
            for s, p in zip(syms, P)
        )
        fixed = _vr.repair_terminal_groups(block, tol=float(tol_deg))
        if fixed == block:
            return 0
        f_syms, f_P, _ = _hp_read(fixed)
        if len(f_syms) != len(syms) or f_P.shape[0] != P.shape[0]:
            return 0
        moved = 0
        for i, s in enumerate(syms):
            if _el.normalise(s) != "H":
                continue                   # Schweratome bleiben, wo sie sind
            if float(np.linalg.norm(f_P[i] - P[i])) <= 1e-9:
                continue
            P[i] = f_P[i]
            moved += 1
        return moved
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# STUFE B -- die Laenge (xh_stretch, xh_orphan).
# ---------------------------------------------------------------------------
def _hp_stage_length(syms: List[str], P: np.ndarray, parents: Dict[int, int],
                     adj: Sequence[Sequence[int]]) -> int:
    """X-H radial auf die Sollaenge setzen.  Richtung bleibt, Winkel bleiben."""
    moved = 0
    for h in sorted(parents):
        p = parents[h]
        v = P[h] - P[p]
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            continue                      # entartet: es gibt keine Richtung
        target = _hp_xh_target(syms[p])
        if target <= 1e-6:
            continue
        if _XH_FIRE_LO <= d / target <= _XH_FIRE_HI:
            continue
        skip = _hp_skip(syms, adj, h, p)
        before = _hp_clearance(syms, P, h, p, skip)
        before_m = _hp_metal_clear(syms, P, h, p)
        old = P[h].copy()
        P[h] = P[p] + v / d * target
        after = _hp_clearance(syms, P, h, p, skip)
        after_m = _hp_metal_clear(syms, P, h, p)
        # ZWEI Rueckrollungen, nicht eine.  Die erste schuetzt den engsten
        # Kontakt ueberhaupt, die zweite die Metallachse als eigene Groesse --
        # sonst verdeckt das Minimum sie (siehe _hp_metal_clear).
        if (after < min(1.0, before) - _EPS_GAIN
                or after_m < min(1.0, before_m) - _EPS_GAIN):
            P[h] = old
            continue
        moved += 1
    return moved


# ---------------------------------------------------------------------------
# STUFE C -- der Rotor (xh_hh_clash, h_axis_H_proximal_via_donor).
# ---------------------------------------------------------------------------
def rotor_groups(syms: Sequence[str], parents: Dict[int, int],
                 adj: Sequence[Sequence[int]]
                 ) -> List[Tuple[int, int, List[int]]]:
    """(Zentrum, schwerer Nachbar, H-Liste) fuer Zentren mit GENAU EINEM schweren
    Nachbarn und mindestens einem H.

    Diese Bedingung ist zugleich der Stereobeweis: alle uebrigen Substituenten
    des Zentrums sind Wasserstoffe, also konstitutionell gleich -- ein solches
    Zentrum kann kein Stereozentrum sein.
    """
    h_of: Dict[int, List[int]] = {}
    for h, p in parents.items():
        h_of.setdefault(p, []).append(h)
    out: List[Tuple[int, int, List[int]]] = []
    for c in sorted(h_of):
        heavy_nb = [j for j in adj[c] if _el.normalise(syms[j]) != "H"]
        if len(heavy_nb) != 1:
            continue
        out.append((c, heavy_nb[0], sorted(h_of[c])))
    return out


def _hp_rotate(pts: np.ndarray, origin: np.ndarray, axis: np.ndarray,
               ang: float) -> np.ndarray:
    """Rodrigues -- starre Drehung; Laengen und Winkel bleiben exakt erhalten."""
    k = axis / float(np.linalg.norm(axis))
    v = pts - origin
    c = float(np.cos(ang))
    s = float(np.sin(ang))
    return (origin + v * c + np.cross(k, v) * s
            + np.outer(v.dot(k), k) * (1.0 - c))


def _hp_stage_rotor(syms: List[str], P: np.ndarray, parents: Dict[int, int],
                    adj: Sequence[Sequence[int]]) -> int:
    """Rotor-H starr um die Zentrum-Nachbar-Achse auf den besten Rasterwinkel."""
    moved = 0
    for centre, nb, hs in rotor_groups(syms, parents, adj):
        axis = P[centre] - P[nb]
        if float(np.linalg.norm(axis)) < 1e-6:
            continue
        skips = {h: _hp_skip(syms, adj, h, centre,
                             extra=list(hs) + [centre, nb]) for h in hs}

        def _group_clear() -> float:
            return min(_hp_clearance(syms, P, h, centre, skips[h]) for h in hs)

        def _group_metal() -> float:
            return min(_hp_metal_clear(syms, P, h, centre) for h in hs)

        base = _group_clear()
        base_m = _group_metal()
        if base >= 1.0:
            continue                       # nichts verletzt -> nichts anfassen
        orig = P[hs].copy()
        best_ang = 0.0
        best_val = base
        for step in range(1, _ROTOR_STEPS):
            ang = 2.0 * np.pi * step / _ROTOR_STEPS
            P[hs] = _hp_rotate(orig, P[centre], axis, ang)
            if _group_metal() < min(1.0, base_m) - _EPS_GAIN:
                continue                   # Metallachse eigenstaendig geschuetzt
            val = _group_clear()
            if val > best_val + _EPS_GAIN:
                best_val = val
                best_ang = ang
        if best_ang == 0.0:
            P[hs] = orig
            continue
        P[hs] = _hp_rotate(orig, P[centre], axis, best_ang)
        moved += len(hs)
    return moved


# ---------------------------------------------------------------------------
# Stereochemie -- Vorzeichenvolumen, damit der Beweis eine Messung hat.
# ---------------------------------------------------------------------------
def stereo_signature(syms: Sequence[str], P: np.ndarray
                     ) -> List[Tuple[int, int, float, int, Tuple[int, ...]]]:
    """(Zentrum, Vorzeichen, |Volumen|, Anzahl H, Nachbarn) je 4-Nachbar-Zentrum.

    Die drei Nachbarn mit den kleinsten Indizes spannen die Determinante auf.
    Da dieses Modul die Atomreihenfolge nie aendert, sind vorher und nachher
    direkt vergleichbar.

    ⚠ DIE ANZAHL H IST NICHT SCHMUCK, SIE IST DIE HALBE MESSUNG.  Ein Zentrum
    mit ZWEI oder mehr Wasserstoffen hat zwei konstitutionell gleiche
    Substituenten und ist damit KEIN Stereozentrum -- sein Determinanten-
    vorzeichen bedeutet nichts, und ein Methyl umzuklappen kippt es
    zwangslaeufig.  Erste Messung auf 50 Archivsystemen: 12 Vorzeichenwechsel
    ueber 2204 Zentren, und die Frage "an welchen" beantwortet nur dieses Feld.
    Wer nur die Gesamtzahl liest, liest ein Artefakt.
    """
    adj = _hp_graph(syms, P)
    # Dieselbe Elternnachtragung wie im Reparateur -- sonst haette ein Zentrum
    # mit gestrecktem X-H vorher DREI und nachher VIER Nachbarn, und der
    # Vorher/Nachher-Vergleich betraefe verschiedene Mengen.
    _hp_link_parents(adj, parents_of_h(syms, P))
    out: List[Tuple[int, int, float]] = []
    for c in range(len(syms)):
        sc = _el.normalise(syms[c])
        if sc == "H" or _hp_metal(sc):
            continue
        nbrs = sorted(adj[c])
        if len(nbrs) != 4:
            continue
        a = P[nbrs[0]] - P[c]
        b = P[nbrs[1]] - P[c]
        d3 = P[nbrs[2]] - P[c]
        vol = float(np.dot(np.cross(a, b), d3))
        n_h = sum(1 for j in nbrs if _el.normalise(syms[j]) == "H")
        out.append((c, (1 if vol > 0 else (-1 if vol < 0 else 0)),
                    abs(vol), n_h, tuple(nbrs)))
    return out


def stereo_delta(before, after) -> Tuple[int, int, int, int, int, int]:
    """(gemeinsam, Wechsel, platt, Kandidaten, Wechsel_Kandidaten, Nachbarwechsel).

    "Kandidat" = Zentrum mit HOECHSTENS EINEM Wasserstoff.  Nur dort kann ein
    Vorzeichen ueberhaupt Chemie bedeuten; ab zwei H sind zwei Substituenten
    gleich, und die Determinante ist eine Nummerierungsfrage.
    "plattgedrueckt" = |Volumen| faellt unter 10 % seines Ausgangswerts -- die
    Signatur der Ebenenprojektion, die am 19.08. 61 Zentren zerstoert hat.

    ⚠ EIN ZENTRUM MIT GEAENDERTER NACHBARSCHAFT WIRD NICHT ALS WECHSEL GEZAEHLT,
    sondern eigens gemeldet.  Gemessen 19.08.: ein einziger scheinbarer Wechsel
    unter 1077 Kandidaten, und er kam daher, dass ein gequetschtes H seinen
    Elternteil wechselte (Te 1,249 A gegen C 1,281 A -- die Regel "naechstes
    Schweratom" entscheidet dort mit 0,03 A).  Damit spannen vorher und nachher
    VERSCHIEDENE Dreibeine auf; das Vorzeichen zu vergleichen waere sinnlos.
    Die Zahl verschwindet dadurch nicht, sie steht nur unter dem richtigen Namen.
    """
    b = {t[0]: t for t in before}
    common = flips = flat = cand = cand_flips = nb_changed = 0
    for c, s, v, nh, nbrs in after:
        if c not in b:
            continue
        common += 1
        _, s0, v0, nh0, nbrs0 = b[c]
        if nbrs0 != nbrs:
            nb_changed += 1
            continue
        flipped = (s0 != 0 and s != 0 and s0 != s)
        if flipped:
            flips += 1
        if v0 > 1e-6 and v < 0.10 * v0:
            flat += 1
        if nh0 <= 1 and nh <= 1:
            cand += 1
            if flipped:
                cand_flips += 1
    return common, flips, flat, cand, cand_flips, nb_changed


# ---------------------------------------------------------------------------
# Oeffentliche Schnittstelle.
# ---------------------------------------------------------------------------
def h_placement_enabled() -> bool:
    """DIE EINE Lesestelle von DELFIN_FFFREE_H_PLACEMENT (Vorgabe 0)."""
    return os.environ.get(FLAG, "0") == "1"


def repair_xyz(xyz: str, *, stats: Optional[dict] = None) -> str:
    """Ungattert reparieren -- fuer Selbsttest und Messung, NICHT im Baupfad.

    Gibt das EINGABEOBJEKT unveraendert zurueck, wenn kein H bewegt wurde, und
    ebenso, wenn die Schlusspruefung ein bewegtes Schweratom findet: dann ist
    der Lauf ungueltig, und ein ungueltiger Lauf darf nichts kosten.
    """
    if not xyz:
        return xyz
    try:
        syms, P, lines = _hp_read(xyz)
        if P.shape[0] == 0 or not _hp_hydrogens(syms):
            return xyz
        P = P.astype(float).copy()
        frozen = P.copy()
        sig_before = stereo_signature(syms, frozen)
        n_a = _hp_stage_umbrella(syms, P)
        adj = _hp_graph(syms, P)
        parents = parents_of_h(syms, P)
        _hp_link_parents(adj, parents)
        n_b = _hp_stage_length(syms, P, parents, adj)
        n_c = _hp_stage_rotor(syms, P, parents, adj)
        # Der Anspruch dieses Moduls, geprueft statt behauptet.
        for i, s in enumerate(syms):
            if _el.normalise(s) != "H" and float(
                    np.linalg.norm(P[i] - frozen[i])) > 1e-9:
                if stats is not None:
                    stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                                  "moved": 0, "aborted_heavy_moved": 1,
                                  "aborted_stereo": 0})
                return xyz
        moved = n_a + n_b + n_c
        if moved == 0:
            if stats is not None:
                stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                              "moved": 0, "aborted_heavy_moved": 0,
                              "aborted_stereo": 0})
            return xyz
        # DAS STEREOTOR.  Kein Verbot vorab, ein Urteil danach: kippt in diesem
        # Frame ein ECHTER Stereokandidat (hoechstens ein H am Zentrum) sein
        # Vorzeichen, oder wird ein Zentrum plattgedrueckt, faellt der GANZE
        # Frame zurueck.  Siehe den Block ueber _hp_within: der billigere Weg
        # (solche H gar nicht erst anfassen) kostet ~90 % der Wirkung auf der
        # haertesten Stufe und wurde deshalb verworfen.
        _, _, _flat, _, _cand_flips, _ = stereo_delta(
            sig_before, stereo_signature(syms, P))
        if _cand_flips or _flat:
            if stats is not None:
                stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                              "moved": 0, "aborted_heavy_moved": 0,
                              "aborted_stereo": 1})
            return xyz
        if stats is not None:
            stats.update({"umbrella": n_a, "length": n_b, "rotor": n_c,
                          "moved": moved, "aborted_heavy_moved": 0,
                          "aborted_stereo": 0})
        return _hp_write(lines, syms, P)
    except Exception:
        return xyz


def apply_xyz(xyz: str) -> str:
    """Gattert reparieren.  AUS -> das Eingabeobjekt, unveraendert."""
    if not h_placement_enabled():
        return xyz
    return repair_xyz(xyz)


def apply_to_results(results):
    """Frameliste (xyz, label, ...) -> Frameliste.  Gattert, ausfallsicher.

    Form wie ``_h_vsepr_realism.correct_results`` und ``_me_bond_snap``, damit
    die Aufrufstelle am FF-freien Ausgang wie ihre Nachbarn aussieht.  Die
    Frameanzahl aendert sich nie: dies ist ein Korrektor, kein Enumerator.
    """
    if not results or not h_placement_enabled():
        return results
    out = []
    for entry in results:
        try:
            out.append((repair_xyz(entry[0]),) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Zensus im Modul -- Richtungsanzeige, KEIN Ersatz fuer das Auge.
# ---------------------------------------------------------------------------
def _hp_frame_census(syms: Sequence[str], P: np.ndarray) -> Dict[str, int]:
    """Die Defektfamilien nach den VEROEFFENTLICHTEN Schwellen des Auges.

    ⚠ Dies ist eine NACHBILDUNG der Schwellen aus ``find_xh_integrity``,
    ``full_verdict.adapt_methyl_quality`` und ``metric_h_axis`` -- keine zweite
    Meinung des Auges.  Sie zeigt die RICHTUNG der Wirkung im Modul selbst; das
    Urteil faellt der echte Detektor auf den geschriebenen Framepaaren.
    """
    n = len(syms)
    res = {"xh_stretch": 0, "xh_collision": 0, "xh_orphan": 0,
           "hh_clash": 0, "methyl_broken": 0, "h_prox_donor": 0}
    hs = _hp_hydrogens(syms)
    metals = _hp_metals(syms)
    for h in hs:
        best_j, best_d = -1, 1e9
        for j in range(n):
            if j == h or _el.normalise(syms[j]) == "H":
                continue
            d = float(np.linalg.norm(P[h] - P[j]))
            if d < best_d:
                best_d, best_j = d, j
        if best_j < 0 or _hp_metal(syms[best_j]):
            continue                       # M-H ist ein eigener Befund
        if best_d > _DET_ORPHAN_MAX:
            res["xh_orphan"] += 1
            continue
        target = _hp_xh_target(syms[best_j])
        if best_d < _DET_COLL_FRAC * target or best_d < _DET_COLL_ABS:
            res["xh_collision"] += 1
        elif best_d > _DET_STRETCH_FRAC * target or best_d > _DET_STRETCH_ABS:
            res["xh_stretch"] += 1
    for a in range(len(hs)):
        for b in range(a + 1, len(hs)):
            if float(np.linalg.norm(P[hs[a]] - P[hs[b]])) < _HH_FLOOR:
                res["hh_clash"] += 1
    parents = parents_of_h(syms, P)
    h_of: Dict[int, List[int]] = {}
    for h, p in parents.items():
        h_of.setdefault(p, []).append(h)
    for c, hlist in sorted(h_of.items()):
        if len(hlist) != 3:
            continue
        worst = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                u = P[hlist[i]] - P[c]
                v = P[hlist[j]] - P[c]
                nn = float(np.linalg.norm(u) * np.linalg.norm(v))
                if nn < 1e-9:
                    continue
                ang = float(np.degrees(np.arccos(
                    max(-1.0, min(1.0, float(np.dot(u, v)) / nn)))))
                worst = max(worst, abs(ang - _TETRA_DEG))
        if worst > _DET_METHYL_DEG:
            res["methyl_broken"] += 1
    for h, p in sorted(parents.items()):
        for m in metals:
            d_mh = float(np.linalg.norm(P[h] - P[m]))
            if d_mh >= _PROX_H_M_MAX:
                continue
            if d_mh <= float(np.linalg.norm(P[p] - P[m])) - _PROX_DELTA:
                res["h_prox_donor"] += 1
                break
    return res


def split_frames(text: str) -> List[str]:
    """Mehrframe-XYZ (aneinandergehaengt) in Einzelframes zerlegen."""
    lines = text.splitlines()
    frames: List[str] = []
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if len(head) == 1 and head[0].isdigit():
            k = int(head[0])
            block = lines[i:i + 2 + k]
            if len(block) == 2 + k:
                frames.append("\n".join(block) + "\n")
            i += 2 + k
        else:
            i += 1
    return frames


# ---------------------------------------------------------------------------
# Selbsttest.
# ---------------------------------------------------------------------------
def _hp_selftest() -> int:
    fails = 0

    def _expect(name: str, ok: bool, detail: str = "") -> None:
        nonlocal fails
        if not ok:
            fails += 1
        mark = "ok  " if ok else "FAIL"
        print(f"  [{mark}] {name}" + (f" -- {detail}" if detail else ""))

    print("== Schalter ==")
    os.environ.pop(FLAG, None)
    _expect("Vorgabe AUS", h_placement_enabled() is False)
    stretched = ("5\ntest\n"
                 "C       0.000000     0.000000     0.000000\n"
                 "C       1.520000     0.000000     0.000000\n"
                 "H      -0.500000     1.520000     0.000000\n"
                 "H      -0.363000    -0.520000     0.901000\n"
                 "H      -0.363000    -0.520000    -0.901000\n")
    _expect("AUS ist byte-identisch", apply_xyz(stretched) == stretched)
    os.environ[FLAG] = "1"
    _expect("AN wird gelesen", h_placement_enabled() is True)

    print("== Stufe B: Laenge ==")
    st: Dict[str, int] = {}
    out = repair_xyz(stretched, stats=st)
    _s0, P0, _ = _hp_read(stretched)
    _s1, P1, _ = _hp_read(out)
    d_before = float(np.linalg.norm(P0[2] - P0[0]))
    d_after = float(np.linalg.norm(P1[2] - P1[0]))
    _expect("gestrecktes C-H wird auf die Sollaenge gesetzt",
            d_before > 1.5 and abs(d_after - 1.07) < 1e-3,
            f"{d_before:.3f} -> {d_after:.3f}")
    _expect("Schweratome unbewegt",
            float(np.linalg.norm(P1[0] - P0[0])) < 1e-9
            and float(np.linalg.norm(P1[1] - P0[1])) < 1e-9)
    v0 = P0[2] - P0[0]
    v1 = P1[2] - P1[0]
    cosang = float(np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1)))
    _expect("Richtung unveraendert (kein Winkel bewegt sich)",
            abs(cosang - 1.0) < 1e-9)
    _expect("Zaehlung meldet Stufe B", st.get("length", 0) >= 1, str(st))

    print("== Stufe B: Waise ==")
    orphan = ("4\ntest\n"
              "N       0.000000     0.000000     0.000000\n"
              "C       1.470000     0.000000     0.000000\n"
              "H       0.000000     1.850000     0.000000\n"
              "H      -0.340000    -0.480000     0.830000\n")
    so, Po, _ = _hp_read(orphan)
    c0 = _hp_frame_census(so, Po)
    sr, Pr, _ = _hp_read(repair_xyz(orphan))
    c1 = _hp_frame_census(sr, Pr)
    _expect("verwaistes H bekommt seinen Elternteil zurueck",
            c0["xh_orphan"] == 1 and c1["xh_orphan"] == 0,
            f"{c0['xh_orphan']} -> {c1['xh_orphan']}")

    print("== Stufe A: Dach ==")
    # C2 haengt an C1, damit C1 KEINE terminale Gruppe ist -- sonst haelt
    # repair_terminal_groups das Zentrum fuer ankerlos und ruehrt es nicht an.
    broken = ("6\ntest\n"
              "C       0.000000     0.000000     0.000000\n"
              "C       1.520000     0.000000     0.000000\n"
              "H      -0.400000     0.980000     0.000000\n"
              "H      -0.400000    -0.490000     0.849000\n"
              "H      -1.070000     0.000000    -0.100000\n"
              "C       2.040000     1.430000     0.000000\n")
    sb, Pb, _ = _hp_read(broken)
    cb0 = _hp_frame_census(sb, Pb)
    sb1, Pb1, _ = _hp_read(repair_xyz(broken))
    cb1 = _hp_frame_census(sb1, Pb1)
    _expect("gebrochenes Methyl wird repariert",
            cb0["methyl_broken"] == 1 and cb1["methyl_broken"] == 0,
            f"{cb0['methyl_broken']} -> {cb1['methyl_broken']}")
    _expect("Methylreparatur laesst die Schweratome stehen",
            float(np.linalg.norm(Pb1[0] - Pb[0])) < 1e-9
            and float(np.linalg.norm(Pb1[1] - Pb[1])) < 1e-9)

    print("== Stufe C: Rotor ==")
    # Methyl an C0-C1-C2, dazu ein N-H, dessen H genau auf einem Methyl-H sitzt
    # (1,09 A, unter dem Kristallboden 1,50).  Alle X-H stehen exakt auf ihrer
    # Sollaenge, damit die Stufen A und B nachweislich nichts beitragen.
    import math as _m
    _rows = ["9", "test",
             f"C    {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}",
             f"C    {1.53:12.6f} {0.0:12.6f} {0.0:12.6f}",
             f"C    {2.05:12.6f} {1.45:12.6f} {0.0:12.6f}"]
    _r = 1.07 * _m.sin(_m.radians(109.471))
    _x = -1.07 * _m.cos(_m.radians(109.471)) * -1.0
    for k in range(3):
        a = 2 * _m.pi * k / 3
        _rows.append(f"H    {_x:12.6f} {_r * _m.cos(a):12.6f} "
                     f"{_r * _m.sin(a):12.6f}")
    _rows.append(f"N    {_x:12.6f} {_r + 2.11:12.6f} {0.0:12.6f}")
    _rows.append(f"H    {_x:12.6f} {_r + 1.09:12.6f} {0.0:12.6f}")
    _rows.append(f"C    {_x + 1.47:12.6f} {_r + 2.11:12.6f} {0.0:12.6f}")
    clashing = "\n".join(_rows) + "\n"
    se, Pe, _ = _hp_read(clashing)
    ce0 = _hp_frame_census(se, Pe)
    _st_c: Dict[str, int] = {}
    se1, Pe1, _ = _hp_read(repair_xyz(clashing, stats=_st_c))
    ce1 = _hp_frame_census(se1, Pe1)
    _expect("H...H-Konflikt wird entdreht",
            ce0["hh_clash"] > 0 and ce1["hh_clash"] < ce0["hh_clash"],
            f"{ce0['hh_clash']} -> {ce1['hh_clash']}")
    _expect("nur der Rotor hat gearbeitet",
            _st_c.get("umbrella", -1) == 0 and _st_c.get("length", -1) == 0
            and _st_c.get("rotor", 0) > 0, str(_st_c))
    _expect("Rotor laesst ALLE Schweratome stehen",
            all(float(np.linalg.norm(Pe1[i] - Pe[i])) < 1e-9
                for i, s in enumerate(se) if s != "H"))
    dh0 = float(np.linalg.norm(Pe[3] - Pe[0]))
    dh1 = float(np.linalg.norm(Pe1[3] - Pe1[0]))
    _expect("Rotor erhaelt die C-H-Laenge exakt", abs(dh0 - dh1) < 1e-6,
            f"{dh0:.6f} -> {dh1:.6f}")

    print("== Stereochemie ==")
    chiral = ("5\ntest\n"
              "C       0.000000     0.000000     0.000000\n"
              "F       1.350000     0.000000     0.000000\n"
              "Cl     -0.560000     1.680000     0.000000\n"
              "Br     -0.640000    -0.900000     1.640000\n"
              "H      -0.300000    -0.420000    -0.760000\n")
    sc0, Pc0, _ = _hp_read(chiral)
    sig0 = stereo_signature(sc0, Pc0)
    sc1, Pc1, _ = _hp_read(repair_xyz(chiral))
    common, flips, flat, cand, cand_flips, _nbch = stereo_delta(
        sig0, stereo_signature(sc1, Pc1))
    _expect("Stereozentrum erkannt", cand >= 1,
            f"gemeinsam={common} Kandidaten={cand}")
    _expect("kein Vorzeichenwechsel", flips == 0 and cand_flips == 0,
            f"Wechsel={flips}/{cand_flips}")
    _expect("nicht plattgedrueckt", flat == 0, f"flach={flat}")
    # Ein Methyl ist KEIN Stereozentrum -- es muss als Nicht-Kandidat gelten,
    # sonst zaehlt der Zensus jedes entdrehte Methyl als Stereoschaden.
    _mc = [t for t in stereo_signature(sb, Pb) if t[3] >= 2]
    _expect("Methylzentrum ist kein Stereokandidat", len(_mc) >= 1,
            f"Zentren mit >=2 H: {len(_mc)}")

    print("== Idempotenz ==")
    once = repair_xyz(stretched)
    _expect("zweiter Lauf aendert nichts mehr", once == repair_xyz(once))

    print("== Metallnaehe ==")
    # naechstes Schweratom IST das Metall -> Hydrid, tabu.
    hydride = ("3\ntest\n"
               "Fe      0.000000     0.000000     0.000000\n"
               "H       1.600000     0.000000     0.000000\n"
               "C       3.500000     0.000000     0.000000\n")
    _expect("terminales Hydrid wird nicht angefasst",
            repair_xyz(hydride) == hydride)
    # Donor-H, das AUF das Metall zeigt: naechstes Schweratom ist der Donor,
    # also KEIN Hydrid -- diese Klasse muss der Reparateur sehen koennen.
    prox = ("5\ntest\n"
            "Fe      0.000000     0.000000     0.000000\n"
            "N       2.150000     0.000000     0.000000\n"
            "C       3.620000     0.000000     0.000000\n"
            "H       1.900000     0.480000     0.000000\n"
            "H       2.500000    -0.900000     0.400000\n")
    sp_, Pp_, _ = _hp_read(prox)
    _expect("Donor-H am Metall bekommt seinen Elternteil (kein Hydrid)",
            3 in parents_of_h(sp_, Pp_) and parents_of_h(sp_, Pp_)[3] == 1,
            str(parents_of_h(sp_, Pp_)))
    _expect("und der Zensus meldet ihn",
            _hp_frame_census(sp_, Pp_)["h_prox_donor"] == 1,
            str(_hp_frame_census(sp_, Pp_)))

    os.environ.pop(FLAG, None)
    _expect("Schalter wieder AUS", h_placement_enabled() is False)
    return fails


def _hp_run_census(paths: List[str], limit: int = 0) -> None:
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    keys = ["xh_stretch", "xh_collision", "xh_orphan", "hh_clash",
            "methyl_broken", "h_prox_donor"]
    tot_b = {k: 0 for k in keys}
    tot_a = {k: 0 for k in keys}
    fr_b = {k: 0 for k in keys}
    fr_a = {k: 0 for k in keys}
    n_frames = n_changed = n_ident = 0
    st_sum = {"umbrella": 0, "length": 0, "rotor": 0,
              "aborted_heavy_moved": 0, "aborted_stereo": 0}
    s_common = s_flips = s_flat = s_cand = s_cand_flips = s_nbch = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        for frame in split_frames(text):
            syms, P, _ = _hp_read(frame)
            if P.shape[0] == 0:
                continue
            n_frames += 1
            cb = _hp_frame_census(syms, P)
            for k in keys:
                tot_b[k] += cb[k]
                fr_b[k] += 1 if cb[k] else 0
            sig0 = stereo_signature(syms, P)
            st: Dict[str, int] = {}
            out = repair_xyz(frame, stats=st)
            for k in st_sum:
                st_sum[k] += st.get(k, 0)
            if out == frame:
                n_ident += 1
                for k in keys:
                    tot_a[k] += cb[k]
                    fr_a[k] += 1 if cb[k] else 0
                continue
            n_changed += 1
            syms2, P2, _ = _hp_read(out)
            ca = _hp_frame_census(syms2, P2)
            for k in keys:
                tot_a[k] += ca[k]
                fr_a[k] += 1 if ca[k] else 0
            c, f, fl, cd, cf, nbc = stereo_delta(
                sig0, stereo_signature(syms2, P2))
            s_common += c
            s_flips += f
            s_flat += fl
            s_cand += cd
            s_cand_flips += cf
            s_nbch += nbc
    print(f"Dateien {len(files)}  Frames {n_frames}  "
          f"veraendert {n_changed}  unveraendert {n_ident}")
    print(f"bewegte H: Dach {st_sum['umbrella']}  Laenge {st_sum['length']}  "
          f"Rotor {st_sum['rotor']}  "
          f"Abbruch-Schweratom {st_sum['aborted_heavy_moved']}  "
          f"Abbruch-Stereotor {st_sum['aborted_stereo']}")
    print(f"{'Befund':<16}{'Treffer vor':>13}{'nach':>9}"
          f"{'Frames vor':>13}{'nach':>9}")
    for k in keys:
        print(f"{k:<16}{tot_b[k]:>13}{tot_a[k]:>9}{fr_b[k]:>13}{fr_a[k]:>9}")
    print(f"Zentren mit 4 Nachbarn verglichen {s_common}  "
          f"Vorzeichenwechsel {s_flips}  plattgedrueckt {s_flat}")
    print(f"davon ECHTE Stereokandidaten (<=1 H) {s_cand}  "
          f"Vorzeichenwechsel {s_cand_flips}")
    print(f"Zentren mit GEAENDERTER Nachbarschaft (nicht vergleichbar) {s_nbch}")


def _hp_byteid(paths: List[str], limit: int = 0) -> int:
    """AUS-Beweis auf ECHTEN Frames: identisches OBJEKT, nicht nur gleicher Text.

    Der Schalter wird hier NICHT gesetzt.  Geprueft wird die gattierte
    Schnittstelle (``apply_xyz`` / ``apply_to_results``) -- also genau das, was
    eine kuenftige Aufrufstelle benutzen wuerde.  ``is``-Vergleich statt ``==``:
    ein gleicher String waere schon gut, dasselbe Objekt ist besser, weil es
    beweist, dass nicht einmal neu formatiert wurde.
    """
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    n = bad = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        frames = split_frames(text)
        for frame in frames:
            n += 1
            if apply_xyz(frame) is not frame:
                bad += 1
        got = apply_to_results([(f, "x") for f in frames])
        for k, f in enumerate(frames):
            if got[k][0] is not f:
                bad += 1
    print(f"Schalter AUS: {n} Frames, ueber apply_xyz UND apply_to_results, "
          f"nicht-identische Rueckgaben: {bad}")
    print("BYTE-IDENTISCH" if bad == 0 else "VERLETZT")
    return bad


def _hp_write_pairs(paths: List[str], out_before: str, out_after: str,
                    limit: int = 0) -> None:
    """Vorher/Nachher als zwei Verzeichnisse -- damit der ECHTE Detektor urteilt."""
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    os.makedirs(out_before, exist_ok=True)
    os.makedirs(out_after, exist_ok=True)
    n = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        stem = os.path.splitext(os.path.basename(fp))[0]
        for i, frame in enumerate(split_frames(text)):
            syms, P, _ = _hp_read(frame)
            if P.shape[0] == 0:
                continue
            with open(os.path.join(out_before, f"{stem}_{i:04d}.xyz"), "w") as fh:
                fh.write(frame)
            with open(os.path.join(out_after, f"{stem}_{i:04d}.xyz"), "w") as fh:
                fh.write(repair_xyz(frame))
            n += 1
    print(f"{n} Framepaare geschrieben: {out_before} / {out_after}")


if __name__ == "__main__":
    import sys

    _argv = sys.argv[1:]

    def _hp_pop_limit(rest: List[str]) -> Tuple[List[str], int]:
        if "--limit" in rest:
            k = rest.index("--limit")
            return rest[:k] + rest[k + 2:], int(rest[k + 1])
        return rest, 0

    if _argv and _argv[0] == "--census":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ[FLAG] = "1"
        _hp_run_census(_rest, limit=_lim)
        raise SystemExit(0)
    if _argv and _argv[0] == "--byteid":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ.pop(FLAG, None)
        raise SystemExit(1 if _hp_byteid(_rest, limit=_lim) else 0)
    if _argv and _argv[0] == "--pairs":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ[FLAG] = "1"
        _hp_write_pairs(_rest[:-2], _rest[-2], _rest[-1], limit=_lim)
        raise SystemExit(0)
    print("== _h_placement Selbsttest ==")
    _n_fail = _hp_selftest()
    print("== " + ("ALLE BESTANDEN" if _n_fail == 0
                   else f"{_n_fail} FEHLER") + " ==")
    raise SystemExit(1 if _n_fail else 0)
