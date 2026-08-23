"""_atropisomer_enum.py — ADDITIVE Vollstaendigkeit fuer AXIALE Chiralitaet (Atropisomere).

DIE LUECKE, GEMESSEN AM 16.08.2026 (nicht vermutet).  Die Achse `ccdc_atropisomer_realized`
stand auf `rows_spy5geom` (965 Systeme) bei **0 von 44** -- kein einziges System realisiert die
axiale Haendigkeit seines Kristalls.  Bei zufaelliger Haendigkeit waeren ~22 Treffer zu
erwarten; null ist systematisch.

DIE URSACHE IST NICHT DIE GEOMETRIE.  Ueber alle 44 Systeme und **145 Achsen** gemessen
(Verdrillung ueber ALLE Frames, nicht nur Frame 0):

    mittlere Verdrillung   42.8 Grad      89 von 145 Achsen >= 30 Grad, nur 22 < 10 Grad
    Vorzeichen im Manifold:  beide 48 (33%) | NUR EINES 59 (41%) | nie verlaesslich 38 (26%)

Der Bauer erzeugt echte, verdrillte Achsen -- `BIRVUW` erreicht 85.2 Grad.  Was fehlt, ist die
GEGENHAENDIGKEIT.  Und die Achse verlangt, dass JEDE Kristallachse ihr Vorzeichen im Manifold
findet: bei ~3.3 Achsen je System und nur einem Drittel Achsen mit beiden Haenden ist ein
Volltreffer rechnerisch fast ausgeschlossen.  **Das erklaert 0/44 quantitativ, und es macht die
Enumeration zur Antwort -- nicht eine Geometriekorrektur.**

MESOMERIE (User, 16.08.).  `BIRVUW` traegt die Achsen N29-C30 und N18-C19 -- **Aryl-Amid**-
Achsen.  Ihre Stereogenitaet entsteht ueberhaupt erst durch den partiellen C-N-Doppelbindungs-
charakter: das freie N-Elektronenpaar konjugiert in die Carbonylgruppe, die Bindung bekommt
Doppelbindungsanteil, und daraus folgen Barriere und bevorzugte Verdrillung.  Wer sie als
Einfachbindung behandelt, sieht dort gar keine Achse.  Der Fall ist weiterhin abgedeckt --
seit 23.08. aber durch den pi-Zweig von `_atrop_is_sp2_center` (ein Amid-N mit C-N <= 1,37 A
faellt von selbst hinein), nicht mehr durch einen eigenen Sonderfall.  Der Unterschied ist
nicht kosmetisch: der Sonderfall pruefte dieselbe Chemie mit einer ANDEREN Grenze als das
Auge, und damit sah der Bauer dort Achsen, die das Auge nicht fuehrte, und umgekehrt.

🔴 WARUM DAS MODUL BIS 23.08.2026 NICHTS BEWIRKT HAT -- GEMESSEN, NICHT VERMUTET.
`harness/atrop_schluessel_vergleich.py` haelt beide Achsenerkennungen auf DENSELBEN Frames
gegeneinander.  Auf 120 Systemen / 3350 Frames:

                              vorher      nachher
    Achsen Bauer                 112          282
    Signatur verschieden     103/103        0/282
    Vorzeichen verschieden        87            0
    Achsen nur im Bauer            9            0

**Auf JEDER gemeinsamen Bindung war die Signatur verschieden, auf 84 % auch das Vorzeichen.**
Damit war jede angehaengte Spiegelung ein Frame im FALSCHEN Eimer -- und ein Frame im
falschen Eimer kostet genau so viel wie ein fehlendes.  Zwei unabhaengige Ursachen:

  * die Elementtabelle bildete alles ausserhalb ihrer 16 Eintraege auf **0** ab (Silber: hier
    0, im Auge 29) -> anderer Flankenrang -> anderes Bezugsatom -> anderes Vorzeichen;
  * `_atrop_dihedral` bildete den ersten Vektor umgekehrt -> der Dieder war um **180 Grad**
    verschoben.  Der gefaltete Betrag bleibt dabei gleich, deshalb fiel es nie auf.

Seit 23.08. ist die Achsendefinition Stufe fuer Stufe die des Auges (Adjazenz, Ringe,
sp2-Test, Einfachbindungsband, Flankenrang, Dieder, Verdrillungsband 20-88 statt 10-80).
`atrop_schluessel_vergleich.py` ist der Waechter dagegen, dass sie wieder auseinanderlaufen.

DIE OPERATION IST EXAKT, KEINE OPTIMIERUNG.  Das Spiegelatropisomer hat denselben BETRAG der
Verdrillung und das umgekehrte VORZEICHEN.  Eine Drehung der einen Seite um die Achse um
`-2*theta` bildet `theta -> -theta` ab: Vorzeichen gekippt, Betrag erhalten.  Eine starre
Drehung um die BINDUNGSACHSE aendert ausschliesslich die Torsion -- alle Bindungslaengen und
alle Bindungswinkel bleiben exakt gleich.  Deshalb braucht dieser Schritt KEINE Relaxation und
kann per Konstruktion keine Bindung BRECHEN.

⚠ ER KANN ABER EINE BINDUNG ERZEUGEN, und das war bis 23.08. unbemerkt.  Der Kollisionsboden
lag bei 1,70 A, die organische Bindungsschwelle des Bauers bei r_i+r_j+0,25 -- fuer C-C
**1,77 A**.  Ein gedrehtes Atom durfte also auf 1,71 A heranruecken, bestand die
Kollisionspruefung und bekam eine Scheinbindung; damit aenderte sich das Flankenprofil und
der Achsenschluessel.  Statt eine Zahl gegen eine andere zu setzen, prueft
`_atrop_topologie_gleich` seither die EIGENSCHAFT selbst: bleibt die Bindungstopologie
unveraendert -- in BEIDEN Adjazenzen, der des Bauers und der des Auges?  Sonst kein Frame.

ADDITIV UND NIE-SCHLECHTER.  Originale bleiben unveraendert; es werden nur FEHLENDE Vorzeichen
angehaengt.  Vorgabe AUS -> byte-identisch.  Gedeckelt und beim Deckel PROTOKOLLIERT (keine
stille Kuerzung).

⚠ WAS DIESES MODUL NICHT TUT: es beurteilt nicht, WELCHES Vorzeichen richtig ist.  Das waere
eine Kristallfrage, und die Konstruktion muss ohne Kristall auskommen.  Es stellt sicher, dass
BEIDE im Manifold stehen -- Vollstaendigkeitsregel, keine Auswahl.

⚠ ABGRENZUNG (geprueft 16.08., bevor gebaut wurde): `_chirality_enumerator.py` macht die
Lambda/Delta-HELIZITAET des Koordinationspolyeders aus Chelatpaaren -- metallzentrierte
Chiralitaet.  Seine "axial"-Stellen meinen axiale POSITIONEN (axial gegen aequatorial), nicht
axiale Chiralitaet.  Es gibt im Bauer keine zweite Atropisomer-Enumeration.

⚠ NAMEN: die privaten Helfer tragen alle das Praefix `_atrop_`, und der Einstieg heisst
`expand_atropisomers` statt `expand_results`.  Der Wachhund `check_exists_first` hatte die
generischen Namen (`_env_int`, `_is_enabled`, `_dihedral`, `expand_results`) als Kollision
gemeldet; statt ihn zu uebergehen sind sie eindeutig gemacht.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Set

import numpy as np

# Dieselben Primitive wie `_stereocenter_enum` -- eine Quelle, keine Drift.
from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

_LOG = logging.getLogger(__name__)

# ===== DIE ACHSENDEFINITION IST DIE DES AUGES -- ZEICHENGENAU (23.08.2026) ==================
#
# GEMESSEN, nicht vermutet (`harness/atrop_schluessel_vergleich.py` auf 120 Systemen,
# 3350 Frames): auf den 103 Bindungen, die BEIDE Seiten fuer eine Achse hielten, war die
# Signatur in **103 von 103 Faellen verschieden** und das VORZEICHEN in **87 von 103**.
# Der Bauer fand 112 Achsen, das Auge 2090.  Damit war jede Spiegelung dieses Moduls ein
# Frame im FALSCHEN Eimer -- und ein Frame im falschen Eimer kostet genau so viel wie ein
# fehlendes.  Das erklaert `atrop44` (16.08., 9x false in BEIDEN Armen) und `atrop10k`
# (21.08., 35 -> 35 bei +316 Frames) ohne Rest.
#
# DAS EINE BEISPIEL, DAS ES ZEIGT (ABOZIB C51-C52, dieselbe Bindung, derselbe Frame):
#     Bauer: ('C', ((6, 6, 6), 0), ...)  und  ('C', ((8, 0), 0), ...)
#     Auge : ('C', ((29, 8), 0), ...)    und  ('C', ((6, 6, 6), 0), ...)
# Die `0` ist SILBER.  Die alte Tabelle kannte 16 Elemente und bildete alles uebrige auf 0
# ab; das Auge faellt auf `round(Kovalenzradius * 20)` zurueck.  Damit kippt der Flankenrang,
# damit das Bezugsatom des Diederwinkels, damit das Vorzeichen.  "Bewusst KLEIN gehalten,
# es geht um eine ORDNUNG" war genau der Denkfehler: eine Ordnung ist nur dann dieselbe,
# wenn beide Seiten dieselbe Ordnung benutzen.
#
# KEIN IMPORT AUS DEM AUGE -- die Konstruktion darf nicht vom Messgeraet abhaengen
# ("Auge = OBERGRENZE der Konstruktion").  Darum steht die Definition hier ZWEITES MAL,
# absichtlich, und `harness/atrop_schluessel_vergleich.py` ist der Waechter dagegen, dass
# die beiden wieder auseinanderlaufen: er muss 0 % SIG_VERSCHIEDEN melden.
#
# Herkunft jeder einzelnen Zahl: `weddell/detectors/atropisomer_sign.py`, Zeilen 101-131.

# Verdrillungsband: unterhalb planar (kein Vorzeichen), oberhalb senkrecht (die beiden Haende
# werden ununterscheidbar).  FRUEHER 10/80 -- unten zu weit (das Auge schreibt unter 20 Grad
# nichts gut, jedes solche Frame war Ausschuss), oben zu eng (`BIRVUW` traegt Achsen bei
# 83,7 und 85,2 Grad, die das Auge VERLANGT und der Bauer nicht sehen konnte).
_ATROP_TWIST_MIN_DEG = 20.0
_ATROP_TWIST_MAX_DEG = 88.0

# Ordnungszahlen fuer den CIP-artigen Flankenrang -- identisch mit `_Z` des Auges.
_ATROP_Z = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "As": 33, "Se": 34, "Br": 35, "Te": 52, "I": 53,
}

# Kovalenzradien -- identisch mit `_COV_R` des Auges.  Sie tragen DREI Rollen: den Rueckfall
# der Ordnungszahl, die Adjazenz und die Achsenlaenge.  Eine Abweichung hier verschiebt alle
# drei auf einmal.
_ATROP_COV_R = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03,
    "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

_ATROP_SP2_ELEMS = {"C", "N", "O", "S", "Se", "B", "P"}   # sp2-faehige Ring-/Zentrumselemente
_ATROP_ARO_BOND_MAX = 1.46        # mittlere Ringbindung darunter = aromatisch/konjugiert
_ATROP_RING_OOP_MAX = 0.35        # RMS-Abweichung aus der Ringebene (A) darunter = planar
_ATROP_SP2_SUM_MIN = 348.0        # Summe der drei Nachbarwinkel darueber = planar sp2
# Achsen-Einfachbindungsband: oberhalb der Doppelbindungsschulter (ein C=C/C=N gehoert zu
# `ez_stereo`, nicht hierher) und unterhalb Marge * Summe der Kovalenzradien.
_ATROP_DOUBLE_CUT = {("C", "C"): 1.42, ("C", "N"): 1.37, ("N", "N"): 1.37, ("C", "O"): 1.36}
_ATROP_LEN_MARGIN = 1.18
_ATROP_ADJ_CUT_FRAC = 1.30        # Adjazenz des Auges: d <= 1,30 * (r_i + r_j)


def _atrop_znum(sym: str) -> int:
    """Identisch mit `_znum` des Auges -- Tabelle, sonst Kovalenzradius * 20."""
    return _ATROP_Z.get(sym, int(round(_ATROP_COV_R.get(sym, 0.9) * 20)))


def _atrop_env_int(name: str, default: int) -> int:
    try:
        return int(str(os.environ.get(name, default)).strip())
    except Exception:
        return default


def _atrop_env_float(name: str, default: float) -> float:
    try:
        return float(str(os.environ.get(name, default)).strip())
    except Exception:
        return default


def _atrop_enabled() -> bool:
    """DIE EINE Lesestelle.  Vorgabe AUS -> byte-identisch."""
    return (os.environ.get("DELFIN_FFFREE_ATROPISOMER_ENUM", "0") == "1"
            or os.environ.get("DELFIN_ATROPISOMER_ENUM", "0") == "1")


# ---------------------------------------------------------------------------------------------
# Achsenerkennung -- geometrisch, ohne Elementtabelle und ohne das Auge
# ---------------------------------------------------------------------------------------------

def _atrop_heavy_nbrs(syms: Sequence[str], nbrs: List[List[int]], i: int) -> List[int]:
    return [j for j in nbrs[i] if syms[j] != "H"]


def _atrop_adjacency(syms: Sequence[str], pts: np.ndarray,
                     cut_frac: float = _ATROP_ADJ_CUT_FRAC) -> List[List[int]]:
    """Adjazenz DES AUGES: d <= cut_frac * (r_i + r_j), schwer UND H.

    ⚠ NICHT `_build_geometric_adjacency` (r_i+r_j+0,25 organisch).  Fuer C-C sind das
    1,98 gegen 1,77 A -- verschiedene Nachbarschaften, verschiedene Flankenprofile,
    verschiedene Signaturen.  Die Achsenerkennung MUSS die Nachbarn sehen, die das Auge
    sieht; die Kollisionspruefung der Spiegelung benutzt weiterhin die engere Bauadjazenz.
    """
    n = len(syms)
    r = np.array([_ATROP_COV_R.get(s, 1.5) for s in syms])
    P = np.asarray(pts, float)
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    cut = cut_frac * (r[:, None] + r[None, :])
    np.fill_diagonal(D, 1e9)
    nbr: List[List[int]] = [[] for _ in range(n)]
    ii, jj = np.where(D <= cut)
    for a, b in zip(ii.tolist(), jj.tolist()):
        if a < b:
            nbr[a].append(b)
            nbr[b].append(a)
    return nbr


def _atrop_aromatic_rings(nbrs: List[List[int]], syms: Sequence[str],
                          pts: np.ndarray) -> List[frozenset]:
    """Aromatisch-artige sp2-Ringe (Groesse 5-6, sp2-faehige Elemente, kurze Ringbindungen,
    planar) -- identisch mit `_aromatic_rings` des Auges.

    ⚠ ERSETZT `_atrop_ring_of` (kleinster Ring aus einer Breitensuche, bis Groesse 8).  Der
    Ring bestimmt, welche Atome beim Flankenprofil AUSGESCHLOSSEN werden; ein anderer Ring
    heisst ein anderes Profil.  Ein Cyclohexanring war fuer die alte Fassung ein Ring und
    fuer das Auge keiner.
    """
    P = np.asarray(pts, float)
    adjA = {i: [k for k in nbrs[i] if syms[k] in _ATROP_SP2_ELEMS]
            for i in range(len(syms)) if syms[i] in _ATROP_SP2_ELEMS}
    found = set()
    for start in adjA:
        stack = [(start, (start,))]
        while stack:
            node, path = stack.pop()
            for nb in adjA.get(node, ()):
                if nb == start and 5 <= len(path) <= 6:
                    found.add(frozenset(path))
                elif nb not in path and len(path) < 6:
                    stack.append((nb, path + (nb,)))
    rings = []
    for r in found:
        idx = list(r)
        bl = [float(np.linalg.norm(P[a] - P[b])) for ai, a in enumerate(idx)
              for b in idx[ai + 1:] if b in nbrs[a]]
        if not bl or (sum(bl) / len(bl)) > _ATROP_ARO_BOND_MAX:
            continue
        Q = P[idx]
        c = Q.mean(0)
        try:
            _u, _s, vt = np.linalg.svd(Q - c)
        except np.linalg.LinAlgError:
            continue
        oop = float(np.sqrt(np.mean(np.dot(Q - c, vt[2]) ** 2)))
        if oop > _ATROP_RING_OOP_MAX:
            continue
        rings.append(r)
    return rings


def _atrop_sum_angles(pts: np.ndarray, c: int, neigh: List[int]) -> float:
    if len(neigh) < 3:
        return 0.0
    s = 0.0
    for a in range(len(neigh)):
        for b in range(a + 1, len(neigh)):
            v1 = pts[neigh[a]] - pts[c]
            v2 = pts[neigh[b]] - pts[c]
            nn = float(np.linalg.norm(v1)) * float(np.linalg.norm(v2))
            if nn < 1e-9:
                continue
            s += math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2)) / nn))))
    return s


def _atrop_is_sp2_center(i: int, nbrs: List[List[int]], syms: Sequence[str],
                         pts: np.ndarray, ring_of: List[List[frozenset]]) -> bool:
    """sp2-Einheit -- identisch mit `_is_sp2_center` des Auges: Ringmitglied ODER
    dreifach koordiniert und planar ODER mit einem Nachbarn auf Doppelbindungsabstand.

    ⚠ DER LETZTE ZWEIG ERSETZT `_atrop_is_mesomeric_amide`.  Ein Amid-N mit C-N <= 1,37 A
    faellt dort von selbst hinein -- ohne Sonderfall, und vor allem: GENAU DANN, wenn das
    Auge es auch tut.  Der Sonderfall hat dieselbe Chemie mit einer anderen Grenze gepruft.
    """
    if syms[i] == "H" or _is_metal_sym(syms[i]) or syms[i] not in _ATROP_SP2_ELEMS:
        return False
    if ring_of[i]:
        return True
    heavy = [k for k in nbrs[i] if syms[k] != "H" and not _is_metal_sym(syms[k])]
    if len(heavy) == 3 and _atrop_sum_angles(pts, i, heavy) >= _ATROP_SP2_SUM_MIN:
        return True
    for k in heavy:
        cut = (_ATROP_DOUBLE_CUT.get((syms[i], syms[k]))
               or _ATROP_DOUBLE_CUT.get((syms[k], syms[i])))
        d = float(np.linalg.norm(pts[i] - pts[k]))
        rr = _ATROP_COV_R.get(syms[i], 1.5) + _ATROP_COV_R.get(syms[k], 1.5)
        if (cut is not None and d <= cut) or d <= 0.93 * rr:
            return True
    return False


def _atrop_axis_is_single(i: int, j: int, syms: Sequence[str], pts: np.ndarray) -> bool:
    """Echte Einfachbindung zwischen zwei Einheiten -- identisch mit `_axis_is_single`."""
    d = float(np.linalg.norm(pts[i] - pts[j]))
    lo = (_ATROP_DOUBLE_CUT.get((syms[i], syms[j]))
          or _ATROP_DOUBLE_CUT.get((syms[j], syms[i])))
    rr = _ATROP_COV_R.get(syms[i], 1.5) + _ATROP_COV_R.get(syms[j], 1.5)
    if lo is not None and d <= lo:
        return False                 # klare Doppelbindung (E/Z-Bereich) -- keine Atropachse
    return d <= _ATROP_LEN_MARGIN * rr


def _atrop_share_ring(i: int, j: int, ring_of: List[List[frozenset]]) -> bool:
    return any(r in ring_of[j] for r in ring_of[i])


def _atrop_side_atoms(syms: Sequence[str], nbrs: List[List[int]],
                      i: int, j: int) -> Optional[Set[int]]:
    """Alle Atome auf der j-Seite der Bindung i-j (ohne i), Metalle als Grenze.

    None, wenn die Bindung IM ORGANISCHEN RING liegt -- dann gibt es keine zwei Seiten, und
    eine Drehung wuerde die Struktur zerreissen.

    ⚠ METALLE WERDEN NICHT DURCHLAUFEN (16.08.2026, aus dem Selbsttest).  Auf `BEBGUL` fand
    dieses Modul KEINE Achse, waehrend das Auge dort `C26-C27` mit 81.2 Grad fuehrt.  Grund:
    bei einem CHELATISIERENDEN Biaryl koordinieren BEIDE Arylringe -- es gibt also einen Weg
    von der einen Seite zur anderen, ueber das METALL.  Das Absuchen hielt die Achse damit
    fuer eine Ringbindung.  Das ist ein METALLACYCLUS, kein organischer Ring; das Auge sieht
    ihn nicht, weil es RDKits Ringerkennung auf dem organischen Graphen benutzt.

    ⚠⚠ DARAUS FOLGT EINE BAUGRENZE, DIE NICHT WEGZUPATCHEN IST.  Haengen BEIDE Seiten am
    Metall, laesst sich die Gegenhaendigkeit NICHT durch starre Drehung erzeugen -- sie
    wuerde die M-D-Bindung zerreissen.  Genau die Eigenschaft, auf die sich dieses Modul
    beruft ("eine Drehung um die Bindungsachse aendert ausschliesslich die Torsion"), gilt
    dort nicht mehr.  `_atrop_find_axes` lehnt solche Achsen AUSDRUECKLICH ab, statt sie
    falsch zu bauen: ihre Gegenhaendigkeit braucht eine NEUSETZUNG, keine Drehung.
    """
    seen = {j}
    stack = [k for k in nbrs[j] if k != i]
    while stack:
        cur = stack.pop()
        if cur == i:
            return None                       # organischer Ringschluss -> keine Drehachse
        if cur in seen:
            continue
        seen.add(cur)
        if _is_metal_sym(syms[cur]):
            continue                          # Metall: Grenze der organischen Seite
        stack.extend(k for k in nbrs[cur] if k not in seen)
    return seen


def _atrop_dihedral(pts: np.ndarray, a: int, i: int, j: int, b: int) -> Optional[float]:
    """Vorzeichenbehafteter Dieder a-i-j-b in (-180, 180] -- identisch mit `_signed_dihedral`
    des Auges.

    🔴 DIE ALTE FASSUNG WAR UM 180 GRAD VERSCHOBEN.  Sie bildete den ersten Vektor als
    `pts[a] - pts[i]`, das Auge als `P[i] - P[a]` -- umgekehrtes Vorzeichen, damit `n1` und
    `m` umgekehrt, damit `atan2(-y, -x) = atan2(y, x) +- pi`.  Der GEFALTETE Betrag bleibt
    dabei gleich (deshalb fiel es nie auf: 1 von 103 Bindungen wich in der Verdrillung ab),
    aber das VORZEICHEN kippt: +60 des Auges war hier -120.  Gemessen: 87 von 103
    gemeinsamen Bindungen trugen entgegengesetzte Haendigkeit.  Ein Enumerator, der das
    fehlende Vorzeichen ergaenzt, ergaenzte damit systematisch das VORHANDENE.
    """
    b1 = pts[i] - pts[a]
    b2 = pts[j] - pts[i]
    b3 = pts[b] - pts[j]
    nb2 = float(np.linalg.norm(b2))
    if nb2 < 1e-9:
        return None
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    if float(np.linalg.norm(n1)) < 1e-9 or float(np.linalg.norm(n2)) < 1e-9:
        return None
    m = np.cross(n1, b2 / nb2)
    return math.degrees(math.atan2(float(np.dot(m, n2)), float(np.dot(n1, n2))))


def _atrop_fold(dih: float) -> float:
    """|Abweichung von planar| in [0, 90] -- identisch mit `_fold` des Auges."""
    a = abs(dih)
    return min(a, 180.0 - a)


# ===== FLANKEN STATT GANZER SEITE (16.08.2026, nach dem `atrop44`-Verdikt) ==================
#
# WARUM UMGESTELLT.  `atrop44` lief additiv und schadensfrei -- `never_worse_ok = true`, alle
# elf `ccdc_*_lost = 0`, 15 Frames auf 9 Systemen ergaenzt -- und die Achse
# `ccdc_atropisomer_realized` blieb trotzdem bei **9x false in BEIDEN Armen**.  Das ist der
# vorregistrierte H0-Fall: ich ergaenze die Gegenhaendigkeit MEINER Achsen, das Auge prueft
# SEINE, und die beiden Signaturen decken sich nicht.
#
# Die alte Fassung verschluesselte auf der GANZEN Seite hinter der Achse (Elementzahl +
# sortierte Elemente).  Das Auge verschluesselt auf den ZWEI FLANKEN am Achsenatom -- und das
# ist die Groesse, die die Haendigkeit ueberhaupt definiert: zwei Achsen mit gleichen Flanken,
# aber verschieden grossen Resten dahinter, sind dieselbe stereogene Situation.
#
# ⚠ UND DER BEZUGSPUNKT DES DIEDERWINKELS MUSS MIT.  Die alte Fassung nahm den schweren
# Nachbarn mit dem KLEINSTEN INDEX.  Das Auge misst ueber die HOCHRANGIGE Flanke.  Mit einem
# anderen Bezugsatom kann derselbe Frame das umgekehrte Vorzeichen tragen -- mein "P" waere
# dann das "M" des Auges, und die Enumeration ergaenzte die bereits vorhandene Haendigkeit.
# Das ist der subtilste Teil der Umstellung und faellt in keinem Selbsttest auf, der nur die
# eigene Ausgabe vergleicht.
#
# KEIN IMPORT AUS DEM AUGE.  Die Konstruktion darf nicht vom Messgeraet abhaengen ("Auge =
# OBERGRENZE der Konstruktion"), darum ist die Ringerkennung hier aus der Adjazenz gebaut und
# nicht aus RDKit uebernommen.


def _atrop_sub_profile(o: int, exclude: frozenset,
                       nbrs: List[List[int]], syms: Sequence[str]) -> tuple:
    """CIP-artiger Rang des Substituenten an Flankenatom `o`, von der Achse WEG gesehen --
    identisch mit `_sub_profile` des Auges.

    Schluessel = (absteigend sortierte schwere Ordnungszahlen bis zwei Bindungen weit, -nH).
    Ohne Atomindex -- ein rein chemischer Rang, ueber Frames hinweg vergleichbar.

    ⚠ `_atrop_znum` statt `_ATROP_Z.get(..., 0)`: ein Element ausserhalb der Tabelle bekommt
    den Kovalenzradius-Rueckfall des Auges, nicht die Null.  Genau daran hing ABOZIB
    (Silber: hier 0, dort 29) -- und damit der Flankenrang, das Bezugsatom und das Vorzeichen.
    """
    heavy_z: List[int] = []
    n_h = 0
    seen = {o}
    shell1 = [k for k in nbrs[o] if k not in exclude]
    for k in shell1:
        if syms[k] == "H":
            n_h += 1
        else:
            heavy_z.append(_atrop_znum(syms[k]))
            seen.add(k)
    for k in [x for x in seen if x != o]:
        for q in nbrs[k]:
            if q in exclude or q in seen or q == o:
                continue
            if syms[q] == "H":
                n_h += 1
            else:
                heavy_z.append(_atrop_znum(syms[q]))
    return (tuple(sorted(heavy_z, reverse=True)), -n_h)


def _atrop_flanks(i: int, partner: int, nbrs: List[List[int]],
                  syms: Sequence[str], ring_of: List[List[frozenset]]):
    """(flankeHoch, flankeNiedrig, schluesselHoch, schluesselNiedrig) oder None --
    identisch mit `_flanks` des Auges.

    Ringatom -> die zwei ORTHO-Ringnachbarn des ERSTEN Rings; sonst die zwei schweren
    Nicht-Partner.  None, wenn es nicht genau zwei UNTERSCHEIDBARE Flanken gibt -- dann ist
    die Achse nicht stereogen (lokale Spiegelebene) und darf gar nicht enumeriert werden.
    """
    rings_i = ring_of[i]
    if rings_i:
        ring = rings_i[0]
        flanks = [k for k in nbrs[i] if k in ring and k != partner]
        excl_base = frozenset(ring)
    else:
        flanks = [k for k in nbrs[i] if k != partner and syms[k] != "H"]
        excl_base = frozenset({i})
    if len(flanks) != 2:
        return None
    f1, f2 = flanks
    ex = excl_base | {i, partner}
    k1 = _atrop_sub_profile(f1, ex, nbrs, syms)
    k2 = _atrop_sub_profile(f2, ex, nbrs, syms)
    if k1 == (tuple(), 0) and k2 == (tuple(), 0):
        return None                      # kein Substituent -> keine Rotationsbarriere
    if k1 == k2:
        return None                      # ununterscheidbar -> achiral
    return (f1, f2, k1, k2) if k1 > k2 else (f2, f1, k2, k1)


def _atrop_axis_sig_flanks(syms: Sequence[str], i: int, j: int,
                           kHiA, kLoA, kHiB, kLoB) -> tuple:
    """Kanonische Achsensignatur aus den FLANKEN -- dieselbe Groesse, die das Auge fuehrt.
    Die Seiten sind vertauschbar, darum sortiert; das Vorzeichen ist unter Seitentausch
    invariant, weil dihedral(a,i,j,b) == dihedral(b,j,i,a)."""
    sideA = (syms[i], kHiA, kLoA)
    sideB = (syms[j], kHiB, kLoB)
    return tuple(sorted((repr(sideA), repr(sideB))))


def _atrop_find_axes(syms: Sequence[str], pts: np.ndarray,
                     nbrs: Optional[List[List[int]]] = None) -> List[dict]:
    """Alle stereogenen Achsen eines Frames -- Stufe fuer Stufe die des Auges (`_axes`),
    zuzueglich der EINEN Groesse, die das Auge nicht braucht und der Bauer nicht entbehren
    kann: `side`, die Atome, die gedreht werden.

    ⚠ `nbrs` wird ENTGEGENGENOMMEN UND VERWORFEN.  Die Aufrufer reichen die Bauadjazenz
    durch (r_i+r_j+0,25); die Achsenerkennung braucht die des Auges (1,30*(r_i+r_j)).  Der
    Parameter bleibt nur, damit kein Aufrufer stillschweigend etwas anderes bekommt als er
    denkt -- er ist ausdruecklich unbenutzt.

    ⚠⚠ DER BAUER FINDET WENIGER ALS DAS AUGE, UND DAS IST KEIN FEHLER: haengen beide Seiten
    einer Achse am Metall, gibt es keine Drehachse (`_atrop_side_atoms` -> None).  Ihre
    Gegenhaendigkeit braucht eine NEUSETZUNG.  Das wird GEZAEHLT und gemeldet, nicht
    verschwiegen -- eine stille Kuerzung sieht von aussen aus wie Vollstaendigkeit.
    """
    del nbrs                                  # siehe Docstring: bewusst verworfen
    P = np.asarray(pts, float)
    n = len(syms)
    if n < 6:
        return []
    nbr = _atrop_adjacency(syms, P)
    rings = _atrop_aromatic_rings(nbr, syms, P)
    ring_of: List[List[frozenset]] = [[] for _ in range(n)]
    for r in rings:
        for a in r:
            ring_of[a].append(r)
    out: List[dict] = []
    n_metallbruecke = 0
    seen_bond = set()
    for i in range(n):
        if syms[i] == "H" or _is_metal_sym(syms[i]):
            continue
        for j in nbr[i]:
            if j <= i or syms[j] == "H" or _is_metal_sym(syms[j]):
                continue
            if (i, j) in seen_bond:
                continue
            seen_bond.add((i, j))
            # Eine Achse ist eine EINFACHBINDUNG ZWISCHEN ZWEI sp2-EINHEITEN.
            if not (_atrop_is_sp2_center(i, nbr, syms, P, ring_of)
                    and _atrop_is_sp2_center(j, nbr, syms, P, ring_of)):
                continue
            if ring_of[i] and ring_of[j] and _atrop_share_ring(i, j, ring_of):
                continue                  # kondensiert / selber Ring -> keine Achse zwischen Einheiten
            if not (ring_of[i] or ring_of[j]):
                continue                  # mindestens eine Aryleinheit noetig
            if not _atrop_axis_is_single(i, j, syms, P):
                continue
            fa = _atrop_flanks(i, j, nbr, syms, ring_of)
            fb = _atrop_flanks(j, i, nbr, syms, ring_of)
            if fa is None or fb is None:
                continue                  # keine zwei unterscheidbaren Flanken -> nicht stereogen
            fHiA, _fLoA, kHiA, kLoA = fa
            fHiB, _fLoB, kHiB, kLoB = fb
            # ⚠ Der Dieder geht ueber die HOCHRANGIGEN Flanken -- dieselbe Wahl wie im Auge.
            dih = _atrop_dihedral(P, fHiA, i, j, fHiB)
            if dih is None:
                continue
            fold = _atrop_fold(dih)
            if not (_ATROP_TWIST_MIN_DEG <= fold <= _ATROP_TWIST_MAX_DEG):
                continue                  # planar oder senkrecht -> Vorzeichen ist Rauschen
            # AB HIER ist die Achse fuer das Auge eine Achse.  Erst jetzt fragt der Bauer,
            # ob er sie ueberhaupt DREHEN kann -- die Reihenfolge ist wichtig, sonst
            # verschwindet eine Baugrenze in der Achsenerkennung.
            #
            # ⚠ NICHT DREHBAR HEISST NICHT UNBAUBAR (23.08.2026).  Frueher fielen diese
            # Achsen hier ganz heraus, und damit war der Manifold auf ihnen still
            # unvollstaendig -- gemessen 89 Stueck auf 120 Systemen.  Sie bleiben jetzt in
            # der Liste, mit `rotatable=False`: die Gegenhaendigkeit braucht dort eine
            # NEUSETZUNG (Gesamtspiegelung), keine Drehung.  Wer sie herausfiltert, misst
            # seine eigene Werkzeuggrenze und haelt sie fuer Chemie.
            side = _atrop_side_atoms(syms, nbr, i, j)      # None = Ring / Metallbruecke
            _rot = bool(side) and len(side) >= 2
            if not _rot:
                n_metallbruecke += 1
            out.append({
                "i": i, "j": j, "dih": dih, "fold": fold,
                "sign": "P" if dih > 0 else "M",
                "side": side if _rot else None,
                "rotatable": _rot,
                "sig": _atrop_axis_sig_flanks(syms, i, j, kHiA, kLoA, kHiB, kLoB),
            })
    if n_metallbruecke:
        _LOG.debug("atropisomer-enum: %d Achse(n) sind fuer das Auge Achsen, aber nicht "
                   "drehbar (Ring oder Metallbruecke) -- sie brauchen eine Neusetzung",
                   n_metallbruecke)
    return sorted(out, key=lambda d: (d["i"], d["j"]))


# ---------------------------------------------------------------------------------------------
# Die Spiegelung
# ---------------------------------------------------------------------------------------------

def _atrop_rotate_side(pts: np.ndarray, ax: dict, angle_deg: float) -> np.ndarray:
    """Starre Drehung der j-Seite um die Achse i->j (Rodrigues).

    Aendert ausschliesslich die Torsion: Bindungslaengen und Bindungswinkel bleiben exakt.
    """
    p = pts.copy()
    o = pts[ax["i"]]
    k = pts[ax["j"]] - o
    nk = float(np.linalg.norm(k))
    if nk < 1e-9:
        return p
    k = k / nk
    t = math.radians(angle_deg)
    ct, st = math.cos(t), math.sin(t)
    idx = sorted(ax["side"])
    v = p[idx] - o
    p[idx] = (v * ct
              + np.cross(np.broadcast_to(k, v.shape), v) * st
              + np.outer(v @ k, k) * (1.0 - ct)) + o
    return p


def _atrop_min_nonbonded(syms: Sequence[str], pts: np.ndarray, nbrs: List[List[int]],
                         moved: Set[int]) -> float:
    """Kleinster Abstand zwischen einem bewegten und einem ruhenden Atom; 1-2 und 1-3 sind
    ausgenommen, weil die sich per Konstruktion nicht aendern."""
    worst = 9.9
    still = [k for k in range(len(syms)) if k not in moved]
    if not still:
        return worst
    for m in sorted(moved):
        excl = set(nbrs[m]) | {m}
        for q in nbrs[m]:
            excl |= set(nbrs[q])
        for s in still:
            if s in excl:
                continue
            d = float(np.linalg.norm(pts[m] - pts[s]))
            if d < worst:
                worst = d
    return worst


def _atrop_topologie_gleich(syms: Sequence[str], alt: np.ndarray, neu: np.ndarray) -> bool:
    """Hat die Drehung die BINDUNGSTOPOLOGIE unveraendert gelassen?

    🔴 DIE ZWEITE WURZEL (23.08.2026).  Der Kollisionsboden allein reicht NICHT: er lag bei
    1,70 A, die organische Bindungsschwelle des Bauers bei r_i+r_j+0,25 -- fuer C-C **1,77 A**,
    fuer C-N 1,72.  Ein gedrehtes Atom durfte also auf 1,71 A heranruecken, bestand die
    Kollisionspruefung und bekam eine SCHEINBINDUNG.  Damit aendert sich `_sub_profile`,
    damit der Achsenschluessel -- und der angehaengte Frame landet in einem NEUEN Eimer,
    statt den zu fuellen, dem das Vorzeichen fehlt.  Genau das zeigt PEHWEH: der ON-Arm hat
    zwei Signaturen fuer dieselbe Bindung, beide weiter mit demselben Vorzeichen.

    Eine Zahl gegen eine andere Zahl zu setzen waere geraten.  Geprueft wird die EIGENSCHAFT,
    auf die sich das Modul beruft: "eine starre Drehung um die Bindungsachse aendert
    ausschliesslich die Torsion".  Aendert sich dabei ein Bindungspartner, ist das nicht mehr
    wahr -- unabhaengig davon, bei welchem Abstand es passiert.  Geprueft wird in BEIDEN
    Adjazenzen: der des Bauers (die spaeter jeder Reparateur sieht) und der des Auges (die
    ueber den Achsenschluessel entscheidet).
    """
    for adj in (lambda P: _build_geometric_adjacency(list(syms), P)[0],
                lambda P: _atrop_adjacency(syms, P)):
        a = [frozenset(x) for x in adj(alt)]
        b = [frozenset(x) for x in adj(neu)]
        if a != b:
            return False
    return True


def _atrop_mirror_frame(xyz: str, ax: dict, clash_min: float) -> Optional[str]:
    """Ein Frame mit gekipptem Vorzeichen dieser EINEN Achse.

    None, wenn es kollidiert ODER wenn die Drehung die Bindungstopologie veraendert --
    ein Frame im falschen Signatur-Eimer kostet genau so viel wie ein fehlendes.
    """
    syms, pts, lines = _parse_xyz(xyz)
    nbrs, _bd = _build_geometric_adjacency(syms, pts)
    new = _atrop_rotate_side(pts, ax, -2.0 * ax["dih"])
    if _atrop_min_nonbonded(syms, new, nbrs, ax["side"]) < clash_min:
        return None
    if not _atrop_topologie_gleich(syms, pts, new):
        _LOG.debug("atropisomer-enum: Drehung um %d-%d aendert die Bindungstopologie "
                   "-- Frame verworfen", ax["i"], ax["j"])
        return None
    # Reihenfolge beachten: _format_xyz(orig_lines, syms, positions) -- siehe
    # _coord_angle_corrector.py:135.  Vertauscht liefert es stillen Unsinn.
    return _format_xyz(lines, syms, new)


def _atrop_reseat_frame(xyz: str) -> Optional[str]:
    """NEUSETZUNG statt Drehung: die Gesamtspiegelung des Frames.

    WOZU.  Haengen BEIDE Seiten einer Achse am Metall -- ein chelatisierendes Biaryl --,
    laesst sich die Gegenhaendigkeit NICHT durch starre Drehung erzeugen: sie wuerde die
    M-D-Bindung zerreissen.  Gemessen am 23.08. auf 120 Systemen: **89 Achsen**, die das
    Auge fuehrt und der Bauer nicht drehen kann.  Bisher fielen sie still heraus.

    Die Gesamtspiegelung ist die richtige Operation, und zwar aus einem Grund, nicht aus
    Bequemlichkeit: sie ist eine ISOMETRIE.  Jeder Abstand bleibt exakt erhalten -- alle
    Bindungslaengen, alle Winkel, auch jede M-D-Bindung.  Sie kann daher per Konstruktion
    nichts zerreissen, und sie kippt das Vorzeichen JEDER stereogenen Achse auf einmal.

    ⚠ SIE KIPPT AUCH ALLES ANDERE: Lambda/Delta am Metall, Rueckgrat-Zentren.  Das ist
    kein Nebenschaden, sondern ein ANDERES Isomer -- und additiv angehaengt heisst das
    Vollstaendigkeit, nicht Ersatz.  Das Original bleibt unberuehrt.

    KEINE ZWEITE SPIEGELUNG.  `_mirror_enum.mirror_frame` macht genau das, samt
    Achiralitaetspruefung (ist der Frame sein eigenes Spiegelbild, gibt es nichts zu
    gewinnen und die Funktion gibt None).  Hier wird sie benutzt, nicht nachgebaut.
    """
    try:
        from delfin.manta._mirror_enum import mirror_frame as _mf
    except Exception as exc:                       # pragma: no cover
        _LOG.debug("atropisomer-enum: Spiegelung nicht verfuegbar: %s", exc)
        return None
    try:
        return _mf(xyz)
    except Exception as exc:
        _LOG.debug("atropisomer-enum: Neusetzung fehlgeschlagen: %s", exc)
        return None


def _atrop_realized(xyz: str, sig: tuple, want: str) -> bool:
    """Traegt dieser Frame die Signatur `sig` mit dem Vorzeichen `want` -- WIRKLICH?

    🔑 DIE FRAGE, DIE DEN GANZEN 23.08. GEFEHLT HAT.  `atrop44` und `atrop10k` haben Frames
    angehaengt und die Achse blieb stehen; niemand hat je nachgesehen, ob das angehaengte
    Frame das Vorzeichen ueberhaupt traegt, das es tragen sollte.  Bei PEHWEH tat es das
    nicht -- es landete in einem NEUEN Signatur-Eimer, mit demselben Vorzeichen.

    Ein Frame, das seinen Zweck nicht erfuellt, ist kein Gewinn: es kostet Bau- und
    Augenzeit, hebt Maxima und traegt nichts bei.  Darum wird hier nach der Operation
    GEPRUEFT, nicht vorher gehofft.
    """
    try:
        A = _atrop_analyze(xyz)
    except Exception:
        return False
    if not A:
        return False
    return any(ax["sig"] == sig and ax["sign"] == want for ax in A["axes"])


def _atrop_analyze(xyz: str) -> Optional[dict]:
    syms, pts, _lines = _parse_xyz(xyz)
    if len(syms) < 4:
        return None
    axes = _atrop_find_axes(syms, pts)
    return {"axes": axes} if axes else None


# ---------------------------------------------------------------------------------------------
# Der additive Durchgang
# ---------------------------------------------------------------------------------------------

def expand_atropisomers(results):
    """ADDITIV: haenge fuer jede stereogene Achse das FEHLENDE Vorzeichen an.

    Durchgang 1 sammelt ueber ALLE Frames, welche (Achsensignatur, Vorzeichen) vorkommen --
    genau der Punkt, an dem `_stereocenter_enum` am 10.08. gescheitert ist: es meldete
    Vollstaendigkeit auf einer GROEBEREN Partition, als es vervollstaendigte.  Hier ist der
    Schluessel die ACHSE selbst, also die Groesse, um die es geht.

    Durchgang 2 baut fuer jede Signatur, der ein Vorzeichen fehlt, EIN Frame aus dem
    Vertreter.  Originale bleiben unangetastet -> nie-schlechter per Konstruktion.
    """
    if not results:
        return results
    max_added = _atrop_env_int("DELFIN_ATROPISOMER_MAX_ADDED", 64)
    clash_min = _atrop_env_float("DELFIN_ATROPISOMER_CLASH_MIN", 1.70)

    present: Set[tuple] = set()                 # (sig, sign) bereits im Manifold
    reps: Dict[tuple, tuple] = {}               # sig -> (order, xyz, label, axis)
    for order, entry in enumerate(results):
        try:
            xyz = entry[0]
            lbl = entry[1] if len(entry) > 1 else ""
        except Exception:
            continue
        try:
            A = _atrop_analyze(xyz)
        except Exception as exc:
            _LOG.debug("atropisomer-enum: Frame %d nicht lesbar: %s", order, exc)
            continue
        if not A:
            continue
        for ax in A["axes"]:
            present.add((ax["sig"], ax["sign"]))
            reps.setdefault(ax["sig"], (order, xyz, lbl, ax))

    if not reps:
        return results

    # NEUSETZUNG fuer nicht drehbare Achsen -- eigener Schalter, Vorgabe AUS, damit der
    # Drehpfad unveraendert messbar bleibt und `atrop2` genau das misst, was vorregistriert
    # wurde.  Ein zweiter Mechanismus im selben Lauf waere ein zweite Achse im A/B.
    _reseat_on = _atrop_env_int("DELFIN_ATROPISOMER_RESEAT", 0) == 1
    added: List[tuple] = []
    capped = 0
    n_unverifiziert = 0
    n_nicht_drehbar = 0
    for sig, (order, xyz, lbl, ax) in sorted(reps.items(), key=lambda kv: kv[1][0]):
        want = "M" if ax["sign"] == "P" else "P"
        if (sig, want) in present:
            continue                             # beide Haende schon da
        if len(added) >= max_added:
            capped += 1
            continue
        mx = None
        how = "atrop"
        if ax.get("rotatable", True):
            try:
                mx = _atrop_mirror_frame(xyz, ax, clash_min)
            except Exception as exc:
                _LOG.debug("atropisomer-enum: Drehung fehlgeschlagen (%s): %s", lbl, exc)
                mx = None
        else:
            # Beide Seiten am Metall: eine Drehung wuerde die M-D-Bindung zerreissen.
            # ⚠ KEIN RUECKFALL VON DER DREHUNG AUF DIE SPIEGELUNG.  Eine Drehung, die an
            # Kollision oder Topologie SCHEITERT, ist eine begruendete Absage -- sie mit
            # einer Operation zu ueberdecken, die nebenbei jedes andere Stereoelement
            # kippt, waere Kosmetik.  Die Neusetzung gilt nur, wo gar nicht gedreht
            # werden KANN.
            n_nicht_drehbar += 1
            if _reseat_on:
                mx = _atrop_reseat_frame(xyz)
                how = "atropR"
        if mx is None:
            continue                             # nicht bauen, nichts verlieren
        # ===== VERIFIKATION: TUT DAS FRAME, WOFUER ES GEBAUT WURDE? =====================
        if not _atrop_realized(mx, sig, want):
            n_unverifiziert += 1
            _LOG.debug("atropisomer-enum: %s traegt %s auf %r NICHT -- verworfen",
                       lbl, want, sig)
            continue
        added.append((mx, f"{lbl}_{how}-{want}"))
        # Die Gesamtspiegelung kippt ALLE Achsen auf einmal.  Was sie mitrealisiert, wird
        # als vorhanden vermerkt -- sonst baut die naechste Runde dasselbe Frame noch
        # einmal fuer eine andere Signatur.
        if how == "atropR":
            try:
                for _a in (_atrop_analyze(mx) or {}).get("axes", []):
                    present.add((_a["sig"], _a["sign"]))
            except Exception:
                pass
        present.add((sig, want))
    if n_unverifiziert:
        _LOG.warning("atropisomer-enum: %d Frame(s) gebaut und VERWORFEN, weil sie das "
                     "fehlende Vorzeichen nicht trugen", n_unverifiziert)
    if n_nicht_drehbar and not _reseat_on:
        _LOG.info("atropisomer-enum: %d Achse(n) nicht drehbar (Metallbruecke); "
                  "DELFIN_ATROPISOMER_RESEAT=1 setzt sie per Neusetzung", n_nicht_drehbar)

    if capped:
        # KEINE STILLE KUERZUNG: was der Deckel wegnimmt, wird gemeldet.
        _LOG.warning("atropisomer-enum: Deckel DELFIN_ATROPISOMER_MAX_ADDED=%d erreicht, "
                     "%d Achsen NICHT vervollstaendigt", max_added, capped)
    if not added:
        return results
    _LOG.info("atropisomer-enum: %d Gegenhaendigkeiten ergaenzt (%d Achsen gesamt)",
              len(added), len(reps))

    out = list(results)
    proto = results[0]
    for mx, lab in added:
        if isinstance(proto, tuple):
            out.append((mx, lab))
        elif isinstance(proto, list):
            out.append([mx, lab])
        else:
            out.append(mx)
    return out


def _atrop_why(syms, pts, nbrs, i, j) -> str:
    """An WELCHER Stufe faellt die Bindung i-j durch?  Ein Selbsttest, der nur 'nichts
    gefunden' meldet, ist wertlos -- er muss die Stufe nennen.

    ⚠ Die Stufen stehen hier in DERSELBEN Reihenfolge wie in `_atrop_find_axes`; laufen die
    beiden auseinander, nennt der Selbsttest eine Stufe, an der es gar nicht scheitert.
    `nbrs` wird auch hier verworfen -- geprueft wird gegen die Adjazenz des Auges.
    """
    del nbrs
    P = np.asarray(pts, float)
    if syms[i] == "H" or syms[j] == "H":
        return "H"
    if _is_metal_sym(syms[i]) or _is_metal_sym(syms[j]):
        return "Metall"
    nbr = _atrop_adjacency(syms, P)
    rings = _atrop_aromatic_rings(nbr, syms, P)
    ring_of: List[List[frozenset]] = [[] for _ in range(len(syms))]
    for r in rings:
        for a in r:
            ring_of[a].append(r)
    sp2i = _atrop_is_sp2_center(i, nbr, syms, P, ring_of)
    sp2j = _atrop_is_sp2_center(j, nbr, syms, P, ring_of)
    if not (sp2i and sp2j):
        return f"nicht sp2 (i={sp2i} deg={len(nbr[i])}, j={sp2j} deg={len(nbr[j])})"
    if ring_of[i] and ring_of[j] and _atrop_share_ring(i, j, ring_of):
        return "kondensiert / selber Ring"
    if not (ring_of[i] or ring_of[j]):
        return "keine Seite ist ein aromatischer Ring"
    if not _atrop_axis_is_single(i, j, syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        return f"keine Einfachbindung (d={d:.2f})"
    fa = _atrop_flanks(i, j, nbr, syms, ring_of)
    fb = _atrop_flanks(j, i, nbr, syms, ring_of)
    if fa is None:
        return "Seite i: keine zwei unterscheidbaren Flanken"
    if fb is None:
        return "Seite j: keine zwei unterscheidbaren Flanken"
    d = _atrop_dihedral(P, fa[0], i, j, fb[0])
    if d is None:
        return "Dieder nicht berechenbar"
    f = _atrop_fold(d)
    if not (_ATROP_TWIST_MIN_DEG <= f <= _ATROP_TWIST_MAX_DEG):
        return f"Verdrillung {f:.1f} ausserhalb {_ATROP_TWIST_MIN_DEG}-{_ATROP_TWIST_MAX_DEG}"
    side = _atrop_side_atoms(syms, nbr, i, j)
    if not side:
        return "Bindung liegt im Ring oder ueberbrueckt das Metall -- nicht drehbar"
    if len(side) < 2:
        return "Seite zu klein"
    return "OK"


if __name__ == "__main__":  # pragma: no cover -- Achsenvergleich gegen das Auge
    # WOZU.  `atrop44` hat am 16.08. gemessen, dass die Enumeration Frames ergaenzt und die
    # Achse `ccdc_atropisomer_realized` trotzdem bei null bleibt: meine Achsensignatur war
    # nicht die des Auges.  Dieser Einstieg druckt die gefundenen Achsen in derselben Form
    # wie `weddell/detectors/atropisomer_sign.py`, damit sich beide Listen NEBENEINANDER
    # lesen lassen -- Atompaar, Verdrillung, Vorzeichen.  Decken sie sich nicht, ist ein
    # weiterer A/B-Lauf sinnlos.
    #
    #     python delfin/manta/_atropisomer_enum.py <frame.xyz> [weitere.xyz ...]
    import sys as _sys
    for _p in _sys.argv[1:]:
        try:
            with open(_p) as _fh:
                _txt = _fh.read()
        except Exception as _e:
            print(f"  {_p}: nicht lesbar ({_e})")
            continue
        # Mehrframe-XYZ: nur der erste Frame, wie die Augen-Kommandozeile es auch tut.
        _lines = _txt.splitlines()
        try:
            _n = int(_lines[0].split()[0])
            _one = "\n".join(_lines[:_n + 2])
        except Exception:
            _one = _txt
        _A = _atrop_analyze(_one)
        print(f"\n=== {_p} ===")
        if not _A:
            print("  keine stereogene Achse gefunden -- Ablehnungsgruende je Kandidat:")
            _s, _P, _l = _parse_xyz(_one)
            _nb, _bd = _build_geometric_adjacency(_s, _P)
            _seen_reason = {}
            for _i in range(len(_s)):
                for _j in _nb[_i]:
                    if _j <= _i:
                        continue
                    _r = _atrop_why(_s, _P, _nb, _i, _j)
                    if _r in ("H", "Metall"):
                        continue
                    _seen_reason.setdefault(_r, []).append(f"{_s[_i]}{_i}-{_s[_j]}{_j}")
            for _r, _ex in sorted(_seen_reason.items(), key=lambda kv: -len(kv[1])):
                print(f"    {len(_ex):4d}x  {_r}   z.B. {' '.join(_ex[:3])}")
            continue
        for _ax in _A["axes"]:
            print(f"  {_ax['i']}-{_ax['j']}  twist {_ax['fold']:5.1f}  "
                  f"sign {_ax['sign']}  sig {_ax['sig']}")


