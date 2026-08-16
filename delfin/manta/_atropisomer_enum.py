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
Einfachbindung behandelt, sieht dort gar keine Achse.  Darum ist der Amidfall hier
AUSDRUECKLICH mitaufgenommen (`_atrop_is_mesomeric_amide`), nicht nur der Biaryl-Fall.

DIE OPERATION IST EXAKT, KEINE OPTIMIERUNG.  Das Spiegelatropisomer hat denselben BETRAG der
Verdrillung und das umgekehrte VORZEICHEN.  Eine Drehung der einen Seite um die Achse um
`-2*theta` bildet `theta -> -theta` ab: Vorzeichen gekippt, Betrag erhalten.  Eine starre
Drehung um die BINDUNGSACHSE aendert ausschliesslich die Torsion -- alle Bindungslaengen und
alle Bindungswinkel bleiben exakt gleich.  Deshalb braucht dieser Schritt KEINE Relaxation und
kann per Konstruktion keine Bindung brechen; er kann nur eine Kollision erzeugen, und die wird
geprueft.

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

# Verdrillungsband, in dem eine Achse ueberhaupt stereogen ist.  Unterhalb planar (kein
# Vorzeichen), oberhalb senkrecht (die beiden Haende werden ununterscheidbar).  Bewusst weiter
# gefasst als das Auge: hier wird ENUMERIERT, nicht geurteilt -- ein Frame zu viel kostet
# wenig, ein fehlendes kostet die ganze Achse.
_ATROP_TWIST_MIN_DEG = 10.0
_ATROP_TWIST_MAX_DEG = 80.0

# Ordnungszahlen fuer den CIP-artigen Flankenrang.  Bewusst KLEIN gehalten: es geht um eine
# ORDNUNG zwischen zwei Flanken, nicht um Chemie -- ein unbekanntes Element bekommt 0 und
# sortiert damit nach unten, statt den Rang unbrauchbar zu machen.
_ATROP_Z = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "As": 33, "Se": 34, "Br": 35, "Sb": 51, "Te": 52, "I": 53,
}

_ATROP_SP2_DEG = 3                # sp2: genau drei Nachbarn (schwer + H)
_ATROP_PLANAR_TOL_DEG = 25.0      # Winkelsumme am Zentrum nahe 360 -> planar


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


def _atrop_is_planar_sp2(pts: np.ndarray, nbrs: List[List[int]], c: int) -> bool:
    """Drei Nachbarn, Winkelsumme nahe 360 Grad -> planares sp2-Zentrum."""
    nb = nbrs[c]
    if len(nb) != _ATROP_SP2_DEG:
        return False
    tot = 0.0
    for a in range(3):
        for b in range(a + 1, 3):
            u = pts[nb[a]] - pts[c]
            v = pts[nb[b]] - pts[c]
            nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
            if nu < 1e-9 or nv < 1e-9:
                return False
            cos = max(-1.0, min(1.0, float(np.dot(u, v)) / (nu * nv)))
            tot += math.degrees(math.acos(cos))
    return abs(tot - 360.0) <= _ATROP_PLANAR_TOL_DEG


def _atrop_is_mesomeric_amide(syms: Sequence[str], nbrs: List[List[int]],
                              i: int, j: int) -> bool:
    """MESOMERIE: eine C-N-Bindung, deren C eine Carbonylgruppe traegt.

    Das freie Elektronenpaar am N konjugiert in C=O, die C-N-Bindung bekommt Doppelbindungs-
    anteil -- daraus entstehen Rotationsbarriere und bevorzugte Verdrillung.  Genau das macht
    Aryl-Amide (BIRVUW N29-C30, N18-C19) ueberhaupt erst atropisomer.  Ein rein geometrischer
    sp2-Test sieht dort keine Achse, weil das N pyramidal gebaut sein kann.
    """
    for c, n in ((i, j), (j, i)):
        if syms[c] != "C" or syms[n] != "N":
            continue
        for k in nbrs[c]:
            if k != n and syms[k] == "O" and len(_atrop_heavy_nbrs(syms, nbrs, k)) == 1:
                return True          # terminales O am selben C = Carbonyl
    return False


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
    b0 = pts[a] - pts[i]
    b1 = pts[j] - pts[i]
    b2 = pts[b] - pts[j]
    nb1 = float(np.linalg.norm(b1))
    if nb1 < 1e-9:
        return None
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    if float(np.linalg.norm(n1)) < 1e-9 or float(np.linalg.norm(n2)) < 1e-9:
        return None
    m = np.cross(n1, b1 / nb1)
    return math.degrees(math.atan2(float(np.dot(m, n2)), float(np.dot(n1, n2))))


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


def _atrop_ring_of(nbrs: List[List[int]], i: int, max_size: int = 8) -> Optional[frozenset]:
    """Kleinster Ring, der `i` enthaelt -- oder None.  Aus der Adjazenz, ohne RDKit.

    `i` liegt in einem Ring, wenn zwei seiner Nachbarn ohne Umweg ueber `i` verbunden sind.
    Gesucht wird der KUERZESTE solche Weg (Breitensuche), damit ein Atom in einem
    kondensierten System seinen kleinsten Ring bekommt und nicht den umschliessenden.
    """
    best = None
    nb = list(nbrs[i])
    for a in range(len(nb)):
        for b in range(a + 1, len(nb)):
            src, dst = nb[a], nb[b]
            prev = {src: None}
            queue = [src]
            while queue:
                cur = queue.pop(0)
                if cur == dst:
                    break
                for k in nbrs[cur]:
                    if k == i or k in prev:
                        continue
                    prev[k] = cur
                    queue.append(k)
            if dst not in prev:
                continue
            path, cur = [], dst
            while cur is not None:
                path.append(cur)
                cur = prev[cur]
            ring = frozenset(path) | {i}
            if len(ring) <= max_size and (best is None or len(ring) < len(best)):
                best = ring
    return best


def _atrop_sub_profile(o: int, exclude: frozenset,
                       nbrs: List[List[int]], syms: Sequence[str]) -> tuple:
    """CIP-artiger Rang des Substituenten an Flankenatom `o`, von der Achse WEG gesehen.

    Schluessel = (absteigend sortierte schwere Ordnungszahlen bis zwei Bindungen weit, -nH).
    Ohne Atomindex -- ein rein chemischer Rang, ueber Frames hinweg vergleichbar.
    """
    heavy_z: List[int] = []
    n_h = 0
    seen = {o}
    shell1 = [k for k in nbrs[o] if k not in exclude]
    for k in shell1:
        if syms[k] == "H":
            n_h += 1
        else:
            heavy_z.append(_ATROP_Z.get(syms[k], 0))
            seen.add(k)
    for k in [x for x in seen if x != o]:
        for q in nbrs[k]:
            if q in exclude or q in seen:
                continue
            if syms[q] == "H":
                n_h += 1
            else:
                heavy_z.append(_ATROP_Z.get(syms[q], 0))
    return (tuple(sorted(heavy_z, reverse=True)), -n_h)


def _atrop_flanks(i: int, partner: int, nbrs: List[List[int]],
                  syms: Sequence[str], ring_i: Optional[frozenset]):
    """(flankeHoch, flankeNiedrig, schluesselHoch, schluesselNiedrig) oder None.

    Ringatom -> die zwei ORTHO-Ringnachbarn; sonst die zwei schweren Nicht-Partner.
    None, wenn es nicht genau zwei UNTERSCHEIDBARE Flanken gibt -- dann ist die Achse nicht
    stereogen (lokale Spiegelebene) und darf gar nicht enumeriert werden.
    """
    if ring_i:
        flanks = [k for k in nbrs[i] if k in ring_i and k != partner]
        excl_base = ring_i
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
                     nbrs: List[List[int]]) -> List[dict]:
    """Alle stereogenen Achsen eines Frames.  Deterministisch sortiert."""
    out: List[dict] = []
    for i in range(len(syms)):
        if syms[i] == "H" or _is_metal_sym(syms[i]):
            continue
        for j in nbrs[i]:
            if j <= i or syms[j] == "H" or _is_metal_sym(syms[j]):
                continue
            # Eine Achse ist eine EINFACHBINDUNG ZWISCHEN ZWEI KONJUGIERTEN ZENTREN.
            sp2 = (_atrop_is_planar_sp2(pts, nbrs, i)
                   and _atrop_is_planar_sp2(pts, nbrs, j))
            if not (sp2 or _atrop_is_mesomeric_amide(syms, nbrs, i, j)):
                continue
            side = _atrop_side_atoms(syms, nbrs, i, j)      # None = im Ring
            if not side or len(side) < 2:
                continue
            # RINGKRITERIUM (16.08.): mindestens EINE Seite muss ein Ringsystem sein -- eine
            # Aryleinheit.  Das fehlte in der ersten Fassung ganz, waehrend das Auge es
            # verlangt; damit nahm ich jede sp2-sp2-Einfachbindung fuer eine Achse.
            ring_i = _atrop_ring_of(nbrs, i)
            ring_j = _atrop_ring_of(nbrs, j)
            if not (ring_i or ring_j):
                continue
            fa = _atrop_flanks(i, j, nbrs, syms, ring_i)
            fb = _atrop_flanks(j, i, nbrs, syms, ring_j)
            if fa is None or fb is None:
                continue                  # keine zwei unterscheidbaren Flanken -> nicht stereogen
            fHiA, _fLoA, kHiA, kLoA = fa
            fHiB, _fLoB, kHiB, kLoB = fb
            # ⚠ Der Dieder wird ueber die HOCHRANGIGEN Flanken gemessen -- dieselbe Wahl wie im
            # Auge.  Mit einem anderen Bezugsatom traegt derselbe Frame das umgekehrte
            # Vorzeichen, und die Enumeration ergaenzte die bereits vorhandene Haendigkeit.
            dih = _atrop_dihedral(pts, fHiA, i, j, fHiB)
            if dih is None:
                continue
            fold = abs(dih)
            if fold > 90.0:
                fold = 180.0 - fold
            if not (_ATROP_TWIST_MIN_DEG <= fold <= _ATROP_TWIST_MAX_DEG):
                continue                              # planar oder senkrecht -> nicht stereogen
            out.append({
                "i": i, "j": j, "dih": dih, "fold": fold,
                "sign": "P" if dih > 0 else "M",
                "side": side,
                "sig": _atrop_axis_sig_flanks(syms, i, j, kHiA, kLoA, kHiB, kLoB),
            })
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


def _atrop_mirror_frame(xyz: str, ax: dict, clash_min: float) -> Optional[str]:
    """Ein Frame mit gekipptem Vorzeichen dieser EINEN Achse.  None, wenn es kollidiert."""
    syms, pts, lines = _parse_xyz(xyz)
    nbrs, _bd = _build_geometric_adjacency(syms, pts)
    new = _atrop_rotate_side(pts, ax, -2.0 * ax["dih"])
    if _atrop_min_nonbonded(syms, new, nbrs, ax["side"]) < clash_min:
        return None
    # Reihenfolge beachten: _format_xyz(orig_lines, syms, positions) -- siehe
    # _coord_angle_corrector.py:135.  Vertauscht liefert es stillen Unsinn.
    return _format_xyz(lines, syms, new)


def _atrop_analyze(xyz: str) -> Optional[dict]:
    syms, pts, _lines = _parse_xyz(xyz)
    if len(syms) < 4:
        return None
    nbrs, _bd = _build_geometric_adjacency(syms, pts)
    axes = _atrop_find_axes(syms, pts, nbrs)
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

    added: List[tuple] = []
    capped = 0
    for sig, (order, xyz, lbl, ax) in sorted(reps.items(), key=lambda kv: kv[1][0]):
        want = "M" if ax["sign"] == "P" else "P"
        if (sig, want) in present:
            continue                             # beide Haende schon da
        if len(added) >= max_added:
            capped += 1
            continue
        try:
            mx = _atrop_mirror_frame(xyz, ax, clash_min)
        except Exception as exc:
            _LOG.debug("atropisomer-enum: Spiegelung fehlgeschlagen (%s): %s", lbl, exc)
            continue
        if mx is None:
            continue                             # Kollision -> nicht bauen, nichts verlieren
        added.append((mx, f"{lbl}_atrop-{want}"))
        present.add((sig, want))

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
    gefunden' meldet, ist wertlos -- er muss die Stufe nennen."""
    if syms[i] == "H" or syms[j] == "H":
        return "H"
    if _is_metal_sym(syms[i]) or _is_metal_sym(syms[j]):
        return "Metall"
    sp2i = _atrop_is_planar_sp2(pts, nbrs, i)
    sp2j = _atrop_is_planar_sp2(pts, nbrs, j)
    if not ((sp2i and sp2j) or _atrop_is_mesomeric_amide(syms, nbrs, i, j)):
        return f"nicht sp2/Amid (sp2 i={sp2i} deg={len(nbrs[i])}, j={sp2j} deg={len(nbrs[j])})"
    side = _atrop_side_atoms(syms, nbrs, i, j)
    if not side:
        return "Bindung liegt im Ring"
    if len(side) < 2:
        return "Seite zu klein"
    ri = _atrop_ring_of(nbrs, i)
    rj = _atrop_ring_of(nbrs, j)
    if not (ri or rj):
        return "keine Seite ist ein Ring"
    fa = _atrop_flanks(i, j, nbrs, syms, ri)
    fb = _atrop_flanks(j, i, nbrs, syms, rj)
    if fa is None:
        return "Seite i: keine zwei unterscheidbaren Flanken"
    if fb is None:
        return "Seite j: keine zwei unterscheidbaren Flanken"
    d = _atrop_dihedral(pts, fa[0], i, j, fb[0])
    if d is None:
        return "Dieder nicht berechenbar"
    f = abs(d)
    if f > 90.0:
        f = 180.0 - f
    if not (_ATROP_TWIST_MIN_DEG <= f <= _ATROP_TWIST_MAX_DEG):
        return f"Verdrillung {f:.1f} ausserhalb {_ATROP_TWIST_MIN_DEG}-{_ATROP_TWIST_MAX_DEG}"
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


