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


def _atrop_side_atoms(nbrs: List[List[int]], i: int, j: int) -> Optional[Set[int]]:
    """Alle Atome auf der j-Seite der Bindung i-j (ohne i).

    None, wenn die Bindung IM RING liegt -- dann gibt es keine zwei Seiten, und eine Drehung
    wuerde die Struktur zerreissen.  Das ist zugleich der Ringtest: wird i beim Absuchen der
    j-Seite erreicht, schliesst sich ein Ring ueber die Bindung.
    """
    seen = {j}
    stack = [k for k in nbrs[j] if k != i]
    while stack:
        cur = stack.pop()
        if cur == i:
            return None                       # Ringschluss -> keine Drehachse
        if cur in seen:
            continue
        seen.add(cur)
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


def _atrop_ref_partner(syms: Sequence[str], nbrs: List[List[int]],
                       c: int, other: int) -> Optional[int]:
    """Deterministischer Bezugsnachbar: schweres Atom, kleinster Index.  KEINE CIP-Regel --
    das Vorzeichen dient hier nur der UNTERSCHEIDUNG zweier Frames, nicht ihrer Benennung.
    Deshalb darf es auch nicht als M/P im chemischen Sinn gelesen werden."""
    h = sorted(k for k in _atrop_heavy_nbrs(syms, nbrs, c) if k != other)
    return h[0] if h else None


def _atrop_axis_sig(syms: Sequence[str], nbrs: List[List[int]], i: int, j: int) -> tuple:
    """Kanonische Achsensignatur -- unabhaengig von der Atomnummerierung, damit dieselbe Achse
    ueber verschiedene Frames hinweg wiedererkannt wird."""
    def side_key(c: int, other: int) -> tuple:
        s = _atrop_side_atoms(nbrs, other, c)
        if not s:
            return (syms[c], -1, ())
        return (syms[c], len(s), tuple(sorted(syms[k] for k in s)))
    return tuple(sorted((side_key(i, j), side_key(j, i)), key=repr))


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
            side = _atrop_side_atoms(nbrs, i, j)      # None = im Ring
            if not side or len(side) < 2:
                continue
            a = _atrop_ref_partner(syms, nbrs, i, j)
            b = _atrop_ref_partner(syms, nbrs, j, i)
            if a is None or b is None:
                continue                              # ohne zwei Bezugsatome kein Dieder
            dih = _atrop_dihedral(pts, a, i, j, b)
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
                "sig": _atrop_axis_sig(syms, nbrs, i, j),
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
