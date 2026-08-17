"""delfin.manta._mirror_enum — DER SPIEGELABSCHLUSS des Manifolds.

DER BEFUND (gemessen 17.08.2026 auf `rows_census50k`, 38 556 Saetze)
-------------------------------------------------------------------
**Das Korpus traegt keine Stereochemie.**  9 von 129 314 SMILES haben eine
Chiralitaetsmarke `@` (0,007 %), 2 eine Doppelbindungsmarke.  Geprueft auf den 2269
Fehlschlaegen (0 mit `@`), den 588 reinen C-Faellen (0) UND den 5726 Erfolgen (0) --
gleiche Grundrate, also eine Eigenschaft der Grundgesamtheit, kein Merkmal der Fehler.

=> Der Bauer bekommt an KEINEM Zentrum eine Vorgabe.  Jede Haendigkeit ist freie Wahl.
Die Vollstaendigkeitsregel ("immer + und - bauen") ist damit nicht eine Zusatzforderung,
sondern der EINZIGE Mechanismus, der die Kristallhaendigkeit je treffen kann.

Die Fehlrate ist FLACH ueber alle Zentrenarten -- C[sp3] 21,4 % · N[sp3] 19,0 % ·
P 24,9 % · Cu 16,7 % · Zn 24,9 % · Co 21,2 %; Summe sp3 20,7 % gegen Metall 22,4 %.
Flach heisst: es ist die allgemeine Erzeugung, keine Regel eines Moduls.  Und die
Ursache ist ABDECKUNG, nicht Auswahl: von 2269 Fehlschlaegen sind **1946 (85,8 %)
"Haendigkeit NIE gebaut"**, nur 323 sind `joint`.

WARUM GERADE DIE GESAMTSPIEGELUNG
---------------------------------
Von den 1946 fehlen bei **1326 (68,1 %) ALLE** Zentren des Systems -- dort ist das
Spiegelbild EXAKT der fehlende Isomer.  1166 (59,9 %) haben ueberhaupt nur ein Zentrum.
**889 der 1326 (67 %) bauen <= 5 Frames**, dort ist ein Zusatzframe am billigsten.

Und es ist die EINZIGE Operation, die fuer Kohlenstoff, Stickstoff und Metall gleich
funktioniert -- sie fragt nicht, was das Zentrum IST.  Das zaehlt, weil die drei Klassen
sonst drei Module braechten: `_stereocenter_enum` deckt per Elementzeile (`:191`) nur
N/P/As/Sb/Bi ab und verlangt Metallbindung (`:195`); Kohlenstoff (41 % der fehlenden
Zentren) und Metalle (40 %) haben GAR KEINEN Enumerator.

Eine Spiegelung ist eine ISOMETRIE: jede Bindungslaenge, jeder Winkel, jeder
M-D-Abstand, jeder Torsionsbetrag bleibt exakt erhalten -- nur die Vorzeichen kippen.
Keine Relaxation noetig, kein Clash moeglich, kein Qualitaetsrisiko.  Ein Frame je
Struktur.

VORPRUEFUNG (17.08., im Quelltext belegt): BEIDE Entdopplungen sind chiralitaetssicher
GEBAUT und lassen ein Spiegelbild stehen --
  * `permute_dedup._kabsch_rmsd_perm:197` "proper rotation only ... reflections
    FORBIDDEN ... enantiomeric frames never align and are KEPT" (und ohnehin AUS)
  * `assemble_complex._kabsch_rot:63` "determinant-corrected to forbid reflection"

GRENZEN, AUSDRUECKLICH
----------------------
* Die **620 Teilfaelle** (im Mittel 3,74 Zentren, davon 1,67 falsch; haeufigstes Muster
  **2 Zentren / 1 fehlend**, 196x) trifft dieser Pass NICHT.  Die Spiegelung verbindet
  RR<->SS, nicht RR<->RS.  Diastereomere brauchen zentrumsweise Inversion -- ein anderer
  Mechanismus, hier bewusst NICHT mitgebaut.
* Ein ACHIRALES Molekuel ist auf sein Spiegelbild abbildbar; der Zusatzframe waere eine
  Dublette.  Der Test dagegen ist ein Kabsch mit VERBOTENER Spiegelung in FESTER
  Atomreihenfolge (`assemble_complex._kabsch_rot`).  ⚠ Er UNTERerkennt Achiralitaet,
  wenn die Deckabbildung eine Atompermutation braucht -- dann bleibt eine Dublette
  stehen, die die nachgelagerte Entdopplung faengt.  Das ist die sichere Richtung:
  nie ein fehlender Frame, gelegentlich ein ueberfluessiger.
* REIHENFOLGE: dieser Pass gehoert ZULETZT.  `_stereocenter_enum` liest `present`, BEVOR
  es ergaenzt -- genau so hat `trans208` (29 gemischte trans-Anordnungen) die
  Stereozentren-Falten verdraengt und ein CCDC-Isomer gekostet, obwohl beide Paesse
  additiv sind.  Der FF-freie Anschluss sitzt darum im `_ffree_shared_tail` (Aufruf
  :32429), also NACH der Faltenenumeration (:32361).

Additiv (Originale bleiben unberuehrt), deterministisch (kein RNG, feste Reihenfolge),
FF-frei (reine Koordinatenoperation, kein Kraftfeld).  Vorgabe AUS -> byte-identisch.
"""
from __future__ import annotations

import logging
import os
from typing import List, Optional, Tuple

import numpy as np

from delfin.manta._coord_angle_corrector import _format_xyz, _parse_xyz

_LOG = logging.getLogger(__name__)

# Spiegelung an der xy-Ebene.  JEDE uneigentliche Operation taugt (det = -1); diese ist
# die einfachste und braucht keine Zentrierung, weil eine Ebenenspiegelung durch den
# Ursprung alle paarweisen Abstaende ohnehin erhaelt.
_MIRROR = np.diag([1.0, 1.0, -1.0])


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _is_enabled() -> bool:
    return (_env_int("DELFIN_MIRROR_ENUM", 0) == 1
            or _env_int("DELFIN_FFFREE_MIRROR_ENUM", 0) == 1)


def _self_mirror_rmsd(syms: List[str], P: np.ndarray, Pm: np.ndarray) -> Optional[float]:
    """Schwer-Atom-RMSD zwischen Frame und Spiegelbild unter EIGENTLICHER Rotation.

    Nahe null => das Molekuel ist (in dieser Atomreihenfolge) achiral, der Spiegelframe
    waere eine Dublette.  Nutzt `assemble_complex._kabsch_rot`, das Spiegelungen per
    Determinantenkorrektur VERBIETET -- ohne das waere jeder Frame trivial auf sein
    Spiegelbild abbildbar und der Test wertlos.
    """
    try:
        from delfin.manta.assemble_complex import _kabsch_rot
    except Exception:
        return None
    heavy = [i for i, s in enumerate(syms) if s != "H"] or list(range(len(syms)))
    A = P[heavy]
    B = Pm[heavy]
    A = A - A.mean(axis=0)
    B = B - B.mean(axis=0)
    try:
        R = _kabsch_rot(A, B)
        return float(np.sqrt(((A @ R.T - B) ** 2).sum(axis=1).mean()))
    except Exception:
        return None


def mirror_frame(xyz: str) -> Optional[str]:
    """Das Spiegelbild EINES Frames, oder None wenn es keines gibt / achiral ist."""
    if not xyz:
        return None
    try:
        syms, P, lines = _parse_xyz(xyz)
    except Exception:
        return None
    if not syms or P is None or len(syms) < 4:
        return None
    P = np.asarray(P, dtype=float)
    Pm = P @ _MIRROR
    if not np.all(np.isfinite(Pm)):
        return None
    min_rmsd = _env_float("DELFIN_MIRROR_MIN_RMSD", 0.10)
    r = _self_mirror_rmsd(syms, P, Pm)
    if r is not None and r < min_rmsd:
        return None                       # achiral in dieser Reihenfolge -> kein Zugewinn
    try:
        return _format_xyz(lines, syms, Pm)
    except Exception:
        return None


def expand_results(results):
    """ADDITIV: haengt je Frame sein Spiegelbild an.  Originale bleiben unberuehrt.

    ``results`` ist die Liste von ``(xyz, label)`` des FF-freien Pfads.  Bit-genauer
    No-op, wenn der Schalter aus ist oder kein Frame ein Spiegelbild hat.
    """
    if not results or not _is_enabled():
        return results
    max_added = _env_int("DELFIN_MIRROR_MAX_ADDED", 128)
    added: List[Tuple[str, str]] = []
    n_achiral = 0
    n_failed = 0
    for (xyz, label) in results:
        if len(added) >= max_added:
            break
        m = mirror_frame(xyz)
        if m is None:
            n_achiral += 1
            continue
        if m == xyz:
            n_failed += 1
            continue                      # sollte nicht vorkommen; nie eine Dublette anhaengen
        added.append((m, f"{label}_mirror"))
    if not added:
        return results
    if len(added) >= max_added:
        # KEINE STILLE KUERZUNG.  Ein Deckel, der nicht meldet, liest sich hinterher wie
        # "mehr gab es nicht" -- derselbe Fehler wie bei den Faltungen.
        _LOG.warning("mirror_enum: bei %d Spiegelframes gedeckelt "
                     "(DELFIN_MIRROR_MAX_ADDED); %d Frames nicht mehr geprueft",
                     max_added, max(0, len(results) - max_added))
    _LOG.debug("mirror_enum: %d Spiegelframes ergaenzt (%d achiral/uebersprungen, %d ohne)",
               len(added), n_achiral, n_failed)
    return list(results) + added


# ---------------------------------------------------------------------------
# Selbsttest:  python delfin/manta/_mirror_enum.py
# ---------------------------------------------------------------------------
def _self_test() -> int:
    def _xyz(rows):
        out = [str(len(rows)), "test"]
        for s, x, y, z in rows:
            out.append(f"{s:<2}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
        return "\n".join(out) + "\n"

    def _dmat(xyz):
        s, p, _ = _parse_xyz(xyz)
        return np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)

    def _chirality(xyz):
        s, p, _ = _parse_xyz(xyz)
        return float(np.dot(p[1] - p[0], np.cross(p[2] - p[0], p[3] - p[0])))

    fails = 0
    os.environ["DELFIN_MIRROR_ENUM"] = "1"

    # Ein CHIRALES Zentrum: C mit vier verschiedenen Substituenten.
    chiral = _xyz([("C", 0.0, 0.0, 0.0),
                   ("F", 1.10, 0.0, 0.30),
                   ("Cl", -0.55, 0.95, 0.30),
                   ("Br", -0.55, -0.95, 0.30),
                   ("H", 0.0, 0.0, -1.10)])
    m = mirror_frame(chiral)

    ok = m is not None
    print(f"1 chirales Zentrum hat ein Spiegelbild: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    if m:
        d0, d1 = _dmat(chiral), _dmat(m)
        ok = bool(np.max(np.abs(d0 - d1)) < 1e-9)
        print(f"2 ISOMETRIE -- Abstandsmatrix identisch (max Delta "
              f"{np.max(np.abs(d0 - d1)):.2e}): {'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1

        c0, c1 = _chirality(chiral), _chirality(m)
        ok = (c0 * c1 < 0) and abs(abs(c0) - abs(c1)) < 1e-9
        print(f"3 Haendigkeit gekippt ({c0:+.4f} -> {c1:+.4f}), Betrag gleich: "
              f"{'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1

    # Ein ACHIRALES Molekuel (planar, quadratisch): Spiegelbild ist deckungsgleich.
    planar = _xyz([("Pt", 0.0, 0.0, 0.0),
                   ("Cl", 2.30, 0.0, 0.0),
                   ("Cl", -2.30, 0.0, 0.0),
                   ("N", 0.0, 2.05, 0.0),
                   ("N", 0.0, -2.05, 0.0)])
    ok = mirror_frame(planar) is None
    print(f"4 achiral (planar) wird uebersprungen: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # ADDITIV: Originale bleiben, genau ein Spiegelframe kommt dazu.
    res = [(chiral, "iso0")]
    out = expand_results(res)
    ok = (len(out) == 2 and out[0] == res[0] and out[1][1] == "iso0_mirror")
    print(f"5 additiv, Original unberuehrt, Label 'iso0_mirror': {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # AUSGESCHALTET -> byte-identisch (dieselbe Liste, unveraendert).
    os.environ["DELFIN_MIRROR_ENUM"] = "0"
    ok = (expand_results(res) == res)
    print(f"6 ausgeschaltet byte-identisch: {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    # Zweimal spiegeln ergibt wieder das Original (Involution).
    os.environ["DELFIN_MIRROR_ENUM"] = "1"
    back = mirror_frame(m) if m else None
    ok = back is not None and np.max(np.abs(_dmat(back) - _dmat(chiral))) < 1e-9 \
        and _chirality(back) * _chirality(chiral) > 0
    print(f"7 zweimal gespiegelt = Original (Involution): {'OK' if ok else 'FEHLER'}")
    fails += 0 if ok else 1

    print(f"\n{7 - fails}/7 bestanden")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
