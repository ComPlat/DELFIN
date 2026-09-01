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


def _hat_stereozentrum(xyz) -> Optional[bool]:
    """Traegt das MOLEKUEL mindestens ein Stereozentrum?  None = nicht lesbar.

    ⚠ WARUM DIESE FRAGE UEBERHAUPT NOETIG IST.  `_self_mirror_rmsd` prueft, ob ein
    Frame in FESTER ATOMREIHENFOLGE auf sein Spiegelbild passt -- der Docstring
    dort sagt es selbst ("in dieser Atomreihenfolge").  Ein weicher Konformer tut
    das praktisch nie, egal ob das Molekuel chiral ist.  Der Test misst also
    KONFORMER-Haendigkeit, nicht MOLEKUEL-Chiralitaet.

    GEMESSEN 27.08. auf 4069 legacy gebauten Systemen: `expand_results` haengt bei
    4069 von 4069 = 100 % etwas an.  Dubletten erklaeren das nicht (gegen alle
    uebrigen Frames desselben Systems bleiben 99,8 % neu; unter `permute_dedup`
    mit echten Automorphismen ueberleben 900 von 988).  Aufgeteilt nach dem
    MOLEKUEL:

        >=1 Stereozentrum   1325 / 4069 = 32,6 %   Spiegel = ECHTES neues Isomer
        kein Stereozentrum  2744 / 4069 = 67,4 %   nur ein zweiter Konformer

    Auf einem zweiten Archiv bestaetigt: 32,8 %, Abweichung 0,2 pp.
    ⇒ 68 % der angehaengten Frames bringen NULL Isomergewinn.

    Gelesen wird ueber die Bindungsperzeption aus dem Frame selbst -- ein SMILES
    liegt an dieser Stelle nicht vor.  Nicht lesbar -> None -> das Tor laesst
    durch (nie WENIGER bauen, wenn die Messung fehlt).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return None
    try:
        syms, P, _lines = _parse_xyz(xyz)
    except Exception:
        return None
    if not syms:
        return None
    try:
        block = "%d\n\n" % len(syms) + "\n".join(
            "%s %.6f %.6f %.6f" % (s, p[0], p[1], p[2])
            for s, p in zip(syms, P))
        mol = Chem.MolFromXYZBlock(block)
        if mol is None:
            return None
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.AssignStereochemistryFrom3D(mol)
        zentren = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
        return bool(zentren)
    except Exception:
        return None


def expand_results(results):
    """ADDITIV: haengt je Frame sein Spiegelbild an.  Originale bleiben unberuehrt.

    ``results`` ist die Liste von ``(xyz, label)`` des FF-freien Pfads.  Bit-genauer
    No-op, wenn der Schalter aus ist oder kein Frame ein Spiegelbild hat.

    ZWEI TORE, beide EINZELN schaltbar, beide Vorgabe AUS -> byte-identisch:

      DELFIN_MIRROR_STEREO_GATE=1   nur Systeme mit >=1 Stereozentrum spiegeln.
          Streicht 67,4 % der Systeme und damit 68 % des Preises, bei NULL
          Isomerverlust.

      DELFIN_MIRROR_ONE_PER_SYSTEM=1   EIN Repraesentant statt je Konformer.
          Fuer die ISOMERabdeckung genuegt ein Frame mit gekippter Haendigkeit;
          jeden Konformer zu spiegeln verdoppelt das Archiv, ohne ein Isomer mehr
          zu treffen.  Gemessen: 1325 statt 39 254 Zusatzframes = +1,0 % statt
          +93,7 %.

      DELFIN_MIRROR_QUALITY_GATE=1   nur HEILE Frames spiegeln.

          WARUM (01.09.2026, gemessen an mirrleg6k, 1571 Systeme).  Der Pass hat
          bis heute KEINE Qualitaetspruefung: er spiegelt und haengt an, ohne je
          zu fragen, ob die Vorlage heil ist.  Ergebnis: von 2645 angehaengten
          Frames tragen 1727 einen HARTEN Befund (65,3 %).

          🔑 UND DAS IST REINE VERERBUNG, KEIN NEUER SCHADEN.  `mirror_frame`
          ist eine REFLEXION (`P @ _MIRROR`) und damit eine ISOMETRIE: alle
          Abstaende und Winkel sind exakt invariant, nur Torsionsvorzeichen
          kippen.  Ein Spiegel kann also weder eine Kollision noch eine
          Bindungslaenge verschlechtern -- er ist genau dann hart, wenn seine
          VORLAGE hart war.  Die Zahlen bestaetigen es: 65,3 % der Spiegel gegen
          68,7 % im Bestand.

          ⇒ Darum prueft dieses Tor die VORLAGE, nicht den Spiegel.  Das ist
          nicht nur billiger, es ist die einzig richtige Stelle: `_rg_score` ist
          reflexionsinvariant, am Spiegel gemessen kaeme dasselbe heraus.

          WAS ES KOSTET.  Der Spiegel eines kaputten Frames ist ein zweiter
          kaputter Frame -- er traegt kein Isomer bei, das zaehlt (Nutzerregel:
          "nur mit sehr schlechter Geometrie erreichbare Isomere zaehlen NICHT").
          GEMESSEN auf mirrleg6k: auf allen 8 sperrenden Systemen sind die
          angehaengten Frames ausnahmslos hart; 1032 von 1566 Systemen bekommen
          AUSSCHLIESSLICH harte Spiegel.

    ⚠ DIE TORE SIND GETRENNT, weil sie VERSCHIEDENE Fragen beantworten -- das
      erste "welche Systeme", das zweite "wie viele Frames je System", das dritte
      "welche Vorlagen ueberhaupt".  Sie zu buendeln machte jedes Verdikt
      unzuordenbar.
    """
    if not results or not _is_enabled():
        return results
    max_added = _env_int("DELFIN_MIRROR_MAX_ADDED", 128)
    _stereo_tor = _env_int("DELFIN_MIRROR_STEREO_GATE", 0) == 1
    _einer = _env_int("DELFIN_MIRROR_ONE_PER_SYSTEM", 0) == 1
    _qual_tor = _env_int("DELFIN_MIRROR_QUALITY_GATE", 0) == 1
    _rg = None
    if _qual_tor:
        # SPAET importiert: `_refine_gate` zieht `_h_placement` nach, und ein
        # Modulimport auf Dateiebene waere ein Zyklus.  Faellt der Import aus,
        # ist das Tor AUS -- eine fehlende Abhaengigkeit darf nie stillschweigend
        # Frames streichen.
        try:
            from delfin.manta._refine_gate import _rg_score as _rg
        except Exception as _e:          # pragma: no cover - Verdrahtungsschutz
            _LOG.warning("mirror_enum: QUALITAETSTOR angefordert, aber _rg_score "
                         "nicht importierbar (%s) -- Tor bleibt AUS, es wird "
                         "NICHTS gestrichen", type(_e).__name__)
            _qual_tor = False

    if _stereo_tor:
        # EINMAL je System fragen, nicht je Frame: die Stereozentren des MOLEKUELS
        # aendern sich zwischen Konformeren nicht.  Der erste lesbare Frame
        # entscheidet; ist keiner lesbar, laesst das Tor durch.
        _hat = None
        for (xyz, _lab) in results:
            _hat = _hat_stereozentrum(xyz)
            if _hat is not None:
                break
        if _hat is False:
            _LOG.debug("mirror_enum: STEREO-TOR -- kein Stereozentrum im Molekuel, "
                       "der Spiegel waere nur ein zweiter Konformer (0 Frames)")
            return results
    added: List[Tuple[str, str]] = []
    n_achiral = 0
    n_failed = 0
    n_already = 0
    n_kaputt = 0
    # ── EIN REPRAESENTANT HEISST EINER JE SYSTEM, NICHT EINER JE AUFRUF ─────────────
    # GEMESSEN (01.09.2026, `LOOP_FIRE_TRACE` auf ABUSAU und JEJROI, beide Arme):
    # der Legacy-Aufruf `smiles_converter.py:35932` laeuft ZWEIMAL je System, beide
    # Male im Wiedereintritt der Konformer-Vollstaendigkeit
    # (`outermost = not _CONF_COMPLETE_ACTIVE.value`).
    #
    # Die Idempotenz-Pruefung unten (`endswith("_mirror")`) verhindert nur, dass ein
    # SPIEGEL gespiegelt wird.  Sie verhindert NICHT, dass der zweite Aufruf die
    # naechste unberuehrte Vorlage spiegelt -- unter `_einer` bricht die Schleife
    # nach dem ersten Anhaengen ab, also liefert jeder Aufruf einen weiteren Spiegel.
    # Im Archiv sichtbar als ZWEI `..._mirror`-Etiketten je System, und ABUSAU
    # kommt auf 58 + 2 - 1 = 59 Frames: der zweite Spiegel kostet einen Basis-Frame
    # (`...Δ-conf4_stereo-u` verschwindet).
    #
    # ⚠ NICHT der Aufrufer wurde geschuetzt.  Ein `outermost`-Tor dort schaltet die
    #   Achse auf dem Legacy-Pfad GANZ ab (gemessen an `addroot10`: 58 -> 58, null
    #   Spiegel) -- ein stiller Faehigkeitsverlust, zurueckgenommen als b424e6cd.
    #   Der Vertrag gehoert dorthin, wo er formuliert ist: EIN Repraesentant JE
    #   SYSTEM.  Liegt schon einer vor, ist dieser Pass fertig.
    if _einer and any(str(_l).endswith("_mirror") for _x, _l in results):
        _LOG.debug("mirror_enum: EIN-REPRAESENTANT -- es liegt bereits ein Spiegel "
                   "vor, dieser Aufruf haengt nichts an (Wiedereintritt)")
        return results
    for (xyz, label) in results:
        if len(added) >= max_added:
            break
        # ===== IDEMPOTENZ (18.08.2026) ==========================================
        # Der Pass war NICHT idempotent: das Spiegelbild eines Spiegelbildes ist
        # wieder das Original, und die Dublettenpruefung ``m == xyz`` sieht das
        # nicht, weil sie gegen die EINGABE vergleicht, nicht gegen die Menge.
        # Zweimaliges Anwenden haette also jedes Original ein zweites Mal
        # angehaengt.  Das war bisher folgenlos, weil es genau EINE Aufrufstelle
        # gab -- und genau das aendert sich mit der zweiten, die den Pass vom
        # Schalter DELFIN_FFFREE_SHARED_TAIL unabhaengig macht.
        # Eine Bedingung, die nur unter der heutigen Verdrahtung stimmt, ist eine
        # Falle fuer die naechste.
        if str(label).endswith("_mirror"):
            n_already += 1
            continue
        if _qual_tor:
            # DIE VORLAGE ENTSCHEIDET, nicht der Spiegel (Isometrie, s. Docstring).
            # ⚠ `(-1, ...)` heisst UNLESBAR, nicht "kaputt" -- ein unlesbarer Frame
            #   wird DURCHGELASSEN.  Wer Unlesbarkeit als Defekt zaehlt, streicht
            #   auf einer Nichtmessung, und das ist genau die Bauform, die hier
            #   schon dreimal eine stille Null erzeugt hat.
            try:
                _sc = _rg(xyz)
            except Exception:
                _sc = None
            # ⚠ NUR die KOLLISION vetoiert, NICHT `n_bond_out`.
            #
            # GEMESSEN 01.09. im Selbsttest, und es hat den ersten Entwurf gekippt:
            #     "heile"  Testvorlage -> _rg_score = (0, 4)
            #     "kaputte" Testvorlage -> _rg_score = (0, 2)
            # Die heile scort SCHLECHTER.  `n_bond_out` zaehlt jede Bindung ausser
            # halb des Zielbands und ist auf handgebauten wie auf echten Frames
            # dicht besetzt -- als ABSOLUTE Schwelle ist es unbrauchbar.
            #
            # 🔑 `_rg_score` ist ein VERGLEICHSMASS ("wurde es schlechter?", so
            #    benutzt es `keep_better`), keine Schwelle.  Wer es absolut liest,
            #    liest einen Detektornamen statt einer Messung.
            #    `n_clash` dagegen ist eine Zaehlung echter Ueberlappungen und
            #    braucht keine Kalibrierung: 0 heisst keine, >0 heisst welche.
            if _sc is not None and len(_sc) >= 1 and _sc[0] > 0:
                n_kaputt += 1
                continue
        m = mirror_frame(xyz)
        if m is None:
            n_achiral += 1
            continue
        if m == xyz:
            n_failed += 1
            continue                      # sollte nicht vorkommen; nie eine Dublette anhaengen
        added.append((m, f"{label}_mirror"))
        if _einer:
            # EIN REPRAESENTANT.  Fuer die Isomerabdeckung genuegt ein Frame mit
            # gekippter Haendigkeit -- weitere Spiegel sind Konformere DESSELBEN
            # Spiegelisomers und treffen kein Isomer mehr.
            # ⚠ Der Abbruch steht NACH dem Anhaengen, nicht davor: sonst bricht er
            #   auch dann ab, wenn `mirror_frame` gerade None geliefert hat, und
            #   das System bekaeme GAR keinen Spiegel.  Genau die Bauform, die
            #   heute schon dreimal eine stille Nullmessung erzeugt hat.
            _LOG.debug("mirror_enum: EIN-REPRAESENTANT -- 1 Spiegelframe statt %d",
                       len(results))
            break
    if _qual_tor and n_kaputt:
        # KEINE STILLE STREICHUNG.  Wer nicht sagt, wie viel er weggelassen hat,
        # liest sich hinterher wie "mehr gab es nicht" -- dieselbe Falle wie beim
        # Deckel unten und bei den Faltungen.
        _LOG.info("mirror_enum: QUALITAETSTOR -- %d von %d Vorlagen nicht gespiegelt "
                  "(kaputt laut _rg_score); %d Spiegel angehaengt",
                  n_kaputt, len(results), len(added))
    if not added:
        if _qual_tor and n_kaputt:
            _LOG.warning("mirror_enum: QUALITAETSTOR hat ALLE %d Vorlagen gestrichen "
                         "-- dieses System bekommt KEINEN Spiegel", n_kaputt)
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

    # ===== QUALITAETSTOR (01.09.2026) =========================================
    # ⚠ ALS SKRIPT gestartet liegt `delfin` NICHT im Pfad -- `_refine_gate` waere
    #   dann nicht importierbar und das Tor schaltete sich (korrekt) selbst ab.
    #   Genau das ist beim ersten Lauf passiert: der Selbsttest haette das Tor
    #   fuer kaputt gehalten, obwohl die Fail-safe arbeitete.  Im Paketbetrieb
    #   gibt es das Problem nicht; hier wird die Wurzel nachgetragen.
    import sys as _sys, os.path as _op
    _root = _op.dirname(_op.dirname(_op.dirname(_op.abspath(__file__))))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    # Zwei Proben, und die zweite ist die wichtigere: ein Tor, das nur streicht,
    # ist kein Tor -- es muss eine HEILE Vorlage auch durchlassen.
    os.environ["DELFIN_MIRROR_ENUM"] = "1"
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "1"
    # (a) HEILE Vorlage -> wird gespiegelt.  `chiral` ist die Probe aus Test 1.
    try:
        from delfin.manta._refine_gate import _rg_score as _dbg0
        _sc_gut = _dbg0(chiral)
    except Exception as _e:
        _sc_gut = f"IMPORT-FEHLER {type(_e).__name__}"
    r_gut = expand_results([(chiral, "iso0")])
    ok = len(r_gut) == 2 and r_gut[1][1] == "iso0_mirror"
    print(f"8 QUALITAETSTOR laesst heile Vorlage durch: {'OK' if ok else 'FEHLER'}"
          f"   [_rg_score={_sc_gut}]")
    fails += 0 if ok else 1

    # (b) KAPUTTE Vorlage -> wird NICHT gespiegelt.  Zwei Kohlenstoffe auf 0,40 A
    #     sind eine Kollision, die `_rg_score` sicher sieht.
    # ECHTE Kollision: zwei SUBSTITUENTEN uebereinander.  Cl und Br haengen beide
    # am C, sind untereinander NICHT gebunden und liegen 0,12 A auseinander -- das
    # ist eine Ueberlappung, keine kurze Bindung.  (Der erste Entwurf setzte zwei
    # Kohlenstoffe auf 0,40 A; der Graph machte daraus eine BINDUNG und n_clash
    # blieb null.  Eine Probe, die den Detektor nicht ausloest, prueft nichts.)
    kaputt = _xyz([("C", 0.0, 0.0, 0.0), ("H", 0.0, 1.09, 0.0),
                   ("F", 1.03, -0.36, 0.0), ("Cl", -0.51, -0.36, 1.55),
                   ("Br", -0.51, -0.36, 1.67)])
    try:
        from delfin.manta._refine_gate import _rg_score as _dbg
        _sc_dbg = _dbg(kaputt)
    except Exception as _e:
        _sc_dbg = f"IMPORT-FEHLER {type(_e).__name__}: {_e}"
    # ⚠ EHRLICHKEITSPRUEFUNG VOR DER PROBE.  Ist die Vorlage gar nicht spiegelbar,
    #   haengt `expand_results` auch OHNE Tor nichts an -- ein "zurueckgehalten"
    #   waere dann ein Fehlschluss.  Genau das ist beim zweiten Entwurf passiert.
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "0"
    _spiegelbar = len(expand_results([(kaputt, "iso0")])) == 2
    _hat_clash = isinstance(_sc_dbg, tuple) and len(_sc_dbg) >= 1 and _sc_dbg[0] > 0
    os.environ["DELFIN_MIRROR_QUALITY_GATE"] = "1"
    r_bad = expand_results([(kaputt, "iso0")])
    if not (_spiegelbar and _hat_clash):
        print(f"9 QUALITAETSTOR gegen echte Kollision: UNGEPRUEFT -- die Probe ist "
              f"{'nicht spiegelbar' if not _spiegelbar else 'kollisionsfrei'} "
              f"[_rg_score={_sc_dbg}].  Eine von Hand gebaute Probe, die zugleich "
              f"CHIRAL und KOLLIDIEREND ist, ist mir nicht gelungen; die Kalibrierung "
              f"gehoert auf echte Archivframes, nicht hierher.")
    else:
        ok = len(r_bad) == 1
        print(f"9 QUALITAETSTOR haelt kollidierende Vorlage zurueck: "
              f"{'OK' if ok else 'FEHLER'}   [_rg_score={_sc_dbg}]")
        fails += 0 if ok else 1

    # ── 10  EIN REPRAESENTANT JE SYSTEM, AUCH BEI ZWEI AUFRUFEN ────────────────────
    # Der Fall, der `ABUSAU` einen Frame gekostet hat: der Legacy-Pfad ruft den Pass
    # ZWEIMAL (gemessen mit LOOP_FIRE_TRACE), und ohne diese Probe haengt der zweite
    # Aufruf einen weiteren Spiegel an.
    _alt_einer = os.environ.get("DELFIN_MIRROR_ONE_PER_SYSTEM")
    try:
        os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = "1"
        _rows = [("C", 0.0, 0.0, 0.0), ("N", 1.5, 0.0, 0.0),
                 ("O", 0.0, 1.5, 0.0), ("F", 0.0, 0.0, 1.5)]
        _a = _xyz(_rows)
        _rows2 = [("C", 0.1, 0.0, 0.0), ("N", 1.6, 0.0, 0.0),
                  ("O", 0.0, 1.6, 0.0), ("F", 0.0, 0.0, 1.6)]
        _b = _xyz(_rows2)
        _erst = expand_results([(_a, "iso1"), (_b, "iso2")])
        _zweit = expand_results(_erst)
        ok = (len(_erst) == 3 and len(_zweit) == 3)
        print(f"10 ZWEITER Aufruf haengt NICHTS an: erst {len(_erst)}, dann "
              f"{len(_zweit)} Frames  {'OK' if ok else 'FEHLER'}")
        fails += 0 if ok else 1
        # Gegenprobe: OHNE den Ein-Repraesentant-Schalter darf er sehr wohl erneut
        # anhaengen -- sonst waere die Probe oben nur ein abgeschalteter Pass.
        os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = "0"
        _drei = expand_results(_erst)
        ok2 = len(_drei) > len(_erst)
        print(f"   GEGENPROBE ohne ONE_PER_SYSTEM haengt weiter an: "
              f"{len(_erst)} -> {len(_drei)}  {'OK' if ok2 else 'FEHLER'}")
        fails += 0 if ok2 else 1
    finally:
        if _alt_einer is None:
            os.environ.pop("DELFIN_MIRROR_ONE_PER_SYSTEM", None)
        else:
            os.environ["DELFIN_MIRROR_ONE_PER_SYSTEM"] = _alt_einer

    print(f"\n{11 - fails}/11 bestanden (Probe 9 nur wenn sie den Detektor ausloest)")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
