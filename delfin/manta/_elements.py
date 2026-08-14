"""_elements -- DIE EINE Quelle fuer "ist dieses Element ein Metall?".

WARUM ES DAS GIBT (Inventur 14.08.2026).  Diese Frage wurde im Doppelbaum an rund
95 Stellen eigenstaendig beantwortet -- 22 in DELFIN, ~73 im Auge.  Der vorhandene
Pruefer `weddell/tools/metal_predicate_audit.py` nennt selbst die richtige Regel:

    "the question is not 'how many copies' but 'on WHICH ELEMENTS do they differ',
     and that is a measurement, not an opinion."

Gemessen wurde: bei URAN sagen 4 Praedikate Metall und 5 sagen nein.  Dieselbe
Spaltung bei allen Lanthanoiden und Actinoiden, beim s-Block und bei Tl.

ZWEI FEHLERKLASSEN, die dabei sichtbar wurden und die kein A/B je finden koennte:

  (a) NEUN Module widersprechen SICH SELBST.  Ihr `_METAL_Z_RANGES` enthaelt
      `range(57, 81)`, deckt also Ce (58) bis Lu (71) ab -- ihre Symboltabelle
      springt aber von `"La": 57` direkt auf `"Hf": 72`.  Der Nachschlag liefert
      `None`, das Praedikat gibt False.  Alle Lanthanoide ausser Lanthan sind dort
      keine Metalle, obwohl der eigene Zahlenbereich sie einschliesst.
      (_h_vsepr_realism, _cp_piano_stool, _pi_h_projector, _fix_sp2n_planarize,
       _fix_sp3_n_pyramidality, _fix_sp3_h_tetrahedrality, _coord_angle_corrector,
       _fix_bridging_anion, _fix_wuxqak_sp3_c_linear)

  (b) DREI Module schalten bei einem Importfehler die Metallerkennung GANZ ab:
      `return sym in _METAL_SET if _METAL_SET else False` -- bei leerer Menge fuer
      JEDES Element False, ohne Fehlermeldung.  Das ist keine reduzierte Liste,
      das ist Blindheit.  (_post_optimizer, _energy_terms, _fragment_archetypes)

DIESES MODUL HAT KEINE ABHAENGIGKEITEN ausser `os`.  Das ist der Punkt: der Kanon
stand bisher in `smiles_converter` (36k Zeilen), weshalb jeder Importeur einen
try/except-Fallback baute -- und GENAU diese Fallbacks sind die divergierenden
Kopien.  Ein Modul, das nichts zieht, kann von ueberall geholt werden, auch vom
Auge, ohne Zirkel und ohne Fallback.

=====  DIE ENTSCHEIDUNG, ELEMENT FUER ELEMENT  =====

DRIN (68) -- alles, was im Korpus als Koordinationszentrum auftritt:
    s-Block      Li Na K Rb Cs | Be Mg Ca Sr Ba
    d-Block      Sc..Zn | Y..Cd | Hf..Hg
    f-Block      La..Lu | Ac Th Pa U Np Pu
    p-Block      Al Ga In Tl Sn Pb Bi Po

DRAUSSEN, und warum:
    Ge As Sb Te   METALLOIDE.  Sie sind im Bau bereits als DONOREN gefuehrt
                  (`decompose._METALLOID_DONORS`, `smiles_converter._METALLOID_MD_DONORS`,
                  dieselben Elemente).  Sie hier zusaetzlich als Metall zu fuehren
                  hiesse, dasselbe Atom gleichzeitig als Zentrum und als Donor zu
                  behandeln.  ~30 Fundstellen taten das bisher.
    Fr Ra         kommen im CCDC-Korpus als Koordinationszentrum praktisch nicht vor;
                  11 Fundstellen fuehrten sie trotzdem.
    Am..Lr        dito, 22 Fundstellen fuehrten sie.

WER DIESE MENGE AENDERT, aendert die Bedeutung jeder Kollisions-, Kollaps- und
Koordinationsmessung gleichzeitig.  Vorher `weddell/tools/metal_predicate_audit.py`
laufen lassen: er sagt, WELCHE Elemente betroffen sind, nicht nur wie viele.
"""
from __future__ import annotations

import os

# Der Kanon.  Wortgleich mit smiles_converter._METALS (dort historisch entstanden);
# ab jetzt liegt die Wahrheit HIER und smiles_converter liest sie.
METALS = frozenset("""
    Li Na K Rb Cs  Be Mg Ca Sr Ba
    Sc Ti V Cr Mn Fe Co Ni Cu Zn
    Y Zr Nb Mo Tc Ru Rh Pd Ag Cd
    La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
    Hf Ta W Re Os Ir Pt Au Hg
    Al Ga In Tl Sn Pb Bi Po
    Ac Th Pa U Np Pu
""".split())

# Die Metalloide, ausdruecklich als EIGENE Menge statt als Streitfall in METALS.
# Wer sie braucht, holt sie hier -- statt sie still in die eigene Metallliste zu
# schreiben, wie es ~30 Fundstellen taten.
METALLOIDS = frozenset("Ge As Sb Te Se B Si".split())

# Isotopenschreibweisen, die sonst als unbekanntes Schweratom durchgehen.
_ISO = {"D": "H", "T": "H"}


def normalise(sym: str) -> str:
    """Symbol auf sein Element abbilden: 'D'/'T' -> 'H', 'Fe2+' -> 'Fe'.

    GROSSSCHREIBUNG IST NICHT EINDEUTIG.  Ein CCDC-Label 'CL1' kann Chlor oder ein
    Kohlenstoff-Label sein; `^[A-Z][a-z]?` liest daraus 'C'.  Diese Funktion loest das
    NICHT -- sie normalisiert nur, was eindeutig ist, und laesst den Rest unveraendert,
    damit der Aufrufer den Zweifelsfall sieht statt ihn geschenkt zu bekommen.
    """
    s = (sym or "").split(".")[0].strip()
    if s in _ISO:
        return _ISO[s]
    if len(s) >= 2 and s[0].isupper() and s[1].islower():
        return s[:2]
    return s[:1].upper() if s else s


def is_metal(sym: str) -> bool:
    """Die eine Metallfrage.  Nimmt rohe Symbole entgegen, auch 'Fe2+' oder 'D'."""
    return normalise(sym) in METALS


def is_metalloid(sym: str) -> bool:
    return normalise(sym) in METALLOIDS


def unified_enabled() -> bool:
    """Ob migrierte Module die EINE Quelle benutzen sollen.

    Vorgabe AUS -> jedes migrierte Modul verhaelt sich byte-identisch wie zuvor.
    Das ist Absicht: die Vereinheitlichung aendert das Urteil ueber Lanthanoide,
    Actinoide und den s-Block gleichzeitig in vielen Modulen, und so etwas wird
    gemessen und nicht geglaubt.
    """
    return os.environ.get("DELFIN_FFFREE_METAL_UNIFIED", "0") == "1"


# ===== Z-BASIERTE VARIANTE =====
# Zwei Praedikate im Bauer nehmen die ORDNUNGSZAHL statt des Symbols
# (`_system_classifier._is_metal`, `_rotamer_diversity._is_metal`).  Gleicher Name,
# andere Signatur -- der Pruefer musste das eigens abfangen ("SAME NAME, DIFFERENT
# SIGNATURE").  Damit die Vereinheitlichung auch sie erreicht, steht die Umrechnung
# HIER und nicht in jedem Modul noch einmal.
_Z_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La "
    "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po "
    "At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr").split()
SYMBOL_BY_Z = {i + 1: s for i, s in enumerate(_Z_SYMBOLS)}
METAL_Z = frozenset(z for z, s in SYMBOL_BY_Z.items() if s in METALS)


def is_metal_z(z) -> bool:
    """Die Metallfrage ueber die Ordnungszahl -- dieselbe Menge wie `is_metal`."""
    try:
        return int(z) in METAL_Z
    except (TypeError, ValueError):
        return False
