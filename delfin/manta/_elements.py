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


# ===== DER KOVALENZRADIUS -- ZWEITE GROESSE IN DERSELBEN QUELLE (16.08.2026) =====
#
# INVENTUR: 100 Radientabellen im Doppelbaum (30 im Bauer, 70 im Auge).  Fuer EISEN
# stehen darin  0.90 · 0.95 · 1.16 · 1.25 · 1.30 · 1.32 · 1.42 · 1.50 · 1.52  --
# Spanne 0.62 A.  Mit dem Bindungsfaktor 1.30 heisst "Fe-N gebunden" je nach Tabelle
# unter 2.09 oder unter 2.90 A.  Das ist keine Toleranz, das ist Beliebigkeit.
#
# DIE SCHLIMMSTE STELLE IST DER KANON DES BAUERS SELBST: `_bond_decollapse._COV`
# fuehrt 15 Elemente und KEIN EINZIGES METALL.  `_ideal_bond("Fe","N")` faellt auf
# den Default und liefert 0.90 + 0.71 = 1.61 A -- eine erfundene Zahl, und sie traegt
# das SELBSTGATE.  Sieben Augen-Detektoren lesen dieselbe Funktion und kompensieren
# das Loch durch aufgeblaehte Faktoren (1.40 / 1.45 / 1.65); `metric_md_short_collapse`
# schreibt den Grund selbst hin: "bd._COV is missing all TM radii".
#
# WARUM CORDERO 2008 UND NICHT PYYKKO: rund 90 % des Baums fuehrt faktisch schon
# Cordero-Zahlen -- nur an zwei Stellen falsch als "Pyykkoe 2009" etikettiert
# (`smiles_converter.py:139`, `find_md_break.py:98`; Pyykkoes Fe waere 1.16, nicht
# 1.32).  Auf Pyykko umzustellen wuerde JEDE kalibrierte Schwelle im Auge gleichzeitig
# verschieben -- eine Zahl aendern und alles neu eichen ist teurer als eine Etikette
# korrigieren.  Hier wird die im Baum HERRSCHENDE Konvention explizit gemacht, keine
# neue eingefuehrt: die Werte sind byte-genau `smiles_converter._COVALENT_RADII`.
#
# SPIN: Cordero fuehrt Cr/Mn/Fe/Co zweiwertig (low/high spin).  Diese Tabelle ist
# durchgehend LOW SPIN (Mn 1.39, Fe 1.32, Co 1.26).  Beide Spinzustaende in EINER
# Tabelle zu mischen ist ein eigener Fehler und steht heute in
# `find_coord_geometry_realism.py:99` (Mn high spin neben Ni low spin).
#
# ⚠ DIES IST HEUTE DIE DRITTE IMPLEMENTIERUNG, NICHT DIE EINE QUELLE.  Es gibt
# bereits `find_metal_atom_overlap.cov_radius` (Cordero + mendeleev-Fallback,
# Default 1.0, D-Behandlung) und `find_ligand_specific.cov_radius` (eigene Tabelle,
# eine Zeile).  Beide sind MIGRATIONSZIELE, keine Konkurrenz -- "eine Quelle" ist
# dieses Modul erst, wenn sie hierher delegieren.  Solange das nicht geschehen ist,
# ist der Satz "wir haben das vereinheitlicht" FALSCH und darf nicht behauptet werden.
#
# Quelle: B. Cordero et al., Dalton Trans. 2008, 2832-2838.
COV_R = {
    # Hauptgruppe (in fast allen Kopien identisch -- unstrittig)
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39,
    "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    # s-Block
    "Li": 1.28, "Na": 1.66, "K": 2.03, "Rb": 2.20, "Cs": 2.44,
    "Be": 0.96, "Mg": 1.41, "Ca": 1.76, "Sr": 1.95, "Ba": 2.15,
    "Al": 1.21, "Ga": 1.22, "In": 1.42, "Tl": 1.45,
    # 3d -- LOW SPIN fuer Cr/Mn/Fe/Co
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    # 4d
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47,
    "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    # 5d
    "La": 2.07, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51,
    "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    # Lanthanoide
    "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98,
    "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90,
    "Yb": 1.87, "Lu": 1.87,
    # Actinoide
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
}

COV_R_DEFAULT = 1.50        # dieselbe Vorgabe wie weddell/detectors/_bond_criterion.py


def covalent_radius(sym: str) -> float:
    """Kovalenzradius in Angstroem -- MIT Metallen, mit Isotopennormierung (D/T -> H).

    Die Vorgabe 1.50 ist bewusst GROSSZUEGIG: ein unbekanntes Element ist eher schwer
    als leicht, und ein zu KLEINER Radius laesst eine echte Bindung VERSCHWINDEN --
    genau der Fehler, der heute mit 0.90 fuer jedes Metall im Bauer steht.  Ein zu
    grosser fasst sie nur zu weit.
    """
    return COV_R.get(normalise(sym), COV_R_DEFAULT)


def cov_radii_enabled() -> bool:
    """DIE EINE Lesestelle fuer DELFIN_FFFREE_COV_METALS (Vorgabe AUS -> byte-identisch).

    Getrennt von `unified_enabled()`, weil es zwei verschiedene Fragen sind: WELCHE
    Symbole sind Metalle (Praedikat) und WIE GROSS sind sie (Radius).  Ein gemeinsamer
    Schalter waere zwei Aenderungen in EINER Achse -- und dann sagt ein Verdikt nicht
    mehr, welche von beiden gewirkt hat.
    """
    return os.environ.get("DELFIN_FFFREE_COV_METALS", "0") == "1"
