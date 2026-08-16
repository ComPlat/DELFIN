"""_oxidation_state.py — Oxidationsstufe und d-Elektronenzahl aus der LIGANDENBILANZ.

WOZU (16.08.2026, User).  Drei gemessene Baustellen haengen an derselben fehlenden Zahl:

  1. `smiles_converter._PREFERRED_CN4_GEOMETRY` fuehrt `'Cu': 'TH'` -- Tetraeder.  Das gilt
     fuer Cu(I) d10.  **Cu(II) ist d9 und quadratisch-planar bis JT-gestreckt-oktaedrisch**,
     und Cu(II) ist der weit haeufigere Fall.  Dieselbe Zeile fuer beide Stufen.  Ebenso
     `'Au': 'SQ'` -- richtig fuer Au(III) d8, falsch fuer Au(I) d10, das LINEAR ist.
     Der Bauer denkt bereits in d-Zahlen (die Kommentare dort sagen "d8 metals",
     "d6 low-spin", "d0-d5, d7, d10") -- aber als fest verdrahtete ELEMENTLISTE.
  2. JAHN-TELLER ist genau bei ungleicher e_g-Besetzung stark: **d4 HS, d7 LS, d9**.  Die
     e_g-Orbitale zeigen auf die Liganden, die Verzerrung liegt bei 0.2-0.4 A -- messbar.
     (Die t2g-Faelle d1/d2/d4 LS/d5 LS/d6 HS/d7 HS liegen bei ~0.05 A und verschwinden im
     Rauschen.)  Gemessen am 16.08. auf `_ml_prior.json`: **Cu-N Schiefe 3.7x bei n=12339,
     Cu-O 3.3x bei n=7987** -- der lange obere Schwanz der "zwei langen" Bindungen.  Mn und
     Cr zeigen ihn NICHT (0.9x / 1.2x), weil sie in der CSD ueberwiegend als Mn(II) d5 und
     Cr(III) d3 vorkommen -- **nicht JT-aktiv**.  Am ELEMENT ist das nicht zu entscheiden,
     nur an der STUFE.
  3. Der M-D-Laengenschluessel: HS und LS unterscheiden sich um 0.1-0.2 A (Co(II) HS ~2.1,
     Co(III) LS ~1.9).  Die Bauer-Tabelle fuehrt EINEN Wert je (Metall, Donor).

⚠ WARUM NICHT AUS DER FORMALLADUNG.  DELFIN schreibt Formalladungen, die der SMILES-Schreiber
zur VALENZERFUELLUNG setzt: `[Cu-2]`, `[Re-2]`, `[Ti-2]`, `[O-][Re]`.  Das sind keine
physikalischen Oxidationsstufen -- das Metall traegt dort eine NEGATIVE Formalladung, waehrend
es real positiv ist.  Die SUMME aller Formalladungen ist dagegen korrekt: sie ist die
Gesamtladung des Komplexes.

DIE RECHNUNG ist die klassische ionische Zaehlung:

    OS(Metall)  =  Gesamtladung(Komplex)  -  Summe(Ligandenladungen)

    Beispiel [Cu(acac)2], neutral:  0 - (-1 -1) = +2  ->  Cu(II), d9   ✓
    Beispiel [Re(=O)Cl3(PPh3)2]:    0 - (-2 -3 +0) = +5 -> Re(V), d2   ✓

    d-Zahl = Gruppennummer - OS

⚠⚠ DIE WICHTIGSTE REGEL: IM ZWEIFEL NICHTS.  Ist auch nur EIN Donor nicht klassifizierbar,
liefert dieses Modul `None` -- keine geschaetzte Stufe.  Eine FALSCHE Oxidationsstufe ist
schlechter als keine, weil sie die Geometriewahl aktiv verdirbt: sie wuerde ein Cu(I) d10 als
d9 behandeln und ihm eine JT-Streckung aufzwingen, die es nicht hat.  Der Aufrufer faellt bei
`None` auf das heutige Verhalten zurueck, und das ist per Konstruktion byte-identisch.

Ebenso `None` bei MEHREREN Metallen: die Gesamtladung laesst sich dann nicht eindeutig
verteilen (ein gemischtvalentes Fe(II)/Fe(III) waere sonst zwei Mal Fe(2.5)).

Geometrie-frei, deterministisch, ohne Kristall, ohne Elementtabelle mit Radien.
"""
from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

_LOG = logging.getLogger(__name__)

# Gruppennummer der Uebergangsmetalle -- d-Zahl = Gruppe - OS.
_OX_GROUP: Dict[str, int] = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}

# JT-STARK: ungleiche e_g-Besetzung.  (d-Zahl, low_spin) -> True.
# d4 HS (t2g3 eg1) | d7 LS (t2g6 eg1) | d9 (t2g6 eg3, ein Loch).
# d9 ist spinunabhaengig -- es gibt nur eine Besetzung.
_OX_JT_STRONG = {(4, False), (7, True), (9, True), (9, False)}


def _ox_donor_charge(atom, metal_idx: int) -> Optional[int]:
    """Ionische Ladung EINES Donoratoms.  None = nicht klassifizierbar -> Abbruch.

    Klassifiziert wird ueber das Element, die schweren Nicht-Metall-Nachbarn und die
    H-Zahl -- nicht ueber die Formalladung, die hier gerade unbrauchbar ist.
    """
    sym = atom.GetSymbol()
    nh = atom.GetTotalNumHs()
    heavy = [n for n in atom.GetNeighbors()
             if n.GetSymbol() != "H" and n.GetIdx() != metal_idx]
    nheavy = len(heavy)

    if sym in ("F", "Cl", "Br", "I"):
        return -1 if (nheavy == 0 and nh == 0) else None      # Halogenid; sonst organisch gebunden

    if sym == "O":
        if nheavy == 0 and nh == 0:
            return -2                                          # terminales Oxo
        if nheavy == 0 and nh == 1:
            return -1                                          # Hydroxid
        if nheavy == 0 and nh == 2:
            return 0                                           # Aqua
        if nheavy == 1 and nh == 0:
            return -1                                          # Alkoxid / Carboxylat-O / Phenolat
        if nheavy == 2 and nh == 0:
            return 0                                           # Ether
        return None

    if sym == "S":
        if nheavy == 0 and nh == 0:
            return -2                                          # Sulfido
        if nheavy == 1 and nh == 0:
            return -1                                          # Thiolat
        if nheavy == 2 and nh == 0:
            return 0                                           # Thioether
        return None

    if sym == "N":
        if nheavy == 0 and nh == 0:
            return -3                                          # Nitrido
        if atom.GetIsAromatic():
            return 0                                           # Pyridin / Imidazol -- L-Typ
        # ⚠ AMID GEGEN AMMIN -- der Selbsttest hat es aufgedeckt (Cisplatin las sich als
        # Pt(IV) d6 statt Pt(II) d8).  Ein metallgebundenes NH3 und ein NR2-Amid haben im
        # Graphen BEIDE "zwei Reste neben dem Metall"; die Bindung zum Metall verbraucht
        # eine Valenz.  Unterschieden wird am WASSERSTOFF: traegt das N ein H, ist es in
        # diesen SMILES praktisch immer ein neutrales Amin/Ammin -- ein deprotoniertes Amid
        # wuerde mit expliziter Ladung geschrieben.  Nur ein N mit ZWEI schweren Resten und
        # KEINEM H ist ein Amid.
        if nh > 0:
            return 0                                           # Ammin / Amin / Imin -- L-Typ
        if nheavy == 1:
            return -2                                          # Imido
        if nheavy == 2:
            return -1                                          # Amid
        if nheavy == 3:
            return 0                                           # tertiaeres Amin -- L-Typ
        return None

    if sym == "P":
        if nheavy + nh == 3:
            return 0                                           # Phosphan -- L-Typ
        if nheavy + nh == 2:
            return -1                                          # Phosphid
        return None

    if sym == "C":
        # Carbonyl (C mit terminalem O) und Isonitril sind L-Typ; Alkyl/Aryl sind X-Typ.
        for n in heavy:
            if n.GetSymbol() == "O" and len([q for q in n.GetNeighbors()
                                             if q.GetSymbol() != "H"]) == 1:
                return 0                                       # CO
        if atom.GetIsAromatic() or nheavy + nh >= 1:
            return -1                                          # Alkyl / Aryl / Cyclopentadienyl-C
        return None

    return None                                                # unbekannter Donor -> Abbruch


def oxidation_state(mol, metal_idx: int) -> Optional[Tuple[int, int]]:
    """(Oxidationsstufe, d-Elektronenzahl) oder None.

    None heisst AUSDRUECKLICH "nicht bestimmbar", nicht "null" -- der Aufrufer muss dann auf
    sein heutiges Verhalten zurueckfallen.  Gruende fuer None: mehrere Metalle, unbekannter
    Donor, unbekanntes Metall, unphysikalisches Ergebnis.
    """
    if mol is None:
        return None
    try:
        m = mol.GetAtomWithIdx(int(metal_idx))
        msym = m.GetSymbol()
        if msym not in _OX_GROUP:
            return None
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP]
        if len(metals) != 1:
            return None                       # Verteilung auf mehrere Zentren ist nicht eindeutig
        total = int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
        lig = 0
        for n in m.GetNeighbors():
            q = _ox_donor_charge(n, int(metal_idx))
            if q is None:
                return None                   # EIN unbekannter Donor genuegt -- nichts raten
            lig += q
        os_ = total - lig
        if os_ < 0 or os_ > 8:
            return None                       # unphysikalisch -> die Bilanz stimmt nicht
        d = _OX_GROUP[msym] - os_
        if d < 0 or d > 10:
            return None
        return os_, d
    except Exception as exc:
        _LOG.debug("oxidation-state: nicht bestimmbar (%s)", exc)
        return None




# ⚠ EINE JT-ABFRAGE STEHT HIER BEWUSST NOCH NICHT.  Sie haette heute keinen Aufrufer -- die
# JT-Enumeration (drei Streckachsen je Oktaeder) ist nicht gebaut.  Der Reichweiten-Waechter
# hat genau das gemeldet, und er hat recht: unerreichbarer Code kann nichts aendern, und ihn
# trotzdem zu schreiben ist derselbe Fehler, den dieses Projekt heute dreimal gefunden hat.
# `_OX_JT_STRONG` oben traegt das Wissen als DATEN; die Abfrage kommt mit ihrem Verbraucher.


if __name__ == "__main__":  # pragma: no cover -- Selbsttest an bekannten Faellen
    # Jeder Fall ist von Hand nachgerechnet; die erwartete Stufe steht aus der Chemie fest,
    # nicht aus einem frueheren Lauf dieses Moduls.  Ein Selbsttest, der nur die eigene
    # Ausgabe von gestern bestaetigt, prueft nichts.
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    CASES = [
        ("CC(=O)C=C(C)[O][Cu][O]C(C)=CC(C)=O", 2, 9, "Cu(acac)2 -- d9, JT-aktiv"),
        ("Cl[Pt](Cl)(N)N", 2, 8, "cis-Platin -- d8"),
        ("[O][Re](Cl)(Cl)Cl", 5, 2, "Re(V) oxo-trichlorid"),
        ("[Ni](C#[O])(C#[O])(C#[O])C#[O]", 0, 10, "Ni(CO)4 -- d10, alle L-Typ"),
        ("Cl[Fe](Cl)Cl", 3, 5, "FeCl3 -- d5"),
    ]
    ok = bad = 0
    for smi, exp_os, exp_d, note in CASES:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                except Exception:
                    pass
        if mol is None:
            print(f"  ?? {smi}  -- SMILES nicht lesbar ({note})"); bad += 1; continue
        mi = next((a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP), None)
        r = oxidation_state(mol, mi) if mi is not None else None
        if r is None:
            print(f"  -- {smi}\n     None (nicht bestimmbar) -- erwartet OS {exp_os}, d{exp_d}  [{note}]")
            bad += 1
            continue
        got_os, got_d = r
        mark = "OK" if (got_os, got_d) == (exp_os, exp_d) else "!!"
        ok += (mark == "OK"); bad += (mark != "OK")
        print(f"  {mark} {smi}\n     OS {got_os} d{got_d}   erwartet OS {exp_os} d{exp_d}   [{note}]")
    print(f"\n  {ok} richtig, {bad} nicht.  "
          f"'None' ist KEIN Fehler des Verfahrens -- es ist die Weigerung zu raten.")
