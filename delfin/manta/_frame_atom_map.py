"""_frame_atom_map.py — Atomzuordnung FRAME -> MOL fuer den FF-freien Pfad.

DAS PROBLEM, GEMESSEN AM 16.08.2026.  Die post-hoc-Korrektoren des Bauers brauchen
`mol`-ATOMINDIZES und wenden sie auf XYZ-KOORDINATEN an.  Auf dem legacy-Pfad stimmt das,
weil das XYZ aus demselben `mol` stammt.  Auf dem FF-FREIEN Pfad nicht: dort liegt das
Metall auf 0, danach ein `AddHs`-Block je Ligand in BAUREIHENFOLGE
(`converter_backend._heteroleptic_block_offsets`, `_config_block_offsets`) -- eine voellig
andere Ordnung als RDKits.

WAS DAS ANRICHTET.  `_ffree_shared_tail` warnt woertlich davor; derselbe Bruch hat dort schon
den Ring-Pucker-Emitter zum Nullhebel gemacht (185 von 187 Systemen byte-identisch).  Und die
Wache der Planarisierer prueft nur die ANZAHL der Atome: bei gleicher Anzahl und vertauschter
Ordnung haetten sie die FALSCHEN Atome verflacht -- keine stille Wirkungslosigkeit, sondern
stille Zerstoerung.  (Seit 16.08. prueft dort ein Riegel die Reihenfolge, `864cba3f`.)

⚠ WARUM NICHT DIE ZUORDNUNG MITFUEHREN.  Naheliegend waere, die zur BAUZEIT bekannten Offsets
mit dem Frame zu transportieren.  Zwei Gruende dagegen: die XYZ-Kommentarzeile traegt bereits
das LABEL ("BEBGUL frame0 alt-bind-C"), und Labels werden stromabwaerts geparst
(`_arrangement_key`, Anordnungsfamilien) -- dort etwas anzuhaengen gefaehrdet sie.  Und die
Tupelbreite `(xyz, label)` ist an vielen Stellen angenommen.

DER WEG HIER: REKONSTRUIEREN, NICHT RATEN.  Aus der Geometrie des Frames wird ein Graph
gebaut (dieselbe Adjazenz, die die Korrektoren ohnehin benutzen), daraus ein RDKit-Molekuel
mit GENERISCHEN Bindungen, und dann sucht RDKit die Untergraph-Uebereinstimmung.  Ein
VOLLSTAENDIGER Substruktur-Treffer ueber alle Atome IST die Isomorphie -- er erhaelt Elemente
und Konnektivitaet per Definition.  Das ist keine Heuristik mit Rueckfallkette.

⚠⚠ GENERISCHE BINDUNGEN SIND PFLICHT.  Der Frame-Graph kennt nur "gebunden ja/nein" (er kommt
aus Abstaenden), `mol` kennt Einfach/Doppel/aromatisch/dativ.  Ein Abgleich mit Bindungsordnung
wuerde IMMER scheitern -- und zwar still, als "keine Zuordnung".  `makeBondsGeneric` hebt das
auf.

IM ZWEIFEL NICHTS.  Kein Treffer, mehrdeutiger Treffer, ungleiche Atomzahl oder fehlendes
RDKit -> `None`.  Der Aufrufer faellt dann auf sein heutiges Verhalten zurueck, und das ist
byte-identisch.  Eine FALSCHE Zuordnung waere schlimmer als keine: sie liesse die Korrektoren
genau die falschen Atome anfassen -- der Fehler, den der Riegel vom 16.08. abfaengt.
"""
from __future__ import annotations

import logging
from typing import List, Optional, Sequence

_LOG = logging.getLogger(__name__)


def frame_to_mol_map(mol, syms: Sequence[str], nbrs: List[List[int]]) -> Optional[List[int]]:
    """Zuordnung FRAME-Index -> MOL-Index, oder None.

    ``nbrs`` ist die Adjazenz des Frames (aus `_build_geometric_adjacency`), also genau der
    Graph, mit dem die Korrektoren ohnehin arbeiten -- keine zweite Bindungsdefinition.

    Rueckgabe: ``m[frame_idx] = mol_idx``.  None heisst AUSDRUECKLICH "nicht bestimmbar".
    """
    if mol is None:
        return None
    try:
        from rdkit import Chem
    except Exception:
        return None
    n = len(syms)
    if n == 0 or mol.GetNumAtoms() != n:
        return None
    # Frame-Graph als Molekuel nachbauen -- Elemente und Kanten, keine Ordnungen.
    try:
        rw = Chem.RWMol()
        for s in syms:
            a = Chem.Atom(s)
            a.SetNoImplicit(True)
            rw.AddAtom(a)
        seen = set()
        for i in range(n):
            for j in nbrs[i]:
                if j <= i or (i, j) in seen:
                    continue
                seen.add((i, j))
                rw.AddBond(i, j, Chem.BondType.SINGLE)
        frame = rw.GetMol()
    except Exception as exc:
        _LOG.debug("frame-map: Frame-Graph nicht baubar: %s", exc)
        return None
    # Generische Bindungen auf BEIDEN Seiten: der Frame kennt keine Ordnungen, `mol` schon.
    # Ohne das scheitert der Abgleich immer -- und zwar still.
    try:
        params = Chem.AdjustQueryParameters.NoAdjustments()
        params.makeBondsGeneric = True
        params.makeDummiesQueries = False
        q = Chem.AdjustQueryProperties(frame, params)
        target = Chem.AdjustQueryProperties(Chem.Mol(mol), params)
    except Exception as exc:
        _LOG.debug("frame-map: Bindungen nicht generisch zu machen: %s", exc)
        return None
    try:
        matches = target.GetSubstructMatches(q, useChirality=False, uniquify=False,
                                             maxMatches=2)
    except Exception as exc:
        _LOG.debug("frame-map: Substruktursuche fehlgeschlagen: %s", exc)
        return None
    if not matches:
        return None
    m = list(matches[0])
    if len(m) != n:
        return None
    # Elementgleichheit ist durch die Suche garantiert; hier wird sie trotzdem geprueft.  Eine
    # Zuordnung, die man nicht nachrechnet, ist eine Behauptung -- und genau daran ist heute
    # schon einmal eine Messung gescheitert.
    try:
        for fi, mi in enumerate(m):
            if mol.GetAtomWithIdx(int(mi)).GetSymbol() != syms[fi]:
                return None
    except Exception:
        return None
    return m


if __name__ == "__main__":  # pragma: no cover -- Selbsttest
    # ZWEI FAELLE, und der zweite ist der eigentliche:
    #   (a) gleiche Reihenfolge  -> die Zuordnung MUSS die Identitaet sein
    #   (b) VERTAUSCHTE Reihenfolge -> sie muss die Vertauschung zurueckgeben
    # Faellt (b) durch, ist das Modul wertlos: genau dafuer existiert es.
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from delfin.manta._coord_angle_corrector import _build_geometric_adjacency
    import numpy as _np

    smi = "CC(=O)C=C(C)[O][Cu][O]C(C)=CC(C)=O"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print("  SMILES nicht lesbar"); raise SystemExit(1)
    molH = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    if AllChem.EmbedMolecule(molH, randomSeed=1) != 0:
        print("  Einbettung fehlgeschlagen"); raise SystemExit(1)
    conf = molH.GetConformer()
    syms = [a.GetSymbol() for a in molH.GetAtoms()]
    pts = _np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                      conf.GetAtomPosition(i).z] for i in range(molH.GetNumAtoms())])
    nbrs, _ = _build_geometric_adjacency(syms, pts)

    m = frame_to_mol_map(molH, syms, nbrs)
    ident = (m == list(range(len(syms)))) if m else False
    print(f"  (a) gleiche Ordnung  -> {'IDENTITAET' if ident else m if m else 'None'}")

    # (b) Frame umsortieren: Metall nach vorn, wie es der FF-freie Bauer tut.
    mi = next(i for i, s in enumerate(syms) if s == "Cu")
    order = [mi] + [i for i in range(len(syms)) if i != mi]
    inv = {old: new for new, old in enumerate(order)}
    syms2 = [syms[o] for o in order]
    nbrs2 = [[inv[k] for k in nbrs[o]] for o in order]
    m2 = frame_to_mol_map(molH, syms2, nbrs2)
    ok2 = bool(m2) and all(molH.GetAtomWithIdx(m2[f]).GetSymbol() == syms2[f]
                           for f in range(len(syms2)))
    print(f"  (b) Metall-auf-0    -> {'ZUORDNUNG GEFUNDEN' if ok2 else (m2 if m2 else 'None')}")
    if m2:
        print(f"      Frame 0 ({syms2[0]}) -> mol {m2[0]} ({molH.GetAtomWithIdx(m2[0]).GetSymbol()})")
