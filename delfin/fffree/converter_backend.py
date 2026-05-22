"""delfin.fffree.converter_backend — adapter wiring the metal-FF-free foundation
into the converter's smiles_to_xyz_isomers contract.

_fffree_isomers(smiles) -> [(xyz_string, label), ...]  or  None (-> legacy fallback).

v1 handles mononuclear Werner complexes with explicit metal-donor bonds, all
monodentate donors, CN 4/5/6 (decompose.py).  Anything else returns None so the
caller falls through to the existing pipeline.  Deterministic.
"""
from __future__ import annotations
from collections import Counter
from typing import List, Optional, Tuple
from rdkit import Chem

from delfin.fffree import decompose as DEC
from delfin.fffree import polya_isomer_count as PIC
from delfin.fffree import assemble_complex as AC

_GEOM_TO_POLYA = {
    "OC-6 octahedron": "octahedron",
    "SP-4 square planar": "square_planar",
    "T-4 tetrahedron": "tetrahedron",
    "TBP-5 trigonal bipyramid": "trigonal_bipyramid",
    "SPY-5 square pyramid": "square_pyramid",
}


# antipodal vertex pairs per geometry (polya vertex ordering) — for universal
# cis/trans/fac/mer classification from the coloring.
_ANTIPODE = {
    "octahedron": {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4},
    "square_planar": {0: 2, 1: 3, 2: 0, 3: 1},
}


def _classify_coloring(geom_key, vertex_elems) -> str:
    """Universal scientific isomer name from per-vertex DONOR ELEMENTS + the
    polyhedron antipode structure (element-based trans-pair analysis, matching the
    project's _classify_isomer_label scheme — no system-specific rules):
      MA4B2 -> cis/trans · MA3B3 -> fac/mer · MA2B2C2 -> all-cis/all-trans/El-trans
      (square-planar MA2B2 -> cis/trans).  Returns '' for single-isomer cases."""
    from collections import Counter
    anti = _ANTIPODE.get(geom_key)
    if anti is None:
        return ""
    cnt = Counter(vertex_elems)
    n = len(vertex_elems)

    def is_trans(el):
        v = [i for i, e in enumerate(vertex_elems) if e == el]
        return any(anti.get(v[a]) == v[b] for a in range(len(v)) for b in range(a + 1, len(v)))

    pairs2 = [el for el, c in cnt.items() if c == 2]
    threes = [el for el, c in cnt.items() if c == 3]
    if n == 6:
        if threes:                                    # MA3B3 / MA3B2C
            return "mer" if is_trans(threes[0]) else "fac"
        if len(pairs2) == 3:                           # MA2B2C2
            trans_els = sorted(el for el in pairs2 if is_trans(el))
            if not trans_els:
                return "all-cis"
            if len(trans_els) == 3:
                return "all-trans"
            return "-".join(f"{e}trans" for e in trans_els)
        if len(pairs2) == 1:                           # MA4B2
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 4:                                       # SP-4 / T-4 MA2B2
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    return ""


def _xyz(syms, P) -> str:
    # HEADER-LESS atom block in the CANONICAL converter format ("{sym:4s}
    # {x:12.6f}..."), byte-identical to every other pool so viewers (Avogadro,
    # etc.) treat fffree output exactly like all other archives.  The pool
    # evaluator prepends "{count}\n{comment}".
    return "\n".join(f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}"
                     for s, (x, y, z) in zip(syms, P))


def _fffree_isomers(smiles: str, max_isomers: int = 50
                    ) -> Optional[List[Tuple[str, str]]]:
    d = DEC.decompose(smiles)
    if d is None:
        return None
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return None
    # ligand identity = canonical SMILES of each fragment; group by it
    lig_label, lig_ref, lab_elem = [], {}, {}
    for lg in d["ligands"]:
        try:
            lab = Chem.MolToSmiles(lg["mol"])
        except Exception:
            return None
        lig_label.append(lab)
        lig_ref.setdefault(lab, (lg["mol"], lg["donor_local_idx"]))
        lab_elem[lab] = lg["donor_elem"]
    spec = dict(Counter(lig_label))
    try:
        colorings = PIC.enumerate_isomers(geom_key, spec)
    except Exception:
        return None
    if not colorings:
        return None
    results: List[Tuple[str, str]] = []
    for k, coloring in enumerate(colorings[:max_isomers]):
        vertex_specs = [lig_ref[lab] for lab in coloring]
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], d["geometry"], vertex_specs)
        except Exception:
            return None
        if built is None:
            return None
        syms, P = built
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        geom_tag = d["geometry"].split()[0]
        label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        results.append((_xyz(syms, P), label))
    # generate-gate-floor: never return zero isomers if the decomposition succeeded
    return results or None


if __name__ == "__main__":
    for label, smi in [("cisplatin", "N[Pt](N)(Cl)Cl"),
                       ("[CoCl3(NH3)3]", "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"),
                       ("hexammineCo", "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]")]:
        r = _fffree_isomers(smi)
        if r is None:
            print(f"{label:<16} -> None (legacy)")
        else:
            print(f"{label:<16} -> {len(r)} isomers: {[lab for _, lab in r]}")
            print("   first xyz head:", r[0][0].splitlines()[0], r[0][0].splitlines()[1])
