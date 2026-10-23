"""delfin.fffree.aromatic_bond_targets — delocalised aromatic bond-length targets.

Aromatic (bond order "aromatic") bond lengths per element pair, in angstrom, for
the fully delocalised (mesomerism-equalised) ring bond — the length at which a
benzene C-C, a pyridine C-N, a furan C-O, a thiophene C-S or a pyrazole N-N
actually sits in real crystal structures, NOT the single-bond covalent sum (~1.52
for C-C) that a covalent-radius prior would give and NOT a Kekule-localised
1.34/1.50 alternation.

PROVENANCE.  The five values are medians over the Crystallography Open Database
(COD, CC0): the SMILES graph mapped onto each crystal structure, median per
element pair of the bonds of bond order "aromatic", over 26 566 organic COD
structures (C-C n=335 129, C-N 22 914, C-O 2 820, C-S 2 576, N-N 2 053).  See
delfin/manta/PROVENANCE.md.

They are deliberately VENDORED here (a plain dict, no external dependency) so the
public DELFIN construction path can seat aromatic rings at the right length with
NO import from any other repository.  Do NOT replace them with a Pyykko
single<->double interpolation (1.393 / 1.333 / ...) or any other covalent-radius
estimate — those optimise the ring to the wrong point.

Universality note: the table is intentionally small (the common aromatic ring-bond
element pairs).  Any pair NOT listed returns ``None`` from :func:`aromatic_target`,
and every caller falls back to its existing generic covalent-sum ideal for that
bond — so an unlisted heteroaromatic pair is simply left at its prior length
(never shortened to a guessed value).
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

# Keyed by tuple(sorted((element_1, element_2))) -> delocalised aromatic length (A).
AROMATIC_BOND_TARGETS: Dict[Tuple[str, str], float] = {
    ("C", "C"): 1.387,
    ("C", "N"): 1.360,
    ("C", "O"): 1.374,
    ("C", "S"): 1.729,
    ("N", "N"): 1.352,
}


def aromatic_target(e1: str, e2: str) -> Optional[float]:
    """Return the delocalised aromatic bond-length target (A) for the unordered
    element pair ``(e1, e2)``, or ``None`` when the pair is not tabulated (the
    caller must then fall back to its generic covalent-sum ideal).  Order of the
    two symbols does not matter."""
    return AROMATIC_BOND_TARGETS.get((min(e1, e2), max(e1, e2)))
