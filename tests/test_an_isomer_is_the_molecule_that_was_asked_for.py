"""Acceptance criterion for the SMILES → 3D construction path.

Every isomer that ``smiles_to_xyz_isomers`` returns is supposed to be the
same molecule as the SMILES it was given, laid out in space a different way.
Today, for complexes whose ligands are small — water, hydroxide, ammine,
chloride — it is not:

    [Fe+2](O)(O)(O)(O)(O)(O)O   4 isomers: three of 8 atoms, one of 15
    [Zr+4](O)(O)...O           57 isomers: 26 of 9 atoms, 31 of 17
    [Fe+2]1(N)(N)(Cl[Fe+2]...) 80 isomers: 41 of 9 atoms, 39 of 17

The short ones are the bare coordination skeleton — ``Fe O O O O O O O`` —
with every hydrogen gone. That is not a conformer of the input; it is a
different chemical species with a different electron count.

The two sizes are not sorted by label. For Fe(OH)7 the three stripped ones
are the ones carrying the chemistry (``capped-octahedral cap-O/O-O``) and
the whole one is an anonymous ``Isomer 1``; for Zr and Fe2 both kinds appear
at both sizes. So a caller cannot tell from the label which it has, and
picking the informative entry is not a way around it.

Nothing downstream repairs it. The dashboard passes the XYZ through and only
counts atoms (input_processing.py), the MANTA CLI writes it to a file, and
GUPPY takes the coordinate lines verbatim as starting geometries. What the
converter returns is what a user gets.

The full-size structures have their own problem: their atom ORDER differs
from the template's. Much of the pipeline indexes coordinates by template
atom index — ``post_optimize_geometry`` states it in its own signature, "a
single conformer matching xyz atom-for-atom" — and
``_verify_topology_from_graph`` does the same. Given the wrong order it
compares a metal-donor bond against two oxygens 4 Å apart.

Between them these two also explain why none of this was visible: the
verifier returns True when the atom count does not match, because it cannot
check. Not-checkable is reported as fine. So the structures missing seven
hydrogens pass the gate, and the one complete enough to be checked is the
only one that fails.

Ir(ppy)2(acac) is deliberately not here: 55 atoms out of 55, order correct.
Its topology failure is a real geometry violation and a separate question.

These tests are xfail because the construction path is being rewritten. They
are the acceptance criterion for that work: when they XPASS, the rewrite has
fixed what this describes, and these markers come off.
"""

from __future__ import annotations

import pytest

from rdkit import Chem

from delfin.smiles_converter import (
    smiles_to_xyz_isomers,
    _normalize_metal_smiles,
)

_REASON = (
    "SMILES->3D construction returns isomers of the bare coordination "
    "skeleton, and full-size ones in a different atom order than the "
    "template; being rewritten"
)

# (name, smiles, seconds for one build measured 2026-08-14)
_CASES = [
    pytest.param(
        "Fe(OH)7 CN=7", "[Fe+2](O)(O)(O)(O)(O)(O)O",
        id="Fe(OH)7 CN=7", marks=pytest.mark.xfail(reason=_REASON, strict=False),
    ),
    pytest.param(
        "Zr(OH)8 CN=8", "[Zr+4](O)(O)(O)(O)(O)(O)(O)O",
        id="Zr(OH)8 CN=8", marks=pytest.mark.xfail(reason=_REASON, strict=False),
    ),
    # 31 s: the nightly lane, not the merge gate.
    pytest.param(
        "Fe2 (mu-Cl)2 bimetallic", "[Fe+2]1(N)(N)(Cl[Fe+2](N)(N)Cl1)Cl",
        id="Fe2 (mu-Cl)2 bimetallic",
        marks=[pytest.mark.slow,
               pytest.mark.xfail(reason=_REASON, strict=False)],
    ),
]


def _template(smi: str):
    """The molecule the SMILES describes, hydrogens included."""
    norm = _normalize_metal_smiles(smi) or smi
    mol = Chem.MolFromSmiles(norm)
    if mol is None:
        pytest.skip(f"{smi!r}: SMILES failed to parse")
    return Chem.AddHs(mol)


def _elements(xyz: str) -> list[str]:
    return [l.split()[0] for l in xyz.strip().splitlines() if l.strip()]


_CACHE: dict = {}


def _build(smi: str):
    """One build per SMILES, shared by the tests below.

    All three read the same structures and none of them mutates anything, so
    building once turns 19.6 s into about 7 for the two gate cases.
    """
    if smi not in _CACHE:
        res, err = smiles_to_xyz_isomers(
            smi, apply_uff=True, deterministic=True,
            collapse_label_variants=True,
        )
        assert err is None, f"smiles_to_xyz_isomers returned error: {err}"
        assert res, "no isomers returned"
        _CACHE[smi] = res
    return _CACHE[smi]


# --------------------------------------------------------------------------
# the molecule
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,smiles", _CASES)
def test_every_isomer_is_the_whole_molecule(name, smiles):
    """A conformer search may move atoms. It may not remove them."""
    tmpl = _template(smiles)
    expected = tmpl.GetNumAtoms()

    short = {}
    for xyz, label in _build(smiles):
        n = len(_elements(xyz))
        if n != expected:
            short.setdefault(n, []).append(label or "<unlabelled>")

    assert not short, (
        f"{name}: isomers came back with atom counts {sorted(short)} instead "
        f"of {expected}. Examples: "
        + "; ".join(f"{n} atoms: {labels[:3]}" for n, labels in short.items())
    )


@pytest.mark.parametrize("name,smiles", _CASES)
def test_no_isomer_loses_its_hydrogens(name, smiles):
    """Said separately from the count because it is the specific thing that
    happens and the specific thing that makes the result a different
    species: Fe(OH)7 arrives as Fe with seven bare oxygens."""
    tmpl = _template(smiles)
    expected_h = sum(1 for a in tmpl.GetAtoms() if a.GetSymbol() == "H")
    if expected_h == 0:
        pytest.skip("nothing to lose")

    offenders = [
        (label or "<unlabelled>", _elements(xyz).count("H"))
        for xyz, label in _build(smiles)
        if _elements(xyz).count("H") != expected_h
    ]

    assert not offenders, (
        f"{name}: expected {expected_h} hydrogens in every isomer, got "
        f"{offenders[:5]}"
    )


# --------------------------------------------------------------------------
# the order
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,smiles", _CASES)
def test_every_isomer_lists_its_atoms_in_template_order(name, smiles):
    """Coordinates are indexed by template atom index across the pipeline —
    post_optimize_geometry says "matching xyz atom-for-atom" in its own
    signature, and the topology gate reads coords[i] for template atom i. A
    permuted order does not make those checks fail loudly; it makes them
    measure the wrong pair."""
    tmpl = _template(smiles)
    expected = [a.GetSymbol() for a in tmpl.GetAtoms()]

    wrong = []
    for xyz, label in _build(smiles):
        got = _elements(xyz)
        if len(got) == len(expected) and got != expected:
            first = next(i for i, (a, b) in enumerate(zip(got, expected))
                         if a != b)
            wrong.append((label or "<unlabelled>", first, got[first],
                          expected[first]))

    assert not wrong, (
        f"{name}: {len(wrong)} isomer(s) in a different atom order than the "
        f"template. First divergence, (label, index, got, expected): "
        f"{wrong[:3]}"
    )
