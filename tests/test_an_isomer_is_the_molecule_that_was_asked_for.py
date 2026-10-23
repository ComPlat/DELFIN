"""Every isomer is the molecule that was asked for.

An isomer is the same molecule laid out differently: same atoms, same count,
same order as the template the SMILES describes. The pipeline depends on
that last part — ``post_optimize_geometry`` asks for "a single conformer
matching xyz atom-for-atom" in its own signature, and the topology gate
reads ``coords[i]`` for template atom i. A permuted order does not make
those checks fail loudly, it makes them measure the wrong pair.

Measured on the production runtime, Python 3.11 / RDKit 2025.09.6:

    [Fe+2](O)(O)(O)(O)(O)(O)O     4 isomers, all 15 atoms, order correct
    [Zr+4](O)(O)...O             58 isomers, all 17 atoms, order correct
    [Fe+2]1(N)(N)(Cl[Fe+2]...)   71 isomers, all 17 atoms, order correct

These tests were written on 2026-08-14 as xfail, describing a defect: three
of four Fe(OH)7 isomers came back as the bare coordination skeleton with
every hydrogen gone, and the full-size ones in a different atom order. That
measurement was taken on Python 3.13, which is this developer's shell and
not what DELFIN runs on. On 3.11 none of it reproduces. The defect was an
artifact of the measuring environment and the report of it was wrong.

So the markers are off and these are ordinary tests now, which is what the
repository's own convention asks for — the neighbouring xfails in
test_decompose_kekulize_split carry the note "py3.11 CI status unverified —
un-xfail once confirmed green there", and it is now confirmed.

They will fail on 3.13. That is worth knowing when a run is red locally and
green in CI, and it is not a reason to weaken them: the contract they state
is the one the pipeline relies on, and it holds where the code runs.
"""


from __future__ import annotations

import pytest

from rdkit import Chem

from delfin.smiles_converter import (
    smiles_to_xyz_isomers,
    _normalize_metal_smiles,
)

# (name, smiles, seconds for one build measured 2026-08-14)
_CASES = [
    pytest.param(
        "Fe(OH)7 CN=7", "[Fe+2](O)(O)(O)(O)(O)(O)O",
        id="Fe(OH)7 CN=7",
    ),
    pytest.param(
        "Zr(OH)8 CN=8", "[Zr+4](O)(O)(O)(O)(O)(O)(O)O",
        id="Zr(OH)8 CN=8",
    ),
    # 31 s: the nightly lane, not the merge gate.
    pytest.param(
        "Fe2 (mu-Cl)2 bimetallic", "[Fe+2]1(N)(N)(Cl[Fe+2](N)(N)Cl1)Cl",
        id="Fe2 (mu-Cl)2 bimetallic",
        marks=pytest.mark.slow,
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
