"""A SMILES that describes two systems should come back as two systems.

``A.B`` means two molecules that are not bonded to each other.  Handed to a
converter whole they are embedded into one coordinate set and land inside one
another, which is a picture nobody can work in.
"""

from __future__ import annotations

import math

import pytest

from delfin.dashboard import separate_systems as separate

TWO = ("c1ccc(-c2c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)"
       "c2-c2ccccc2)cc1.c1ccccc1")


def test_the_dot_that_separates_molecules_is_the_one_at_the_top():
    """A dot inside brackets is part of an atom and one inside parentheses is
    still the same molecule, so the nesting is counted rather than the string
    split."""
    assert separate.split_smiles("CCO") == ["CCO"]
    assert separate.split_smiles(TWO) == [
        "c1ccc(-c2c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)"
        "c2-c2ccccc2)cc1",
        "c1ccccc1",
    ]
    assert separate.split_smiles("[Na+].[Cl-]") == ["[Na+]", "[Cl-]"]
    assert len(separate.split_smiles("[Cd-3].CC.O")) == 3
    assert separate.split_smiles("C(C.C)C") == ["C(C.C)C"], (
        "a dot inside parentheses does not separate anything"
    )
    assert separate.split_smiles("") == []

    assert separate.has_separate_systems(TWO) is True
    assert separate.has_separate_systems("CCO") is False


def test_they_are_set_beside_each_other_rather_than_inside():
    """Two spheres that just contain them, set apart by the gap."""
    first = "3\na\nC 0 0 0\nC 1.5 0 0\nC 0 1.5 0\n"
    second = "2\nb\nO 0 0 0\nO 1.2 0 0\n"

    row = separate.place_beside([first, second])

    assert row.splitlines()[0] == "5"
    assert separate.closest_approach([first, second]) < 1.0, (
        "on top of each other before"
    )
    # after placing, measured between the two groups of the row
    lines = [l for l in row.splitlines()[2:] if l.strip()]
    assert separate.closest_approach(
        ["\n".join(lines[:3]), "\n".join(lines[3:])]) >= separate.GAP - 0.01


def test_a_single_system_is_left_exactly_as_the_converter_made_it():
    """There is nothing to place beside anything, and rebuilding it would only
    move it."""
    def convert(smiles):
        return "1\nx\nHe 1.0 2.0 3.0\n", 1, "quick", None

    out = separate.convert_each("[He]", convert)

    assert out["parts"] == 1
    assert out["xyz"] == "1\nx\nHe 1.0 2.0 3.0\n"


def test_a_part_that_cannot_be_built_says_which_one():
    def convert(smiles):
        if smiles == "bad":
            return None, 0, None, "no such element"
        return "1\nx\nHe 0 0 0\n", 1, "quick", None

    out = separate.convert_each("[He].bad", convert)

    assert out["ok"] is False
    assert "part 2 of 2" in out["error"]
    assert "no such element" in out["error"]


def test_one_arrangement_beside_many_does_not_multiply_into_many():
    """A metal complex may come back as twelve arrangements and the counter-ion
    beside it as one.  Every pairing would be twelve frames of the same
    counter-ion, which is what a combinatorial product buys here: nothing."""
    many = [(f"1\niso{i}\nH {i} 0 0\n", f"label{i}") for i in range(12)]
    one = [("1\nonly\nHe 0 0 0\n", "counter-ion")]

    frames = separate.combine_isomers([many, one])

    assert len(frames) == 12, "the part with the most arrangements leads"
    assert [f[1] for f in frames] == [f"label{i}" for i in range(12)], (
        "and the labels come from it, not from the one that stands still"
    )
    for frame in frames:
        rows = [l for l in frame[0].splitlines()[2:] if l.strip()]
        assert len(rows) == 2, "both systems are in every frame"


# ---------------------------------------------------------------------------
# the molecule the whole thing started from
# ---------------------------------------------------------------------------
@pytest.mark.skipif(pytest.importorskip("rdkit") is None, reason="needs rdkit")
def test_hexaphenylbenzene_and_benzene_stop_being_inside_each_other():
    """Measured before: the benzene came out inside the other molecule, 0.877 A
    at the closest and the two centres 2.22 A apart."""
    from delfin.dashboard.input_processing import smiles_to_xyz_quick

    together, atoms, _method, error = smiles_to_xyz_quick(TWO)
    assert error is None and atoms == 84

    built = separate.convert_each(TWO, smiles_to_xyz_quick)

    assert built["ok"] is True, built["error"]
    assert built["parts"] == 2
    assert built["atoms"] == 84, "no atom is lost by building them apart"

    rows = [l for l in built["xyz"].splitlines()[2:] if l.strip()]
    big, small = "\n".join(rows[:72]), "\n".join(rows[72:])
    assert separate.closest_approach([big, small]) > 3.0, (
        "the two systems have to be far enough apart to be two systems"
    )


def test_both_conversion_paths_build_the_parts_apart():
    """Quick and the isomer manifold alike -- CONVERT, CONVERT + UFF and MANTA
    all come through the second one."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    worker = source.split("def _start_smiles_conversion")[1].split("\n    def ")[0]
    assert "_separate.has_separate_systems(cleaned_data)" in worker
    assert "_separate.convert_each(" in worker, "the quick path"
    assert "_separate.combine_isomers(per_part)" in worker, "the manifold path"
    assert "if not error and isomers and not separate:" in worker, (
        "a hapticity preview made from the whole string would describe a "
        "molecule that is not any of the frames"
    )
