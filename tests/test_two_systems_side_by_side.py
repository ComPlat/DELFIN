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

    lines = [l for l in row.splitlines() if l.strip()]
    assert len(lines) == 5, "coordinate lines only, no header of its own"
    assert separate.closest_approach([first, second]) < 1.0, (
        "on top of each other before"
    )
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
    # (coordinates, how many atoms, what to call it) -- what every converter
    # here answers with
    many = [(f"H {i} 0 0\n", 1, f"label{i}") for i in range(12)]
    one = [("He 0 0 0\n", 1, "counter-ion")]

    frames = separate.combine_isomers([many, one])

    assert len(frames) == 12, "the part with the most arrangements leads"
    assert all(f[2].startswith("2 systems") for f in frames), (
        "the label says how many systems are in the picture"
    )
    assert [f[2].split(", ", 1)[1] for f in frames] == [
        f"label{i}" for i in range(12)], (
        "and the rest of it comes from the part that leads, not from the one "
        "that stands still"
    )
    for frame in frames:
        rows = [l for l in frame[0].splitlines() if l.strip()]
        assert len(rows) == 2, "both systems are in every frame"
        assert frame[1] == 2, "and the count is the row's own"


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

    rows = [l for l in built["xyz"].splitlines() if l.strip()]
    assert len(rows) == 84, "the block is coordinates and nothing else"
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
    assert "_separate.combine_isomers(per_part)" in worker, (
        "both paths build the parts apart and then set them side by side"
    )
    assert worker.count("_separate.combine_isomers(per_part)") == 2, (
        "the quick path and the manifold path, one each"
    )
    assert "if not error and isomers and not separate:" in worker, (
        "a hapticity preview made from the whole string would describe a "
        "molecule that is not any of the frames"
    )


def test_a_hapto_ligand_is_never_written_with_a_dot():
    """Splitting on the dot cannot take a complex apart.

    Of the 185,549 structures in the batch file DELFIN and MANTA are built on,
    not one contains a dot: every complex is written as a single connected
    SMILES, hapto ligands included, through the ring closures and the
    charge-separated dative bonds.  So a dot means what it says.
    """
    for connected in (
        "CC1=[N+]2NC(C3=CC=CC=N3)=[O+][Cd-3]2([Br])([Br])[N+]2=CC=CC=C12",
        "[I][Cd-3]123[N+](=CC4=CC=CC=[N+]41)CCC[N+]2=CC1=CC=CC=[N+]13",
        "[Cl][Cd-3]12([Cl])[N+]3=C(C=CC=C3C3=[N+]1C1=CC=CC=C1N3C1=CC=CC=C1)"
        "C1=[N+]2C2=CC=CC=C2N1C1=CC=CC=C1",
    ):
        assert "." not in connected
        assert separate.has_separate_systems(connected) is False, (
            "a complex must never be split into pieces"
        )

    source = open(separate.__file__, encoding="utf-8").read()
    assert "185,549" in source, "the count that makes this safe belongs here"


def test_each_part_gets_the_hapticity_previews_of_its_own_ligands(editor=None):
    """They are the ways a ligand can sit on its metal -- a question about the
    part that has the metal, and about no other part."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    worker = source.split("def _start_smiles_conversion")[1].split("\n    def ")[0]

    quick = worker.split("if quick and separate:")[1].split("elif quick:")[0]
    assert "smiles_to_xyz_quick_with_previews(part)" in quick, (
        "the quick path makes them per part"
    )
    assert "combine_isomers(per_part)" in quick

    manifold = worker.split("if separate:")[1].split("else:")[0]
    assert "append_hapto_previews_to_isomers(\n" in manifold
    assert "made, part, include_quick=apply_uff" in manifold, (
        "from the part's own SMILES, not from the whole string"
    )


def test_the_block_carries_no_header_of_its_own():
    """The converters here answer with coordinate lines and the tab writes the
    count and the comment itself.  A whole file handed back instead put a
    second header inside the block, and whoever read the count took those two
    lines for atoms and lost the last two -- two hydrogens off the end of a
    benzene, in a structure that said 24 and had 24.
    """
    first = "2\na\nH 0 0 0\nH 1 0 0\n"
    second = "1\nb\nHe 0 0 0\n"

    row = separate.place_beside([first, second])

    lines = [line for line in row.splitlines() if line.strip()]
    assert len(lines) == 3, "three atoms, three lines, nothing else"
    for line in lines:
        assert len(line.split()) == 4, f"not a coordinate line: {line!r}"
    assert not lines[0].strip().isdigit()

    out = separate.convert_each("[H][H].[He]",
                                lambda s: (first if "H][H" in s else second,
                                           0, "quick", None))
    assert out["atoms"] == 3
    assert len([l for l in out["xyz"].splitlines() if l.strip()]) == 3, (
        "the count and the block have to agree, or atoms fall off the end"
    )


def test_three_systems_work_as_well_as_two():
    """Nothing in this is about there being two of them."""
    blocks = ["1\na\nH 0 0 0\n", "1\nb\nHe 0 0 0\n", "1\nc\nLi 0 0 0\n"]

    row = separate.place_beside(blocks)
    lines = [l for l in row.splitlines() if l.strip()]

    assert len(lines) == 3
    xs = sorted(float(l.split()[1]) for l in lines)
    assert xs[1] - xs[0] >= separate.GAP - 0.01
    assert xs[2] - xs[1] >= separate.GAP - 0.01


@pytest.mark.skipif(pytest.importorskip("rdkit") is None, reason="needs rdkit")
def test_three_benzenes_come_out_as_three_benzenes():
    from delfin.dashboard.input_processing import smiles_to_xyz_quick

    out = separate.convert_each("c1ccccc1.c1ccccc1.c1ccccc1",
                                smiles_to_xyz_quick)

    assert out["ok"] is True and out["parts"] == 3
    assert out["atoms"] == 36, "twelve atoms each, and none lost"
    lines = [l for l in out["xyz"].splitlines() if l.strip()]
    assert len(lines) == 36
    groups = ["\n".join(lines[i * 12:(i + 1) * 12]) for i in range(3)]
    assert separate.closest_approach(groups) > 3.0
