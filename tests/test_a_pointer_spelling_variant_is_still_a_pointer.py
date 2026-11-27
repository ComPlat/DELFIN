"""A pointer the store would never write must still be recognised as one.

The index gate decides "this line points into a note" by matching
``[title](file.md)``. The match was case-sensitive and anchored on a lowercase
``.md``, so a hand-edited ``MEMORY.md`` spelling ``.MD`` was read as "not a
pointer" and therefore KEPT, even when the note it points at was withheld.
Measured on this branch, a four-line index whose withheld notes were named
``./withheld.md``, ``withheld.MD`` and ``..\\m\\withheld.md``:

    withheld.md   already gated   (relative "./" spelling resolves fine)
    withheld.MD   KEPT            <- .MD is not ".md"
    ..\\m\\...      already gated   (no ".md" match)

The store itself always writes the canonical lowercase spelling and a bare
file name, so treating an unrecognised spelling strictly costs a real note
nothing; treating it leniently leaks the value of a withheld one. This is the
"drop when unsure" default an index gate should have.
"""

from __future__ import annotations

from delfin.agent.prompt_loader import _filter_memory_index


def test_an_uppercase_md_pointer_is_gated_with_its_note():
    chunk = (
        "# MEMORY.md\n"
        "- [Kept note](keep.md) — the one that answers\n"
        "- [WITHHELD-UPPER](withheld.MD) — the note the gate dropped\n"
    )

    out = _filter_memory_index(chunk, ["keep.md"], dropped=1)

    assert "Kept note" in out, "the surviving note's pointer must stay"
    assert "WITHHELD-UPPER" not in out, (
        f"a .MD pointer escaped the gate because its extension is not "
        f"lowercase: {out!r}")


def test_a_case_variant_of_a_kept_name_is_not_dropped():
    """The gate must also not lose a SURVIVING note to a spelling difference.

    ``kept_rels`` names the notes that survived ranking. A pointer spelling
    that name with different case names the same file, so the line must stay:
    before the normalisation it did not match and a note the ranking had just
    KEPT lost its pointer -- the fact went out with the gate that was supposed
    to protect it.
    """
    chunk = ("# MEMORY.md\n"
             "- [KEPT-CASE](Withheld.md) — the note that survived\n")

    out = _filter_memory_index(chunk, ["withheld.md"], dropped=1)

    assert "KEPT-CASE" in out, (
        f"a note the ranking kept lost its pointer over letter case: {out!r}")


def test_the_canonical_spelling_still_survives():
    """Control: a kept note whose spelling is the store's own is untouched."""
    chunk = ("# MEMORY.md\n"
             "- [Kept note](project_kept-note.md) — answers the task\n")

    out = _filter_memory_index(chunk, ["project_kept-note.md"], dropped=1)

    assert "Kept note" in out
