"""An index trimmed twice must not carry two different counts.

``_filter_memory_index`` appends a count line ("... and N more memories") when
it withholds pointers, and ``_fit_memory_index`` appends one when the index is
too long for its allowance. When both run on the same chunk -- the ordinary
case for a store larger than the budget -- the second one classifies the FIRST
count line as a "loose" fact (a line with no pointer, which the docstring says
must be kept) and emits its own beside it. Measured on this branch:

    allowance=140   count lines: 2
      "... and 9 more memories, ..."     <- stale, and it sorts FIRST
      "- [Note two](n2.md) — two"
      "... and 2 more memories, ..."     <- the live count

Two numbers for one store, and the stale one is placed ahead of a real
pointer, because loose lines are emitted before the ranked pointers. The model
therefore reads a count of 9 for a store that is showing 2 withheld.

This is a docstring/behaviour mismatch as well: the count line is not a fact
the user typed into MEMORY.md, which is what the "loose" branch exists to
protect.
"""

from __future__ import annotations

from delfin.agent.prompt_loader import _fit_memory_index

WITH_COUNT = "\n".join([
    "# MEMORY.md",
    "- [Note one](n1.md) — one",
    "- [Note two](n2.md) — two",
    "- [Note three](n3.md) — three",
    "... and 9 more memories, not listed here — ask about a subject to have "
    "its memory recalled",
])


def test_a_re_trimmed_index_carries_one_count_not_two():
    """At most one count line, and never the STALE number.

    Both allowances are checked because the honest outcome differs: at 80 the
    trim still has to hide a pointer and writes its own live count; at 140
    dropping the stale line frees enough room that nothing is hidden, so the
    right answer is NO count line at all. What is never right is the earlier
    trim's number sitting in the prompt.
    """
    for allowance in (80, 140):
        out = _fit_memory_index(WITH_COUNT, ["n2.md", "n1.md"], allowance)
        counts = [ln for ln in out.splitlines() if "more memories" in ln]
        assert len(counts) <= 1, (
            f"allowance {allowance}: {len(counts)} count lines with different "
            f"numbers: {[c.strip() for c in counts]}")
        assert "and 9 more" not in out, (
            f"allowance {allowance}: the STALE count survived: {out!r}")


def test_a_stale_count_does_not_outrank_a_real_pointer():
    """The live count belongs at the end, after the pointers it counts."""
    out = _fit_memory_index(WITH_COUNT, ["n2.md", "n1.md"], 80)
    lines = [ln for ln in out.splitlines() if ln.strip()]

    first_pointer = next((i for i, ln in enumerate(lines)
                          if "(" in ln and ".md" in ln), None)
    first_count = next((i for i, ln in enumerate(lines)
                        if "more memories" in ln), None)
    assert first_pointer is not None, "a ranked pointer must still be shown"
    assert first_count is not None, "the count line must still say what is hidden"
    assert first_count > first_pointer, (
        "the count line is emitted ahead of the pointers it counts: "
        f"{[ln.strip()[:40] for ln in lines]}")


def test_a_typed_in_fact_is_still_kept():
    """Control: the loose-line rule exists to protect a fact typed into the index."""
    chunk = "\n".join([
        "# MEMORY.md",
        "- [Note one](n1.md) — one",
        "- [Note two](n2.md) — two",
        "- [Note three](n3.md) — three",
        "The deploy window is Tuesday nights.",
    ])
    out = _fit_memory_index(chunk, ["n2.md"], 150)

    assert "The deploy window is Tuesday nights." in out, (
        "a fact typed straight into MEMORY.md must survive the trim")
