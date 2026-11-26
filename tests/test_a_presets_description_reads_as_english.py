"""A preset's description is read twice, and both readings must be grammatical.

Two surfaces carry a preset's own words:

  * the routing table in the built prompt, which renders ONE cue per preset
    (``subagents._situation_cue``, cut at ``_CUE_LIMIT``);
  * the ``/subagents`` listing, which shows the description WHOLE.

A reworded description can miss the first surface entirely and still matter on
the second: ``method-reviewer``'s cue is
``"Use before comparing results from two or more calculations …"`` on both
sides of the change, byte for byte, because the edited words fall past the
64-character cut. Measured 2026-11-26, and independently observed in the real
runs, where both sides quoted that truncated row and both said the role's
*what* was only inferable from the preset name.

What this file pins, therefore, is the surface the change actually reaches --
the whole description -- and the two defects it repaired:

  * ``method-reviewer`` ended on a dangling participle ("ranked or
    interpreted" with nothing they attach to);
  * ``verifier`` spent its second sentence restating its first.

Both are red on the merge base and green after the rewording, which is why the
assertions name the wording rather than a property of the renderer.
"""
from __future__ import annotations

from delfin.agent import subagents as SA


def _description(name: str) -> str:
    return SA.SUBAGENT_PRESETS[name].description


def test_a_presets_dangling_clause_names_its_own_subject():
    """`ranked or interpreted` attaches to nothing; `ranking them` does.

    The clause is what a reader of ``/subagents`` sees at the end of the
    sentence, so an unattached participle there is a sentence that stops
    making sense halfway.
    """
    description = _description("method-reviewer")
    assert "ranking them or interpreting them" in description, description
    assert "ranked or interpreted" not in description, description


def test_a_presets_second_sentence_does_not_repeat_the_first():
    """`verifier` said "refutes a claim by the same check" twice in one cell.

    The table is charged on every turn, so a duplicate sentence spends
    characters for no new information -- and hides the restriction that
    actually discriminates it ("needs refuting rather than confirming").
    """
    description = _description("verifier")
    assert "needs refuting rather than confirming" in description, description
    first, *rest = [s.strip() for s in description.split(".") if s.strip()]
    repeated = [s for s in rest if s.lower() == first.lower()]
    assert not repeated, f"the description repeats its first sentence: {repeated}"
