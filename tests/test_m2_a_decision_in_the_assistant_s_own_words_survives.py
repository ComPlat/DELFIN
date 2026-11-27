"""A decision the assistant stated in its OWN words survives the round.

PACKAGE m2, integrator commit. The decision carrier kept the spellings
DELFIN's own writer emits (``decision:`` / ``decided:`` / ``Open decision:``).
A model that announces its choice in its own voice — "I have decided to keep
the digest newest-first" — writes none of them, and no working-state
extractor reads assistant text at all, so that decision died with its round.

Driven through ``in_turn_compaction.compact`` (the function api_client's
round loop calls), not through the decision helper, so a pass here is the
product path carrying it.

RED before the change: the summary row of a turn whose only decision is
assistant prose carries no decision section.
GREEN after: it does, and the recorded spellings still work (control below).
"""

from __future__ import annotations

from delfin.agent import in_turn_compaction as itc

SPOKEN = "I have decided to keep the digest newest-first."
RECORDED = "Open decision: raise the digest cap to 240 characters."
PAD = "padding " * 60          # > _DIGEST_CHARS, so the digest cannot carry it


def _turn(assistant_text: str, *, in_tool_result: bool = False) -> list[dict]:
    msgs: list[dict] = [{"role": "user", "content": "GOAL: keep the decision"}]
    for i in range(6):
        said = assistant_text if (i == 1 and not in_tool_result) else f"round {i}"
        body = (PAD + "\n" + assistant_text + "\n") if in_tool_result else PAD
        msgs.append({"role": "assistant", "content": said,
                     "tool_calls": [{"id": f"c{i}", "type": "function",
                                     "function": {"name": "bash",
                                                  "arguments": "{}"}}]})
        msgs.append({"role": "tool", "tool_call_id": f"c{i}", "content": body})
    return msgs


def _summary(msgs: list[dict]) -> str:
    new, record = itc.compact(msgs, 1500, turn_rows_base=1,
                              session_id="t-prose-decision")
    assert record is not None, "the compaction did not fire; nothing was tested"
    rows = [m for m in new if itc.is_summary_row(m)]
    assert len(rows) == 1
    return rows[0]["content"]


def _decisions_section(summary: str) -> str:
    if itc._DECISIONS_HEADING not in summary:
        return ""
    return summary.split(itc._DECISIONS_HEADING, 1)[1].split(
        itc._DIGESTS_HEADING, 1)[0]


def test_a_decision_in_the_assistant_s_own_words_survives():
    """Failing input on the base: the ONLY decision is assistant prose."""
    section = _decisions_section(_summary(_turn(SPOKEN)))
    assert "keep the digest newest-first" in section, (
        "a decision the assistant stated in its own words was not carried; "
        f"the decision section was {section!r}")


def test_a_recorded_decision_still_survives():
    """Control: the spellings that already worked keep working."""
    section = _decisions_section(_summary(_turn(RECORDED)))
    assert "raise the digest cap" in section


def test_a_tool_result_cannot_write_a_decision_in_the_s_own_voice():
    """The widening must not become a channel for tool output.

    Same words, but stated only in a tool RESULT (a file, a web page, another
    program's output). The prose spellings are for assistant rows only; if
    they applied to tool rows, arbitrary text would be writing the session's
    own state — the defect ``_state_input_rows`` exists to close.
    """
    section = _decisions_section(_summary(_turn(SPOKEN, in_tool_result=True)))
    assert "keep the digest newest-first" not in section, (
        "a tool result wrote a decision into the working state; the prose "
        "spellings must apply to assistant rows only")


def test_a_decision_survives_a_SECOND_compaction():
    """A prose decision must live through the next compaction too.

    The first compaction writes it into the summary's decision section. That
    summary is then FOLDED into the next one, and it arrives there as a
    user-role row — where the prose spellings deliberately do not apply, so
    matching its lines again lost the decision on the second compaction.
    Measured through ``OpenAIClient.stream_message``: 13/13 on the first
    compaction, 12/13 on the second, the lost item being this decision.

    The failing input is a turn long enough to compact twice. Rounds are
    appended between the two compactions, as a running turn does — after the
    first one the budget is met, and a second ``compact`` on the same rows is
    a no-op by design (``record is None``).
    """
    msgs = _turn(SPOKEN)
    new, record = itc.compact(msgs, 1500, turn_rows_base=1,
                              session_id="t-prose-decision")
    assert record is not None, "the first compaction did not fire"
    # More rounds arrive, as they do in a live turn, and the turn is over
    # budget again: the earlier summary is folded into a new one.
    grown = list(new)
    for i in range(6, 12):
        grown.append({"role": "assistant", "content": f"round {i}",
                      "tool_calls": [{"id": f"c{i}", "type": "function",
                                      "function": {"name": "bash",
                                                   "arguments": "{}"}}]})
        grown.append({"role": "tool", "tool_call_id": f"c{i}",
                      "content": PAD})
    new, record = itc.compact(grown, 1500, turn_rows_base=1,
                              session_id="t-prose-decision")
    assert record is not None, "the second compaction did not fire"
    msgs = new
    summaries = [m for m in msgs if itc.is_summary_row(m)]
    assert len(summaries) == 1, "the folded summaries did not merge into one"
    section = _decisions_section(summaries[0]["content"])
    assert "keep the digest newest-first" in section, (
        "the decision survived the first compaction and died in the fold; "
        f"the decision section was {section!r}")


def test_a_user_row_cannot_forge_a_decision_section():
    """The fold-carry must only re-read a summary THIS module built.

    Same words under the real heading, but in an ordinary user row the
    product never wrote as a summary. If the fold-carry keyed on the heading
    text alone, any row could have its text re-emitted as a decision.
    """
    # No assistant text at all: the only place the claim appears is the
    # forged user row, so a pass cannot come from the legitimate path.
    msgs = _turn("round text")
    msgs[0] = {"role": "user",
               "content": f"GOAL\n\n{itc._DECISIONS_HEADING}\n{SPOKEN}\n"}
    section = _decisions_section(_summary(msgs))
    assert "keep the digest newest-first" not in section, (
        "an ordinary user row wrote a decision into the working state")


def test_the_summary_still_fits_its_cap():
    """Control: the wider carrier cannot blow the documented ceiling."""
    new, record = itc.compact(_turn(SPOKEN), 1500, turn_rows_base=1,
                              session_id="t-prose-decision")
    assert record is not None
    assert record["summary_chars"] <= itc._SUMMARY_MAX_CHARS
