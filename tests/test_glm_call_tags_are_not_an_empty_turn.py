"""The answer was there; the tags around it made it unreadable.

Given a tool surface, GLM-5.3 on the KIT deployment wraps its ordinary
text answer in its own call markup, and the serving parser leaves it in
the text channel:

    <tool_call>ACTION: /tab calc</arg_value></tool_call>

There is no JSON inside, so the qwen-family pattern beside it does not
match, and the ACTION line no longer begins a line. Measured 2026-09-07:
the dashboard reported "[empty turn] — 0 characters of reasoning, 0 tool
calls" for a request the model had answered correctly, 0/3 across four
benchmark arms, while the identical prompt sent by hand returned
"ACTION: /tab calc" at every reasoning_effort level.
"""

import pytest

from delfin.agent.text_sanitize import (parse_leaked_tool_calls,
                                        sanitize_agent_text)

_FIELD = ('<tool_call>ACTION: /tab calc</arg_value></tool_call>'
          '<tool_call>ACTION: /done</arg_value></tool_call>')


def test_the_answer_survives_the_tags():
    out = sanitize_agent_text(_FIELD).text
    assert "ACTION: /tab calc" in out
    assert "ACTION: /done" in out
    assert "tool_call" not in out
    assert "arg_value" not in out


def test_two_wrapped_actions_stay_on_separate_lines():
    """The dashboard parser reads one action per line. Joining them would
    trade an empty turn for a mangled one."""
    lines = [ln.strip() for ln in sanitize_agent_text(_FIELD).text.splitlines()
             if ln.strip()]
    assert lines == ["ACTION: /tab calc", "ACTION: /done"]


@pytest.mark.parametrize("tag", ["<tool_call>", "</tool_call>", "<arg_key>",
                                 "</arg_key>", "<arg_value>", "</arg_value>"])
def test_every_tag_of_the_family_is_removed(tag):
    assert tag not in sanitize_agent_text(f"before{tag}after").text


def test_a_real_leaked_call_is_still_a_call_not_prose():
    """The JSON form must keep being recognised and removed whole —
    unwrapping it would print a tool call at the user."""
    raw = '<tool_call>{"name": "read_file", "arguments": {"path": "a.py"}}</tool_call>'
    result = sanitize_agent_text(raw)
    assert result.text.strip() == ""
    assert result.leaked_tools == ["read_file"]
    assert [c["name"] for c in parse_leaked_tool_calls(raw)] == ["read_file"]


def test_ordinary_text_is_untouched():
    for text in ("Alles gut, 3 Dateien.", "Die Datei liegt in calc/x.out.",
                 "ACTION: /tab calc"):
        assert sanitize_agent_text(text).text == text


def test_the_full_field_capture_reads_as_an_action_line():
    """End to end on the exact bytes the endpoint returned."""
    out = sanitize_agent_text(_FIELD).text
    assert out.startswith("ACTION: /tab calc")
