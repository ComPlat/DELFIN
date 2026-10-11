"""A leaked tool call that follows its own announcement is still a leak.

Wave 14b: deepseek-v4-flash wrote "Let me commit." and then the call as
text. The strict rule (nothing but the block) saw prose and let the turn
end, so three sessions stood still with their work uncommitted.
"""

import pytest

from delfin.agent.text_sanitize import leaked_tool_call_dominates

_CALL = ('<invoke name="bash">\n<parameter name="command">git status'
         '</parameter>\n</invoke>')


def test_an_announced_call_at_the_end_is_asked_for_again():
    got = leaked_tool_call_dominates("Let me commit the two files.\n" + _CALL)
    assert got == ("bash", {"command": "git status"})


def test_two_calls_back_to_back_after_the_announcement():
    got = leaked_tool_call_dominates("Now both.\n" + _CALL + "\n" + _CALL + "\n")
    assert got and got[0] == "bash"


@pytest.mark.parametrize("text", [
    "ok " + _CALL,                                  # same line: an aside
    "Maybe: " + _CALL,
    _CALL + " done",                                # text after the call
    "The model writes:\n" + _CALL + "\nwhich never runs.",
    "Example:\n```\n" + _CALL + "\n```",            # fenced
    "Example:\n```xml\n" + _CALL,                   # fence still open
    "A half <invoke name=\"x\"> and then\n" + _CALL,
])
def test_a_quoted_or_explained_block_stays_prose(text):
    assert leaked_tool_call_dominates(text) is None
