"""A file the agent names has to arrive with its name intact.

The sanitiser removes harmony special-token leftovers -- the words
``json_schema``, ``json`` and ``constrain`` that a mis-decoded tool
channel drops into the visible text. It recognised them by "the token,
followed by anything that is not a space or a letter", and a file
extension is exactly that: ``settings.json`` inside backticks, before a
period, or before a comma matched the rule and lost four characters.

Observed 2026-09-08 in a live answer that named the file the user had
asked about eleven times and wrote ``.delfin/settings. `` every time.
The user cannot open that path, and nothing in the answer says why it
looks like that.

The direction of the trade is the point: a leftover the rule no longer
strips is one stray word in prose, and an extension it strips is a path
that does not exist. So the follow-set is narrowed to what a leftover
actually sits against -- a payload brace, a special-token bar, or the
glitch characters it arrived with -- and never to ordinary punctuation.
"""
import pytest

from delfin.agent.text_sanitize import sanitize_agent_text


# The shape the rule exists for, from the report that produced it:
# a decoded special token glued to the JSON payload of a tool call.
_LEFTOVER = 'json_schema{"doc_id":"orca","section_id":"6.15"}'


@pytest.mark.parametrize("written", [
    "`settings.json`",
    "`.delfin/settings.json`",
    "Trag es in .delfin/settings.json ein.",
    "Dateien: bookmarks.json, package.json und tsconfig.json.",
    "siehe `package.json`",
    "(bookmarks.json)",
    '"settings.json"',
])
def test_a_named_json_file_keeps_its_extension(written):
    out = sanitize_agent_text(written).text
    for name in ("settings.json", "bookmarks.json", "package.json",
                 "tsconfig.json"):
        if name in written:
            assert name in out, (
                f"{name} lost its extension in {written!r}: {out!r}")


def test_the_leftover_the_rule_exists_for_is_still_removed():
    out = sanitize_agent_text(f"Antwort {_LEFTOVER} Ende").text
    assert "json_schema" not in out


def test_a_leftover_before_its_payload_is_still_removed():
    """The other decoded token, in the same position: against the brace
    that opens the arguments it was announcing."""
    out = sanitize_agent_text('Antwort constrain{"a":1} Ende').text
    assert "constrain" not in out


def test_the_word_in_a_sentence_is_not_machinery():
    """Under-stripping is the safe side of this rule, and prose is where
    it shows: a user asking for JSON gets the word back."""
    out = sanitize_agent_text("Gib mir das als json, bitte.").text
    assert "json" in out
