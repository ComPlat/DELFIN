"""A turn that answers only in tool calls must not leave an empty bubble.

When the model emits tool_calls and no text, the chat received an
assistant entry whose entire content was the role name -- an empty box
between the request and the tool line it triggered, repeated for every
tool call in a long run.

What must NOT change: an answer that has text, however short, is an
answer. The rule keys on prose being absent, never on the turn having
called tools, so a reply that both explains and calls a tool keeps its
bubble.
"""

from __future__ import annotations

import pathlib
import re

from delfin.dashboard.tab_agent import _assistant_turn_is_wordless


def test_a_tool_only_turn_is_wordless():
    assert _assistant_turn_is_wordless({"role": "assistant", "content": ""})
    assert _assistant_turn_is_wordless({"role": "assistant", "content": "   \n\t "})
    assert _assistant_turn_is_wordless({"role": "assistant"})
    assert _assistant_turn_is_wordless({"role": "assistant", "content": None})


def test_an_answer_with_text_is_never_wordless():
    for content in (
        "Done.",
        "0",
        "…",
        "Die Reisekostenanträge wurden erstellt.",
    ):
        assert not _assistant_turn_is_wordless(
            {"role": "assistant", "content": content}), (
            f'{content!r} is an answer and must keep its bubble')


def test_a_non_string_body_does_not_raise():
    """Content is not always a str -- structured bodies must not crash the
    renderer, which runs on every token during streaming."""
    assert not _assistant_turn_is_wordless(
        {"role": "assistant", "content": [{"type": "text", "text": "hi"}]})


# ---------------------------------------------------------------------------
# The three places that have to agree
# ---------------------------------------------------------------------------

_SOURCE = pathlib.Path(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
).read_text(encoding="utf-8")


def test_every_site_asks_the_one_definition():
    """Renderer, search rebuild and updater must share the helper.

    Two of them drifting apart is how the empty box survived: the entry was
    stored by one place and drawn by another, each with its own idea of
    what counts as empty.
    """
    uses = _SOURCE.count("_assistant_turn_is_wordless")
    assert uses >= 4, (
        f'only {uses} references to the helper (1 definition + 3 call sites '
        'expected) -- a site is deciding emptiness on its own again')


def test_no_site_hand_rolls_the_emptiness_check():
    """No open-coded variant of the same question next to a role label."""
    hand_rolled = re.findall(
        r'not \(msg\.get\("content"\)[^)]*\)\.strip\(\)', _SOURCE)
    assert not hand_rolled, (
        f'{len(hand_rolled)} hand-rolled emptiness checks remain: {hand_rolled}')


def test_the_updater_keeps_it_out_of_the_stored_history():
    """Skipping it at render time alone would still save it, and a resumed
    session would carry the empty entries back into the conversation."""
    assert "elif not blank:" in _SOURCE, (
        'the updater appends an assistant entry unconditionally again')
    assert re.search(r"if finalize and blank and not had_text:\s*\n\s*#", _SOURCE), (
        'a finalized empty entry is no longer dropped from the history')


def test_an_answer_that_already_has_text_is_never_blanked():
    """A late empty finalize must not wipe a bubble that holds a real answer
    -- with the render rule in place that would make it vanish, not just
    look empty."""
    assert "if not (blank and had_text):" in _SOURCE, (
        'existing assistant text can be overwritten with an empty body')
