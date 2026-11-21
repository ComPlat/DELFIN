"""TB1 adversarial: only the EXACT label is a choice — look-alikes are not.

The refusal classifies TEXT: the two labels are the whole allow-list, and one
of them ("Start a new session from it") is the destructive branch (it forks
the held conversation). So the interesting direction is not "does the exact
label work" (three cases cover that) but "does anything that merely RESEMBLES
a label stay refused".

An allow-list, not a deny-list: instead of naming bad strings, this states the
rule the patch implements -- `user_text.lower() == choice.lower()`, whole
message, and the person gets the refusal line back with the text untouched.
Each candidate below differs from a real label in exactly one way a person's
keyboard can produce:

    "wait."                     trailing punctuation
    "Wait!"                     trailing exclamation
    "W\u0430it"                  CYRILLIC a in place of the latin one
    "\uff37ait"                  FULLWIDTH latin W
    "Start a new session..."    the destructive label with a polite tail
    "start a new session"       the destructive label without "from it"

Red on main: no refused state exists, so every one of them is an ordinary
message and the tab mints a session for it -- the text lands in a stranger.
"""
from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

from test_tb1_adversarial_a_refused_session_is_not_replaced_silently import (  # noqa: E402
    SID,
    _chat_text,
    _held_by_a_live_process,
    _open_held_conversation,
    _session_files,
    _tab,
    _type_and_send,
    _write_session,
    env,                                    # noqa: F401 -- the fixture below
)

# (what a person typed, what makes it a look-alike rather than a label)
LOOK_ALIKES = [
    ("wait.", "a trailing period"),
    ("Wait!", "a trailing exclamation"),
    ("W\u0430it", "a CYRILLIC 'a' instead of the latin one"),
    ("\uff37ait", "a FULLWIDTH latin 'W'"),
    ("Start a new session from it, please", "a polite tail"),
    ("start a new session", "the label without 'from it'"),
]


@pytest.mark.parametrize("typed,why", LOOK_ALIKES)
def test_a_look_alike_is_not_a_choice_and_starts_no_session(
        env, typed, why):                    # noqa: F811 -- the imported fixture
    _write_session(env, SID, "held job", "held-turn")
    proc = _held_by_a_live_process(env)
    try:
        before = set(_session_files(env))
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            _open_held_conversation(tab, env)
            assert "held by a live session" in _chat_text(state), (
                "the fixture must reach the refusal; otherwise this case "
                "measures the ordinary send path and proves nothing")

            _type_and_send(tab, typed)

            assert env["created"] == [], (
                f"{why} is NOT the choice label, so no session may be minted "
                f"for it; typed {typed!r}")
            assert set(_session_files(env)) == before, (
                f"{why} must write no session; typed {typed!r}")
            assert state.get("active_session_id") != SID, (
                "a look-alike must not open the held conversation")
            assert tab["input"].value == typed, (
                "a look-alike is refused, so the text stays in the box "
                f"exactly as typed; got {tab['input'].value!r}")
            assert "held by a live session" in _chat_text(state), (
                "the refusal must be shown again so the person can see their "
                f"text was not taken as a choice; chat said: "
                f"{_chat_text(state)[-300:]!r}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=30)
