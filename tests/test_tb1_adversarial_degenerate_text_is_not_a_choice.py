"""TB1 adversarial: degenerate input after a refusal must not mint a session.

The task's checklist names the edge cases explicitly -- empty, None, huge,
unicode, concurrent -- and my first ten cases cover none of them. Every one of
those cases sends a *plausible* message. A person, a paste or a stray click
sends something else, and the refusal state is exactly where a degenerate send
is most likely to slip through to the "mint a session" branch.

The invariant under attack is the same one the package exists for: after the
tab was refused a held conversation, NOTHING the person does to the send button
may create a session or move their text out of the input box, unless they
chose one of the two offered labels in so many words.

Red on unchanged `delfin/dashboard/tab_agent.py`: there the tab holds no
session after the refusal, so every send below mints one and the typed text
leaves the box -- the package's defect, seen from the degenerate end.
"""
import pytest

pytest.importorskip("ipywidgets")

from test_tb1_adversarial_a_refused_session_is_not_replaced_silently import (  # noqa: E402
    SID,
    _chat_text,
    _held_by_a_live_process,
    _open_held_conversation,
    _session_files,
    _sessions_holding,
    _tab,
    _type_and_send,
    _write_session,
    env,                                    # noqa: F401 -- the fixture below
)

# (a short id for the report, the text a send carries, what makes it degenerate)
# The id is separate on purpose: a 200 000-character parameter would otherwise
# become the test ID, and pytest prints it in full for every failure.
DEGENERATE = [
    ("empty", "", "an empty string"),
    ("spaces", "   ", "spaces only"),
    ("newlines", "\n\t \n", "whitespace and newlines only"),
    ("huge", "x" * 200_000, "200 000 characters"),
    ("nbsp", "\u00a0wait\u00a0",
     "the Wait label wrapped in non-breaking spaces"),
    ("nul", "Wait\x00", "the Wait label with a NUL appended"),
]


@pytest.mark.parametrize(
    "typed,why", [pytest.param(t, w, id=i) for i, t, w in DEGENERATE])
def test_a_degenerate_send_after_a_refusal_mints_nothing(
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
                f"{why} is not a choice label, so no session may be minted "
                f"for it; typed {len(typed)} chars")
            assert set(_session_files(env)) == before, (
                f"{why} must write no session file")
            assert _sessions_holding(env, "held job") == [f"{SID}.json"], (
                "the held conversation's own file must be the only one "
                f"carrying its text; found {_sessions_holding(env, 'held job')}")
            assert state.get("active_session_id") != SID, (
                f"{why} must not open the held conversation either")
            assert tab["input"].value == typed, (
                f"{why} is refused, so the text stays in the box exactly as "
                f"typed; got {len(tab['input'].value or '')} chars")
            assert "held by a live session" in _chat_text(state), (
                "the refusal must still be on screen so the person can see "
                f"their text was not taken as a choice; chat said: "
                f"{_chat_text(state)[-300:]!r}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(10)
