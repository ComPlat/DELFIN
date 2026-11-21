"""TB1 adversarial: the refused state is cleared by the person's own choice.

Target 2 of the package names three clearing routes: choosing a session in the
list, "New session", and the holder ending. My ten earlier cases cover the
second and third. This file covers the first, and the failure a patch can
produce while fixing the defect: a refused state that outlives the refusal and
swallows every later send.

Direction, stated honestly per case:

* ``test_choosing_another_session_clears_the_refusal_so_a_send_lands`` and
  ``test_a_second_load_of_the_same_held_session_still_refuses`` are GUARDS: they
  pass on unchanged ``tab_agent.py``, because no refused state exists yet --
  which is exactly the defect. They earn their place by failing a patch that
  gets the state wrong in the other direction: one that keeps it forever (the
  tab then refuses to send anything again) or one that drops it on every list
  interaction (the second load of a held session opens it). A guard that passes
  both ways measures nothing about the change and is NOT counted as a control.
* ``test_the_clearing_choice_is_the_person_s_not_a_background_tick`` IS a
  control: it is red on unchanged code, because the base defect (a send after
  a refusal mints a session) is present. Its own claim -- that only the person
  clears the state, never a timer -- cannot be measured before a refused state
  exists; after the patch this same case must stay green, which is what makes
  the timer question testable at all.
"""
from __future__ import annotations

import pytest

# The helpers and the fixture live in the main adversarial file; importing them
# keeps one fixture definition instead of a drifting copy.
from test_tb1_adversarial_a_refused_session_is_not_replaced_silently import (  # noqa: E402,F401
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

# A second, unheld session: what the person picks instead.
OTHER = "tb1a000000000002"


def _write_other(env, turn="someone else's turn"):            # noqa: F811
    """A session nobody holds, offered in the same list as the held one."""
    _write_session(env, OTHER, "other job", turn)


def test_choosing_another_session_clears_the_refusal_so_a_send_lands(
        env):                                # noqa: F811 -- the imported fixture
    """Picking a different session must end the refusal, not freeze the tab.

    GUARD (passes on unchanged code; see the module docstring). A patch that
    remembers "this tab was refused" without clearing it on the person's next
    choice turns every later send into a refused one -- the person would be
    unable to use the tab at all. The outcome asserted is the person's, not the
    patch's internals: after choosing another session, typing and sending puts
    the text in a session that is not the held one.
    """
    _write_session(env, SID, "held job", "held-turn")
    _write_other(env)
    proc = _held_by_a_live_process(env)
    try:
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            # The refusal happens first: without it this measures nothing.
            _open_held_conversation(tab, env)
            assert "held by a live session" in _chat_text(state), (
                "the fixture must reach the refusal; otherwise the rest "
                "measures the ordinary load path")

            # The person changes their mind and picks the other session.
            offered = [v for _label, v in tab["session"].options]
            assert OTHER in offered, (
                "the unheld session must be offered in the list the person "
                f"picks from; it offers {offered}")
            tab["session"].value = OTHER
            tab["load"].click()

            # Now a send is the ordinary path again and must reach a session.
            typed = "this belongs in the other session"
            _type_and_send(tab, typed)

            assert _sessions_holding(env, typed), (
                "after choosing another session the send was swallowed: the "
                "refused state outlived the refusal, so the person cannot "
                "write into the tab at all")
            # The held file carries "held-turn" by construction, so the
            # question is not whether that text exists but whether the NEW
            # text landed in the HELD session instead of the chosen one.
            landed = _sessions_holding(env, typed)
            assert f"{SID}.json" not in landed, (
                "the typed text was written into the HELD session, which is "
                f"still held by a live process: {landed}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(timeout=30)


def test_a_second_load_of_the_same_held_session_still_refuses(
        env):                                # noqa: F811 -- the imported fixture
    """Clearing on a choice must not clear the refusal itself.

    GUARD, of the opposite over-reach: if the patch drops the refused state on
    any list interaction, the person who picks the held session twice gets a
    session on the second pick and the original defect returns silently. Both
    loads must refuse, and neither may mint a session.
    """
    _write_session(env, SID, "held job", "held-turn")
    proc = _held_by_a_live_process(env)
    try:
        before = set(_session_files(env))
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            for attempt in ("first", "second"):
                _open_held_conversation(tab, env)
                said = _chat_text(state)
                assert "held by a live session" in said, (
                    f"the {attempt} load of a held session was NOT refused: "
                    "the lock was ignored and the tab opened a conversation "
                    f"another process is writing to. Chat says:\n{said}")
                assert env["created"] == [], (
                    f"the {attempt} load of a held session minted an engine: "
                    "a held conversation was opened behind the refusal")

            after = set(_session_files(env))
            minted = sorted(after - before)
            assert not minted, (
                "loading a held session twice wrote new session files, so the "
                f"person's next send has a stranger to land in: {minted}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(timeout=30)


def test_the_clearing_choice_is_the_person_s_not_a_background_tick(
        env):                                # noqa: F811 -- the imported fixture
    """Nothing clears the refusal on its own while the holder is alive.

    CONTROL: red on unchanged code, where a send after a refusal mints a
    session (the base defect). The claim unique to this case -- that only the
    person clears the state and never a timer -- becomes measurable only once
    a refused state exists, so after the patch this case must stay green.
    A patch that clears the state on a timer would pass a "the choice line
    appeared" test; here time passes with no person acting and the refusal
    must still be in force.
    """
    _write_session(env, SID, "held job", "held-turn")
    proc = _held_by_a_live_process(env)
    try:
        before = set(_session_files(env))
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            _open_held_conversation(tab, env)
            assert "held by a live session" in _chat_text(state)

            # Let whatever the patch polls run. The fixture shortens the
            # interval to 0.2 s, so this is many polls, not a machine-speed
            # race: the bound is on polls elapsing, not on wall-clock luck.
            import time
            deadline = time.monotonic() + 2.0
            while time.monotonic() < deadline:
                time.sleep(0.05)

            typed = "typed while the holder is still alive"
            _type_and_send(tab, typed)

            assert env["created"] == [], (
                "with the holder still alive the tab minted a session after "
                "time passed on its own: the refusal was cleared by a timer "
                "instead of by the person")
            assert not _sessions_holding(env, typed), (
                "the typed text was sent into a session while the held "
                "conversation was still held by a live process")
            after = set(_session_files(env))
            assert not (after - before), (
                "a session file appeared with no person acting: "
                f"{sorted(after - before)}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(timeout=30)
