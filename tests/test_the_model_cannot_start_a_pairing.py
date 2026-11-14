"""Pairing starts only when the person asks: /pair or the Phone button.

Seen 2026-11-12: the transcript showed "Asking the messenger gateway for
a one-time invitation ..." twice, without the person asking to pair. The
dashboard's auto-executor runs every ``ACTION: /command`` line of model
output as a slash command -- and ``/pair`` was not on the person-only
list, so one generated line started a pairing (and its later answers
printed "Pairing refused: the asking window is gone."). Pairing is the
person's act, like /perms or /share: never something the model can
start, never a line in the transcript.
"""
from __future__ import annotations

from delfin.dashboard.tab_agent import _slash_is_user_only


def test_pair_from_the_model_is_ignored():
    # The ACTION prefix is stripped by the executor's parser before the
    # check runs (tab_agent.py:_action_cmd), so the guard sees the bare
    # command; the prefixed form is what the model writes, the bare one
    # what the guard receives.
    for text in ("/pair", "/pair cancel", "/pair revoke all", "/pair code"):
        assert _slash_is_user_only(text), text


def test_pair_from_the_person_is_not_affected():
    # The user's own sends reach _handle_slash_command directly (the send
    # path at tab_agent.py:21289 and the Phone button at :11691) -- the
    # user-only check guards only the model's ACTION lines, so these stay
    # callable by the person.
    from delfin.dashboard.tab_agent import _USER_ONLY_SLASH_COMMANDS
    assert "/pair" in _USER_ONLY_SLASH_COMMANDS
    assert "/sync" in _USER_ONLY_SLASH_COMMANDS  # the existing precedent


def test_the_real_executor_chain_drops_an_action_pair_line():
    # The exact pair of functions the auto-executor uses: extraction from
    # the model text (tab_agent.py:_extract_action_commands) then the
    # person-only guard per command. An ACTION: /pair line must extract
    # and then be dropped -- the model writing it starts nothing.
    from delfin.dashboard.tab_agent import _extract_action_commands
    cmds = _extract_action_commands(
        "I will connect your phone.\n\nACTION: /pair\n")
    assert cmds == ["/pair"]
    assert all(_slash_is_user_only(c) for c in cmds)
