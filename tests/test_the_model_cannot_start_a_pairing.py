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
    for text in ("/pair", "ACTION: /pair", "/pair cancel", "/pair revoke all"):
        assert _slash_is_user_only(text), text


def test_pair_from_the_person_is_not_affected():
    # The user's own sends reach _handle_slash_command directly (the send
    # path at tab_agent.py:21289 and the Phone button at :11691) -- the
    # user-only check guards only the model's ACTION lines, so these stay
    # callable by the person.
    from delfin.dashboard.tab_agent import _USER_ONLY_SLASH_COMMANDS
    assert "/pair" in _USER_ONLY_SLASH_COMMANDS
    assert "/sync" in _USER_ONLY_SLASH_COMMANDS  # the existing precedent
