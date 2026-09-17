"""«Push your branch. Do not push to main.» is a request for a push.

A push needs the user to have asked for one. The reading of that was
asymmetric: an ask counted where it stood, but a refusal counted for the
whole message -- so a task that said "push your branch" and then, three
sentences later, "do not push to main" granted nothing at all. In the
field on 2026-09-17 three sessions were told exactly that; one of them
gave up and handed the user a git command to run by hand.

Both are now read at the push word they belong to:

  ask survives a later refusal   the branch push was asked for
  a refused push stays refused   "nicht pushen" grants nothing
  a report is not a request      "der push ist fehlgeschlagen"
  German and English alike       the matcher reads both
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _asks_for_push

#: The task text that granted nothing in the field.
FIELD_TASK = (
    "Push your branch so your work exists outside this machine:\n"
    "  git push -u origin HEAD\n"
    "Report the branch name and the commit hash. Do not push to main — "
    "that is refused, and it is not your step."
)


def test_the_task_that_granted_nothing_now_asks_for_a_push():
    assert _asks_for_push(FIELD_TASK) is True


@pytest.mark.parametrize("text", [
    "push it",
    "push the branch, but do not push to main",
    "nicht auf main pushen, sondern den branch pushen",
    "bitte den branch pushen und dann den PR stellen",
    "Push the feature branch. Never push to the default branch.",
])
def test_a_message_that_asks_for_one_grants_one(text):
    assert _asks_for_push(text) is True


@pytest.mark.parametrize("text", [
    "bitte nicht pushen",
    "do not push",
    "don't push anything for now",
    "nie pushen ohne mich zu fragen",
    "no push today",
])
def test_a_message_that_refuses_grants_nothing(text):
    assert _asks_for_push(text) is False


@pytest.mark.parametrize("text", [
    "der push ist fehlgeschlagen",
    "the push went through",
    "mein letzter push war kaputt",
    "that push has failed",
])
def test_talking_about_a_push_is_not_asking_for_one(text):
    assert _asks_for_push(text) is False


def test_an_empty_message_asks_for_nothing():
    assert _asks_for_push("") is False
    assert _asks_for_push(None) is False
