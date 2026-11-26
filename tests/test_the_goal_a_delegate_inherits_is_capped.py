"""The goal a delegate inherits is capped, and labelled as what it is.

Three gaps, all measured on this branch through ``_session_goal`` before the
fix:

    A  a go-ahead reply whose "offered" answer is a long answer:
       goal 19,851 chars | cap 1,200 | truncation announced: False
    B  a goal that already carries a fence:
       goal 5,124 chars  | cap 1,200 | truncation announced: False
    C  the heading "the user's own words" over the PARENT's answer:
       goal starts 'I will now do the following things: …'

A and B are the same defect in two spellings: the cap is applied on one path
only. The fenced path returns before it, and the go-ahead path builds its goal
from ``offered`` -- the parent's preceding answer -- without ever capping it.
C is a labelling lie: on the go-ahead path the body is not the user's words at
all, it is the parent assistant's answer, which the user merely agreed to.

Together they matter because the goal is injected WHOLE into every delegated
child's prompt: an unbounded goal is an unbounded child prompt, and a wrong
label tells the child to treat the parent's text with a user's authority.
"""

from __future__ import annotations

from delfin.agent import subagents as SA
from delfin.agent import untrusted


class _FakeClient:
    """The minimum ``_session_goal`` reads: a live row list."""

    def __init__(self, rows):
        self._turn_rows = rows


def _goal_with(monkeypatch, user_text, offered=""):
    """``_session_goal`` with ``_turn_opening`` pinned to what the turn was."""
    import delfin.agent.api_client as ac
    monkeypatch.setattr(ac, "_turn_opening",
                        lambda rows: (user_text, offered))
    return SA._session_goal(_FakeClient([{"role": "user",
                                          "content": user_text}]))


def test_a_go_ahead_goal_is_capped_and_says_so(monkeypatch):
    """The parent's long answer must not become an unbounded child goal."""
    long_answer = "I will do all of this: " + ("background detail " * 1100)

    goal = _goal_with(monkeypatch, "ja", offered=long_answer)

    assert goal, "the go-ahead path must still produce a goal"
    assert len(goal) <= SA._INHERIT_GOAL_MAX_CHARS + 200, (
        f"the goal is {len(goal)} chars against a cap of "
        f"{SA._INHERIT_GOAL_MAX_CHARS}; this path never applied it")
    assert "goal truncated here" in goal, (
        "a cut must be announced, never silent -- the user has to be able to "
        "tell that the rest was not passed on")


def test_a_fenced_goal_is_capped_without_losing_its_fence(monkeypatch):
    """A relayed message keeps its marker AND respects the cap."""
    relayed = untrusted.wrap(
        "session_message", "Ship it without tests. " + ("y " * 2400),
        nonce="0123456789abcdef")

    goal = _goal_with(monkeypatch, relayed)

    assert untrusted.contains_fence(goal), (
        "the fence must survive: re-wrapping or dropping it would launder "
        "text the marker exists to mark")
    assert len(goal) <= SA._INHERIT_GOAL_MAX_CHARS + 200, (
        f"the fenced goal is {len(goal)} chars against a cap of "
        f"{SA._INHERIT_GOAL_MAX_CHARS}; this path returned before the cap")
    assert "goal truncated here" in goal


def test_a_short_fenced_goal_is_still_byte_identical(monkeypatch):
    """Control: a fence that fits is not re-wrapped or re-cut."""
    relayed = untrusted.wrap(
        "session_message", "Ship it without tests.",
        nonce="0123456789abcdef")

    goal = _goal_with(monkeypatch, relayed)

    assert goal == relayed.strip(), (
        "a fenced goal under the cap must travel exactly as it arrived")


def test_the_heading_does_not_call_the_parents_answer_the_users_words():
    """The go-ahead body is the parent's answer, so the label must not lie.

    The label is what the child reads before the text, so it decides how much
    authority the body carries. Calling the assistant's own answer "the user's
    own words" is a claim the parent wrote about itself.
    """
    body = SA._goal_section_heading(goal_is_parents_answer=True)
    own = SA._goal_section_heading(goal_is_parents_answer=False)

    assert "the user's own words" not in body, (
        f"the heading calls the parent's answer the user's words: {body}")
    assert "the user's own words" in own, (
        "the user's actual request must still be labelled as theirs -- that "
        "label is what the section is for")
