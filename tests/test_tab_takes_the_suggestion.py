"""Tab fills in what could come next, without sending it.

The suggestions were already there — drawn from the agent's own open
tasks, so they cost nothing and are never a guess — but taking one meant
reading a number and typing it:

    next
      1  run the full suite on the rebased chain
      2  push the branch and dispatch CI

Typing "2" works and stays. Tab is the key a hand is already on.

On an EMPTY line Tab walks the offer: first press takes the first
suggestion, each further press the next, and one press past the last
empties the line again — a cycle you can always leave, which matters
because the alternative is a prompt you have to erase by hand.

On a line with something in it Tab keeps its old job, completing a
``/command``. That is the whole reason the empty line is the trigger:
the two meanings never meet.

And it FILLS, it does not send. A suggestion is an offer; pressing
return stays the user's.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl as R


@pytest.fixture()
def agent():
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    a._next_steps = ["run the full suite", "push the branch", "write it up"]
    a._tab_index = 0
    return a


def test_the_first_press_takes_the_first_suggestion(agent):
    assert agent._tab_suggestion("") == "run the full suite"


def test_each_press_walks_on(agent):
    first = agent._tab_suggestion("")
    second = agent._tab_suggestion(first)
    third = agent._tab_suggestion(second)
    assert [first, second, third] == [
        "run the full suite", "push the branch", "write it up"]


def test_one_past_the_last_gives_the_line_back(agent):
    text = ""
    for _ in range(3):
        text = agent._tab_suggestion(text)
    assert agent._tab_suggestion(text) == "", (
        "a cycle you cannot leave is a line you have to erase by hand")


def test_a_line_with_something_in_it_is_not_touched(agent):
    assert agent._tab_suggestion("/he") is None, (
        "None means: not mine, let the command completer have it"
    )


def test_nothing_on_offer_means_tab_is_not_mine(agent):
    agent._next_steps = []
    assert agent._tab_suggestion("") is None


def test_it_never_raises_on_a_broken_offer(agent):
    agent._next_steps = None
    assert agent._tab_suggestion("") is None


def test_the_offer_is_forgotten_when_the_user_types(agent):
    """Half a typed word is not the second suggestion."""
    agent._tab_suggestion("")
    assert agent._tab_suggestion("run the full suite and more") is None


def test_tab_fills_and_does_not_send():
    """The key handler puts the text in the buffer; submitting is the
    user's. A Tab that sent the turn would make the offer a trap."""
    import inspect
    src = inspect.getsource(R.TerminalAgent.read_boxed)
    branch = src.split("elif kind == rk.COMPLETE:")[1].split(
        "elif kind == rk.")[0]
    assert "decoder.buffer" in branch, branch
    # Submitting happens by setting submit_text; a Tab that set it would
    # send the offer instead of offering it.
    assert "submit_text" not in branch, branch


def test_the_help_says_so():
    from delfin.agent import repl_commands
    import inspect
    src = inspect.getsource(repl_commands)
    assert "suggestion" in src.lower(), (
        "a key that does something new says so where keys are listed")
