"""What could come next, offered as numbers — and not invented.

A suggestion the model produces costs a call and is a guess. DELFIN
already knows what is open: the agent's own task list, which it fills as
it works. Offering the top of that list costs nothing, is never a guess,
and is wrong only when the agent's own bookkeeping is.

    next
      1  run the full suite on the rebased chain
      2  push the branch and dispatch CI

Typing `2` at the prompt sends the second one. Only a bare number, and
only while an offer stands — "2" is a message in its own right in a
conversation about numbers.

Pending before blocked, because a blocked task is waiting on something
and offering it as the next step is an invitation into a wall.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl as R
from delfin.agent import task_ticker as TT


class _Store:
    def __init__(self, tasks):
        self._tasks = tasks

    def list(self, **kw):
        return list(self._tasks)


@pytest.fixture()
def tasks(monkeypatch):
    def _set(items):
        monkeypatch.setattr(TT, "get_store", lambda ws: _Store(items))
        monkeypatch.setattr(TT, "resolve_session_scope", lambda s: s)
    return _set


def _task(subject, status="pending", seq=1):
    return {"subject": subject, "status": status, "seq": seq, "id": seq}


# -- what is offered --------------------------------------------------------

def test_the_open_work_is_what_is_offered(tasks, tmp_path):
    tasks([_task("run the full suite", seq=1),
           _task("push the branch", seq=2)])
    assert TT.next_steps(tmp_path) == ["run the full suite",
                                       "push the branch"]


def test_nothing_open_offers_nothing(tasks, tmp_path):
    tasks([])
    assert TT.next_steps(tmp_path) == []


def test_a_finished_turn_with_all_done_offers_nothing(tasks, tmp_path):
    tasks([_task("already done", status="completed")])
    assert TT.next_steps(tmp_path) == []


def test_blocked_work_comes_last(tasks, tmp_path):
    """A blocked task is waiting on something; offering it first is an
    invitation into a wall."""
    tasks([_task("waiting on CI", status="blocked", seq=1),
           _task("write the test", status="pending", seq=2)])
    assert TT.next_steps(tmp_path)[0] == "write the test"


def test_it_offers_at_most_three(tasks, tmp_path):
    tasks([_task(f"task {i}", seq=i) for i in range(10)])
    assert len(TT.next_steps(tmp_path)) == 3


def test_a_broken_store_offers_nothing(monkeypatch, tmp_path):
    def _boom(ws):
        raise RuntimeError("no store")
    monkeypatch.setattr(TT, "get_store", _boom)
    assert TT.next_steps(tmp_path) == []


# -- taking one -------------------------------------------------------------

@pytest.fixture()
def agent():
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    said: list[str] = []

    class _T:
        class theme:
            @staticmethod
            def dim(x):
                return x
        @staticmethod
        def chrome(text):
            said.append(text)
    a.transcript = _T()
    a._said = said
    return a


def test_a_number_becomes_the_suggestion(agent):
    agent._next_steps = ["run the full suite", "push the branch"]
    assert agent._expand_next_step("2") == "push the branch"


def test_a_number_outside_the_offer_is_a_message(agent):
    agent._next_steps = ["one thing"]
    assert agent._expand_next_step("7") == "7"


def test_a_number_with_no_offer_is_a_message(agent):
    agent._next_steps = []
    assert agent._expand_next_step("2") == "2", (
        "a bare number is a message in a conversation about numbers")


def test_a_sentence_is_never_expanded(agent):
    agent._next_steps = ["run the full suite"]
    assert agent._expand_next_step("2 files changed") == "2 files changed"


def test_taking_one_says_which(agent):
    agent._next_steps = ["run the full suite"]
    agent._expand_next_step("1")
    assert any("run the full suite" in s for s in agent._said)
