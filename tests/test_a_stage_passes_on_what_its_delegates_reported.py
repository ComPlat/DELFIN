"""A delegate's finding could not reach the others, and nothing relayed it.

`subagent_message` gave a delegate a way to reach the session running it.
Two gaps were left.

During an ORCHESTRATION the session is inside one tool call, so it has no
round boundary of its own to drain at: a message raised in stage one sat
in the queue until the whole orchestration returned, and then arrived in
the session's next round with no stage to belong to. Stages now drain at
their own barrier, and a later stage can read the result with
``{{messages:NAME}}`` exactly as it reads ``{{stage:NAME}}``.

And a session that learned something every delegate needed had to say it
one call at a time -- four delegates, four of its rounds, for one
sentence. ``to=all`` sends it once.

**What is deliberately NOT built: a channel between siblings.** Within a
stage the calls are independent by construction -- that independence is
what makes the fan-out safe, and what lets writers run in separate
worktrees. A side channel would let them create an ordering dependency
the spec does not declare. The failure mode is on record in this
installation's own field reports: three sessions that could talk directly
assigned each other roles and froze interfaces nobody had asked for.
Where one delegate's finding must reach another, that is a stage boundary
or a relay, and in both the session decides what crosses.

So the information flow between delegates is one-way, declared, and
visible: up to the session, then out again if the session sends it.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import _doc_executor


class _Perms:
    def __init__(self, subagent_id=""):
        if subagent_id:
            self.subagent_id = subagent_id


@pytest.fixture(autouse=True)
def _clean():
    SA._NOTES.clear()
    SA._REPLIES.clear()
    yield
    SA._NOTES.clear()
    SA._REPLIES.clear()


def _call(args, perms):
    return json.loads(_doc_executor._execute_subagent_message(args, perms))


# ---------------------------------------------------------------------------
# A later stage reads what the earlier one raised
# ---------------------------------------------------------------------------

def test_a_stage_message_reaches_the_next_stage():
    out = SA._substitute_stage_refs(
        "earlier said: {{messages:one}}",
        {"one": [{"result": 1}]},
        {"one": [{"from": "aa11", "text": "the interface does not exist"}]})
    assert "the interface does not exist" in out
    assert "aa11" in out


def test_results_and_messages_are_separate_tokens():
    """A result is what a delegate was asked for; a message is what it
    raised on its own. Merging them would also change what the
    verification votes run over."""
    out = SA._substitute_stage_refs(
        "R={{stage:one}} M={{messages:one}}",
        {"one": [{"result": 1}]},
        {"one": [{"from": "aa11", "text": "hi"}]})
    assert 'R=[{"result":1}]' in out
    assert '"text":"hi"' in out.split("M=")[1]


def test_an_absent_message_token_is_left_alone():
    """Same barrier semantics as the results token: a name that has not
    run is not substituted."""
    out = SA._substitute_stage_refs("{{messages:later}}", {"one": []},
                                    {"one": []})
    assert out == "{{messages:later}}"


def test_a_spec_that_asks_for_neither_is_untouched():
    text = "plain prompt with no placeholders"
    assert SA._substitute_stage_refs(text, {"one": [{"r": 1}]},
                                     {"one": [{"text": "x"}]}) == text


def test_messages_are_optional_for_older_callers():
    """The third argument defaults, so a caller that predates it still
    substitutes results."""
    assert SA._substitute_stage_refs("{{stage:one}}", {"one": [1]}) == "[1]"


def test_unrenderable_payloads_do_not_break_a_prompt():
    out = SA._substitute_stage_refs("{{messages:one}}", {},
                                    {"one": [{"text": object()}]})
    assert out == "[]"


# ---------------------------------------------------------------------------
# The stage drains at its own barrier
# ---------------------------------------------------------------------------

def test_the_stage_loop_drains_at_the_barrier():
    """It has to be the stage, not the session's round: during an
    orchestration the session is inside one tool call."""
    import inspect

    src = inspect.getsource(SA.run_orchestration)
    i_run = src.index("stage_results[name] = _run_stage_calls")
    i_drain = src.index("take_replies()")
    i_append = src.index("stage_names.append(name)")
    assert i_run < i_drain < i_append, (
        "messages must be drained after the stage runs and before the "
        "stage is marked complete")


def test_messages_are_not_mixed_into_the_results():
    """stage_results is what the verification votes run over and what a
    caller indexes per call."""
    import inspect

    src = inspect.getsource(SA.run_orchestration)
    assert "stage_messages[name] =" in src
    i = src.index("stage_messages[name] =")
    assert "stage_results[name]" not in src[i:i + 200]


def test_the_report_carries_them_even_with_one_stage():
    """A single-stage spec has no later stage to read the token, and a
    message only a template could have surfaced would be lost."""
    import inspect

    src = inspect.getsource(SA.run_orchestration)
    i = src.index('"stages": stage_results,')
    assert '"messages": stage_messages,' in src[i:i + 500]


# ---------------------------------------------------------------------------
# One call reaches every delegate
# ---------------------------------------------------------------------------

@pytest.fixture
def running(monkeypatch):
    entries = {"aa11": {"type": "explore", "description": "one"},
               "bb22": {"type": "plan", "description": "two"},
               "cc33": {"type": "verifier", "description": "three"}}
    monkeypatch.setattr(SA, "read_running", lambda **kw: dict(entries))
    return entries


def test_a_session_reaches_every_delegate_in_one_call(running):
    out = _call({"to": "all", "message": "the schema changed"}, _Perms())
    assert out["status"] == "sent"
    assert sorted(out["to"]) == ["aa11", "bb22", "cc33"]
    for sa_id in running:
        assert SA.notes_waiting(sa_id) == 1


@pytest.mark.parametrize("word", ["all", "ALL", "*", "everyone"])
def test_the_spellings_a_model_reaches_for(running, word):
    assert _call({"to": word, "message": "x"}, _Perms())["status"] == "sent"


def test_a_delegate_already_at_its_limit_is_named_not_hidden(running):
    for i in range(SA._MAX_NOTES):
        SA.post_note("bb22", f"n{i}")
    out = _call({"to": "all", "message": "one more"}, _Perms())
    assert sorted(out["to"]) == ["aa11", "cc33"]
    assert out["not_delivered"] == ["bb22"]


def test_a_broadcast_with_nothing_to_say_sends_nothing(running):
    assert "error" in _call({"to": "all", "message": "  "}, _Perms())
    assert all(SA.notes_waiting(s) == 0 for s in running)


def test_a_broadcast_with_no_delegates_says_so(monkeypatch):
    monkeypatch.setattr(SA, "read_running", lambda **kw: {})
    assert "error" in _call({"to": "all", "message": "x"}, _Perms())


def test_a_delegate_cannot_broadcast(running):
    """`to` is ignored for a delegate, so "all" is not a way around the
    rule that it has exactly one correspondent."""
    out = _call({"to": "all", "message": "listen everyone"}, _Perms("aa11"))
    assert out["status"] == "sent"
    assert SA.take_replies() == [("aa11", "listen everyone")]
    assert all(SA.notes_waiting(s) == 0 for s in running)


# ---------------------------------------------------------------------------
# The delegate is told the route
# ---------------------------------------------------------------------------

def test_the_rule_says_it_cannot_address_a_sibling():
    rule = SA._CHANNEL_RULE.lower()
    assert "cannot address another delegate" in rule


def test_the_rule_names_the_route_that_does_work():
    """A prohibition with no alternative gets worked around."""
    rule = SA._CHANNEL_RULE.lower()
    assert "tell the session" in rule
    assert "decides what crosses" in rule


def test_the_reason_is_recorded_where_the_flow_is_defined():
    """Not in a commit message only: the next person to want a sibling
    channel reads this function."""
    doc = SA._substitute_stage_refs.__doc__ or ""
    assert "BETWEEN siblings" in doc
    assert "ordering dependency" in doc
