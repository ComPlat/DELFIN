"""Delegating held the session's turn, so the user could not reach it.

The capability was already there and unused. `subagent(background=true)`
returns at once -- the schema has said so for a while -- and the
dashboard keeps the message box live during a turn and drains the queue
when the turn ends. But nothing told the agent to use it, so a
delegation blocked, the turn ran for as long as the delegate did, and
everything the user typed waited. With the caps now at an hour that is a
very long silence.

So this is a prompt change, and what it is worth depends entirely on
being IN the prompt: a rule nobody reads is not a rule. Asserted on the
shipped file, which is what gets injected, and on the composed prompt.

Also asserted: the queue the rule depends on. The box stays enabled
during a turn and the queue is drained afterwards -- if either changed,
the rule would be advice that does not pay off.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import _DOC_TOOLS_OPENAI
from delfin.agent.prompt_loader import PromptLoader

_PACK = Path(SA.__file__).resolve().parent / "pack" / "agents"


def _role_prompt() -> str:
    return (_PACK / "solo_agent.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# The rule is in the prompt
# ---------------------------------------------------------------------------

def test_the_prompt_tells_the_agent_to_delegate_in_the_background():
    text = _role_prompt()
    assert "background=true" in text
    assert "END the turn" in text or "end your turn" in text.lower()


def test_it_says_why_blocking_costs_the_user():
    """Without the reason the rule reads as a style preference, and the
    model trades it away the first time blocking looks simpler."""
    text = _role_prompt()
    i = text.index("delegate in the BACKGROUND")
    arm = text[i:i + 600]
    assert "cannot reach" in arm
    assert "waits" in arm


def test_it_names_how_to_collect_and_how_to_talk():
    text = _role_prompt()
    i = text.index("delegate in the BACKGROUND")
    arm = text[i:i + 600]
    assert "subagent_result" in arm
    assert "subagent_message" in arm


def test_it_still_allows_blocking_where_blocking_is_right():
    """A rule with no exception gets broken rather than applied: a short
    run whose answer the next step needs should block."""
    text = _role_prompt()
    i = text.index("delegate in the BACKGROUND")
    arm = text[i:i + 600]
    assert "Block only" in arm


def test_the_rule_survives_into_the_composed_prompt(tmp_path):
    """The role prompt is stripped of lazy modules before it is sent. A
    rule inside a module that the task text does not trigger would be
    absent exactly when it is needed."""
    prompt = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="solo",
        route=["solo_agent"], role_index=0,
        task_text="look at this repository")
    assert "background=true" in prompt


# ---------------------------------------------------------------------------
# The capability the rule leans on
# ---------------------------------------------------------------------------

def test_the_tool_can_return_at_once():
    schema = next(t["function"] for t in _DOC_TOOLS_OPENAI
                  if t["function"]["name"] == "subagent")
    assert "background" in schema["parameters"]["properties"]


def test_the_message_box_stays_live_during_a_turn():
    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    i = src.index("def _update_button_states")
    body = src[i:i + 900]
    assert "input_textarea.disabled = False" in body
    assert "send_btn.disabled = False" in body


def test_what_the_user_typed_is_processed_when_the_turn_ends():
    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    assert "_process_queue()" in src
    i = src.index("def _process_queue")
    body = src[i:i + 500]
    assert 'state["message_queue"]' in body
    assert 'not state["streaming"]' in body


def test_a_queued_message_is_visible_to_the_user():
    """Queued and silent would be indistinguishable from ignored."""
    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    i = src.index("def _update_queue_display")
    assert "queued" in src[i:src.index("def _refresh_context_bar", i)]


# ---------------------------------------------------------------------------
# And the parallel guidance still holds
# ---------------------------------------------------------------------------

def test_the_pool_width_is_still_stated():
    """It was trimmed to pay for the rule above; the number that changes
    what a fan-out costs had to survive that."""
    text = _role_prompt()
    assert "4 workers" in text


def test_parallel_delegation_is_still_one_message():
    text = _role_prompt()
    assert "ONE assistant message" in text


def test_a_queued_message_enters_the_chat_once_when_it_is_sent():
    """Queued: listed at the foot of the chat with its text, not put in the
    transcript -- the send that drains the queue shows it there, and doing
    both showed every queued line twice."""
    from pathlib import Path
    src = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    i = src.index('state["message_queue"].append(user_text)')
    branch = src[i:src.index("return", i)]
    assert "_append_chat_message" not in branch
    assert "_update_queue_display()" in branch
    j = src.index("def _update_queue_display")
    body = src[j:src.index("def _refresh_context_bar", j)]
    assert "delfin-agent-queue-row" in body and "_html.escape(" in body
