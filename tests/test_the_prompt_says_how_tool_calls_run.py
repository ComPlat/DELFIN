"""The solo prompt says how a round's tool calls are run, and nothing more.

Input: the built system prompt of the solo agent.

The OpenAI-compatible loop (api_client) runs the tool calls of one
assistant message one after another, in their order; only `subagent`
calls are fanned out to a pool of four workers (`_fan_out_subagents`).
Bundling independent calls still saves model rounds, which is what the
prompt asks for. Two of its sentences said otherwise:
- "multiple tool_use blocks" is the Anthropic API's term; this loop
  receives `tool_calls`;
- "the runtime runs same-turn calls concurrently" stood in the subagent
  paragraph but read as true of every call.
"""
from __future__ import annotations

import re

from delfin.agent.prompt_loader import PromptLoader


def _prompt() -> str:
    return re.sub(r"\s+", " ", PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo"))


def test_the_prompt_uses_this_loop_s_term_for_a_tool_call():
    assert "tool_use block" not in _prompt()


def test_only_subagent_calls_are_said_to_run_concurrently():
    p = _prompt()
    assert "the runtime runs same-turn calls concurrently" not in p
    assert "runs same-turn subagent calls concurrently" in p
