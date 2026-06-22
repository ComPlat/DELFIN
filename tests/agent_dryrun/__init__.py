"""Dry-run capability harness for the dashboard agent.

The harness spawns the actual ``claude`` CLI subprocess in plan-mode
(no mutations), feeds it the dashboard-agent system prompt, sends a
realistic user query via the stream-json input, and captures every
tool-use block the model produces. Tests then assert on which tool
the agent picked and with what arguments — without running anything.

Why subprocess instead of the Anthropic SDK directly?

1. The user's actual setup uses the CLI via OAuth (no
   ``ANTHROPIC_API_KEY`` env). The CLI is the one source of truth for
   how the dashboard agent will behave in production.
2. The CLI applies our MCP-tool catalog the same way the dashboard
   does, so we test the real surface area.
3. Plan-mode (``--permission-mode plan``) hard-blocks every mutating
   tool — exactly the dry-run semantics the user asked for.
"""

from .runner import (
    AgentDryRunResult,
    ToolCall,
    run_agent_dryrun,
    has_tool_call,
    tool_call_args,
)

__all__ = [
    "AgentDryRunResult",
    "ToolCall",
    "run_agent_dryrun",
    "has_tool_call",
    "tool_call_args",
]
