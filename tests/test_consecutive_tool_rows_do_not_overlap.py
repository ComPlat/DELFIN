"""Tool rows in the agent chat sit close together but never on top of each other.

A negative top margin on a tool row that follows another one pulled a
multi-line command (a heredoc, for instance) up into the row above it.
"""
from __future__ import annotations

import re

from delfin.dashboard.tab_agent import _AGENT_CSS


def _rules(selector: str) -> list[str]:
    pattern = re.escape(selector) + r"\s*\{([^}]*)\}"
    return re.findall(pattern, _AGENT_CSS)


def test_a_tool_row_after_a_tool_row_has_no_negative_margin():
    bodies = _rules(".delfin-chat-tool + .delfin-chat-tool")
    assert bodies, "the rule for consecutive tool rows is gone"
    for body in bodies:
        assert not re.search(r"margin[\w-]*\s*:\s*-", body), body
