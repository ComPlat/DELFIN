"""Which role may run which tool, at the layer that can actually refuse.

``engine._ROLE_TOOL_WHITELIST`` reads like enforcement — "If a role emits
a tool_use event for a tool NOT in its whitelist, the engine silently
blocks it. This prevents prompt-injection or model hallucination from
bypassing role restrictions." It cannot do that. The client yields the
``tool_use`` event and then executes the tool inside the same generator
(``api_client`` ~16431); by the time the engine sees the event the call
has been made. ``continue`` there does not stop anything — it skips the
UI callback, the turn's tool counter, the execution ledger the
functional-claim guard reads for evidence, the trace, and the stray-write
check. A call that RAN, hidden from every record of it.

Every coding tool on the OpenAI-compatible backends arrives as
``mcp__kit-coding__…`` and took the namespace exemption anyway, so on the
backend the KIT models are served from the whitelist did not even reach
that far.

The layer that can refuse is ``_tool_denied_for_role``, checked inside
``_DocToolExecutor.execute`` before the tool runs, and it covers the two
roles someone thought to add. This file adds the reviewer roles — the
ones whose own declared tool set has neither Edit nor Write, which is a
statement of intent that nothing was keeping.

A DENY list, not an allow list. The declared sets are written in the CLI
backend's vocabulary (Read/Grep/Glob/Bash) and the executor's surface is
sixty tools wide; translating them would refuse most of what these roles
legitimately do. Naming the writes refuses exactly what was already
declared out of bounds, and leaves everything else working.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _tool_denied_for_role
from delfin.agent.engine import _ROLE_TOOL_WHITELIST

# Roles whose declared set contains neither Edit nor Write.
_READ_ONLY_ROLES = ("critic_agent", "reviewer_agent", "chief_agent",
                    "session_manager", "runtime_agent", "research_agent")

# The executor's own spellings. The CLI backend's names (Write, Edit, …)
# are deliberately absent: that backend runs its tools in its own
# subprocess and DELFIN sees the call only after it happened, so no gate
# on this side could refuse one. Naming them would be a promise this
# layer cannot keep.
_WRITE_TOOLS = ("write_file", "edit_file", "multi_edit", "apply_patch",
                "notebook_edit")


@pytest.mark.parametrize("role", _READ_ONLY_ROLES)
def test_the_declared_set_really_has_no_write(role):
    """The premise, read from the engine rather than assumed: these roles
    were always meant to be unable to write."""
    declared = _ROLE_TOOL_WHITELIST[role]
    assert "Write" not in declared and "Edit" not in declared


@pytest.mark.parametrize("role", _READ_ONLY_ROLES)
@pytest.mark.parametrize("tool", _WRITE_TOOLS)
def test_a_reviewer_cannot_write(role, tool):
    assert _tool_denied_for_role(role, tool), f"{role} may call {tool}"


@pytest.mark.parametrize("role", _READ_ONLY_ROLES)
@pytest.mark.parametrize("tool", _WRITE_TOOLS)
def test_the_namespace_is_not_the_way_around_it(role, tool):
    """The whole reason the engine-side check missed: on the KIT and
    Ollama backends every coding tool arrives namespaced."""
    assert _tool_denied_for_role(role, f"mcp__kit-coding__{tool}")


@pytest.mark.parametrize("role", _READ_ONLY_ROLES)
@pytest.mark.parametrize("tool", ["read_file", "grep_file", "list_files",
                                  "bash", "search_docs", "report_verdict"])
def test_reading_and_reporting_are_untouched(role, tool):
    """A reviewer that cannot read is not a narrower reviewer, it is a
    broken one. The deny list names writes and nothing else."""
    assert not _tool_denied_for_role(role, tool), f"{role} lost {tool}"


@pytest.mark.parametrize("role", ["builder_agent", "test_agent",
                                  "solo_agent"])
@pytest.mark.parametrize("tool", _WRITE_TOOLS)
def test_the_roles_that_build_still_build(role, tool):
    assert not _tool_denied_for_role(role, tool), f"{role} lost {tool}"


def test_an_unknown_role_is_not_denied_anything():
    """Deny-by-default belongs to roles that declared a list. A role with
    neither list is unrestricted here and gated by permissions alone —
    unchanged, and the reason this file adds a deny list rather than an
    allow list."""
    assert not _tool_denied_for_role("", "write_file")
    assert not _tool_denied_for_role("some_future_agent", "write_file")
