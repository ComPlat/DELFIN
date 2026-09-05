"""Dead MCP servers added up until a turn would not start.

Discovery asks every configured server three questions — tools,
resources, prompts — one after another, and each can sit on the per-call
deadline of 15 s. A server that is dead still costs the full wait, and
``initialize`` sends its own request first, so one unreachable entry can
hold the pass for twice the deadline before anything moves on. Five of
them, serially, across three questions, is minutes.

Every failure was swallowed (``except Exception: pass``), so the only
symptom was a turn that took a very long time to begin, with nothing
saying why. It is the shape this project keeps removing: the worse the
configuration, the longer the silence.

The per-call deadline cannot see this — it is doing its job each time.
What was missing is a ceiling on the WHOLE pass, and a report of what
the pass did not reach. A discovery that quietly returned fewer tools
than the config declares is indistinguishable from a config with fewer
tools in it.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import mcp_client as MC


class _Slow:
    """A server that takes its time and answers."""

    def __init__(self, name, seconds=0.05, tools=1):
        self.name = name
        self.seconds = seconds
        self.last_error = ""
        self._tools = tools
        self.calls = 0

    def list_tools(self):
        self.calls += 1
        time.sleep(self.seconds)
        return [MC.MCPTool(server=self.name, name=f"t{i}", description="",
                           schema={}) for i in range(self._tools)]

    def list_resources(self):
        time.sleep(self.seconds)
        return []

    def list_prompts(self):
        time.sleep(self.seconds)
        return []


class _Dead:
    def __init__(self, name):
        self.name = name
        self.last_error = "No such file or directory"
        self.calls = 0

    def list_tools(self):
        self.calls += 1
        raise OSError("No such file or directory")

    def list_resources(self):
        raise OSError("No such file or directory")

    def list_prompts(self):
        raise OSError("No such file or directory")


def _registry(**servers):
    return MC.MCPRegistry(servers=dict(servers), workspace=None, loaded=True)


# ---------------------------------------------------------------------------
# The pass has a ceiling
# ---------------------------------------------------------------------------

def test_a_spent_budget_stops_asking():
    reg = _registry(a=_Slow("a", 0.05), b=_Slow("b", 0.05), c=_Slow("c", 0.05))
    budget = MC._DiscoveryBudget(total_s=0.0)   # already spent
    tools = reg.discover_all(budget)
    assert tools == []
    assert reg.servers["a"].calls == 0
    assert len(budget.skipped) == 3


def test_the_servers_before_the_ceiling_still_answer():
    reg = _registry(a=_Slow("a", 0.0, tools=2), b=_Slow("b", 0.0, tools=1))
    tools = reg.discover_all()
    assert len(tools) == 3
    assert reg.last_discovery["skipped"] == []


def test_the_budget_is_shared_across_the_three_questions():
    """Three separate ceilings would let a broken config cost three
    times the ceiling, which is what the ceiling exists to prevent."""
    reg = _registry(a=_Slow("a", 0.0))
    t, r, p = reg.discover_everything()
    assert t and r == [] and p == []
    assert reg.last_discovery["skipped"] == []


# ---------------------------------------------------------------------------
# ...and it says what it did not reach
# ---------------------------------------------------------------------------

def test_a_skipped_server_is_named():
    reg = _registry(a=_Slow("a"), b=_Slow("b"))
    budget = MC._DiscoveryBudget(total_s=0.0)
    reg.discover_all(budget)
    assert budget.skipped == ["a", "b"]


def test_the_owning_call_writes_the_report():
    """A pass that was handed a shared budget must not overwrite the
    report of the pass that owns it — otherwise the three questions each
    clobber the previous one's account of what was skipped."""
    reg = _registry(dead=_Dead("dead"))
    shared = MC._DiscoveryBudget()
    reg.last_discovery = {"marker": "kept"}
    reg.discover_all(shared)
    assert reg.last_discovery == {"marker": "kept"}
    reg.discover_all()
    assert "failed" in reg.last_discovery


def test_a_failing_server_is_named_with_its_reason():
    reg = _registry(dead=_Dead("dead"), ok=_Slow("ok", 0.0))
    tools = reg.discover_all()
    assert len(tools) == 1
    failed = reg.last_discovery["failed"]
    assert failed and "dead" in failed[0]
    assert "No such file" in failed[0]


def test_the_notice_is_empty_when_everything_answered():
    reg = _registry(ok=_Slow("ok", 0.0))
    reg.discover_all()
    assert reg.discovery_notice() == ""


def test_the_notice_names_the_cost_and_the_servers():
    reg = _registry(dead=_Dead("dead"))
    reg.discover_all()
    text = reg.discovery_notice()
    assert "did not complete" in text
    assert "dead" in text


def test_a_skipped_pass_says_the_budget_was_spent():
    reg = _registry(a=_Slow("a"), b=_Slow("b"))
    budget = MC._DiscoveryBudget(total_s=0.0)
    reg.discover_all(budget)
    reg.last_discovery = budget.report()
    text = reg.discovery_notice()
    assert "were not asked" in text
    assert "discovery budget" in text


def test_an_empty_registry_is_not_an_error():
    reg = _registry()
    assert reg.discover_all() == []
    assert reg.discovery_notice() == ""


# ---------------------------------------------------------------------------
# The default ceiling is well under the cost it replaces
# ---------------------------------------------------------------------------

def test_the_ceiling_is_smaller_than_two_dead_servers_cost():
    """One dead stdio server can burn two RPC deadlines — the initialize
    and the call. The whole pass must cost less than that twice over, or
    the ceiling changes nothing for the case it exists for."""
    assert MC._DISCOVERY_BUDGET_S < 4 * MC._RPC_TIMEOUT_S
