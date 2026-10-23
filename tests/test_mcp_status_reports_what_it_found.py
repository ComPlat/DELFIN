"""The MCP status echoed the config back and called it live.

``_hot_reload_mcp`` reported

    f"{n_servers} server(s), {n_tools} tool(s) live."

``n_tools`` came from ``len(getattr(new_reg, "tools", {}) or {})``.
``MCPRegistry`` is a dataclass with the fields ``servers``, ``workspace``
and ``loaded`` — it has no ``tools`` attribute at all, and
``hasattr(MCPRegistry, "tools")`` is False. The default fired every time,
so the count was **zero unconditionally**: for a healthy server exactly
as much as for a dead one. Tools live on each ``MCPServer`` and are
filled only by ``list_tools()``, which this path never calls.

``n_servers`` counts what ``load()`` read out of the config file. No
subprocess is spawned, no ``initialize`` is sent. So

    /mcp add x /does/not/exist

answered "✅ MCP server 'x' added" and then "2 server(s), 0 tool(s)
live." The word "live" was the only feedback the user got, and it was an
echo of what they had just typed.

``MCPServer.start()`` records the real cause in ``last_error`` and
nothing surfaces it, so a server that failed to start looked exactly
like one that worked — and because the tool count was zero either way,
there was no signal anywhere that could tell them apart.

The count now comes from asking the servers. A server that cannot be
reached says so, with the reason it recorded.
"""

from __future__ import annotations

import pytest

from delfin.agent import mcp_client as MC


def test_the_registry_still_has_no_tools_attribute():
    """The premise of the defect, asserted so the fix cannot silently be
    undone by adding one and leaving the caller as it was."""
    import dataclasses
    fields = {f.name for f in dataclasses.fields(MC.MCPRegistry)}
    assert "tools" not in fields


def test_counting_tools_asks_the_servers(monkeypatch):
    reg = MC.MCPRegistry(servers={}, workspace=None, loaded=True)

    class _Srv:
        name = "ok"
        last_error = ""

        def list_tools(self):
            return [{"name": "a"}, {"name": "b"}]

    reg.servers = {"ok": _Srv()}
    assert MC.count_live_tools(reg) == 2


def test_a_server_that_cannot_be_reached_contributes_nothing():
    reg = MC.MCPRegistry(servers={}, workspace=None, loaded=True)

    class _Dead:
        name = "dead"
        last_error = "No such file or directory"

        def list_tools(self):
            raise OSError("No such file or directory")

    reg.servers = {"dead": _Dead()}
    assert MC.count_live_tools(reg) == 0


def test_the_reason_a_server_failed_is_reported():
    reg = MC.MCPRegistry(servers={}, workspace=None, loaded=True)

    class _Dead:
        name = "dead"
        last_error = "No such file or directory"

        def list_tools(self):
            raise OSError("boom")

    reg.servers = {"dead": _Dead()}
    problems = MC.unreachable_servers(reg)
    assert problems and "dead" in problems[0]
    assert "No such file" in problems[0]


def test_a_healthy_registry_reports_no_problems():
    reg = MC.MCPRegistry(servers={}, workspace=None, loaded=True)

    class _Srv:
        name = "ok"
        last_error = ""

        def list_tools(self):
            return [{"name": "a"}]

    reg.servers = {"ok": _Srv()}
    assert MC.unreachable_servers(reg) == []


def test_a_registry_with_no_servers_is_not_an_error():
    reg = MC.MCPRegistry(servers={}, workspace=None, loaded=True)
    assert MC.count_live_tools(reg) == 0
    assert MC.unreachable_servers(reg) == []


def test_counting_never_raises():
    """It runs on a slash command; a broken server must not take the
    dashboard down with it."""
    class _Wild:
        name = "wild"

        @property
        def last_error(self):
            raise RuntimeError("even this")

        def list_tools(self):
            raise RuntimeError("and this")

    reg = MC.MCPRegistry(servers={"wild": _Wild()}, workspace=None,
                         loaded=True)
    assert MC.count_live_tools(reg) == 0
    assert isinstance(MC.unreachable_servers(reg), list)


def test_the_dashboard_uses_the_real_count():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin"
           / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = src.index("def _hot_reload_mcp(")
    block = src[i:i + 1400]
    # Assert on the ASSIGNMENT, not on any occurrence: the comment that
    # explains the defect quotes the old expression, and a plain
    # substring check trips over it.
    assert 'n_tools = len(getattr(new_reg, "tools"' not in block, (
        "the count still reads an attribute the registry does not have")
    assert "n_tools = _mcp.count_live_tools(new_reg)" in block
