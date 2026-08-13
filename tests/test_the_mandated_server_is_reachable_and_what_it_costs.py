"""The mandated ops tools existed on a server nothing registered.

The solo agent's chemistry module says the FIRST tool call for any ORCA
question MUST be ``mcp__delfin-ops__extract_*`` or ``parse_orca_output``.
Only ``delfin-tools`` was in the native registry's built-in set. The config
the dashboard auto-injects for the ops server goes to
``~/.delfin/mcp_ops_config.json`` under the key ``mcpServers``, while this
registry reads ``~/.delfin/mcp_servers.json`` under the key ``servers``, and
the injected file reaches only the legacy CLI subprocess. So the mandated
first move answered ``unknown MCP server: 'delfin-ops'``, and the fallback
the same page forbids was the only thing that worked.

The decision taken: REGISTER the ops server. Its 67 typed tools exist
nowhere else, and the prompt, the skills and the dashboard all point at
them. ``delfin-docs`` was deliberately left out — all seven of its tools
are already advertised natively under the same names against the same
index, so registering it would buy nothing and be paid for on every
request.

The second half of these tests measures what registration costs, because
the answer is not free and should not be discovered by someone reading a
token bill. The advertising loop caps MCP schemas at
``_mcp_schema_budget_chars()`` and REPORTS what it drops, but 93 discovered
schemas do not fit in 12,000 characters. That ceiling belongs to the
advertising loop, not to the registry; these tests state the shortfall in
numbers so the trade-off is a decision someone makes rather than a surprise.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import mcp_client as MC


_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"
_OPS = Path(__file__).resolve().parents[1] / "delfin" / "ops_server" / "server.py"


def _ops_tool_names() -> set[str]:
    return set(re.findall(r"mcp\.tool\(name=\"([a-z_][a-z0-9_]*)\"",
                          _OPS.read_text(encoding="utf-8")))


# ---------------------------------------------------------------------------
# 1. The name resolves
# ---------------------------------------------------------------------------

def test_the_ops_server_is_in_the_builtin_set():
    assert "delfin-ops" in MC._BUILTIN_SERVERS


def test_dispatch_no_longer_answers_unknown_mcp_server():
    reg = MC.MCPRegistry()
    reg.load(None)
    out = reg.call("mcp__delfin-ops__parse_orca_output", {})
    # It may fail for any number of runtime reasons; it must not fail
    # because the registry has never heard of the server.
    assert "unknown MCP server" not in out, out
    reg.shutdown()


def test_a_pack_named_ops_tool_exists_on_the_server_that_is_registered():
    """The prompt names a family with a wildcard; the family has members."""
    names = _ops_tool_names()
    assert any(n.startswith("extract_") for n in names)
    assert "parse_orca_output" in names


def test_the_builtin_spawns_the_running_interpreter():
    """``"python"`` is a name a host may not have; the interpreter running
    DELFIN is the one that can import it. Before reachability was
    reported, that spawn failure looked like a healthy empty server."""
    import sys
    for name, cfg in MC._BUILTIN_SERVERS.items():
        assert cfg["command"] == (sys.executable or "python3"), name


# ---------------------------------------------------------------------------
# 2. What it costs, measured
# ---------------------------------------------------------------------------

def _schema_cost(name: str, description: str, schema: dict) -> int:
    return len(json.dumps({
        "type": "function",
        "function": {"name": name, "description": description,
                     "parameters": schema or {"type": "object"}},
    }, ensure_ascii=False))


@pytest.mark.slow
def test_the_registered_servers_answer_and_what_they_advertise_is_measured():
    """Spawns both built-ins for real: the only way to know what the model
    would actually be offered."""
    reg = MC.MCPRegistry()
    reg.load(None)
    try:
        tools = reg.discover_all()
        assert MC.unreachable_servers(reg) == [], "a built-in did not answer"
    finally:
        reg.shutdown()

    per_server: dict[str, int] = {}
    for t in tools:
        per_server[t.server] = per_server.get(t.server, 0) + 1
    assert per_server.get("delfin-ops", 0) >= 60
    assert per_server.get("delfin-tools", 0) >= 20

    total = sum(_schema_cost(t.namespaced_name, t.description or t.name,
                             t.schema) for t in tools)
    budget = A._mcp_schema_budget_chars()
    # The measurement this test exists to hold. If a future change makes
    # the whole surface fit, that is good news and this assertion should
    # be replaced by one that says so.
    assert total > budget, (
        "the MCP surface now fits in the budget; drop the shortfall note "
        f"from the pack (total={total}, budget={budget})")


def test_the_pack_tells_the_model_what_to_do_when_a_typed_tool_is_absent():
    """The budget drops schemas and says which. A page that mandates a
    tool the surface may not carry has to say what to do then, or the
    model silently takes the fallback the same page forbids."""
    text = (_PACK / "agents" / "solo_agent.md").read_text(encoding="utf-8")
    assert "list_tools" in text
    # The instruction must name the recovery, not just the mandate.
    assert "not in your tool list" in text


def test_the_dropped_surface_is_reported_not_silently_cut():
    """The advertising loop records what it could not afford. Asserted on
    the source because the branch needs a live model turn to run."""
    src = (Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert "_mcp_dropped" in src
    assert '"mcp_budget"' in src
    from delfin.agent import security_events
    assert "mcp_budget" in security_events.known_kinds()
