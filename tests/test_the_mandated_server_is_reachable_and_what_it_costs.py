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
        # The REASON, not a fixed sentence. When mcp 2.0.0 removed
        # ``mcp.server.fastmcp`` this assertion fired reading "a built-in did
        # not answer", which is the one thing already known; the traceback
        # that said which module had moved went to DEVNULL.
        problems = MC.unreachable_servers(reg)
        assert problems == [], f"a built-in did not answer: {problems}"
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


# ---------------------------------------------------------------------------
# 3. It survives the SDK's major version, which is where it died
# ---------------------------------------------------------------------------
#
# ``mcp`` 2.0.0 (2026-07-28) removed ``mcp.server.fastmcp``; the replacement
# is ``mcp.server.mcpserver.MCPServer``. All three DELFIN servers imported
# the old path unconditionally and ``pyproject.toml`` declared ``mcp>=1.0``,
# so a fresh install inside DELFIN's own declared range resolved 2.0.0 and
# every built-in exited 1 before answering a single request — the whole
# tools platform and all 67 ops tools, gone.
#
# The box these tests usually run on has an mcp 1.x, on which the defect is
# invisible. So the 2.x line is SIMULATED here: the 1.x module is hidden and
# a stand-in is put at the 2.x path, which is the only way this box can
# assert on the behaviour that broke.

def test_every_delfin_server_reaches_the_sdk_through_the_compat_loader():
    """The import that broke was correct on the day it was written.

    What makes it break again is a fourth server copying the old line, so
    the guard is on the source: outside ``mcp_compat`` itself, nothing may
    name a version-specific SDK server path.
    """
    root = Path(__file__).resolve().parents[1] / "delfin"
    # An IMPORT of a version-specific path, not a mention of one: these
    # module names are named in prose where the defect is explained, and
    # explaining it is the opposite of repeating it.
    moved = re.compile(r"^\s*(?:from|import)\s+mcp\.server\."
                       r"(?:fastmcp|mcpserver)\b", re.MULTILINE)
    offenders = []
    for path in root.rglob("*.py"):
        if path.name == "mcp_compat.py":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for hit in moved.findall(text):
            offenders.append(f"{path.relative_to(root)}: {hit.strip()}")
    assert offenders == [], (
        "these reach a server class that moved between SDK major versions; "
        f"go through delfin.mcp_compat.load_server_class instead: {offenders}")


def test_the_loader_finds_a_server_class_on_the_installed_sdk():
    """Whatever line is installed here, the loader must return the class the
    three servers are built on — constructible, and carrying the three API
    points they use."""
    from delfin.mcp_compat import load_server_class

    cls = load_server_class()
    server = cls("probe", instructions="probe")
    assert callable(server.tool)
    assert callable(server.run)


def test_the_two_1x_names_and_the_2x_name_are_both_tried():
    from delfin import mcp_compat

    assert ("mcp.server.fastmcp", "FastMCP") in mcp_compat._SERVER_CLASSES
    assert ("mcp.server.mcpserver", "MCPServer") in mcp_compat._SERVER_CLASSES


def test_the_loader_falls_through_to_the_2x_class_when_1x_is_gone(monkeypatch):
    """The defect itself, reproduced on a box that has an mcp 1.x.

    ``mcp.server.fastmcp`` is made unimportable and a stand-in is placed at
    the 2.0.0 path. Before the fallback existed this raised
    ``ModuleNotFoundError`` and the server process exited 1.
    """
    import sys
    from delfin import mcp_compat

    class _Stand_in:
        pass

    module = type(sys)("mcp.server.mcpserver")
    module.MCPServer = _Stand_in
    monkeypatch.setitem(sys.modules, "mcp.server.mcpserver", module)

    # Blocked at the meta path, and evicted from the module cache it would
    # otherwise be served from. Patching ``builtins.__import__`` does not
    # do it: ``importlib.import_module`` goes to ``_bootstrap._gcd_import``
    # and never touches the builtin.
    class _NoFastMCP:
        def find_spec(self, name, path=None, target=None):
            if name == "mcp.server.fastmcp":
                raise ModuleNotFoundError(
                    "No module named 'mcp.server.fastmcp'", name=name)
            return None

    monkeypatch.delitem(sys.modules, "mcp.server.fastmcp", raising=False)
    monkeypatch.setattr(sys, "meta_path", [_NoFastMCP(), *sys.meta_path])
    assert mcp_compat.load_server_class() is _Stand_in


def test_an_sdk_with_neither_class_says_so_in_full(monkeypatch):
    """The message is what a user gets when a built-in will not start, so it
    has to name what was tried and what is installed rather than repeat the
    last ``ModuleNotFoundError`` on its own."""
    from delfin import mcp_compat

    monkeypatch.setattr(mcp_compat, "_SERVER_CLASSES",
                        (("delfin_no_such_sdk_a", "A"),
                         ("delfin_no_such_sdk_b", "B")))
    with pytest.raises(ImportError) as excinfo:
        mcp_compat.load_server_class()
    text = str(excinfo.value)
    assert "delfin_no_such_sdk_a" in text
    assert "delfin_no_such_sdk_b" in text
    assert "installed mcp:" in text
