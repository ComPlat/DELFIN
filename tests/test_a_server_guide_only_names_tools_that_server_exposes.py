"""The tools server's own guide named a tool it does not register.

``get_guide`` is the mandated first call in the ``pipeline-build`` skill,
and its ``discover`` step listed ``compatible_successors`` beside
``list_capabilities`` and ``describe_capability``. The function exists and
is exported from ``delfin.tools``, but ``delfin/tools/mcp_server.py`` never
wrapped it in ``@mcp.tool()`` — so the first thing a model was told to read
sent it to a tool the server does not answer to.

The drift test that covers the prompt pack could not see this: its registry
knew ``delfin-docs`` and ``delfin-ops`` and omitted the tools server
entirely, and the guide is not a pack file anyway — it is data the server
returns about itself.

The rule here is per server: every tool name the server's own guide names
has to be a tool that server registers.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


def _registered(module_relpath: str, ) -> set[str]:
    root = Path(__file__).resolve().parents[1] / "delfin"
    text = (root / module_relpath).read_text(encoding="utf-8")
    return (
        set(re.findall(r"@mcp\.tool\(\)\s*\n\s*def\s+([a-z_][a-z0-9_]*)", text))
        | set(re.findall(r"mcp\.tool\(name=\"([a-z_][a-z0-9_]*)\"", text))
    )


def _guide_tool_names() -> set[str]:
    from delfin.tools.manifest import build_guide

    names: set[str] = set()
    for step in build_guide().get("workflow", []):
        for tool in step.get("tools", []):
            # "pipeline_spec schema" is prose naming a schema, not a call.
            if re.fullmatch(r"[a-z_][a-z0-9_]*", tool):
                names.add(tool)
    return names


def test_the_guide_names_only_tools_the_tools_server_registers():
    registered = _registered("tools/mcp_server.py")
    assert registered, "the tools server registers nothing — the scan broke"
    missing = sorted(_guide_tool_names() - registered)
    assert not missing, (
        "the tools server's guide sends the model to tools it does not "
        f"expose: {missing}. Register them, or stop naming them.")


def test_compatible_successors_is_the_one_that_was_missing():
    """Named explicitly so a silent removal is a failure, not a pass."""
    assert "compatible_successors" in _guide_tool_names()
    assert "compatible_successors" in _registered("tools/mcp_server.py")


def test_the_registered_tool_answers():
    from delfin.tools.mcp_server import h_compatible_successors
    import json

    out = json.loads(h_compatible_successors("orca_freq"))
    assert isinstance(out, list) and out, out
    bad = json.loads(h_compatible_successors("no_such_capability"))
    assert "error" in bad


@pytest.mark.parametrize("skill_name", ["pipeline-build"])
def test_a_skill_only_names_tools_of_the_server_it_names(skill_name):
    """The skill that tells the model to work "only through delfin-tools"
    must not, in the same page, name tools that server does not have."""
    from delfin.agent import skills as S

    sk = S.get_skill(skill_name)
    assert sk is not None
    registered = _registered("tools/mcp_server.py")
    # Backticked identifiers that are also known tools-server tool names
    # confirm the page is talking about that server; anything backticked
    # that LOOKS like a tools-server call must then exist there.
    named = set(re.findall(r"`([a-z_][a-z0-9_]*)\(", sk.body))
    unknown = sorted(n for n in named
                     if n not in registered and n in _guide_tool_names())
    assert not unknown, (
        f"{skill_name} calls tools-server functions that are not "
        f"registered: {unknown}")
