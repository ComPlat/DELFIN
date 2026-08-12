"""The prompt pack named tools the model is never given.

Measured on the shipped pack: `mcp__delfin-docs__search` (the doc server's
tool is `search_docs`), `mcp__kit-coding__remember_permission` and
`…_bundle` (native tools, advertised BARE — the namespace is a display
detail the executor strips), `glob_files` (in no schema at all), the
CamelCase `Read` / `Grep` / `Bash` / `WebSearch` / `Glob` (generic harness
names; every tool here is snake_case), and `ExitPlanMode` for
`exit_plan_mode`.

A name the model cannot call is worse than silence: it either spends a
round on a refusal or quietly does the thing the instruction was meant to
prevent. Nothing checked these against the schema, so they drifted.

The three rules below are mechanical on purpose — they catch the shapes a
tool name has in this codebase without needing a hand-maintained list of
the hundreds of backticked identifiers in the pack that are NOT tools
(ORCA keywords, config keys, Python symbols).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"


def _advertised() -> set[str]:
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    return {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}


def _mcp_registry() -> dict[str, set[str]]:
    """Tool names per MCP server this repository itself defines.

    Read out of the servers' registration calls rather than restated here,
    so a renamed tool fails this test instead of drifting past it.
    """
    root = Path(__file__).resolve().parents[1] / "delfin"
    registry: dict[str, set[str]] = {}
    docs = (root / "doc_server" / "server.py").read_text(encoding="utf-8")
    registry["delfin-docs"] = set(
        re.findall(r"@mcp\.tool\(\)\s*\n\s*def\s+([a-z_][a-z0-9_]*)", docs))
    ops = (root / "ops_server" / "server.py").read_text(encoding="utf-8")
    registry["delfin-ops"] = set(
        re.findall(r"mcp\.tool\(name=\"([a-z_][a-z0-9_]*)\"\)", ops))
    return registry


def _pack_files() -> list[Path]:
    return sorted(_PACK.rglob("*.md"))


# ---------------------------------------------------------------------------
# 1. A namespaced name must resolve
# ---------------------------------------------------------------------------

_MCP_RE = re.compile(r"mcp__([A-Za-z0-9_-]+)__([A-Za-z0-9_]+\*?|\*)")


def test_every_namespaced_tool_name_resolves_to_a_real_server_tool():
    registry = _mcp_registry()
    bad: list[str] = []
    for path in _pack_files():
        for match in _MCP_RE.finditer(path.read_text(encoding="utf-8")):
            server, tool = match.group(1), match.group(2)
            if tool == "*":
                # A namespace, not a tool — "when your tool list includes
                # mcp__kit-coding__*" says something true about a backend.
                continue
            known = registry.get(server)
            if known is None:
                bad.append(
                    f"{path.name}: {match.group(0)} — no MCP server named "
                    f"{server!r} is defined here; a native tool is "
                    f"advertised BARE, so write it bare")
                continue
            if tool.endswith("*"):
                if not any(n.startswith(tool[:-1]) for n in known):
                    bad.append(f"{path.name}: {match.group(0)} matches no "
                               f"tool on {server}")
                continue
            if tool not in known:
                bad.append(
                    f"{path.name}: {match.group(0)} — {server} has no tool "
                    f"{tool!r} (it has {sorted(known)})")
    assert not bad, "namespaced tool names that do not exist:\n" + "\n".join(bad)


# ---------------------------------------------------------------------------
# 2. A tool LIST must contain only tools
# ---------------------------------------------------------------------------
#
# A run of adjacent backticked identifiers glued by ", " / "/" / "or" is
# how the pack writes a tool list. When two or more members of such a run
# are known tools, the run IS a tool list and every member has to be one.
# Two is the bar because a single known name next to a parameter value
# ("`remember_permission`, `extra_dir`, or decline") is not a list of tools.

_TOKEN_RE = re.compile(r"`([A-Za-z_][A-Za-z0-9_]*)`")
_GLUE_RE = re.compile(r"[\s,/+]*(?:or|and|,)?[\s,/+]*")


def _runs(text: str) -> list[list[str]]:
    tokens = list(_TOKEN_RE.finditer(text))
    out: list[list[str]] = []
    i = 0
    while i < len(tokens):
        run = [tokens[i]]
        j = i + 1
        while j < len(tokens):
            glue = text[tokens[j - 1].end():tokens[j].start()]
            if len(glue) <= 6 and _GLUE_RE.fullmatch(glue):
                run.append(tokens[j])
                j += 1
            else:
                break
        if len(run) >= 2:
            out.append([m.group(1) for m in run])
        i = max(j, i + 1)
    return out


def test_a_list_of_tools_contains_only_tools_that_exist():
    known = _advertised()
    bad: list[str] = []
    for path in _pack_files():
        for names in _runs(path.read_text(encoding="utf-8")):
            if sum(1 for n in names if n in known) < 2:
                continue
            for name in names:
                if name not in known:
                    bad.append(f"{path.name}: {name!r} listed among tools "
                               f"{[n for n in names if n in known]}")
    assert not bad, "names listed as tools but in no schema:\n" + "\n".join(bad)


# ---------------------------------------------------------------------------
# 3. No generic CamelCase harness tool names
# ---------------------------------------------------------------------------
#
# Every tool in this codebase is snake_case. A capitalised generic tool
# name is therefore always some other harness's spelling, and the model has
# nothing by that name to call.

_FOREIGN_TOOL_NAMES = (
    "Read", "Write", "Edit", "MultiEdit", "Bash", "Grep", "Glob", "LS",
    "Task", "WebSearch", "WebFetch", "NotebookRead", "NotebookEdit",
    "TodoWrite", "ExitPlanMode",
)


@pytest.mark.parametrize("name", _FOREIGN_TOOL_NAMES)
def test_no_foreign_camelcase_tool_name_is_presented_as_callable(name):
    hits = []
    for path in _pack_files():
        text = path.read_text(encoding="utf-8")
        for pattern in (f"`{name}`", f"``{name}``", f"`{name}(", f"{name}("):
            if pattern in text:
                hits.append(f"{path.name}: {pattern}")
    assert not hits, (
        f"{name} is not a tool here (tools are snake_case): " + ", ".join(hits))
