"""The MCP SDK's server class, under whichever name the installed line uses.

``mcp`` 2.0.0 (published 2026-07-28) removed ``mcp.server.fastmcp``. The
class that replaces ``FastMCP`` is ``mcp.server.mcpserver.MCPServer``.
DELFIN's three stdio servers each imported the old path unconditionally
and ``pyproject.toml`` declares ``mcp>=1.0`` with no ceiling, so anyone
installing DELFIN inside its own declared range resolved 2.0.0 and got
``ModuleNotFoundError: No module named 'mcp.server.fastmcp'``, exit 1, on
every one of them: the 67 ops tools and the whole tools platform gone,
for a module that had merely moved.

A major version is rarely only a rename, so the surface DELFIN actually
uses was measured on both lines rather than assumed. It is three points,
and they are signature-identical in 1.27.1 and 2.0.0: the constructor's
``(name, instructions=...)``, the ``tool(name=..., description=...)``
decorator, and ``run(transport="stdio")``. Both classes were then started
for real over stdio and answered ``initialize`` and ``tools/list`` with
the same tool sets under the same ``2025-06-18`` handshake the agent's
client sends. Nothing else of the SDK is imported anywhere in DELFIN —
the agent's own client speaks JSON-RPC over a pipe by hand and never
loads the package — so this one name is the entire crossing.

Hence a shim rather than an ``mcp<2`` ceiling: pinning would have cost
the same three lines and left DELFIN frozen out of the current SDK line
with someone obliged to come back to it.

The two lines are disjoint — 1.29.0, the last of the 1.x line, carries
only ``fastmcp``, and 2.0.0 only ``mcpserver``, so no release exists in
which the order of these attempts could change the answer. The old name
is tried first regardless, so that the version DELFIN runs on today takes
precisely the path it took before this module existed.
"""

from __future__ import annotations

import importlib
from typing import Any

# Oldest first, so an installed 1.x resolves on its own name. Each entry is
# ``(module, attribute)``.
_SERVER_CLASSES: tuple[tuple[str, str], ...] = (
    ("mcp.server.fastmcp", "FastMCP"),
    ("mcp.server.mcpserver", "MCPServer"),
)


def _installed_version() -> str:
    """The MCP SDK version actually present, for the failure message."""
    try:
        from importlib.metadata import PackageNotFoundError, version
        return version("mcp")
    except PackageNotFoundError:
        return "not installed"
    except Exception:               # pragma: no cover - metadata unreadable
        return "unknown"


def load_server_class() -> Any:
    """Return the SDK server class to build a DELFIN stdio server on.

    Raises ``ImportError`` naming every path that was tried, the reason
    each failed and the SDK version actually installed. That text is the
    whole of what a user gets when a built-in server will not start — it
    leaves as the child's stderr and is surfaced by the agent's client as
    the server's ``last_error`` — so it has to read as an answer rather
    than as a bare complaint about a missing module.
    """
    problems: list[str] = []
    for module, attr in _SERVER_CLASSES:
        try:
            return getattr(importlib.import_module(module), attr)
        except (ImportError, AttributeError) as exc:
            problems.append(f"{module}.{attr}: {exc}")
    raise ImportError(
        "the installed MCP SDK offers no server class DELFIN can use — "
        + "; ".join(problems)
        + f" (installed mcp: {_installed_version()}). DELFIN supports the "
        "1.x line (mcp.server.fastmcp.FastMCP) and the 2.x line "
        "(mcp.server.mcpserver.MCPServer)."
    )
