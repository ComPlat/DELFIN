"""MCP server inventory helpers for the Agent Activity tab.

Shows which MCP servers (delfin-docs, delfin-ops, plus any user-injected
ones) are configured for the agent and what tools they expose.  Pure
functions only — no widgets — so they can be tested without ipywidgets.
"""
from __future__ import annotations

import html as _html
import json
from pathlib import Path


# Known tool catalogs for the DELFIN-shipped servers.  These are defined
# in code rather than introspected at runtime so the inventory works even
# when the MCP package is not installed.
_KNOWN_SERVER_TOOLS: dict[str, list[str]] = {
    "delfin-docs": [
        "search_docs",
        "read_section",
        "list_docs",
        "list_sections",
        "search_calcs",
        "get_calc_info",
        "calc_summary",
    ],
    "delfin-ops": [
        # Read-only
        "qm_check",
        "csp_check",
        "mlp_check",
        "analysis_check",
        "stop_dry_run",
        # Mutating (require allow_mutate=True)
        "cleanup",
        "stop",
        "pipeline_run",
        "pipeline_prepare",
        "run_orca_input",
        "co2",
        "tadf_xtb",
        "hyperpol",
    ],
}


def _default_config_paths() -> list[Path]:
    """Return likely MCP config locations."""
    base = Path.home() / ".delfin"
    return [
        base / "mcp_ops_config.json",
        base / "mcp_docs_config.json",
    ]


def discover_mcp_servers(
    config_paths: list[Path] | None = None,
) -> list[dict]:
    """Inventory MCP servers from one or more config files.

    Each returned entry: ``{"server_id", "command", "args", "config_path"}``.
    Duplicate ``server_id`` values across files are de-duplicated (first wins).

    Returns an empty list if no config exists or all are unreadable.
    """
    paths = config_paths if config_paths is not None else _default_config_paths()

    seen: dict[str, dict] = {}
    for path in paths:
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        servers = data.get("mcpServers", {}) or {}
        for server_id, info in servers.items():
            if server_id in seen:
                continue
            if not isinstance(info, dict):
                continue
            seen[server_id] = {
                "server_id": server_id,
                "command": str(info.get("command", "")),
                "args": list(info.get("args", []) or []),
                "config_path": str(path),
            }
    return list(seen.values())


def known_tools(server_id: str) -> list[str]:
    """Return the curated tool list for a known DELFIN MCP server.

    Unknown servers return an empty list (UI shows "tools not introspected").
    """
    return list(_KNOWN_SERVER_TOOLS.get(server_id, []))


def render_overview_html(servers: list[dict]) -> str:
    """Render the MCP server inventory as a compact HTML block.

    Empty input → empty string so the caller can hide the widget.
    """
    if not servers:
        return ""

    cards: list[str] = []
    for entry in servers:
        sid = entry.get("server_id", "")
        config_path = entry.get("config_path", "")
        tools = known_tools(sid)
        if tools:
            tool_chips = "".join(
                f'<span style="background:#e0e7ff;color:#3730a3;'
                f'font-size:10px;padding:2px 7px;border-radius:8px;'
                f'margin:2px 4px 2px 0;display:inline-block;">'
                f'{_html.escape(t)}</span>'
                for t in tools
            )
            tool_section = (
                f'<div style="margin-top:6px;">{tool_chips}</div>'
                f'<div style="font-size:10px;color:#9ca3af;margin-top:2px;">'
                f'{len(tools)} tools</div>'
            )
        else:
            tool_section = (
                '<div style="font-size:11px;color:#9ca3af;margin-top:4px;">'
                'tools not introspected</div>'
            )

        args_str = " ".join(_html.escape(str(a)) for a in entry.get("args", []))
        cards.append(
            '<div style="border:1px solid #e5e7eb;border-radius:8px;'
            'background:#fff;padding:10px 12px;margin-bottom:8px;">'
            f'<div style="display:flex;align-items:center;gap:8px;'
            f'margin-bottom:4px;">'
            f'<span style="font-size:13px;font-weight:700;color:#111827;">'
            f'{_html.escape(sid)}</span>'
            f'<span style="background:#10b981;color:#fff;font-size:10px;'
            f'padding:1px 6px;border-radius:8px;font-weight:700;">ACTIVE</span>'
            '</div>'
            f'<div style="font-size:11px;color:#6b7280;font-family:monospace;">'
            f'{_html.escape(entry.get("command", ""))} {args_str}</div>'
            f'<div style="font-size:10px;color:#9ca3af;margin-top:2px;">'
            f'config: {_html.escape(config_path)}</div>'
            + tool_section
            + '</div>'
        )

    title = (
        '<h4 style="margin:12px 0 6px 0;color:#111827;'
        'font-size:13px;font-weight:700;">'
        f'MCP Servers ({len(servers)})</h4>'
        '<p style="margin:0 0 8px 0;color:#6b7280;font-size:11px;">'
        'Configured Model Context Protocol servers that the agent can call.'
        '</p>'
    )
    return title + "".join(cards)


__all__ = [
    "discover_mcp_servers",
    "known_tools",
    "render_overview_html",
]
