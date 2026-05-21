"""MCP configuration helpers for the operations server."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def generate_mcp_config(workspace: str | Path | None = None) -> dict:
    """Generate an MCP server configuration dict for the ops server.

    Parameters
    ----------
    workspace : str or Path, optional
        Default workspace passed to the server (used for tools that take
        a workspace argument, e.g. ``cleanup``, ``stop``).  If omitted,
        the server uses ``os.getcwd()`` at startup.
    """
    args = ["-m", "delfin.ops_server"]
    if workspace is not None:
        args.extend(["--workspace", str(workspace)])
    return {
        "mcpServers": {
            "delfin-ops": {
                "command": sys.executable,
                "args": args,
            }
        }
    }


def _merge_mcp_configs(base: dict, overlay: dict) -> dict:
    """Merge two MCP config dicts, combining their server entries."""
    merged = dict(base)
    base_servers = merged.get("mcpServers", {})
    overlay_servers = overlay.get("mcpServers", {})
    merged["mcpServers"] = {**base_servers, **overlay_servers}
    return merged


def ensure_mcp_config(
    existing_mcp_config: str = "",
    workspace: str | Path | None = None,
) -> str:
    """Ensure the ops server is included in the MCP configuration.

    If ``existing_mcp_config`` points to an existing config file, the ops
    server entry is merged into it.  Otherwise a new file is created.

    Returns the path to the (possibly updated) MCP config file.
    """
    ops_config = generate_mcp_config(workspace=workspace)

    config_dir = Path.home() / ".delfin"
    config_dir.mkdir(parents=True, exist_ok=True)
    output_path = config_dir / "mcp_ops_config.json"

    if existing_mcp_config and Path(existing_mcp_config).exists():
        try:
            base = json.loads(Path(existing_mcp_config).read_text(encoding="utf-8"))
            merged = _merge_mcp_configs(base, ops_config)
        except (json.JSONDecodeError, OSError):
            merged = ops_config
    else:
        merged = ops_config

    output_path.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    return str(output_path)
