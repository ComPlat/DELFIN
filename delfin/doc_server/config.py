"""MCP configuration helpers for the documentation server."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from .indexer import get_default_index_path


def generate_mcp_config(index_path: str | Path) -> dict:
    """Generate an MCP server configuration dict for the doc server.

    Parameters
    ----------
    index_path : str or Path
        Path to the JSON index file built by ``delfin-docs-index``.

    Returns
    -------
    dict
        MCP configuration ready to be written as JSON.
    """
    return {
        "mcpServers": {
            "delfin-docs": {
                "command": sys.executable,
                "args": ["-m", "delfin.doc_server", "--index", str(index_path)],
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


def ensure_mcp_config(existing_mcp_config: str = "", index_path: str | Path | None = None) -> str:
    """Ensure the doc server is included in the MCP configuration.

    If ``existing_mcp_config`` points to an existing config file, the doc
    server entry is merged into it.  Otherwise a new config file is created.

    Parameters
    ----------
    existing_mcp_config : str
        Path to an existing MCP config file, or empty string.
    index_path : str or Path, optional
        Path to the doc index.  Defaults to ``~/.delfin/doc_index.json``.

    Returns
    -------
    str
        Path to the (possibly updated) MCP config file.
    """
    if index_path is None:
        index_path = get_default_index_path()
    index_path = Path(index_path)

    # If the index doesn't exist, don't inject the doc server
    if not index_path.exists():
        return existing_mcp_config

    doc_config = generate_mcp_config(index_path)

    # Determine output path
    config_dir = Path.home() / ".delfin"
    config_dir.mkdir(parents=True, exist_ok=True)
    output_path = config_dir / "mcp_docs_config.json"

    if existing_mcp_config and Path(existing_mcp_config).exists():
        try:
            base = json.loads(Path(existing_mcp_config).read_text(encoding="utf-8"))
            merged = _merge_mcp_configs(base, doc_config)
        except (json.JSONDecodeError, OSError):
            merged = doc_config
    else:
        merged = doc_config

    output_path.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    return str(output_path)
