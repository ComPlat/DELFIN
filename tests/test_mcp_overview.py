"""Tests for delfin.dashboard.mcp_overview helpers."""
from __future__ import annotations

import json
from pathlib import Path

from delfin.dashboard import mcp_overview


# ---------------------------------------------------------------------------
# discover_mcp_servers
# ---------------------------------------------------------------------------

def test_discover_no_paths(tmp_path):
    """Empty path list → empty result, no errors."""
    result = mcp_overview.discover_mcp_servers(config_paths=[])
    assert result == []


def test_discover_missing_files(tmp_path):
    """Missing files are silently skipped."""
    paths = [tmp_path / "missing1.json", tmp_path / "missing2.json"]
    result = mcp_overview.discover_mcp_servers(config_paths=paths)
    assert result == []


def test_discover_reads_single_config(tmp_path):
    cfg = tmp_path / "mcp.json"
    cfg.write_text(json.dumps({
        "mcpServers": {
            "delfin-ops": {"command": "/usr/bin/python", "args": ["-m", "delfin.ops_server"]}
        }
    }))
    result = mcp_overview.discover_mcp_servers(config_paths=[cfg])
    assert len(result) == 1
    entry = result[0]
    assert entry["server_id"] == "delfin-ops"
    assert entry["command"] == "/usr/bin/python"
    assert entry["args"] == ["-m", "delfin.ops_server"]
    assert entry["config_path"] == str(cfg)


def test_discover_merges_two_configs(tmp_path):
    cfg_docs = tmp_path / "docs.json"
    cfg_ops = tmp_path / "ops.json"
    cfg_docs.write_text(json.dumps({
        "mcpServers": {
            "delfin-docs": {"command": "py", "args": ["-m", "delfin.doc_server"]}
        }
    }))
    cfg_ops.write_text(json.dumps({
        "mcpServers": {
            "delfin-ops": {"command": "py", "args": ["-m", "delfin.ops_server"]}
        }
    }))
    result = mcp_overview.discover_mcp_servers(config_paths=[cfg_docs, cfg_ops])
    ids = {e["server_id"] for e in result}
    assert ids == {"delfin-docs", "delfin-ops"}


def test_discover_dedups_first_wins(tmp_path):
    """Same server_id in two files → first occurrence is kept."""
    cfg1 = tmp_path / "first.json"
    cfg2 = tmp_path / "second.json"
    cfg1.write_text(json.dumps({
        "mcpServers": {"delfin-ops": {"command": "first", "args": []}}
    }))
    cfg2.write_text(json.dumps({
        "mcpServers": {"delfin-ops": {"command": "second", "args": []}}
    }))
    result = mcp_overview.discover_mcp_servers(config_paths=[cfg1, cfg2])
    assert len(result) == 1
    assert result[0]["command"] == "first"


def test_discover_skips_corrupt_json(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{ not json")
    good = tmp_path / "good.json"
    good.write_text(json.dumps({
        "mcpServers": {"delfin-ops": {"command": "py", "args": []}}
    }))
    result = mcp_overview.discover_mcp_servers(config_paths=[bad, good])
    assert len(result) == 1
    assert result[0]["server_id"] == "delfin-ops"


def test_discover_skips_non_dict_entries(tmp_path):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({
        "mcpServers": {
            "broken": "not a dict",
            "good": {"command": "py", "args": []},
        }
    }))
    result = mcp_overview.discover_mcp_servers(config_paths=[cfg])
    assert {e["server_id"] for e in result} == {"good"}


# ---------------------------------------------------------------------------
# known_tools
# ---------------------------------------------------------------------------

def test_known_tools_delfin_docs():
    tools = mcp_overview.known_tools("delfin-docs")
    assert "search_docs" in tools
    assert "read_section" in tools
    assert "search_calcs" in tools


def test_known_tools_delfin_ops_includes_mutating():
    tools = mcp_overview.known_tools("delfin-ops")
    # Read-only
    assert "qm_check" in tools
    assert "stop_dry_run" in tools
    # Mutating
    assert "pipeline_run" in tools
    assert "co2" in tools


def test_known_tools_unknown_server_returns_empty():
    assert mcp_overview.known_tools("custom-server") == []


# ---------------------------------------------------------------------------
# render_overview_html
# ---------------------------------------------------------------------------

def test_render_empty_returns_empty_string():
    assert mcp_overview.render_overview_html([]) == ""


def test_render_includes_server_id_and_active_badge():
    html = mcp_overview.render_overview_html([{
        "server_id": "delfin-ops",
        "command": "/usr/bin/python",
        "args": ["-m", "delfin.ops_server"],
        "config_path": "/tmp/x.json",
    }])
    assert "delfin-ops" in html
    assert "ACTIVE" in html
    assert "delfin.ops_server" in html


def test_render_includes_known_tool_chips():
    html = mcp_overview.render_overview_html([{
        "server_id": "delfin-ops", "command": "py", "args": [], "config_path": "x",
    }])
    # Known tools are rendered as chips
    assert "qm_check" in html
    assert "pipeline_run" in html


def test_render_unknown_server_omits_tool_chips():
    html = mcp_overview.render_overview_html([{
        "server_id": "custom-thing", "command": "py", "args": [], "config_path": "x",
    }])
    assert "custom-thing" in html
    assert "tools not introspected" in html


def test_render_escapes_user_text():
    html = mcp_overview.render_overview_html([{
        "server_id": "<script>",
        "command": "<bad>",
        "args": ["<arg>"],
        "config_path": "<path>",
    }])
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "&lt;bad&gt;" in html


def test_render_count_in_title():
    html = mcp_overview.render_overview_html([
        {"server_id": "a", "command": "py", "args": [], "config_path": ""},
        {"server_id": "b", "command": "py", "args": [], "config_path": ""},
    ])
    assert "(2)" in html
