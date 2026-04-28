"""Tests for delfin.ops_server tool functions and config helpers.

Tool functions are tested at module level (no FastMCP needed).  The MCP
stdio bootstrap (``run_server``) is not exercised here — it requires the
optional ``mcp`` package and a live JSON-RPC client.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from delfin.ops_server import server as ops_server
from delfin.ops_server import config as ops_config


# ---------------------------------------------------------------------------
# Fake CLI recorder fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def cli_recorder(monkeypatch):
    captured: list[list[str]] = []

    def fake_main(argv):
        captured.append(list(argv))
        print("ops-stdout")
        return 0

    import delfin.cli
    monkeypatch.setattr(delfin.cli, "main", fake_main)
    return captured


# ---------------------------------------------------------------------------
# _format_result / _refuse_mutation
# ---------------------------------------------------------------------------

def test_format_result_basic():
    from delfin.api import CommandResult
    rc = CommandResult(returncode=0, stdout="hi", stderr="", argv=["cleanup"])
    out = ops_server._format_result(rc, action="cleanup")
    payload = json.loads(out)
    assert payload["action"] == "cleanup"
    assert payload["ok"] is True
    assert payload["stdout"] == "hi"
    assert payload["argv"] == ["cleanup"]


def test_format_result_truncates_long_stdout():
    from delfin.api import CommandResult
    long_text = "x" * 25_000
    rc = CommandResult(stdout=long_text)
    out = ops_server._format_result(rc, action="qm_check")
    payload = json.loads(out)
    assert "truncated" in payload["stdout"]
    assert len(payload["stdout"]) < 25_000


def test_refuse_mutation_payload():
    out = ops_server._refuse_mutation("pipeline_run")
    payload = json.loads(out)
    assert payload["ok"] is False
    assert payload["error"] == "mutation_blocked"
    assert "pipeline_run" in payload["message"]


# ---------------------------------------------------------------------------
# Read-only tools
# ---------------------------------------------------------------------------

def test_tool_qm_check_no_filter(cli_recorder):
    out = ops_server.tool_qm_check()
    assert cli_recorder[0] == ["qm_check"]
    payload = json.loads(out)
    assert payload["action"] == "qm_check"
    assert payload["ok"] is True


def test_tool_qm_check_with_filter(cli_recorder):
    ops_server.tool_qm_check(tools="xtb, crest")
    assert cli_recorder[0] == ["qm_check", "xtb", "crest"]


def test_tool_csp_check(cli_recorder):
    ops_server.tool_csp_check()
    assert cli_recorder[0] == ["csp_check"]


def test_tool_mlp_check(cli_recorder):
    ops_server.tool_mlp_check()
    assert cli_recorder[0] == ["mlp_check"]


def test_tool_analysis_check(cli_recorder):
    ops_server.tool_analysis_check()
    assert cli_recorder[0] == ["analysis_check"]


def test_tool_stop_dry_run(cli_recorder):
    ops_server.tool_stop_dry_run("/tmp/ws")
    argv = cli_recorder[0]
    assert "stop" == argv[0]
    assert "/tmp/ws" in argv
    assert "--dry-run" in argv


# ---------------------------------------------------------------------------
# Mutation guard
# ---------------------------------------------------------------------------

def test_pipeline_run_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_pipeline_run()
    payload = json.loads(out)
    assert payload["ok"] is False
    assert payload["error"] == "mutation_blocked"
    assert cli_recorder == []  # CLI never called


def test_pipeline_prepare_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_pipeline_prepare()
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_run_orca_input_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_run_orca_input(input_file="job.inp")
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_co2_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_co2(define=True)
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_tadf_xtb_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_tadf_xtb()
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_hyperpol_refuses_without_allow_mutate(cli_recorder):
    out = ops_server.tool_hyperpol()
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


# ---------------------------------------------------------------------------
# Cleanup / stop with mutation guard semantics
# ---------------------------------------------------------------------------

def test_cleanup_dry_run_default_runs_cli(cli_recorder):
    """dry_run=True is allowed without allow_mutate."""
    out = ops_server.tool_cleanup(workspace="/tmp/x")
    assert cli_recorder[0][0] == "cleanup"
    assert "--dry-run" in cli_recorder[0]
    assert json.loads(out)["dry_run"] is True


def test_cleanup_real_refused_without_allow_mutate(cli_recorder):
    out = ops_server.tool_cleanup(dry_run=False, workspace="/tmp/x")
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_cleanup_real_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_cleanup(dry_run=False, allow_mutate=True, workspace="/tmp/x")
    argv = cli_recorder[0]
    assert "--dry-run" not in argv
    assert "--workspace" in argv and "/tmp/x" in argv


def test_stop_dry_run_default_runs_cli(cli_recorder):
    ops_server.tool_stop(workspace="/tmp/y")
    assert cli_recorder[0][0] == "stop"
    assert "--dry-run" in cli_recorder[0]


def test_stop_real_refused_without_allow_mutate(cli_recorder):
    out = ops_server.tool_stop(dry_run=False, workspace="/tmp/y")
    assert json.loads(out)["error"] == "mutation_blocked"
    assert cli_recorder == []


def test_stop_real_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_stop(dry_run=False, allow_mutate=True, workspace="/tmp/y",
                         signal_name="TERM")
    argv = cli_recorder[0]
    assert "--dry-run" not in argv
    assert "TERM" in argv


# ---------------------------------------------------------------------------
# Mutating tools when allow_mutate=True
# ---------------------------------------------------------------------------

def test_pipeline_run_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_pipeline_run(allow_mutate=True, recalc=True)
    argv = cli_recorder[0]
    assert "--recalc" in argv


def test_pipeline_prepare_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_pipeline_prepare(control_file="X.txt", overwrite=True,
                                     allow_mutate=True)
    assert cli_recorder[0] == ["--define", "X.txt", "--overwrite"]


def test_run_orca_input_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_run_orca_input(input_file="job.inp", allow_mutate=True)
    assert cli_recorder[0] == ["run_orca", "job.inp"]


def test_co2_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_co2(define=True, charge=2, allow_mutate=True)
    argv = cli_recorder[0]
    assert "--define" in argv
    assert "--charge" in argv and "2" in argv


def test_tadf_xtb_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_tadf_xtb(extra_args="--foo bar", allow_mutate=True)
    assert cli_recorder[0] == ["tadf_xtb", "--foo", "bar"]


def test_hyperpol_executes_with_allow_mutate(cli_recorder):
    ops_server.tool_hyperpol(extra_args="--baz", allow_mutate=True)
    assert cli_recorder[0] == ["hyperpol", "--baz"]


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def test_generate_mcp_config_default():
    cfg = ops_config.generate_mcp_config()
    servers = cfg["mcpServers"]
    assert "delfin-ops" in servers
    entry = servers["delfin-ops"]
    assert entry["command"] == sys.executable
    assert entry["args"] == ["-m", "delfin.ops_server"]


def test_generate_mcp_config_with_workspace():
    cfg = ops_config.generate_mcp_config(workspace="/tmp/ws")
    args = cfg["mcpServers"]["delfin-ops"]["args"]
    assert args == ["-m", "delfin.ops_server", "--workspace", "/tmp/ws"]


def test_ensure_mcp_config_writes_new_file(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    out_path = ops_config.ensure_mcp_config(workspace=str(tmp_path))
    assert Path(out_path).exists()
    data = json.loads(Path(out_path).read_text())
    assert "delfin-ops" in data["mcpServers"]


def test_ensure_mcp_config_merges_existing(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    base_path = tmp_path / "existing.json"
    base_path.write_text(json.dumps({
        "mcpServers": {"other-server": {"command": "x"}}
    }))
    out_path = ops_config.ensure_mcp_config(
        existing_mcp_config=str(base_path), workspace=str(tmp_path),
    )
    data = json.loads(Path(out_path).read_text())
    assert "other-server" in data["mcpServers"]
    assert "delfin-ops" in data["mcpServers"]


def test_ensure_mcp_config_recovers_from_corrupt_existing(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    bad_path = tmp_path / "bad.json"
    bad_path.write_text("{ not json")
    out_path = ops_config.ensure_mcp_config(existing_mcp_config=str(bad_path))
    data = json.loads(Path(out_path).read_text())
    # falls back to a fresh ops-only config
    assert list(data["mcpServers"].keys()) == ["delfin-ops"]


# ---------------------------------------------------------------------------
# On-demand operational-pattern lookup tools
# ---------------------------------------------------------------------------

def test_tool_list_dashboard_patterns_returns_names():
    out = ops_server.tool_list_dashboard_patterns()
    assert "batch" in out
    assert "smart_recalc" in out
    assert "Available dashboard pattern names" in out


def test_tool_get_dashboard_pattern_batch():
    out = ops_server.tool_get_dashboard_pattern("batch")
    assert "/batch from-calc" in out
    assert "Never" in out


def test_tool_get_dashboard_pattern_unknown_returns_hint():
    out = ops_server.tool_get_dashboard_pattern("notreal")
    assert "Unknown pattern" in out
    assert "batch" in out
