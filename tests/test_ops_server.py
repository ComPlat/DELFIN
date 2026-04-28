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

from delfin import api as delfin_api
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


# ---------------------------------------------------------------------------
# P1 — output parsing tools
# ---------------------------------------------------------------------------


_FAKE_ORCA_OUT = """\
                                  * O   R   C   A *

|  1> ! PBE0 def2-TZVP D4 RIJCOSX def2/J Opt Freq

Number of atoms                            ...     12

  *** SCF CONVERGED AFTER  17 CYCLES ***

FINAL SINGLE POINT ENERGY      -123.456789012345

  *** OPTIMIZATION RUN DONE ***

THERMOCHEMISTRY AT 298.15K

Temperature        ...    298.15 K
Pressure           ...      1.00 atm
Zero point energy                ...      0.123456 Eh   77.42 kcal/mol
Total thermal correction              0.012345 Eh    7.75 kcal/mol
Total Enthalpy                   ...   -123.444444 Eh
Final entropy term               ...     -0.022222 Eh  -13.94 kcal/mol
G-E(el)                          ...      0.090000 Eh   56.48 kcal/mol
Final Gibbs free energy          ...   -123.366000 Eh

Number of imaginary frequencies   ...     0

TOTAL RUN TIME: 0 days 0 hours 12 minutes 34 seconds 567 msec
"""


@pytest.fixture
def fake_calc_folder(tmp_path):
    """Create a folder with a single fake ORCA .out file."""
    out = tmp_path / "out_main.out"
    out.write_text(_FAKE_ORCA_OUT, encoding="utf-8")
    return tmp_path


def test_tool_parse_orca_output_returns_json_with_energies(fake_calc_folder):
    out_path = fake_calc_folder / "out_main.out"
    txt = ops_server.tool_parse_orca_output(str(out_path))
    data = json.loads(txt)
    assert data["final_single_point"] == pytest.approx(-123.456789012345)
    assert data["gibbs_free_energy"] == pytest.approx(-123.366000)
    assert data["zpe"] == pytest.approx(0.123456)
    assert data["scf_converged"] is True
    assert data["opt_converged"] is True
    assert data["n_atoms"] == 12
    assert data["functional"] == "PBE0"
    assert data["basis"] == "def2-TZVP"


def test_tool_parse_orca_output_missing_file():
    txt = ops_server.tool_parse_orca_output("/nope/missing.out")
    data = json.loads(txt)
    assert data["error_summary"] == "file not found"


def test_tool_extract_thermochem_returns_block(fake_calc_folder):
    txt = ops_server.tool_extract_thermochem(str(fake_calc_folder))
    data = json.loads(txt)
    assert data["temperature_k"] == 298.15
    assert data["pressure_atm"] == 1.00
    assert data["zpe"] == pytest.approx(0.123456)
    assert data["final_gibbs"] == pytest.approx(-123.366000)


def test_tool_extract_energy_table_default_props(fake_calc_folder):
    txt = ops_server.tool_extract_energy_table(str(fake_calc_folder))
    rows = json.loads(txt)
    assert len(rows) == 1
    assert rows[0]["status"] == "ok"
    assert rows[0]["gibbs"] == pytest.approx(-123.366000)
    assert rows[0]["single_point"] == pytest.approx(-123.456789012345)


def test_tool_extract_energy_table_explicit_props(fake_calc_folder):
    txt = ops_server.tool_extract_energy_table(
        str(fake_calc_folder), properties="gibbs,scf_converged,walltime_s",
    )
    rows = json.loads(txt)
    assert "gibbs" in rows[0]
    assert "scf_converged" in rows[0]
    assert "walltime_s" in rows[0]
    assert rows[0]["scf_converged"] is True
    # 12 minutes 34.567 seconds = 754.567 s
    assert rows[0]["walltime_s"] == pytest.approx(754.567, abs=0.5)


def test_tool_extract_energy_table_missing_folder(tmp_path):
    """Non-existent folders show up as status=missing, not crashes."""
    txt = ops_server.tool_extract_energy_table(str(tmp_path / "nope"))
    rows = json.loads(txt)
    assert rows[0]["status"] == "missing"
    assert rows[0]["gibbs"] is None


def test_tool_find_orca_errors_picks_up_oom(tmp_path):
    bad = tmp_path / "out_oom.out"
    bad.write_text(
        "Some normal lines\n"
        "Detected oom_kill events in StepId=4089862.batch\n"
        "more lines\n",
        encoding="utf-8",
    )
    txt = ops_server.tool_find_orca_errors(str(tmp_path))
    errors = json.loads(txt)
    assert len(errors) == 1
    assert errors[0]["type"] == "oom"
    assert "Increase --mem" in errors[0]["suggestion"]


def test_tool_find_calculation_extreme_picks_lowest(tmp_path):
    """End-to-end: 3 folders with different Gibbs → return lowest 2."""
    for name, gibbs in [("a", -100.5), ("b", -150.0), ("c", -120.7)]:
        folder = tmp_path / name
        folder.mkdir()
        (folder / "calc.out").write_text(
            f"  *** SCF CONVERGED ***\n"
            f"FINAL SINGLE POINT ENERGY    {gibbs - 0.5}\n"
            f"Final Gibbs free energy            ...    {gibbs} Eh\n",
            encoding="utf-8",
        )
    folders_csv = ",".join(str(tmp_path / x) for x in ("a", "b", "c"))
    txt = ops_server.tool_find_calculation_extreme(
        folders_csv, property="gibbs", extreme="min", n=2,
    )
    rows = json.loads(txt)
    assert len(rows) == 2
    # Lowest Gibbs first
    assert rows[0]["gibbs"] == pytest.approx(-150.0)
    assert rows[1]["gibbs"] == pytest.approx(-120.7)


def test_tool_find_calculation_extreme_excludes_unparseable(tmp_path):
    """A folder without ORCA output is silently dropped from the ranking."""
    good = tmp_path / "good"
    good.mkdir()
    (good / "calc.out").write_text(
        "Final Gibbs free energy            ...    -50.0 Eh\n",
        encoding="utf-8",
    )
    bad = tmp_path / "bad"
    bad.mkdir()  # no .out → dropped from ranking
    folders_csv = f"{good},{bad}"
    txt = ops_server.tool_find_calculation_extreme(
        folders_csv, property="gibbs",
    )
    rows = json.loads(txt)
    assert len(rows) == 1
    assert rows[0]["folder"] == str(good)


# ---------------------------------------------------------------------------
# P1 — statistical plot tools (PNG output, auto-displayed in chat)
# ---------------------------------------------------------------------------


@pytest.fixture
def fake_calc_folders(tmp_path):
    """Five fake calc folders with varying Gibbs / FSPE."""
    folders = []
    for name, gibbs, spe in [
        ("a", -100.5, -101.0),
        ("b", -150.0, -151.2),
        ("c", -120.7, -121.3),
        ("d", -110.0, -110.4),
        ("e", -135.5, -136.1),
    ]:
        folder = tmp_path / name
        folder.mkdir()
        (folder / "calc.out").write_text(
            f"  *** SCF CONVERGED ***\n"
            f"FINAL SINGLE POINT ENERGY    {spe}\n"
            f"Final Gibbs free energy            ...    {gibbs} Eh\n",
            encoding="utf-8",
        )
        folders.append(folder)
    return folders


def test_tool_plot_energy_distribution_writes_png(fake_calc_folders, tmp_path,
                                                   monkeypatch):
    """Default histogram plot lands in agent_workspace and reports stats."""
    # Redirect the workspace dir into tmp_path so we don't pollute $HOME.
    monkeypatch.setattr(
        delfin_api, "_default_plot_dir",
        lambda: str(tmp_path / "ws"),
    )
    folders_csv = ",".join(str(f) for f in fake_calc_folders)
    txt = ops_server.tool_plot_energy_distribution(folders_csv)
    data = json.loads(txt)
    assert data["error"] == ""
    assert data["path"].endswith(".png")
    assert Path(data["path"]).exists()
    assert data["n_points"] == 5
    # Statistics for both properties (defaults: gibbs + single_point)
    assert data["statistics"]["gibbs"]["n"] == 5
    assert data["statistics"]["gibbs"]["min"] == pytest.approx(-150.0)
    assert data["statistics"]["gibbs"]["max"] == pytest.approx(-100.5)


def test_tool_plot_energy_distribution_custom_props_and_type(
        fake_calc_folders, tmp_path, monkeypatch):
    monkeypatch.setattr(
        delfin_api, "_default_plot_dir",
        lambda: str(tmp_path / "ws"),
    )
    folders_csv = ",".join(str(f) for f in fake_calc_folders)
    txt = ops_server.tool_plot_energy_distribution(
        folders_csv, properties="gibbs", plot_type="boxplot",
    )
    data = json.loads(txt)
    assert data["error"] == ""
    assert "energy_boxplot" in data["path"]
    assert Path(data["path"]).exists()
    assert data["properties"] == ["gibbs"]


def test_tool_plot_energy_distribution_no_data(tmp_path, monkeypatch):
    """No parseable folders → error message, no PNG written."""
    monkeypatch.setattr(
        delfin_api, "_default_plot_dir",
        lambda: str(tmp_path / "ws"),
    )
    empty = tmp_path / "empty"
    empty.mkdir()
    txt = ops_server.tool_plot_energy_distribution(str(empty))
    data = json.loads(txt)
    assert data["error"] != ""
    assert data["path"] == ""
    assert data["n_points"] == 0


def test_tool_plot_energy_correlation_writes_png(fake_calc_folders, tmp_path,
                                                  monkeypatch):
    """Scatter plot includes Pearson r in the result statistics."""
    monkeypatch.setattr(
        delfin_api, "_default_plot_dir",
        lambda: str(tmp_path / "ws"),
    )
    folders_csv = ",".join(str(f) for f in fake_calc_folders)
    txt = ops_server.tool_plot_energy_correlation(
        folders_csv, x="single_point", y="gibbs",
    )
    data = json.loads(txt)
    assert data["error"] == ""
    assert Path(data["path"]).exists()
    assert data["properties"] == ["single_point", "gibbs"]
    assert data["n_points"] == 5
    # spe and gibbs are nearly proportional → r close to 1
    assert data["statistics"]["pearson_r"] > 0.95


def test_tool_plot_energy_correlation_no_overlap(tmp_path, monkeypatch):
    """If no folder has BOTH x and y set, return error not exception."""
    monkeypatch.setattr(
        delfin_api, "_default_plot_dir",
        lambda: str(tmp_path / "ws"),
    )
    empty = tmp_path / "nada"
    empty.mkdir()
    txt = ops_server.tool_plot_energy_correlation(str(empty))
    data = json.loads(txt)
    assert data["error"] != ""
    assert data["path"] == ""
