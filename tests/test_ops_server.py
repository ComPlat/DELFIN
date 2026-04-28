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
# Imaginary-frequency extraction + functional comparison
# ---------------------------------------------------------------------------


def _write_freq_out(
    folder, *, gibbs=-100.0, spe=-100.5, functional="PBE0",
    basis="def2-TZVP", imag_freqs=None,
):
    """Write a synthetic ORCA freq output. ``imag_freqs`` is a list of cm-1."""
    folder.mkdir(parents=True, exist_ok=True)
    lines = [
        f"|  1> ! {functional} {basis} Opt Freq",
        "Number of atoms                                 ...     6",
        "  *** SCF CONVERGED ***",
        f"FINAL SINGLE POINT ENERGY    {spe}",
        "VIBRATIONAL FREQUENCIES",
        "-----------------------",
    ]
    idx = 0
    for f in (imag_freqs or []):
        lines.append(f"  {idx}:    {f:.2f} cm**-1 ***imaginary mode***")
        idx += 1
    for f in (250.0, 800.0, 1450.0):
        lines.append(f"  {idx}:    {f:.2f} cm**-1")
        idx += 1
    lines += [
        "THERMOCHEMISTRY",
        "Temperature        ...    298.15 K",
        "Pressure           ...      1.00 atm",
        f"Final Gibbs free energy          ...    {gibbs} Eh",
    ]
    if imag_freqs is None:
        lines.append("Number of imaginary frequencies   ...     0")
    else:
        lines.append(
            f"Number of imaginary frequencies   ...     {len(imag_freqs)}"
        )
    (folder / "calc.out").write_text("\n".join(lines), encoding="utf-8")


def test_tool_extract_imaginary_frequencies_minimum(tmp_path):
    """No imaginary modes → is_minimum=True, n_imag=0."""
    _write_freq_out(tmp_path / "calc")
    txt = ops_server.tool_extract_imaginary_frequencies(
        str(tmp_path / "calc"),
    )
    data = json.loads(txt)
    assert data["n_imag"] == 0
    assert data["is_minimum"] is True
    assert data["is_ts"] is False
    assert data["most_negative"] is None
    assert data["modes"] == []
    assert data["error"] is None


def test_tool_extract_imaginary_frequencies_ts(tmp_path):
    """Single imaginary mode → is_ts=True, most_negative reported."""
    _write_freq_out(tmp_path / "calc", imag_freqs=[-450.7])
    txt = ops_server.tool_extract_imaginary_frequencies(
        str(tmp_path / "calc"),
    )
    data = json.loads(txt)
    assert data["n_imag"] == 1
    assert data["is_minimum"] is False
    assert data["is_ts"] is True
    assert data["most_negative"] == pytest.approx(-450.7)
    assert len(data["modes"]) == 1
    assert data["modes"][0]["frequency_cm"] == pytest.approx(-450.7)


def test_tool_extract_imaginary_frequencies_multi(tmp_path):
    """Multiple imag modes → not a minimum nor a TS."""
    _write_freq_out(
        tmp_path / "calc", imag_freqs=[-200.5, -55.1, -22.3],
    )
    data = json.loads(
        ops_server.tool_extract_imaginary_frequencies(str(tmp_path / "calc"))
    )
    assert data["n_imag"] == 3
    assert data["is_minimum"] is False
    assert data["is_ts"] is False
    # most_negative is the smallest (most negative) value
    assert data["most_negative"] == pytest.approx(-200.5)


def test_tool_extract_imaginary_frequencies_no_freq_block(tmp_path):
    """Folder without VIBRATIONAL FREQUENCIES → error string returned."""
    folder = tmp_path / "calc"
    folder.mkdir()
    (folder / "calc.out").write_text(
        "  *** SCF CONVERGED ***\nFINAL SINGLE POINT ENERGY  -1.0\n",
        encoding="utf-8",
    )
    data = json.loads(
        ops_server.tool_extract_imaginary_frequencies(str(folder))
    )
    assert data["error"] is not None
    assert "Freq" in data["error"]
    assert data["n_imag"] == 0


def test_tool_extract_imaginary_frequencies_missing_folder(tmp_path):
    data = json.loads(
        ops_server.tool_extract_imaginary_frequencies(str(tmp_path / "nope"))
    )
    assert data["error"] == "folder missing"


def test_tool_compare_calculations_method_and_delta(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    _write_freq_out(a, gibbs=-100.0, spe=-101.0, functional="BP86")
    _write_freq_out(b, gibbs=-100.5, spe=-101.6, functional="PBE0")
    txt = ops_server.tool_compare_calculations(str(a), str(b))
    data = json.loads(txt)
    assert data["a_functional"] == "BP86"
    assert data["b_functional"] == "PBE0"
    assert data["method_match"] is False
    assert data["basis_match"] is True
    # delta_gibbs_kcal = (B - A) * 627.5 = -0.5 * 627.5 ≈ -313.75 kcal/mol
    assert data["delta_gibbs_kcal"] == pytest.approx(
        (-100.5 - -100.0) * 627.5094740631,
    )
    assert any("BP86" in n and "PBE0" in n for n in data["notes"])


def test_tool_compare_calculations_with_imag_warning(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    _write_freq_out(a)
    _write_freq_out(b, imag_freqs=[-300.0])
    data = json.loads(
        ops_server.tool_compare_calculations(str(a), str(b))
    )
    assert data["b_n_imag"] == 1
    assert any("imag" in n.lower() for n in data["notes"])


def test_tool_compare_across_functionals_sorted_by_gibbs(tmp_path):
    """Three folders with different functionals sorted by Gibbs ascending."""
    folders = []
    for name, func, gibbs in [
        ("bp", "BP86", -100.0),
        ("pbe", "PBE0", -150.5),
        ("b3", "B3LYP", -120.2),
    ]:
        f = tmp_path / name
        _write_freq_out(f, functional=func, gibbs=gibbs)
        folders.append(str(f))
    txt = ops_server.tool_compare_across_functionals(",".join(folders))
    rows = json.loads(txt)
    assert len(rows) == 3
    assert rows[0]["functional"] == "PBE0"
    assert rows[0]["gibbs"] == pytest.approx(-150.5)
    assert rows[-1]["functional"] == "BP86"


def test_tool_compare_across_functionals_skip_imag(tmp_path):
    """include_imag=False → n_imag/is_minimum stay None even for true minima."""
    f = tmp_path / "calc"
    _write_freq_out(f)
    txt = ops_server.tool_compare_across_functionals(
        str(f), include_imag=False,
    )
    rows = json.loads(txt)
    assert rows[0]["n_imag"] is None
    assert rows[0]["is_minimum"] is None
    # gibbs/spe still extracted even when imag is skipped
    assert rows[0]["gibbs"] == pytest.approx(-100.0)


def test_tool_compare_across_functionals_status_codes(tmp_path):
    """Missing folder → status=missing; folder w/o output → no_output."""
    missing = tmp_path / "missing"
    empty = tmp_path / "empty"
    empty.mkdir()
    txt = ops_server.tool_compare_across_functionals(
        f"{missing},{empty}", sort_by="folder",
    )
    rows = json.loads(txt)
    statuses = sorted(r["status"] for r in rows)
    assert statuses == ["missing", "no_output"]


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


# ---------------------------------------------------------------------------
# Tool / widget catalog discovery
# ---------------------------------------------------------------------------


def test_tool_list_tools_returns_full_catalog():
    txt = ops_server.tool_list_tools()
    catalog = json.loads(txt)
    assert isinstance(catalog, list)
    # Sanity: at least the canonical tools are listed
    names = {e["name"] for e in catalog}
    assert {"parse_orca_output", "submit_calculation",
            "validate_orca_input", "list_dashboard_widgets",
            "list_tools", "describe_tool"}.issubset(names)


def test_tool_list_tools_filter_by_category():
    txt = ops_server.tool_list_tools(category="parsing")
    rows = json.loads(txt)
    assert all(r["category"] == "parsing" for r in rows)
    # parsing has parse_orca_output etc.
    assert any(r["name"] == "parse_orca_output" for r in rows)


def test_tool_list_tools_filter_by_query():
    txt = ops_server.tool_list_tools(query="cancel")
    rows = json.loads(txt)
    assert any("cancel" in r["name"].lower() for r in rows)


def test_tool_describe_tool_returns_docstring():
    txt = ops_server.tool_describe_tool("validate_orca_input")
    info = json.loads(txt)
    assert info["name"] == "validate_orca_input"
    assert info["category"] == "validation"
    assert "Sanity-check" in info["docstring"]


def test_tool_describe_tool_unknown_lists_options():
    txt = ops_server.tool_describe_tool("nopepe")
    info = json.loads(txt)
    assert "error" in info
    assert "available" in info
    assert "parse_orca_output" in info["available"]


# ---------------------------------------------------------------------------
# Widget catalog
# ---------------------------------------------------------------------------


def test_tool_list_dashboard_widgets_all_tabs():
    txt = ops_server.tool_list_dashboard_widgets()
    rows = json.loads(txt)
    tabs = {r["tab"] for r in rows}
    assert {"submit", "orca", "calc", "agent"}.issubset(tabs)


def test_tool_list_dashboard_widgets_filter_orca():
    txt = ops_server.tool_list_dashboard_widgets(tab="orca")
    rows = json.loads(txt)
    assert all(r["tab"] == "orca" for r in rows)
    assert any(r["name"] == "orca-method" for r in rows)
    assert any(r["name"] == "orca-preview" for r in rows)


def test_tool_get_widget_options_returns_dropdown_values():
    txt = ops_server.tool_get_widget_options("orca-method")
    options = json.loads(txt)
    assert "PBE0" in options
    assert "BP86" in options


def test_tool_get_widget_options_empty_for_textarea():
    """Non-dropdown widgets return [] — agent shouldn't try /ui options."""
    txt = ops_server.tool_get_widget_options("orca-preview")
    assert json.loads(txt) == []


# ---------------------------------------------------------------------------
# ORCA-input validation
# ---------------------------------------------------------------------------


def test_tool_validate_orca_input_clean_returns_no_errors():
    inp = (
        "! PBE0 def2-TZVP D4 OPT Freq\n"
        "%pal nprocs 12 end\n"
        "%maxcore 4000\n"
        "* xyz 0 1\n"
        "C 0 0 0\n"
        "*\n"
    )
    txt = ops_server.tool_validate_orca_input(inp)
    issues = json.loads(txt)
    severities = {i["severity"] for i in issues}
    assert "error" not in severities


def test_tool_validate_orca_input_empty_text_is_error():
    txt = ops_server.tool_validate_orca_input("")
    issues = json.loads(txt)
    assert any(i["code"] == "empty" for i in issues)


def test_tool_validate_orca_input_missing_xyz_block():
    inp = "! PBE0 def2-TZVP\n%pal nprocs 4 end\n"
    txt = ops_server.tool_validate_orca_input(inp)
    issues = json.loads(txt)
    assert any(i["code"] == "no_xyz_block" for i in issues)


def test_tool_validate_orca_input_pal_inline_form_passes():
    """The S2 OOM bug — inline %pal nprocs 40 end must validate clean."""
    inp = (
        "! PBE0 def2-TZVP\n"
        "%pal nprocs 40 end\n"
        "%maxcore 6000\n"
        "* xyz 0 1\nC 0 0 0\n*\n"
    )
    txt = ops_server.tool_validate_orca_input(inp)
    issues = json.loads(txt)
    assert not any(i["code"] == "pal_no_nprocs" for i in issues)


def test_tool_validate_orca_input_zora_basis_mismatch():
    """ZORA Hamiltonian + plain def2 basis → warning."""
    inp = (
        "! PBE0 def2-TZVP ZORA\n"
        "%pal nprocs 4 end\n"
        "* xyz 0 1\nFe 0 0 0\n*\n"
    )
    txt = ops_server.tool_validate_orca_input(inp)
    issues = json.loads(txt)
    assert any(i["code"] == "zora_basis_mismatch" for i in issues)


def test_tool_validate_orca_input_neb_without_block():
    inp = (
        "! NEB-TS PBE0 def2-TZVP\n"
        "%pal nprocs 4 end\n"
        "* xyz 0 1\nC 0 0 0\n*\n"
    )
    txt = ops_server.tool_validate_orca_input(inp)
    issues = json.loads(txt)
    assert any(i["code"] == "neb_no_block" for i in issues)


# ---------------------------------------------------------------------------
# Job lifecycle (mutation gate + dry-run)
# ---------------------------------------------------------------------------


def test_tool_submit_calculation_dry_run_no_sbatch(tmp_path):
    """Default allow_mutate=False returns a 'would submit' preview."""
    job_dir = tmp_path / "job1"
    job_dir.mkdir()
    txt = ops_server.tool_submit_calculation(
        str(job_dir), pal=12, maxcore=4000, time_limit="06:00:00",
    )
    data = json.loads(txt)
    assert data["submitted"] is False
    assert "Would submit" in data["message"]
    assert "06:00:00" in data["message"]


def test_tool_submit_calculation_missing_folder():
    txt = ops_server.tool_submit_calculation("/nope/missing")
    data = json.loads(txt)
    assert data["submitted"] is False
    assert "not found" in data["error"]


def test_tool_cancel_calculation_dry_run_no_scancel():
    txt = ops_server.tool_cancel_calculation("12345")
    data = json.loads(txt)
    assert data["ok"] is False
    assert data.get("skipped") is True
    assert "Would cancel" in data["message"]


def test_tool_cancel_calculation_empty_id_is_error():
    txt = ops_server.tool_cancel_calculation("")
    data = json.loads(txt)
    assert data["ok"] is False
    assert "no job_id" in data["error"]


def test_tool_list_active_calculations_returns_list():
    """Backend may or may not have jobs — must return a list shape."""
    txt = ops_server.tool_list_active_calculations()
    data = json.loads(txt)
    assert isinstance(data, list)


# ---------------------------------------------------------------------------
# ORCA-manual lookup + literature indexing
# ---------------------------------------------------------------------------

def test_check_orca_manual_indexed_no_index(tmp_path, monkeypatch):
    """No index file → indexed=False with a clear hint."""
    monkeypatch.setattr(
        "delfin.doc_server.indexer.get_default_index_path",
        lambda: tmp_path / "missing.json",
    )
    txt = ops_server.tool_check_orca_manual_indexed()
    data = json.loads(txt)
    assert data["indexed"] is False
    assert "Drop a PDF" in data["hint"] or "drop" in data["hint"].lower()


def test_check_orca_manual_indexed_finds_orca_doc(tmp_path, monkeypatch):
    """An index containing an ORCA-titled doc returns indexed=True."""
    fake_index = tmp_path / "doc_index.json"
    fake_index.write_text(json.dumps({
        "documents": {
            "orca_manual_6_1_1_delfin": {
                "title": "ORCA 6.1.1 Manual",
                "sections": {"ch1": {"text": "..."}},
            },
        },
    }))
    monkeypatch.setattr(
        "delfin.doc_server.indexer.get_default_index_path",
        lambda: fake_index,
    )
    txt = ops_server.tool_check_orca_manual_indexed()
    data = json.loads(txt)
    assert data["indexed"] is True
    assert "orca_manual_6_1_1_delfin" in data["doc_ids"]


def test_check_orca_manual_indexed_empty_index_is_unindexed(tmp_path, monkeypatch):
    """Index file exists but has only non-ORCA docs."""
    fake_index = tmp_path / "doc_index.json"
    fake_index.write_text(json.dumps({
        "documents": {
            "some_random_paper": {"title": "Random Paper"},
        },
    }))
    monkeypatch.setattr(
        "delfin.doc_server.indexer.get_default_index_path",
        lambda: fake_index,
    )
    txt = ops_server.tool_check_orca_manual_indexed()
    data = json.loads(txt)
    assert data["indexed"] is False
    assert "Literature tab" in data["hint"]


def test_index_new_pdf_missing_file():
    """Non-existent path → ok=False, no exception."""
    txt = ops_server.tool_index_new_pdf("/nope/missing.pdf")
    data = json.loads(txt)
    assert data["ok"] is False
    assert "not found" in data["error"]


def test_index_new_pdf_wrong_extension(tmp_path):
    not_pdf = tmp_path / "x.txt"
    not_pdf.write_text("hi")
    txt = ops_server.tool_index_new_pdf(str(not_pdf))
    data = json.loads(txt)
    assert data["ok"] is False
    assert "PDF" in data["error"]


# ---------------------------------------------------------------------------
# DELFIN-feature explainer
# ---------------------------------------------------------------------------

def test_list_delfin_features_returns_known_concepts():
    txt = ops_server.tool_list_delfin_features()
    rows = json.loads(txt)
    names = {r["name"] for r in rows}
    assert {"control_keys", "occupier", "smart_recalc", "modes",
            "guppy", "permissions"}.issubset(names)


def test_list_delfin_features_filter_workflow():
    txt = ops_server.tool_list_delfin_features(category="workflow")
    rows = json.loads(txt)
    assert all(r["category"] == "workflow" for r in rows)


def test_explain_delfin_feature_returns_summary_and_seealso():
    txt = ops_server.tool_explain_delfin_feature("smart_recalc")
    info = json.loads(txt)
    assert info["name"] == "smart_recalc"
    assert "calc-override-btn" in info["summary"]
    assert info["see_also"]  # at least one source pointer


def test_explain_delfin_feature_case_insensitive():
    a = ops_server.tool_explain_delfin_feature("OCCUPIER")
    b = ops_server.tool_explain_delfin_feature("occupier")
    assert json.loads(a)["name"] == "occupier"
    assert a == b


def test_explain_delfin_feature_unknown_returns_candidates():
    txt = ops_server.tool_explain_delfin_feature("notreal")
    info = json.loads(txt)
    assert "error" in info
    assert "available" in info
    assert "occupier" in info["available"]


def test_explain_delfin_feature_fuzzy_match():
    """A unique substring match still resolves cleanly."""
    txt = ops_server.tool_explain_delfin_feature("recalc")
    info = json.loads(txt)
    # smart_recalc is the only entry with 'recalc' in the name
    assert info.get("name") == "smart_recalc"


def test_explain_delfin_feature_normalises_separators():
    a = ops_server.tool_explain_delfin_feature("smart-recalc")
    b = ops_server.tool_explain_delfin_feature("smart recalc")
    assert json.loads(a)["name"] == "smart_recalc"
    assert json.loads(b)["name"] == "smart_recalc"


# ---------------------------------------------------------------------------
# PDF on-demand tools (read / search / extract_section / list_literature)
# ---------------------------------------------------------------------------

# We don't need a real PDF — mock _extract_pdf_text so tests don't depend
# on pypdf installation or PDF rendering. The tool API treats the result
# the same regardless of the real-vs-mocked source.
_FAKE_PDF_PAGES = [
    {"page": 1, "text": (
        "1 Introduction\n"
        "ORCA is a quantum chemistry program package.\n"
        "It supports DFT, CC, and many other methods.\n"
    )},
    {"page": 2, "text": (
        "2 Density Functional Theory\n"
        "Several functionals are available, including PBE0,\n"
        "BP86, and B3LYP. The default basis is def2-TZVP.\n"
    )},
    {"page": 3, "text": (
        "3 Geometry Optimization\n"
        "Use OPT keyword. For tight convergence: TIGHTOPT.\n"
        "Combined with FREQ for transition states: NEB-TS.\n"
    )},
]


@pytest.fixture
def fake_pdf(tmp_path, monkeypatch):
    """Create a non-empty placeholder PDF + monkeypatch _extract_pdf_text."""
    pdf_path = tmp_path / "fake_orca_manual.pdf"
    pdf_path.write_bytes(b"%PDF-1.4\n%fake\n")
    monkeypatch.setattr(
        "delfin.doc_server.indexer._extract_pdf_text",
        lambda p, quiet=True: _FAKE_PDF_PAGES,
    )
    return pdf_path


def test_tool_read_pdf_full_document(fake_pdf):
    txt = ops_server.tool_read_pdf(str(fake_pdf))
    data = json.loads(txt)
    assert data["error"] == ""
    assert data["n_pages_total"] == 3
    assert data["n_pages_read"] == 3
    assert "Density Functional Theory" in data["text"]
    assert "PAGE 1" in data["text"]


def test_tool_read_pdf_page_range(fake_pdf):
    txt = ops_server.tool_read_pdf(str(fake_pdf), pages="2-3")
    data = json.loads(txt)
    assert data["n_pages_read"] == 2
    assert "Introduction" not in data["text"]
    assert "Density Functional" in data["text"]
    assert "Geometry Optimization" in data["text"]


def test_tool_read_pdf_page_csv(fake_pdf):
    txt = ops_server.tool_read_pdf(str(fake_pdf), pages="1,3")
    data = json.loads(txt)
    assert data["n_pages_read"] == 2
    assert "Introduction" in data["text"]
    assert "Density Functional" not in data["text"]
    assert "Geometry Optimization" in data["text"]


def test_tool_read_pdf_truncates(fake_pdf):
    txt = ops_server.tool_read_pdf(str(fake_pdf), max_chars=80)
    data = json.loads(txt)
    assert data["truncated"] is True
    assert "[truncated" in data["text"]


def test_tool_read_pdf_missing_file():
    txt = ops_server.tool_read_pdf("/nope/missing.pdf")
    data = json.loads(txt)
    assert data["error"] == "file not found"
    assert data["text"] == ""


def test_tool_read_pdf_wrong_extension(tmp_path):
    bad = tmp_path / "x.txt"
    bad.write_text("hello")
    txt = ops_server.tool_read_pdf(str(bad))
    data = json.loads(txt)
    assert data["error"] == "not a PDF file"


def test_tool_search_pdf_local_finds_matches(fake_pdf):
    txt = ops_server.tool_search_pdf_local(str(fake_pdf), query="def2")
    data = json.loads(txt)
    assert data["error"] == ""
    assert len(data["hits"]) >= 1
    hit = data["hits"][0]
    assert "def2" in hit["snippet"].lower()
    assert hit["page"] == 2


def test_tool_search_pdf_local_case_insensitive_default(fake_pdf):
    txt = ops_server.tool_search_pdf_local(str(fake_pdf), query="TIGHTOPT")
    data = json.loads(txt)
    assert len(data["hits"]) >= 1


def test_tool_search_pdf_local_case_sensitive_misses(fake_pdf):
    txt = ops_server.tool_search_pdf_local(
        str(fake_pdf), query="density",  # title is "Density"
        case_sensitive=True,
    )
    data = json.loads(txt)
    assert data["hits"] == []


def test_tool_search_pdf_local_max_hits_cap(fake_pdf):
    txt = ops_server.tool_search_pdf_local(
        str(fake_pdf), query="ORCA", max_hits=1,
    )
    data = json.loads(txt)
    assert len(data["hits"]) <= 1


def test_tool_search_pdf_local_empty_query():
    txt = ops_server.tool_search_pdf_local("/nope/missing.pdf", query="")
    data = json.loads(txt)
    assert data["error"] == "empty query"


def test_tool_extract_pdf_section_finds_heading(fake_pdf):
    txt = ops_server.tool_extract_pdf_section(
        str(fake_pdf), heading="Geometry Optimization",
    )
    data = json.loads(txt)
    assert data["heading_found"] is True
    assert "OPT" in data["text"] or "TIGHTOPT" in data["text"]


def test_tool_extract_pdf_section_unknown_heading(fake_pdf):
    txt = ops_server.tool_extract_pdf_section(
        str(fake_pdf), heading="Quantum Tunneling Through Time",
    )
    data = json.loads(txt)
    assert data["heading_found"] is False
    assert "not found" in data["error"]


def test_tool_list_literature_files_filter_pdfs(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "delfin.doc_server.indexer.get_default_literature_dir",
        lambda: tmp_path,
    )
    (tmp_path / "paper.pdf").write_bytes(b"%PDF-1.4\n")
    (tmp_path / "notes.md").write_text("# stuff\n")
    (tmp_path / "image.png").write_bytes(b"\x89PNG\r\n")
    txt = ops_server.tool_list_literature_files()
    rows = json.loads(txt)
    names = sorted(r["name"] for r in rows)
    assert names == ["notes.md", "paper.pdf"]
    # png is filtered out


def test_tool_list_literature_files_recursive(tmp_path, monkeypatch):
    """Sub-folders are walked too."""
    monkeypatch.setattr(
        "delfin.doc_server.indexer.get_default_literature_dir",
        lambda: tmp_path,
    )
    sub = tmp_path / "papers"
    sub.mkdir()
    (sub / "a.pdf").write_bytes(b"%PDF-1.4\n")
    (tmp_path / "b.pdf").write_bytes(b"%PDF-1.4\n")
    txt = ops_server.tool_list_literature_files()
    rows = json.loads(txt)
    assert len(rows) == 2
    paths = {r["path"] for r in rows}
    assert any(p.endswith("/papers/a.pdf") for p in paths)


def test_tool_list_literature_files_missing_folder(tmp_path):
    """Non-existent explicit folder → empty list, no exception."""
    txt = ops_server.tool_list_literature_files(folder=str(tmp_path / "nope"))
    rows = json.loads(txt)
    assert rows == []


# ---------------------------------------------------------------------------
# Phase B — calc folder management (rename / create / move / archive / delete)
# ---------------------------------------------------------------------------


@pytest.fixture
def calc_tree(tmp_path):
    """A minimal project tree: calc/foo, calc/bar, archive/."""
    calc = tmp_path / "calc"
    arch = tmp_path / "archive"
    (calc / "foo").mkdir(parents=True)
    (calc / "bar").mkdir(parents=True)
    arch.mkdir()
    (calc / "foo" / "calc.out").write_text("done\n")
    (calc / "bar" / "calc.out").write_text("done\n")
    return {"root": tmp_path, "calc": calc, "arch": arch}


def test_tool_rename_calc_folder_blocks_without_allow(calc_tree):
    """Default-blocked: returns 'would rename' hint, no actual rename."""
    src = calc_tree["calc"] / "foo"
    txt = ops_server.tool_rename_calc_folder(
        str(src), "renamed",
        allowed_roots=str(calc_tree["calc"]),
    )
    data = json.loads(txt)
    assert data["ok"] is False
    assert data["skipped"] is True
    assert "Would rename" in data["message"]
    assert src.exists()
    assert not (calc_tree["calc"] / "renamed").exists()


def test_tool_rename_calc_folder_executes(calc_tree):
    src = calc_tree["calc"] / "foo"
    txt = ops_server.tool_rename_calc_folder(
        str(src), "renamed",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    )
    data = json.loads(txt)
    assert data["ok"] is True
    assert not src.exists()
    assert (calc_tree["calc"] / "renamed").exists()


def test_tool_rename_calc_folder_rejects_slash(calc_tree):
    src = calc_tree["calc"] / "foo"
    data = json.loads(ops_server.tool_rename_calc_folder(
        str(src), "renamed/sub",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "basename" in data["error"]


def test_tool_rename_calc_folder_rejects_archive(calc_tree):
    """Refuses to operate on anything inside archive/."""
    inside_arc = calc_tree["arch"] / "stale"
    inside_arc.mkdir()
    data = json.loads(ops_server.tool_rename_calc_folder(
        str(inside_arc), "moved",
        # NOTE we deliberately allow the archive root here so the path-
        # check passes, then the destructive-dir check should still kick in.
        allowed_roots=str(calc_tree["arch"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "archive" in data["error"]


def test_tool_rename_calc_folder_rejects_outside_roots(calc_tree, tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    data = json.loads(ops_server.tool_rename_calc_folder(
        str(outside), "x",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "outside allowed roots" in data["error"]


def test_tool_create_calc_folder_executes(calc_tree):
    txt = ops_server.tool_create_calc_folder(
        str(calc_tree["calc"]), "newone",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    )
    data = json.loads(txt)
    assert data["ok"] is True
    assert (calc_tree["calc"] / "newone").is_dir()


def test_tool_create_calc_folder_blocks_dup(calc_tree):
    data = json.loads(ops_server.tool_create_calc_folder(
        str(calc_tree["calc"]), "foo",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "already exists" in data["error"]


def test_tool_move_calc_folder_executes(calc_tree):
    src = calc_tree["calc"] / "foo"
    dst_parent = calc_tree["calc"] / "bar"
    txt = ops_server.tool_move_calc_folder(
        str(src), str(dst_parent),
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    )
    data = json.loads(txt)
    assert data["ok"] is True
    assert not src.exists()
    assert (dst_parent / "foo").is_dir()


def test_tool_move_calc_folder_blocks_dup(calc_tree):
    """Target already exists at destination → refusal."""
    (calc_tree["calc"] / "bar" / "foo").mkdir()
    data = json.loads(ops_server.tool_move_calc_folder(
        str(calc_tree["calc"] / "foo"),
        str(calc_tree["calc"] / "bar"),
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "already exists" in data["error"]


def test_tool_move_to_archive_executes(calc_tree):
    src = calc_tree["calc"] / "foo"
    txt = ops_server.tool_move_to_archive(
        str(src),
        archive_root=str(calc_tree["arch"]),
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    )
    data = json.loads(txt)
    assert data["ok"] is True
    assert not src.exists()
    assert (calc_tree["arch"] / "foo").is_dir()


def test_tool_move_to_archive_blocks_when_already_archived(calc_tree):
    """Refuse to 're-archive' something already inside archive/."""
    src = calc_tree["arch"] / "old"
    src.mkdir()
    data = json.loads(ops_server.tool_move_to_archive(
        str(src),
        archive_root=str(calc_tree["arch"]),
        allowed_roots=str(calc_tree["arch"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "already inside archive" in data["error"]


def test_tool_delete_calc_folder_three_lock(calc_tree):
    """All three locks (path, allow_mutate, confirm_token) must align."""
    target = calc_tree["calc"] / "foo"

    # Lock 1: missing confirm_token
    data = json.loads(ops_server.tool_delete_calc_folder(
        str(target),
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert target.exists()

    # Lock 2: confirm_token wrong
    data = json.loads(ops_server.tool_delete_calc_folder(
        str(target), confirm_token="bar",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert target.exists()

    # Lock 3: confirm_token right but allow_mutate missing
    data = json.loads(ops_server.tool_delete_calc_folder(
        str(target), confirm_token="foo",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=False,
    ))
    assert data["ok"] is False
    assert target.exists()

    # All three: actually deletes
    data = json.loads(ops_server.tool_delete_calc_folder(
        str(target), confirm_token="foo",
        allowed_roots=str(calc_tree["calc"]),
        allow_mutate=True,
    ))
    assert data["ok"] is True
    assert not target.exists()


def test_tool_delete_calc_folder_refuses_archive(calc_tree):
    """Cannot delete anything inside archive/, no matter the locks."""
    target = calc_tree["arch"] / "old"
    target.mkdir()
    data = json.loads(ops_server.tool_delete_calc_folder(
        str(target), confirm_token="old",
        allowed_roots=str(calc_tree["arch"]),
        allow_mutate=True,
    ))
    assert data["ok"] is False
    assert "archive" in data["error"]
    assert target.exists()


# ---------------------------------------------------------------------------
# Phase B-2 — bulk job control + recalc preparation
# ---------------------------------------------------------------------------


def test_tool_kill_all_user_jobs_dry_run(monkeypatch):
    """Default-blocked: lists candidates, no actual cancel calls."""
    fake_jobs = [
        {"job_id": "111", "status": "RUNNING"},
        {"job_id": "222", "status": "PENDING"},
        {"job_id": "333", "status": "RUNNING"},
    ]
    cancelled: list = []
    monkeypatch.setattr(
        "delfin.api.list_active_calculations",
        lambda: list(fake_jobs),
    )
    monkeypatch.setattr(
        "delfin.api.cancel_calculation",
        lambda jid, allow_mutate: cancelled.append(jid)
        or {"ok": True, "job_id": jid},
    )
    txt = ops_server.tool_kill_all_user_jobs()
    data = json.loads(txt)
    assert data["ok"] is False
    assert data["skipped"] is True
    assert sorted(data["candidates"]) == ["111", "222", "333"]
    assert cancelled == []


def test_tool_kill_all_user_jobs_only_running(monkeypatch):
    monkeypatch.setattr(
        "delfin.api.list_active_calculations",
        lambda: [
            {"job_id": "111", "status": "RUNNING"},
            {"job_id": "222", "status": "PENDING"},
        ],
    )
    monkeypatch.setattr(
        "delfin.api.cancel_calculation",
        lambda jid, allow_mutate: {"ok": True, "job_id": jid},
    )
    txt = ops_server.tool_kill_all_user_jobs(
        only_running=True, allow_mutate=True,
    )
    data = json.loads(txt)
    assert data["ok"] is True
    assert data["cancelled"] == ["111"]
    assert data["total"] == 1


def test_tool_kill_all_user_jobs_executes(monkeypatch):
    """allow_mutate=True → cancels every candidate."""
    monkeypatch.setattr(
        "delfin.api.list_active_calculations",
        lambda: [
            {"job_id": "111", "status": "R"},
            {"job_id": "222", "status": "PD"},
        ],
    )
    seen: list = []
    monkeypatch.setattr(
        "delfin.api.cancel_calculation",
        lambda jid, allow_mutate: seen.append((jid, allow_mutate))
        or {"ok": True, "job_id": jid},
    )
    txt = ops_server.tool_kill_all_user_jobs(allow_mutate=True)
    data = json.loads(txt)
    assert data["ok"] is True
    assert sorted(data["cancelled"]) == ["111", "222"]
    assert all(allow for _, allow in seen)


def _write_recalc_inputs(folder, *, pal=8, maxcore=2000):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "CONTROL.txt").write_text(
        f"functional = BP86\nPAL = {pal}\nmaxcore = {maxcore}\n",
        encoding="utf-8",
    )
    (folder / "calc.inp").write_text(
        f"! BP86 def2-SVP Opt\n%pal nprocs {pal} end\n%maxcore {maxcore}\n",
        encoding="utf-8",
    )


def test_tool_prepare_recalc_dry_run(tmp_path):
    """allow_mutate=False → returns plan with PAL/maxcore but no submission."""
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    txt = ops_server.tool_prepare_recalc(str(folder))
    plan = json.loads(txt)
    assert plan["mode"] == "delfin-recalc"
    assert plan["pal"] == 8
    assert plan["maxcore"] == 2000
    assert plan["valid"] is True
    assert plan["submitted_job_id"] is None
    assert "Would submit" in plan["submission_message"]


def test_tool_prepare_recalc_classic_mode(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    plan = json.loads(ops_server.tool_prepare_recalc(
        str(folder), mode="classic",
    ))
    assert plan["mode"] == "delfin-recalc-classic"
    assert plan["valid"] is True


def test_tool_prepare_recalc_override_needs_value(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    plan = json.loads(ops_server.tool_prepare_recalc(
        str(folder), mode="override",
    ))
    assert plan["valid"] is False
    assert any("STAGE=INDEX" in i for i in plan["issues"])


def test_tool_prepare_recalc_unknown_mode(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    plan = json.loads(ops_server.tool_prepare_recalc(
        str(folder), mode="bogus",
    ))
    assert plan["valid"] is False
    assert any("mode must be" in i for i in plan["issues"])


def test_tool_prepare_recalc_missing_pal(tmp_path):
    folder = tmp_path / "calc"
    folder.mkdir()
    (folder / "calc.inp").write_text("! BP86 def2-SVP\n", encoding="utf-8")
    plan = json.loads(ops_server.tool_prepare_recalc(str(folder)))
    assert plan["valid"] is False
    assert any("PAL" in i for i in plan["issues"])


def test_tool_prepare_recalc_bad_time_format(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    plan = json.loads(ops_server.tool_prepare_recalc(
        str(folder), time_limit="oneday",
    ))
    assert plan["valid"] is False
    assert any("HH:MM:SS" in i for i in plan["issues"])


def test_tool_prepare_recalc_executes(monkeypatch, tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)

    class FakeBackend:
        def submit_delfin(self, **kwargs):
            self.called_with = kwargs

            class R:
                returncode = 0
                job_id = "JOB42"
                stderr = ""

            return R()

    fake = FakeBackend()
    monkeypatch.setattr("delfin.api._resolve_backend", lambda: fake)
    plan = json.loads(ops_server.tool_prepare_recalc(
        str(folder), allow_mutate=True,
    ))
    assert plan["valid"] is True
    assert plan["submitted_job_id"] == "JOB42"
    assert "Submitted delfin-recalc" in plan["submission_message"]
    assert fake.called_with["mode"] == "delfin-recalc"
    assert fake.called_with["pal"] == 8


# ---------------------------------------------------------------------------
# Phase C — list_calc_options + run_calc_option
# ---------------------------------------------------------------------------


def test_tool_list_calc_options_control_txt():
    data = json.loads(ops_server.tool_list_calc_options("CONTROL.txt"))
    assert data["file_type"] == "control_txt"
    assert "Recalc" in data["options"]
    assert "Smart Recalc" in data["options"]


def test_tool_list_calc_options_inp_only_recalc():
    data = json.loads(ops_server.tool_list_calc_options("calc.inp"))
    assert data["file_type"] == "inp"
    assert data["options"] == ["Recalc"]


def test_tool_list_calc_options_xyz_full():
    data = json.loads(ops_server.tool_list_calc_options("structure.xyz"))
    assert data["file_type"] == "xyz"
    assert "Build Batch from XYZ" in data["options"]
    assert "hyperpol_xtb" in data["options"]


def test_tool_list_calc_options_out_set():
    data = json.loads(ops_server.tool_list_calc_options("calc.out"))
    assert data["file_type"] == "out"
    assert "MO Plot" in data["options"]
    assert "Print NMR" in data["options"]


def test_tool_list_calc_options_occupier():
    data = json.loads(ops_server.tool_list_calc_options("OCCUPIER.txt"))
    assert data["file_type"] == "occupier_txt"
    assert data["options"] == ["Override"]


def test_tool_list_calc_options_unknown_returns_other():
    data = json.loads(ops_server.tool_list_calc_options("readme.md"))
    assert data["file_type"] == "other"
    assert "RMSD" in data["options"]


def test_tool_run_calc_option_routes_smart_recalc(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    data = json.loads(ops_server.tool_run_calc_option(
        str(folder), "Smart Recalc",
    ))
    assert data["option"] == "Smart Recalc"
    assert data["mode"] == "delfin-recalc"
    assert data["valid"] is True
    assert data["submitted_job_id"] is None


def test_tool_run_calc_option_routes_classic_recalc(tmp_path):
    folder = tmp_path / "calc"
    _write_recalc_inputs(folder)
    data = json.loads(ops_server.tool_run_calc_option(
        str(folder), "Recalc",
    ))
    assert data["mode"] == "delfin-recalc-classic"


def test_tool_run_calc_option_ui_only_returns_hint(tmp_path):
    folder = tmp_path / "calc"
    folder.mkdir()
    data = json.loads(ops_server.tool_run_calc_option(
        str(folder), "MO Plot",
    ))
    assert data["ok"] is False
    assert data["ui_only"] is True
    assert "ACTION:" in data["hint"]


def test_tool_run_calc_option_unknown_option(tmp_path):
    folder = tmp_path / "calc"
    folder.mkdir()
    data = json.loads(ops_server.tool_run_calc_option(
        str(folder), "Cheese",
    ))
    assert data["ok"] is False
    assert "Unknown option" in data["error"]
    assert "list_calc_options" in data["error"]


def test_tool_list_ssh_transfer_jobs_returns_list(monkeypatch):
    fake_jobs = [
        {"job_id": "ssh-001", "status": "running", "updated_at": "2026-04-28T12:00:00Z"},
        {"job_id": "ssh-002", "status": "done", "updated_at": "2026-04-28T11:30:00Z"},
    ]
    monkeypatch.setattr(
        "delfin.ssh_transfer_jobs.list_transfer_jobs",
        lambda *, limit=8: list(fake_jobs[: limit]),
    )
    txt = ops_server.tool_list_ssh_transfer_jobs(limit=10)
    rows = json.loads(txt)
    assert len(rows) == 2
    assert rows[0]["job_id"] == "ssh-001"


def test_tool_list_ssh_transfer_jobs_handles_failure(monkeypatch):
    """Backend exception → single error row, not a crash."""
    def boom(*, limit=8):
        raise RuntimeError("ssh control socket missing")
    monkeypatch.setattr(
        "delfin.ssh_transfer_jobs.list_transfer_jobs", boom,
    )
    rows = json.loads(ops_server.tool_list_ssh_transfer_jobs())
    assert len(rows) == 1
    assert "ssh control socket missing" in rows[0]["error"]
