"""Tests for the typed delfin.api wrappers.

Strategy: monkeypatch ``delfin.cli.main`` with a recording stub.  The api
layer is a thin argv-builder, so verifying that each wrapper produces the
correct argv is the right level of coverage.
"""
from __future__ import annotations

import sys

import pytest

from delfin import api as delfin_api


@pytest.fixture
def cli_recorder(monkeypatch):
    """Replace ``delfin.cli.main`` with a recorder that returns 0.

    Returns a list that captures every argv passed to cli_main.
    """
    captured: list[list[str]] = []

    def fake_main(argv):
        captured.append(list(argv))
        # Print to verify capture works
        print("stdout-marker")
        print("stderr-marker", file=sys.stderr)
        return 0

    import delfin.cli
    monkeypatch.setattr(delfin.cli, "main", fake_main)
    return captured


def test_command_result_defaults():
    rc = delfin_api.CommandResult()
    assert rc.returncode == 0
    assert rc.stdout == ""
    assert rc.stderr == ""
    assert rc.argv == []
    assert rc.ok is True


def test_command_result_ok_false_on_nonzero():
    rc = delfin_api.CommandResult(returncode=2)
    assert rc.ok is False


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def test_pipeline_run_default_argv(cli_recorder):
    rc = delfin_api.pipeline_run()
    assert rc.ok
    # Default control_file = CONTROL.txt → no --control flag emitted
    assert cli_recorder == [[]]


def test_pipeline_run_with_recalc_overwrite(cli_recorder):
    delfin_api.pipeline_run(recalc=True, overwrite=True, cleanup=False)
    argv = cli_recorder[0]
    assert "--recalc" in argv
    assert "--overwrite" in argv
    assert "--no-cleanup" in argv


def test_pipeline_run_with_custom_control(cli_recorder):
    delfin_api.pipeline_run(control_file="CUSTOM.txt")
    assert cli_recorder[0] == ["--control", "CUSTOM.txt"]


def test_pipeline_run_extra_args_passthrough(cli_recorder):
    delfin_api.pipeline_run(extra_args=["--foo", "1"])
    assert cli_recorder[0] == ["--foo", "1"]


def test_pipeline_run_does_not_duplicate_flags(cli_recorder):
    delfin_api.pipeline_run(recalc=True, extra_args=["--recalc"])
    assert cli_recorder[0].count("--recalc") == 1


def test_pipeline_prepare_argv(cli_recorder):
    delfin_api.pipeline_prepare(control_file="X.txt", overwrite=True)
    assert cli_recorder[0] == ["--define", "X.txt", "--overwrite"]


def test_run_legacy_returns_int(cli_recorder):
    """Legacy ``run()`` must keep returning a plain int."""
    result = delfin_api.run()
    assert isinstance(result, int)
    assert result == 0


def test_prepare_legacy_returns_int(cli_recorder):
    result = delfin_api.prepare(overwrite=True)
    assert isinstance(result, int)
    assert result == 0


# ---------------------------------------------------------------------------
# Read-only checks
# ---------------------------------------------------------------------------

def test_qm_check_no_tools(cli_recorder):
    delfin_api.qm_check()
    assert cli_recorder[0] == ["qm_check"]


def test_qm_check_with_tools(cli_recorder):
    delfin_api.qm_check(tools=["xtb", "crest"])
    assert cli_recorder[0] == ["qm_check", "xtb", "crest"]


def test_csp_check_argv(cli_recorder):
    delfin_api.csp_check()
    assert cli_recorder[0] == ["csp_check"]


def test_mlp_check_argv(cli_recorder):
    delfin_api.mlp_check()
    assert cli_recorder[0] == ["mlp_check"]


def test_analysis_check_argv(cli_recorder):
    delfin_api.analysis_check()
    assert cli_recorder[0] == ["analysis_check"]


# ---------------------------------------------------------------------------
# QM tool runner
# ---------------------------------------------------------------------------

def test_qm_run_minimal(cli_recorder):
    delfin_api.qm_run("xtb")
    assert cli_recorder[0] == ["qm_run", "xtb", "--cwd", ".", "--capture"]


def test_qm_run_with_args(cli_recorder):
    delfin_api.qm_run("xtb", ["coord", "--gfn", "2"], cwd="/tmp", capture_tool=False)
    assert cli_recorder[0] == [
        "qm_run", "xtb", "--cwd", "/tmp", "--", "coord", "--gfn", "2",
    ]


# ---------------------------------------------------------------------------
# cleanup / stop
# ---------------------------------------------------------------------------

def test_cleanup_default(cli_recorder):
    delfin_api.cleanup()
    assert cli_recorder[0] == ["cleanup", "--workspace", "."]


def test_cleanup_orca_dryrun(cli_recorder):
    delfin_api.cleanup(orca=True, dry_run=True, workspace="/tmp/ws", scratch="/scr")
    argv = cli_recorder[0]
    assert "--orca" in argv
    assert "--dry-run" in argv
    assert "--workspace" in argv and "/tmp/ws" in argv
    assert "--scratch" in argv and "/scr" in argv


def test_stop_default_signal_INT(cli_recorder):
    delfin_api.stop()
    argv = cli_recorder[0]
    assert argv[:2] == ["stop", "--workspace"]
    assert "--signal" in argv and "INT" in argv


def test_stop_invalid_signal_raises():
    with pytest.raises(ValueError):
        delfin_api.stop(signal_name="HUP")


def test_stop_kill_with_cleanup_after(cli_recorder):
    delfin_api.stop(signal_name="KILL", cleanup_after=True, dry_run=False)
    argv = cli_recorder[0]
    assert "KILL" in argv
    assert "--cleanup" in argv
    assert "--dry-run" not in argv


# ---------------------------------------------------------------------------
# ORCA + workflows
# ---------------------------------------------------------------------------

def test_run_orca_input_no_args(cli_recorder):
    delfin_api.run_orca_input()
    assert cli_recorder[0] == ["run_orca"]


def test_run_orca_input_with_files(cli_recorder):
    delfin_api.run_orca_input("job.inp", "job.out")
    assert cli_recorder[0] == ["run_orca", "job.inp", "--output", "job.out"]


def test_co2_define_force(cli_recorder):
    delfin_api.co2(define=True, force=True, charge=2, multiplicity=3,
                   solvent="DMF", metal="Fe")
    argv = cli_recorder[0]
    assert argv[0] == "co2"
    assert "--define" in argv
    assert "--force" in argv
    assert "--charge" in argv and "2" in argv
    assert "--multiplicity" in argv and "3" in argv
    assert "--solvent" in argv and "DMF" in argv
    assert "--metal" in argv and "Fe" in argv


def test_tadf_xtb_passthrough(cli_recorder):
    delfin_api.tadf_xtb(extra_args=["--foo"])
    assert cli_recorder[0] == ["tadf_xtb", "--foo"]


def test_hyperpol_passthrough(cli_recorder):
    delfin_api.hyperpol(extra_args=["--bar"])
    assert cli_recorder[0] == ["hyperpol", "--bar"]


# ---------------------------------------------------------------------------
# Stream capture
# ---------------------------------------------------------------------------

def test_stdout_stderr_captured(cli_recorder):
    rc = delfin_api.qm_check()
    assert "stdout-marker" in rc.stdout
    assert "stderr-marker" in rc.stderr
    assert rc.returncode == 0
    assert rc.argv == ["qm_check"]


def test_capture_false_passes_streams(cli_recorder, capsys):
    rc = delfin_api.qm_check(capture=False)
    assert rc.stdout == ""
    assert rc.stderr == ""
    captured = capsys.readouterr()
    assert "stdout-marker" in captured.out


def test_systemexit_handled(monkeypatch):
    """``cli_main`` may raise SystemExit; api should convert to returncode."""
    def fake_main(argv):
        raise SystemExit(2)

    import delfin.cli
    monkeypatch.setattr(delfin.cli, "main", fake_main)

    rc = delfin_api.qm_check()
    assert rc.returncode == 2
    assert rc.ok is False


def test_systemexit_none_handled(monkeypatch):
    """``SystemExit(None)`` should map to returncode 0."""
    def fake_main(argv):
        raise SystemExit()

    import delfin.cli
    monkeypatch.setattr(delfin.cli, "main", fake_main)

    rc = delfin_api.qm_check()
    assert rc.returncode == 0
    assert rc.ok is True


# ---------------------------------------------------------------------------
# On-demand operational-pattern lookup
# ---------------------------------------------------------------------------

def test_list_dashboard_patterns_returns_known_names():
    names = delfin_api.list_dashboard_patterns()
    assert isinstance(names, list)
    assert names == sorted(names)
    assert {"batch", "control_edit", "smart_recalc", "submit_orca",
            "analyze", "recalc", "cancel"}.issubset(set(names))


def test_get_dashboard_pattern_returns_batch_recipe():
    out = delfin_api.get_dashboard_pattern("batch")
    assert "/batch from-calc" in out
    assert "Never" in out  # carries the don't-reinvent rule


def test_get_dashboard_pattern_smart_recalc_recipe():
    out = delfin_api.get_dashboard_pattern("smart_recalc")
    assert "calc-options value Smart Recalc" in out
    assert "calc-override-btn" in out


def test_get_dashboard_pattern_case_insensitive():
    a = delfin_api.get_dashboard_pattern("Batch")
    b = delfin_api.get_dashboard_pattern("BATCH")
    c = delfin_api.get_dashboard_pattern("batch")
    assert a == b == c
    assert "/batch from-calc" in a


def test_get_dashboard_pattern_normalises_separators():
    via_hyphen = delfin_api.get_dashboard_pattern("control-edit")
    via_space = delfin_api.get_dashboard_pattern("control edit")
    via_underscore = delfin_api.get_dashboard_pattern("control_edit")
    assert via_hyphen == via_space == via_underscore
    assert "/control key" in via_hyphen


def test_get_dashboard_pattern_unknown_lists_options():
    """An unknown name returns a hint, not an exception."""
    out = delfin_api.get_dashboard_pattern("nonexistent")
    assert "Unknown pattern" in out
    assert "batch" in out and "smart_recalc" in out


def test_get_dashboard_pattern_empty_name_lists_options():
    out = delfin_api.get_dashboard_pattern("")
    assert "No pattern requested" in out
    assert "batch" in out
