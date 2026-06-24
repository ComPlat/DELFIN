"""Auto-verification loop: the harness checks code the agent just edited before
letting the turn finish, and forces a bounded fix round on a problem."""

from __future__ import annotations

from pathlib import Path

from delfin.agent.api_client import _run_auto_verify, _resolve_auto_verify


def test_syntax_check_flags_broken_python(tmp_path):
    bad = tmp_path / "bad.py"
    bad.write_text("def f(:\n  return 1\n")             # syntax error
    out = _run_auto_verify([str(bad)], "syntax", "", tmp_path)
    assert out and ("SyntaxError" in out or "invalid syntax" in out)


def test_syntax_check_passes_valid_python(tmp_path):
    ok = tmp_path / "ok.py"
    ok.write_text("def f():\n    return 1\n")
    assert _run_auto_verify([str(ok)], "syntax", "", tmp_path) == ""


def test_off_and_missing_file_are_safe(tmp_path):
    assert _run_auto_verify([str(tmp_path / "nope.py")], "syntax", "", tmp_path) == ""
    assert _run_auto_verify([], "off", "", tmp_path) == ""


def test_command_mode_reports_failure(tmp_path):
    # A command that exits non-zero must surface as a problem.
    out = _run_auto_verify([], "command", "python -c \"import sys; sys.exit(3)\"",
                           tmp_path)
    assert "failed (exit 3)" in out


def test_command_mode_clean_on_success(tmp_path):
    assert _run_auto_verify([], "command", "python -c \"pass\"", tmp_path) == ""


def test_default_mode_is_syntax():
    from delfin.user_settings import DEFAULT_SETTINGS
    assert DEFAULT_SETTINGS["agent"]["auto_verify"] == "syntax"
    mode, _cmd = _resolve_auto_verify()
    assert mode in ("syntax", "command", "off")


def test_loop_is_wired():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert "_run_auto_verify(" in src
    assert "force a fix round instead of ending" in src
    assert 'fn_name in ("edit_file", "multi_edit", "write_file")' in src
