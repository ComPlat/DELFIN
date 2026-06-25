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


def test_default_mode_is_smart():
    from delfin.user_settings import DEFAULT_SETTINGS
    assert DEFAULT_SETTINGS["agent"]["auto_verify"] == "smart"
    mode, _cmd = _resolve_auto_verify()
    assert mode in ("smart", "syntax", "command", "off")


def test_detect_test_command(tmp_path):
    from delfin.agent.api_client import _detect_test_command
    assert _detect_test_command(tmp_path) == ""          # no test setup
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_x.py").write_text("def test_x():\n    assert True\n")
    assert "pytest" in _detect_test_command(tmp_path)     # detected


def test_smart_runs_tests_and_flags_failure(tmp_path):
    # syntax OK but a real test fails → smart mode surfaces it.
    (tmp_path / "mod.py").write_text("def val():\n    return 1\n")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_mod.py").write_text(
        "from mod import val\ndef test_val():\n    assert val() == 2\n")
    out = _run_auto_verify([str(tmp_path / "mod.py")], "smart", "", tmp_path)
    assert out and "failed" in out


def test_smart_clean_when_tests_pass(tmp_path):
    (tmp_path / "mod.py").write_text("def val():\n    return 2\n")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_mod.py").write_text(
        "from mod import val\ndef test_val():\n    assert val() == 2\n")
    assert _run_auto_verify([str(tmp_path / "mod.py")], "smart", "", tmp_path) == ""


def test_smart_syntax_error_short_circuits(tmp_path):
    bad = tmp_path / "bad.py"
    bad.write_text("def f(:\n  return 1\n")
    out = _run_auto_verify([str(bad)], "smart", "", tmp_path)
    assert out and ("SyntaxError" in out or "invalid syntax" in out)


def test_smart_no_tests_is_syntax_only(tmp_path):
    ok = tmp_path / "ok.py"
    ok.write_text("def f():\n    return 1\n")
    assert _run_auto_verify([str(ok)], "smart", "", tmp_path) == ""


def test_loop_is_wired():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert "_run_auto_verify(" in src
    assert "force a fix round instead of ending" in src
    assert 'fn_name in ("edit_file", "multi_edit", "write_file")' in src
    # Adversarial-review fix: the gate must NOT depend on a per-turn "verified"
    # flag — that let a model which only ACKNOWLEDGES (no new edit) skip
    # re-verification and finish with an unfixed problem.
    assert "_verified_this_turn" not in src
    assert "_edited_py and _av_mode != \"off\" and _verify_attempts < 2" in src


def test_auto_verify_scopes_to_edited_package(tmp_path):
    """A broken test in a SIBLING dir must not fail auto-verify of a clean
    package the agent edited. Bug 2026-06-25: a green spreadsheet_test/ package
    was flagged because auto-verify ran the whole workspace's pytest and hit
    stale tests in a sibling directory."""
    from delfin.agent.api_client import _run_auto_verify
    pkg = tmp_path / "spreadsheet_test"
    (pkg / "sheet").mkdir(parents=True)
    (pkg / "tests").mkdir()
    (pkg / "sheet" / "__init__.py").write_text("")
    (pkg / "sheet" / "m.py").write_text("def add(a, b):\n    return a + b\n")
    (pkg / "tests" / "test_m.py").write_text(
        "from sheet.m import add\n\ndef test_add():\n    assert add(1, 2) == 3\n")
    # unrelated broken test elsewhere in the workspace
    junk = tmp_path / "old_junk"
    junk.mkdir()
    (junk / "test_broken.py").write_text("import does_not_exist_xyz_42\n")

    prob = _run_auto_verify([str(pkg / "sheet" / "m.py")], "smart", "",
                            str(tmp_path))
    assert prob == "", f"clean package falsely flagged: {prob[:200]}"


def test_auto_verify_still_catches_failure_inside_the_package(tmp_path):
    """The scoping must not hide a REAL failure in the edited package."""
    from delfin.agent.api_client import _run_auto_verify
    pkg = tmp_path / "pkg"
    (pkg / "tests").mkdir(parents=True)
    (pkg / "tests" / "test_x.py").write_text(
        "def test_x():\n    assert 1 == 2\n")
    prob = _run_auto_verify([str(pkg / "tests" / "test_x.py")], "smart", "",
                            str(tmp_path))
    assert prob != ""
