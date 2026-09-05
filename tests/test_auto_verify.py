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


# The wiring used to be asserted as five substrings of api_client.py,
# including a NEGATIVE one -- ``"_verified_this_turn" not in src``. A
# negative substring is satisfied forever by choosing a different name for
# the same flag, which is the one change the assertion exists to catch.
# The others pin an exact expression, so a legitimate refactor breaks them
# while the behaviour is untouched.
#
# What is promised: a model that only ACKNOWLEDGES, without editing
# anything, is verified again rather than let off. Two acknowledgement
# rounds are the shape that proves it.

def _loop_over_a_broken_edit(tmp_path, replies_after_the_edit):
    """Drive the real client loop: one edit that does not compile, then
    text-only answers. Returns every request the loop sent."""
    import json
    from unittest.mock import patch

    from delfin.agent import api_client as A
    from delfin.agent import mcp_client as _mcp
    from delfin.agent import model_capabilities as _mc

    class _Delta:
        def __init__(self, content=None, tool_calls=None):
            self.content = content
            self.tool_calls = tool_calls
            self.reasoning_content = None

    class _Choice:
        def __init__(self, delta, finish):
            self.index, self.delta, self.finish_reason = 0, delta, finish

    class _Chunk:
        def __init__(self, choices):
            self.choices, self.usage = choices, None

    class _Stream:
        def __init__(self, chunks):
            self._chunks = chunks

        def __iter__(self):
            return iter(self._chunks)

        def close(self):
            pass

    class _ToolCall:
        def __init__(self, name, args):
            self.index, self.id, self.type = 0, "call_1", "function"
            self.function = type("F", (), {"name": name, "arguments": args})()

    class _NoRegistry:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    caps = _mc.ModelCapabilities(model="m", provider="ollama",
                                 context_window=200_000, supports_tools=True)
    sent: list = []
    with patch.object(_mc, "resolve", lambda *a, **k: caps), \
         patch.object(_mcp, "get_registry", lambda *a, **k: _NoRegistry()), \
         patch.object(A, "_resolve_auto_verify", lambda: ("syntax", "")), \
         patch.object(A.time, "sleep", lambda *a, **k: None):
        client = A.create_client(backend="api", provider="ollama",
                                 model="qwen2.5-coder:32b", cwd=str(tmp_path))
        perms = A.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
        perms.confirm_callback = lambda n, a, prev: True
        client.set_permissions(perms)
        rounds = {"n": 0}

        def _create(**kwargs):
            rounds["n"] += 1
            sent.append(list(kwargs.get("messages") or []))
            if rounds["n"] == 1:
                return _Stream([
                    _Chunk([_Choice(_Delta(tool_calls=[_ToolCall(
                        "write_file",
                        json.dumps({"path": "kaputt.py",
                                    "content": "def f(:\n    return 1\n"}))]),
                        None)]),
                    _Chunk([_Choice(_Delta(), "tool_calls")]),
                ])
            i = min(rounds["n"] - 2, len(replies_after_the_edit) - 1)
            return _Stream([_Chunk([_Choice(
                _Delta(content=replies_after_the_edit[i]), "stop")])])

        client.client.chat.completions.create = _create
        list(client.stream_message("sys", [{"role": "user", "content": "los"}],
                                   max_tokens=64))
    return sent


_VERIFY_FEEDBACK = "Auto-verification of the file(s) you just edited"


def _verification_rounds(sent) -> int:
    """How many times the loop injected the verification feedback."""
    last = sent[-1] if sent else []
    return sum(1 for m in last
               if m.get("role") == "user"
               and _VERIFY_FEEDBACK in str(m.get("content") or ""))


def test_a_broken_edit_forces_a_fix_round(tmp_path):
    """The turn does not end on the model's word that it is done."""
    sent = _loop_over_a_broken_edit(
        tmp_path, ["Fertig, sieht gut aus.", "Ja, erledigt.", "Erledigt."])
    assert _verification_rounds(sent) >= 1, (
        "a file that does not compile ended the turn unremarked")


def test_a_model_that_only_acknowledges_is_verified_again(tmp_path):
    """The whole point of the gate: no new edit is not the same fact as
    no remaining problem. Two acknowledgement rounds, two verifications."""
    sent = _loop_over_a_broken_edit(
        tmp_path, ["Fertig, sieht gut aus.", "Ja, ist erledigt.", "Erledigt."])
    assert _verification_rounds(sent) == 2, (
        "the second acknowledgement was let off: a per-turn 'already "
        "verified' flag is back")


def test_the_fix_rounds_are_bounded(tmp_path):
    """A genuinely unfixable failure must not loop forever."""
    sent = _loop_over_a_broken_edit(
        tmp_path, ["nope"] * 8)
    assert _verification_rounds(sent) == 2
    assert len(sent) <= 5, f"the loop ran {len(sent)} rounds"


def test_a_repaired_file_ends_the_turn(tmp_path):
    """A gate that never lets go would be switched off within a week."""
    import json
    from unittest.mock import patch

    from delfin.agent import api_client as A

    real = A._run_auto_verify
    calls = {"n": 0}

    def _verify(paths, mode, cmd, ws, **kw):
        calls["n"] += 1
        # The model "fixed" it between the first and the second check.
        if calls["n"] == 1:
            return real(paths, mode, cmd, ws, **kw)
        return ""

    with patch.object(A, "_run_auto_verify", _verify):
        sent = _loop_over_a_broken_edit(
            tmp_path, ["Fertig.", "Jetzt wirklich fertig."])
    assert calls["n"] == 2, "the verification did not run twice"
    assert _verification_rounds(sent) == 1, (
        "a passing check still forced another fix round")


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


# --- Honest auto-verify: timeout → scoped fallback, no silent clean ---------


def _ws_with_matching_test(tmp_path):
    (tmp_path / "mod.py").write_text("def val():\n    return 1\n")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_mod.py").write_text(
        "from mod import val\ndef test_val():\n    assert val() == 1\n")
    return tmp_path


def test_timeout_falls_back_to_scoped_run(tmp_path, monkeypatch):
    """A full-suite timeout used to blacklist the workspace and return ''
    (clean) forever after — verification silently OFF. Now it retries scoped
    to the edited modules' tests and remembers that command."""
    from delfin.agent import api_client as A
    ws = _ws_with_matching_test(tmp_path)
    calls = []

    def fake_run(cmd, workspace, timeout):
        calls.append(cmd)
        if "test_mod.py" in cmd:
            return f"`{cmd}` failed (exit 1):\nboom", False
        return "", True                       # full suite: timeout

    monkeypatch.setattr(A, "_run_test_command", fake_run)
    status: dict = {}
    prob = A._run_auto_verify([str(ws / "mod.py")], "smart", "", str(ws),
                              status=status)
    try:
        assert "failed" in prob               # scoped result, NOT silent clean
        assert status.get("scoped") is True
        assert "test_mod.py" in status.get("command", "")
        assert str(ws) in A._SLOW_TEST_WS     # full suite stays blacklisted
        assert "test_mod.py" in A._SLOW_WS_SCOPED_CMD[str(ws)]
        assert len(calls) == 2                # full (timed out) + scoped
    finally:
        A._SLOW_TEST_WS.discard(str(ws))
        A._SLOW_WS_SCOPED_CMD.pop(str(ws), None)


def test_blacklisted_workspace_runs_scoped_not_silent_clean(tmp_path, monkeypatch):
    from delfin.agent import api_client as A
    ws = _ws_with_matching_test(tmp_path)
    A._SLOW_TEST_WS.add(str(ws))
    monkeypatch.setattr(
        A, "_run_test_command",
        lambda cmd, workspace, timeout: (f"`{cmd}` failed (exit 1):\nred", False))
    try:
        status: dict = {}
        prob = A._run_auto_verify([str(ws / "mod.py")], "smart", "", str(ws),
                                  status=status)
        assert prob != ""                     # the old behaviour returned ""
        assert status.get("scoped") is True
    finally:
        A._SLOW_TEST_WS.discard(str(ws))
        A._SLOW_WS_SCOPED_CMD.pop(str(ws), None)


def test_blacklisted_workspace_without_scoped_match_reports_skip(tmp_path):
    from delfin.agent import api_client as A
    (tmp_path / "lonely.py").write_text("x = 1\n")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_other.py").write_text(
        "def test_other():\n    assert True\n")
    A._SLOW_TEST_WS.add(str(tmp_path))
    try:
        status: dict = {}
        prob = A._run_auto_verify([str(tmp_path / "lonely.py")], "smart", "",
                                  str(tmp_path), status=status)
        assert prob == ""
        assert status.get("skipped") is True  # visible, not silent
    finally:
        A._SLOW_TEST_WS.discard(str(tmp_path))


def test_no_test_suite_reports_skip_status(tmp_path):
    from delfin.agent import api_client as A
    ok = tmp_path / "ok.py"
    ok.write_text("def f():\n    return 1\n")
    status: dict = {}
    assert A._run_auto_verify([str(ok)], "smart", "", str(tmp_path),
                              status=status) == ""
    assert status.get("skipped") is True
    assert "no test suite" in status.get("reason", "")


def test_scoped_cmd_resolves_conventional_test_files(tmp_path):
    from delfin.agent.api_client import _scoped_test_cmd_for_edits
    ws = _ws_with_matching_test(tmp_path)
    cmd = _scoped_test_cmd_for_edits([str(ws / "mod.py")], str(ws))
    assert "pytest" in cmd and "test_mod.py" in cmd
    # No matching test anywhere → no scoped command.
    (ws / "orphan.py").write_text("y = 2\n")
    assert _scoped_test_cmd_for_edits([str(ws / "orphan.py")], str(ws)) \
        .count("orphan") == 0
