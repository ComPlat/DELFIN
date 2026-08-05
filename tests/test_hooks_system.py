"""Tests for delfin.agent.hooks (settings.json hook system)."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from delfin.agent import hooks as H


def _write_settings(workspace: Path, hooks_obj: dict) -> None:
    settings_dir = workspace / ".delfin"
    settings_dir.mkdir(exist_ok=True)
    (settings_dir / "settings.json").write_text(
        json.dumps({"hooks": hooks_obj}), encoding="utf-8",
    )


def test_load_hooks_empty_when_no_settings():
    with tempfile.TemporaryDirectory() as d:
        cfg = H.load_hooks(d)
        assert cfg.is_empty()


def test_load_hooks_reads_project_settings():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": "edit_file",
                 "hooks": [{"type": "command", "command": "true"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        assert not cfg.is_empty()
        pre = cfg.for_event("PreToolUse")
        assert len(pre) == 1
        assert pre[0].matcher == "edit_file"
        assert pre[0].command == "true"


def test_run_hooks_fires_matching_hook():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        marker = ws / "hook_fired.txt"
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": "edit_file",
                 "hooks": [{"type": "command",
                            "command": f"touch {marker}"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        results = H.run_hooks(
            "PreToolUse", cfg,
            tool_name="edit_file", arguments={"path": "x"},
            workspace=ws,
        )
        assert len(results) == 1
        assert results[0].exit_code == 0
        assert marker.exists()


def test_run_hooks_skips_non_matching():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": "edit_file",
                 "hooks": [{"type": "command", "command": "true"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        results = H.run_hooks(
            "PreToolUse", cfg,
            tool_name="bash", arguments={},
            workspace=ws,
        )
        assert results == []


def test_blocking_hook_via_exit_code():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": ".*",
                 "hooks": [{"type": "command",
                            "command": "echo nope >&2; exit 1"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        results = H.run_hooks(
            "PreToolUse", cfg,
            tool_name="bash", arguments={},
            workspace=ws,
        )
        blk = H.first_block(results)
        assert blk is not None
        assert blk.exit_code == 1
        assert "nope" in blk.stderr


def test_blocking_hook_via_json_decision():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        decision = {"decision": "block", "reason": "tests are red"}
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": ".*",
                 "hooks": [{"type": "command",
                            "command": f"echo '{json.dumps(decision)}'"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        results = H.run_hooks(
            "PreToolUse", cfg,
            tool_name="bash", arguments={},
            workspace=ws,
        )
        blk = H.first_block(results)
        assert blk is not None
        assert blk.decision == "block"
        assert "tests are red" in blk.reason


def test_template_expansion_in_command():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        out = ws / "out.txt"
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": "edit_file",
                 "hooks": [{"type": "command",
                            "command": f"echo ${{path}} > {out}"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        H.run_hooks(
            "PreToolUse", cfg,
            tool_name="edit_file",
            arguments={"path": "/foo/bar.py"},
            workspace=ws,
        )
        assert out.read_text().strip() == "/foo/bar.py"


def test_user_prompt_submit_hook():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        marker = ws / "submit.txt"
        _write_settings(ws, {
            "UserPromptSubmit": [
                {"matcher": "",
                 "hooks": [{"type": "command",
                            "command": f"echo $DELFIN_USER_PROMPT > {marker}"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        H.run_hooks(
            "UserPromptSubmit", cfg,
            user_prompt="hello world",
            workspace=ws,
        )
        assert marker.read_text().strip() == "hello world"


def test_post_tool_use_runs_after_dispatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        marker = ws / "post.txt"
        _write_settings(ws, {
            "PostToolUse": [
                {"matcher": ".*",
                 "hooks": [{"type": "command",
                            "command": f"echo POST > {marker}"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        H.run_hooks(
            "PostToolUse", cfg,
            tool_name="bash", arguments={},
            workspace=ws,
        )
        assert marker.read_text().strip() == "POST"


def test_local_settings_extends_project_settings():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / ".delfin").mkdir()
        (ws / ".delfin" / "settings.json").write_text(json.dumps({
            "hooks": {"PreToolUse": [
                {"matcher": "a", "hooks": [{"type": "command", "command": "true"}]}
            ]}
        }), encoding="utf-8")
        (ws / ".delfin" / "settings.local.json").write_text(json.dumps({
            "hooks": {"PreToolUse": [
                {"matcher": "b", "hooks": [{"type": "command", "command": "true"}]}
            ]}
        }), encoding="utf-8")
        cfg = H.load_hooks(ws)
        matchers = sorted(h.matcher for h in cfg.for_event("PreToolUse"))
        assert matchers == ["a", "b"]


def test_invalid_event_ignored():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_settings(ws, {
            "GarbageEvent": [
                {"matcher": ".*",
                 "hooks": [{"type": "command", "command": "true"}]}
            ]
        })
        cfg = H.load_hooks(ws)
        assert cfg.is_empty()


def test_hook_timeout_returns_124():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": ".*",
                 "hooks": [{"type": "command",
                            "command": "sleep 5",
                            "timeout": 0.2}]}
            ]
        })
        cfg = H.load_hooks(ws)
        results = H.run_hooks(
            "PreToolUse", cfg,
            tool_name="bash", arguments={}, workspace=ws,
        )
        assert results and results[0].exit_code == 124


def test_env_vars_set():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        out = ws / "env.txt"
        _write_settings(ws, {
            "PreToolUse": [
                {"matcher": ".*",
                 "hooks": [{"type": "command",
                            "command": (
                                f"echo $DELFIN_HOOK_EVENT $DELFIN_TOOL_NAME "
                                f"$DELFIN_WORKSPACE > {out}"
                            )}]}
            ]
        })
        cfg = H.load_hooks(ws)
        H.run_hooks(
            "PreToolUse", cfg,
            tool_name="edit_file", arguments={}, workspace=ws,
        )
        line = out.read_text().strip().split()
        assert line[0] == "PreToolUse"
        assert line[1] == "edit_file"
        assert str(ws.resolve()) == str(Path(line[2]).resolve())


# ---------------------------------------------------------------------------
# A PostToolUse hook can tell the agent what it found
# ---------------------------------------------------------------------------

def test_post_tool_hook_output_reaches_the_model():
    """It went to stderr only, so the most useful hook shape there is --
    run the linter after every edit and tell the agent what it said -- was
    impossible. A user could get it only as a BLOCKING PreToolUse, which
    refuses the edit instead of reporting on it."""
    from delfin.agent.api_client import _post_hook_note

    class _Outcome:
        def __init__(self, stdout):
            self.stdout = stdout

    note = _post_hook_note([_Outcome("ruff: 2 issues in mod.py")])
    assert "ruff: 2 issues" in note
    assert "[hook]" in note, "the model must be able to tell hook from tool"


def test_a_silent_hook_adds_nothing():
    """The common case: no output, no note, no cost."""
    from delfin.agent.api_client import _post_hook_note

    class _Outcome:
        stdout = ""

    assert _post_hook_note([]) == ""
    assert _post_hook_note([_Outcome()]) == ""


def test_a_chatty_hook_cannot_bury_the_tool_result():
    """The original reason for discarding this output was real and still
    holds -- it just did not justify throwing the output away."""
    from delfin.agent.api_client import _post_hook_note

    class _Outcome:
        stdout = "x" * 50_000

    note = _post_hook_note([_Outcome()])
    assert len(note) < 2_500
    assert "omitted" in note, "silent truncation reads as complete output"


def test_a_broken_outcome_does_not_break_the_tool_result():
    from delfin.agent.api_client import _post_hook_note

    assert _post_hook_note(None) == ""
    assert _post_hook_note([object()]) == ""


def test_only_events_that_fire_are_configurable():
    """"Notification" validated in settings and was raised nowhere: a user
    could configure a hook for it and nothing would ever run. A declared
    event that never arrives is worse than an absent one, because it looks
    configured."""
    from delfin.agent.hooks import _VALID_EVENTS

    assert "Notification" not in _VALID_EVENTS
    for event in _VALID_EVENTS:
        assert event in ("PreToolUse", "PostToolUse", "UserPromptSubmit",
                         "Stop"), event
