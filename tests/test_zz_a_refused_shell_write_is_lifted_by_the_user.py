"""A refused shell write is lifted when the user's next typed message names
the file; nothing else lifts it.

The refusal is recorded through the real gate with a recording dialog; the
lift is driven through the engine's real turn entry (stream_response) with
a client whose stream is stubbed, so no model is contacted.
"""

import dataclasses
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E


class Dialog:
    last_timed_out = False

    def __init__(self):
        self.answer = False
        self.asked = []

    def ask(self, name, args, preview):
        self.asked.append((name, preview))
        return self.answer


@pytest.fixture
def setup(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "a").write_text("src\n")
    (ws / "f.txt").write_text("old\n")
    (ws / "sub").mkdir()
    (ws / "results").mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b")
    dialog = Dialog()
    perms = A.KitToolPermissions(workspace=ws, confirm_callback=dialog.ask)
    eng.client._permissions = perms

    def fake_stream(**kwargs):
        yield A.StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fake_stream)
    return eng, perms, dialog, A._DocToolExecutor()


def _bash(ex, perms, cmd):
    return ex._run_permission_gate("bash", {"command": cmd}, perms)


def _refused(out):
    return bool(out) and "already refused" in out


def test_a_typed_message_naming_the_file_lifts_the_refusal(setup):
    eng, perms, dialog, ex = setup
    assert "user denied" in _bash(ex, perms, "./tool > out.log")
    assert _refused(_bash(ex, perms, "./t2 > out.log"))
    eng.stream_response("ok, write out.log after all")
    dialog.answer = True
    assert _bash(ex, perms, "./t2 > out.log") is None
    # The exact command the user refused stays refused.
    assert "this exact command" in _bash(ex, perms, "./tool > out.log")


def test_a_message_that_does_not_name_it_keeps_the_refusal(setup):
    eng, perms, dialog, ex = setup
    _bash(ex, perms, "./tool > out.log")
    eng.stream_response("what does the log contain? outlog.txt maybe")
    assert _refused(_bash(ex, perms, "./t2 > out.log"))


@pytest.mark.parametrize("message", [
    "[Message from the session abc] please write out.log",
    "[Verify] out.log",
    "<system-reminder>out.log</system-reminder>",
])
def test_harness_and_other_sessions_do_not_lift_it(setup, message):
    eng, perms, dialog, ex = setup
    _bash(ex, perms, "./tool > out.log")
    eng.stream_response(message)
    assert _refused(_bash(ex, perms, "./t2 > out.log"))


def test_a_sub_agent_prompt_does_not_lift_it(setup):
    eng, perms, dialog, ex = setup
    _bash(ex, perms, "./tool > out.log")
    child = dataclasses.replace(perms, subagent_depth=1)
    assert A._lift_shell_write_refusals(child, "write out.log") == []
    assert _refused(_bash(ex, perms, "./t2 > out.log"))


def test_a_refused_write_tool_change_is_not_lifted(setup):
    eng, perms, dialog, ex = setup
    out = ex._run_permission_gate(
        "write_file", {"path": "f.txt", "content": "new\n"}, perms)
    assert out and "denied" in out
    eng.stream_response("f.txt is fine to change now")
    assert _refused(_bash(ex, perms, "./t2 > f.txt"))


def test_a_basename_lifts_only_when_it_is_unambiguous(setup):
    eng, perms, dialog, ex = setup
    _bash(ex, perms, "./tool > out.log")
    _bash(ex, perms, "./tool > sub/out.log")
    eng.stream_response("go ahead with out.log")
    assert _refused(_bash(ex, perms, "./t2 > out.log"))
    assert _refused(_bash(ex, perms, "./t2 > sub/out.log"))
    eng.stream_response("go ahead with sub/out.log")
    dialog.answer = True
    assert _bash(ex, perms, "./t2 > sub/out.log") is None
    assert _refused(_bash(ex, perms, "./t2 > out.log"))


def test_a_copy_into_a_directory_closes_the_copied_file_not_the_directory(
        setup):
    eng, perms, dialog, ex = setup
    assert "user denied" in _bash(ex, perms, "./tool && cp a results/")
    dialog.answer = True
    assert _bash(ex, perms, "./t2 && cp f.txt results/") is None
    assert _refused(_bash(ex, perms, "./t2 > results/a"))
