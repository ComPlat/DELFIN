"""At the end of a turn whose last code edit nothing checked, one line.

``agent.verify_nudge`` (off | honest | test; DELFIN_VERIFY_NUDGE wins):
when the model is about to answer and the last code edit of the turn --
any language, write tools and shell writes alike -- has no executed test,
check or run after it that touched it (turn_verification), the loop adds
one harness row at the tail and lets the model answer again:

* honest: "[Verify] No check ran after your last edit. Check it, or
  state clearly in your answer what is unverified."
* test:   "[Verify] Run a test that covers your last edit before you
  answer."

At most once per turn; never in plan mode, never without an edit, never
when the answer already names the change as unverified. The system
prompt is not touched (the prefix cache keeps its head). Driven through
the real client loop with a scripted model.
"""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import channels
from delfin.agent import turn_verification as tv

# The wording is part of the measured design (verifynudge PREREG): pinned.
HONEST = ("[Verify] No check ran after your last edit. Check it, or state "
          "clearly in your answer what is unverified.")
TEST = "[Verify] Run a test that covers your last edit before you answer."
DEFAULT = getattr(A, "_VERIFY_NUDGE_DEFAULT", "off")


def _loop(tmp_path, script, *, nudge="honest", mode="acceptEdits",
          monkeypatch=None):
    """Drive the client loop; ``script`` items are ("tool", name, args) or
    ("text", content). Returns every request's messages."""
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
        def __init__(self, n, name, args):
            self.index, self.id, self.type = 0, f"call_{n}", "function"
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
    env = {"DELFIN_VERIFY_NUDGE": nudge}
    with patch.dict("os.environ", env), \
         patch.object(_mc, "resolve", lambda *a, **k: caps), \
         patch.object(_mcp, "get_registry", lambda *a, **k: _NoRegistry()), \
         patch.object(A, "_resolve_auto_verify", lambda: ("off", "")), \
         patch.object(A.time, "sleep", lambda *a, **k: None):
        client = A.create_client(backend="api", provider="ollama",
                                 model="qwen2.5-coder:32b", cwd=str(tmp_path))
        perms = A.KitToolPermissions(workspace=tmp_path, mode=mode)
        perms.confirm_callback = lambda n, a, prev: True
        client.set_permissions(perms)
        rounds = {"n": 0}

        def _create(**kwargs):
            i = rounds["n"]
            rounds["n"] += 1
            sent.append(list(kwargs.get("messages") or []))
            item = script[min(i, len(script) - 1)]
            if item[0] == "tool":
                return _Stream([
                    _Chunk([_Choice(_Delta(tool_calls=[_ToolCall(
                        i, item[1], json.dumps(item[2]))]), None)]),
                    _Chunk([_Choice(_Delta(), "tool_calls")]),
                ])
            return _Stream([_Chunk([_Choice(_Delta(content=item[1]),
                                            "stop")])])

        client.client.chat.completions.create = _create
        list(client.stream_message("sys", [{"role": "user", "content": "los"}],
                                   max_tokens=64))
    return sent


def _lines(sent, text):
    last = sent[-1] if sent else []
    return [m for m in last if m.get("role") == "user"
            and text in str(m.get("content") or "")]


_EDIT_RS = ("tool", "write_file", {"path": "lib.rs",
                                   "content": "pub fn f() -> i32 { 1 }\n"})


def test_an_unchecked_edit_gets_the_honest_line_once(tmp_path):
    channels.reset()
    sent = _loop(tmp_path, [_EDIT_RS, ("text", "Fertig."),
                            ("text", "Fertig, wirklich.")])
    assert len(_lines(sent, HONEST)) == 1
    assert channels.snapshot().get("turn_end_verify") == 1
    assert channels.snapshot().get("turn_end_verify:honest") == 1
    channels.reset()


def test_the_test_arm_has_its_own_line(tmp_path):
    sent = _loop(tmp_path, [_EDIT_RS, ("text", "Fertig."), ("text", "Ok.")],
                 nudge="test")
    assert len(_lines(sent, TEST)) == 1
    assert not _lines(sent, HONEST)


def test_off_adds_nothing(tmp_path):
    sent = _loop(tmp_path, [_EDIT_RS, ("text", "Fertig.")], nudge="off")
    assert not _lines(sent, "[Verify]")
    assert len(sent) == 2


def test_the_line_is_at_the_tail_and_the_system_prompt_is_unchanged(tmp_path):
    sent = _loop(tmp_path, [_EDIT_RS, ("text", "Fertig."), ("text", "Ok.")])
    last = sent[-1]
    assert HONEST in str(last[-1].get("content"))
    assert last[-2].get("role") == "assistant"
    assert "Fertig." in str(last[-2].get("content"))
    assert sent[0][0] == sent[-1][0]          # the system row, byte for byte
    assert "[Verify]" not in str(sent[-1][0].get("content"))


def test_a_question_without_an_edit_gets_no_line(tmp_path):
    sent = _loop(tmp_path, [("text", "Die Antwort ist 42.")])
    assert not _lines(sent, "[Verify]") and len(sent) == 1


def test_an_answer_that_names_the_change_untested_gets_no_line(tmp_path):
    sent = _loop(tmp_path, [_EDIT_RS,
                            ("text", "lib.rs ist geändert, aber ungetestet.")])
    assert not _lines(sent, "[Verify]") and len(sent) == 2


def test_a_run_that_touches_the_edit_satisfies_it(tmp_path):
    sent = _loop(tmp_path, [
        ("tool", "write_file", {"path": "tool.py", "content": "print(6*7)\n"}),
        ("tool", "bash", {"command": "python3 tool.py"}),
        ("text", "Ausgabe 42.")])
    assert not _lines(sent, "[Verify]")


def test_plan_mode_gets_no_line(tmp_path):
    sent = _loop(tmp_path, [_EDIT_RS, ("text", "Plan: lib.rs ändern.")],
                 mode="plan")
    assert not _lines(sent, "[Verify]")


# --- the account the line reads (shared with the benchmark) -----------------

def _ok(cmd, code=0, out=""):
    return json.dumps({"exit_code": code, "stdout": out, "stderr": "",
                       "command": cmd})


def test_the_account_needs_a_run_that_touches_the_last_edit(tmp_path):
    t = tv.TurnVerification(str(tmp_path))
    t.note("edit_file", {"path": "netcfg/port.go"}, "Edited netcfg/port.go")
    assert t.unverified_edit() == "netcfg/port.go"
    t.note("bash", {"command": "go test ./...", "cwd": "textkit"},
           _ok("go test ./...", 0, "ok"))
    assert t.unverified_edit() == "netcfg/port.go"      # another project
    t.note("bash", {"command": "go test ./...", "cwd": "netcfg"},
           '{"error": "bash: \'go test ./...\' is not on the auto-allow list"}')
    assert t.unverified_edit() == "netcfg/port.go"      # refused
    t.note("bash", {"command": "go test ./...", "cwd": "netcfg"},
           _ok("go test ./...", 1, "--- FAIL: TestX"))
    assert t.unverified_edit() is None                   # a red run checks it


def test_a_shell_write_is_an_edit_and_all_skipped_is_not_a_check(tmp_path):
    t = tv.TurnVerification(str(tmp_path))
    cmd = "sed -i 's/a/b/' pkg/mod.py"
    t.note("bash", {"command": cmd}, _ok(cmd))
    assert t.unverified_edit() == "pkg/mod.py"
    t.note("run_tests", {"target": "pkg"}, json.dumps({
        "status": "ok", "framework": "pytest", "exit_code": 0,
        "summary": {"passed": 0, "failed": 0, "errors": 0, "skipped": 3}}))
    assert t.unverified_edit() == "pkg/mod.py"


def test_the_benchmark_and_the_loop_share_the_account():
    from delfin.agent import benchmark as bm
    assert bm._is_code_edit is tv.is_code_edit
    assert bm._CODE_SUFFIXES is tv.CODE_SUFFIXES


def test_the_benchmark_pins_the_shipped_default():
    from delfin.agent import benchmark_runner as br
    assert br.PINNED_BEHAVIOUR_ENV["DELFIN_VERIFY_NUDGE"] == DEFAULT
    assert A._VERIFY_NUDGE_TEXT == {"honest": HONEST, "test": TEST}


@pytest.mark.parametrize("env,setting,want", [
    ("", "", DEFAULT), ("test", "honest", "test"),
    ("", "honest", "honest"), ("bogus", "", DEFAULT)])
def test_the_setting_resolves_env_then_settings_then_default(
        env, setting, want, monkeypatch):
    from delfin import user_settings
    monkeypatch.setenv("DELFIN_VERIFY_NUDGE", env)
    monkeypatch.setattr(user_settings, "load_settings", lambda *a, **k: {
        "agent": {"verify_nudge": setting}})
    assert A._resolve_verify_nudge() == want
