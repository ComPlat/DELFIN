"""A turn does not end on a call that never ran, a cut answer, or a step it
only announced.

Diagnosis 2026-10-15: 34 of 646 DeepSeek benchmark attempts had a write or
edit call whose arguments were invalid JSON, and 9 turns ENDED there with
"Let me write the file properly"; 4 science_pipeline answers ended
mid-sentence; 2 long tasks ended on "...:" with the step undone. The cause
of most of it was the benchmark's own output limit -- `bench run` sent
max_tokens=1024 per response where a session sends the role's 32768 --
and the turn loop had no answer for the rest:

* a write call cut at the output limit is told to go in parts
  (write_file(..., mode="append")) instead of "re-emit the same call";
* the loop adds one "[Continue]" row and lets the model answer again when
  the last tool call did not run, when the answer was cut at the limit, or
  when it closes by announcing a step after a code edit -- once per turn,
  never in plan mode.
"""

from __future__ import annotations

import json
from unittest.mock import patch

from delfin.agent import api_client as A
from delfin.agent import channels


def _loop(tmp_path, script, *, mode="acceptEdits"):
    """``script`` items: ("tool", name, raw_args, finish) or
    ("text", content, finish). Returns every request's messages."""
    from delfin.agent import mcp_client as _mcp
    from delfin.agent import model_capabilities as _mc

    class _Delta:
        def __init__(self, content=None, tool_calls=None):
            self.content, self.tool_calls = content, tool_calls
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
    with patch.dict("os.environ", {"DELFIN_VERIFY_NUDGE": "off"}), \
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
                        i, item[1], item[2])]), None)]),
                    _Chunk([_Choice(_Delta(), item[3])]),
                ])
            return _Stream([_Chunk([_Choice(_Delta(content=item[1]),
                                            item[2])])])

        client.client.chat.completions.create = _create
        list(client.stream_message("sys", [{"role": "user", "content": "los"}],
                                   max_tokens=64))
    return sent


def _rows(sent, text):
    return [m for m in (sent[-1] if sent else [])
            if text in str(m.get("content") or "")]


# vLLM ends a cut call with finish_reason "tool_calls"; the call is
# dispatched and its arguments do not parse.
_CUT_WRITE = ("tool", "write_file",
              '{"path": "export.py", "content": "import csv\\n\\ndef ex',
              "tool_calls")


def test_a_cut_write_is_told_to_go_in_parts(tmp_path):
    sent = _loop(tmp_path, [_CUT_WRITE, ("text", "Done.", "stop")])
    tool_rows = [m for m in sent[1] if m.get("role") == "tool"]
    assert tool_rows and 'mode=\\"append\\"' in tool_rows[-1]["content"]
    assert "did NOT run" in tool_rows[-1]["content"]


def test_a_response_cut_inside_a_call_is_told_the_call_was_dropped(tmp_path):
    cut = ("tool", "write_file",
           '{"path": "export.py", "content": "import csv', "length")
    sent = _loop(tmp_path, [cut, ("text", "export.py is not written.", "stop")])
    rows = _rows(sent, "[Continue] Your last tool call was cut")
    assert len(rows) == 1 and 'mode="append"' in rows[0]["content"]


def test_a_turn_does_not_end_on_a_call_that_did_not_run(tmp_path):
    channels.reset()
    sent = _loop(tmp_path, [_CUT_WRITE,
                            ("text", "Let me write the file properly.", "stop"),
                            ("text", "export.py is not written.", "stop")])
    assert len(_rows(sent, "[Continue] Your last tool call did not run")) == 1
    assert channels.snapshot().get("turn_end_continue:unexecuted") == 1
    channels.reset()


def test_a_cut_answer_is_continued_once(tmp_path):
    sent = _loop(tmp_path, [("text", "Ich baue jetzt den Schritt. Die Aufgabe ist", "length"),
                            ("text", "explizit, Code zu schreiben.", "length"),
                            ("text", "Ende.", "stop")])
    assert len(_rows(sent, "[Continue] Your answer was cut")) == 1
    assert len(sent) == 2


_EDIT = ("tool", "write_file",
         json.dumps({"path": "lib.rs", "content": "pub fn f() {}\n"}), "tool_calls")


def test_a_first_person_announcement_is_left_to_the_existing_nudge(tmp_path):
    # First-person announcements belong to the announced-action
    # continuation further down the loop; this block adds no second row.
    sent = _loop(tmp_path, [_EDIT, ("text", "Let me write the tests now.", "stop"),
                            ("text", "Let me write the tests now.", "stop")])
    assert not _rows(sent, "[Continue] Your answer announces a step")
    assert len(_rows(sent, "You ended your turn by announcing an action")) == 1


def test_a_closing_colon_after_an_edit_is_an_announcement(tmp_path):
    sent = _loop(tmp_path, [_EDIT, ("text", "**tests/test_store.py** — Persistenz:", "stop"),
                            ("text", "Erledigt.", "stop")])
    assert len(_rows(sent, "[Continue] Your answer announces a step")) == 1


def test_a_finished_answer_is_not_continued(tmp_path):
    sent = _loop(tmp_path, [_EDIT, ("text", "lib.rs ist geschrieben.", "stop")])
    assert not _rows(sent, "[Continue]") and len(sent) == 2


def test_a_closing_colon_without_an_edit_is_not_continued(tmp_path):
    sent = _loop(tmp_path, [("text", "Die Ursachen:", "stop")])
    assert not _rows(sent, "[Continue]") and len(sent) == 1


def test_plan_mode_is_not_continued(tmp_path):
    sent = _loop(tmp_path, [("text", "Der Plan ist", "length")], mode="plan")
    assert not _rows(sent, "[Continue]")


def test_the_benchmark_sends_the_role_limit_and_records_it(tmp_path):
    import inspect

    from delfin.agent import benchmark as bm
    from delfin.agent import benchmark_runner as br
    from delfin.agent import cli
    assert inspect.signature(br.run_task).parameters["max_tokens"].default == 0
    parser = cli._build_parser() if hasattr(cli, "_build_parser") else None
    if parser is not None:
        args = parser.parse_args(["bench", "run", "--model", "m"])
        assert args.max_tokens == 0
    path = bm.write_run([], model="m", runs_dir=tmp_path,
                        run_args={"max_tokens": 0, "repeats": 1})
    assert bm.read_run_header(path)["settings"]["run_args"]["max_tokens"] == 0
