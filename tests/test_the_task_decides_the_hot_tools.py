"""The hot tool set follows the task, and a missing tool becomes a script.

Measured on a field install (solo mode, kit.glm-5.3-flash, effort low): a
German request to split an .xlsx into one sheet per value of a column. The
generalist listed the workspace, fetched read_document and edit_sheet
through describe_native_tool, and ended three turns with no text. 243 of
330 failed turns on that install were empty turns. Two causes:

* deferral (2026-11-23) sent the generalist 34 of 98 schemas, a hot set
  derived from traces that are 96% coding, so read_document, edit_sheet
  and sum_column were deferred on exactly the request that needs them.
  The hot set now also carries the tools of every prompt module the
  request triggers -- from the prompt loader's own trigger table;
* an empty end right after a tool lookup ended the turn. It now gets one
  continuation row naming the move that is always there: a short script
  run through bash, or one sentence saying why the work cannot be done.

Both are proven through ``stream_message``, the path a session takes.
"""

from __future__ import annotations

import json
from unittest.mock import patch

from delfin.agent import api_client as A

SPREADSHEET_REQUEST = (
    "lege für jeden collaborator in Spalte E ein Tabellenblatt an, in dem "
    "alle zugehörigen Einträge zu dem Collaborator erfasst sind")
CODING_REQUEST = "fix the failing test in tests/test_parser.py"


def _drive(tmp_path, script, user_text):
    """Run one turn; ``script`` items are ("tool", name, args) or
    ("text", content). Returns every request's kwargs."""
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
        perms = A.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
        perms.confirm_callback = lambda n, a, prev: True
        client.set_permissions(perms)
        rounds = {"n": 0}

        def _create(**kwargs):
            if not kwargs.get("stream"):
                # The once-per-round non-streamed retry of an empty round:
                # it answers empty as well.
                raise RuntimeError("empty again")
            i = rounds["n"]
            rounds["n"] += 1
            sent.append(dict(kwargs, messages=list(kwargs.get("messages")
                                                   or [])))
            item = script[min(i, len(script) - 1)]
            if item[0] == "tool":
                return _Stream([
                    _Chunk([_Choice(_Delta(tool_calls=[_ToolCall(
                        i, item[1], item[2])]), None)]),
                    _Chunk([_Choice(_Delta(), "tool_calls")]),
                ])
            return _Stream([_Chunk([_Choice(_Delta(content=item[1]),
                                            "stop")])])

        client.client.chat.completions.create = _create
        list(client.stream_message(
            "sys", [{"role": "user", "content": user_text}], max_tokens=64))
    return sent


def _schemas(kwargs) -> dict:
    return {t["function"]["name"]: t for t in kwargs.get("tools") or []}


def test_a_spreadsheet_request_carries_the_document_schemas(tmp_path):
    (tmp_path / "liste.xlsx").write_bytes(b"")
    sent = _drive(tmp_path, [("text", "Fertig.")], SPREADSHEET_REQUEST)
    tools = _schemas(sent[0])
    for name in ("read_document", "edit_sheet", "sum_column"):
        assert name in tools, f"{name} deferred on a spreadsheet request"
        assert tools[name]["function"].get("parameters"), name


def test_a_coding_request_carries_none_of_them(tmp_path):
    sent = _drive(tmp_path, [("text", "Done.")], CODING_REQUEST)
    tools = _schemas(sent[0])
    assert not (set(A._OFFICE_TOOL_NAMES) & set(tools)), sorted(
        set(A._OFFICE_TOOL_NAMES) & set(tools))


def test_the_triggers_are_the_prompt_loaders():
    """One keyword list: the modules come from _MODULE_TRIGGERS."""
    from delfin.agent.prompt_loader import PromptLoader
    assert set(A._MODULE_TOOL_NAMES) <= set(PromptLoader._MODULE_TRIGGERS)
    msgs = [{"role": "user", "content": SPREADSHEET_REQUEST}]
    assert A.task_prompt_modules(msgs) == {"documents"}
    assert A.task_prompt_modules(
        [{"role": "user", "content": CODING_REQUEST}]) == frozenset()
    # A tool result naming an .xlsx is not a request about one.
    assert A.task_prompt_modules(
        [{"role": "user", "content": CODING_REQUEST},
         {"role": "tool", "content": "data.xlsx"}]) == frozenset()


def test_a_triggered_module_costs_exactly_its_tools():
    """The budget ratchet: no module, no growth; a module adds its tools."""
    base = A.ToolSurfaceContext(role="solo_agent", chemistry_workspace=False)
    docs = A.ToolSurfaceContext(role="solo_agent", chemistry_workspace=False,
                                task_modules=frozenset({"documents"}))

    def shipped(ctx):
        return {t["function"]["name"]: t for t in A.advertisable_tools(
            A._DOC_TOOLS_OPENAI, ctx, defer=True)}

    a, b = shipped(base), shipped(docs)
    extra = set(b) - set(a)
    assert not (set(a) - set(b))
    assert extra and extra <= set(A._OFFICE_TOOL_NAMES), sorted(extra)


def test_an_empty_end_after_a_tool_lookup_is_told_to_script_it(tmp_path):
    sent = _drive(tmp_path, [
        ("tool", "describe_native_tool", json.dumps({"name": "edit_sheet"})),
        ("text", ""),
        ("text", "Done."),
    ], SPREADSHEET_REQUEST)
    rows = [m for m in sent[-1]["messages"]
            if "write a short script" in str(m.get("content") or "")]
    assert len(rows) == 1, [m.get("content") for m in sent[-1]["messages"]]
    text = rows[0]["content"]
    assert "Python" in text and "bash" in text and "one sentence" in text


def test_an_empty_end_without_a_lookup_gets_no_such_row(tmp_path):
    sent = _drive(tmp_path, [("text", "")], CODING_REQUEST)
    assert not any("write a short script" in str(m.get("content") or "")
                   for kw in sent for m in kw["messages"])
