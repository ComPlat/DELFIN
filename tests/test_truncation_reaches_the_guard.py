"""What the model could not see has to reach the guard that says so.

Two guards read a tool result: one warns when a count comes out of output
that was cut short, the other when a total is stated over a column the
reader declared undecidable. Both consumed the ``tool_result`` event, and
that event carries a 2000-character HEAD slice of the result.

* The "… truncated, N chars …" marker is written to the CONTEXT copy,
  after the event is built, and bash writes its own marker ~4000 chars
  into stdout. Either way it is not in the slice, so the engine's sniff
  for it could never fire on the large results it exists for.
* ``read_document`` appends its NOTE lines AFTER the grid. A 200-row
  sheet is several kB, so the note about an unreadable column was never
  in the first 2000 characters either.

Both now travel as their own fields. Every existing test injected the
ledgers by hand; these drive the real client loop and the real engine
event handling instead.
"""

from __future__ import annotations

import json
import textwrap
from unittest.mock import MagicMock, patch

import pytest

import delfin.agent.api_client as A
from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc
from delfin.agent.api_client import StreamEvent, _smart_truncate


# ---------------------------------------------------------------------------
# A fake OpenAI-compatible stream that asks for one real tool call
# ---------------------------------------------------------------------------

class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls
        self.reasoning_content = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.index = 0
        self.delta = delta
        self.finish_reason = finish


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _ToolCall:
    def __init__(self, name, arguments):
        self.index = 0
        self.id = "call_1"
        self.type = "function"
        self.function = type("F", (), {"name": name,
                                       "arguments": arguments})()


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


def _run_one_tool_call(tmp_path, tool, arguments):
    """Drive the real client loop through one tool call; return its events."""
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    with patch.object(mc, "resolve", lambda *a, **k: caps), \
         patch.object(M, "get_registry", lambda *a, **k: _FakeRegistry()), \
         patch.object(A.time, "sleep", lambda *a, **k: None):
        client = A.create_client(backend="api", provider="ollama",
                                 model="qwen2.5-coder:7b", cwd=str(tmp_path))
        rounds = {"n": 0}

        def _create(**kwargs):
            rounds["n"] += 1
            if rounds["n"] == 1:
                return _Stream([
                    _Chunk([_Choice(_Delta(tool_calls=[
                        _ToolCall(tool, json.dumps(arguments))]))]),
                    _Chunk([_Choice(_Delta(), finish="tool_calls")]),
                ])
            return _Stream([_Chunk([_Choice(_Delta(content="fertig"),
                                            finish="stop")])])

        client.client.chat.completions.create = _create
        return list(client.stream_message(
            "sys", [{"role": "user", "content": "los"}], max_tokens=64))


def _tool_result(events):
    return next(e for e in events if e.type == "tool_result")


@pytest.fixture
def big_sheet(tmp_path):
    """200 rows whose money column reads two ways — the field shape."""
    rows = ["Posten,Anschaffungswert"]
    rows += [f"Geraet {i},8.986" for i in range(200)]
    (tmp_path / "inventar.csv").write_text("\n".join(rows), encoding="utf-8")
    return tmp_path


# ---------------------------------------------------------------------------
# The client reports the cut and the reader's notes as data
# ---------------------------------------------------------------------------

def test_a_large_result_is_reported_as_truncated(big_sheet):
    event = _tool_result(_run_one_tool_call(
        big_sheet, "read_document", {"path": "inventar.csv"}))
    assert event.output_chars > 2000
    assert event.output_truncated is True


def test_the_ambiguity_note_travels_although_it_is_past_the_slice(big_sheet):
    from delfin.agent import verify_guard as vg

    event = _tool_result(_run_one_tool_call(
        big_sheet, "read_document", {"path": "inventar.csv"}))
    # The head slice the event used to be — this is why it never worked.
    assert vg.extract_ambiguous_columns(event.tool_output) == []
    # ... and the note, carried as its own field.
    assert vg.extract_ambiguous_columns(event.output_notes) == [
        "Anschaffungswert"]


def test_a_small_result_is_not_reported_as_truncated(tmp_path):
    (tmp_path / "klein.csv").write_text("Posten,Anzahl\nStuhl,3\n",
                                        encoding="utf-8")
    event = _tool_result(_run_one_tool_call(
        tmp_path, "read_document", {"path": "klein.csv"}))
    assert event.output_truncated is False
    assert event.output_chars == len(event.tool_output)


def test_the_marker_a_tool_writes_itself_is_past_the_slice():
    """bash caps its own stdout and writes the marker roughly 4000 chars
    in — the reason sniffing the event's 2000-char preview could not see
    it."""
    cut = _smart_truncate("x" * 40_000, 12_000, "stdout")
    assert cut.index("truncated,") > 2000
    assert "truncated," not in cut[:2000]


def test_the_emitted_result_is_still_redacted(tmp_path):
    """The event is what the tool trace is written from, and the trace gets
    bundled into bug reports. Restructuring the emission site must not lose
    that."""
    secret = "sk-ant-api03-THISLOOKSLIKEAKEY0123456789abcdefghijklmnop"
    (tmp_path / "notes.txt").write_text(f"token: {secret}\n", encoding="utf-8")
    event = _tool_result(_run_one_tool_call(
        tmp_path, "read_file", {"path": "notes.txt"}))
    assert secret not in event.tool_output
    assert "[redacted:" in event.tool_output


def test_the_event_defaults_stay_backwards_compatible():
    """Backends that build the event positionally or without the new
    fields must keep working."""
    event = StreamEvent(type="tool_result", tool_name="x", tool_output="y")
    assert event.output_truncated is False
    assert event.output_chars == 0
    assert event.output_notes == ""


# ---------------------------------------------------------------------------
# The engine arms its ledgers from those fields
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


_AMBIGUITY_NOTE = (
    "NOTE: column 'Anschaffungswert': values like '8.986' read as two "
    "different numbers depending on the convention, and nothing in the "
    "column decides it.")


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def _client(replies, tool_events=()):
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        for event in tool_events:
            yield event
        yield StreamEvent(type="text_delta",
                          text=replies[min(i, len(replies) - 1)])
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


_BIG_GRID_EVENT = StreamEvent(
    type="tool_result",
    tool_name="mcp__delfin-docs__read_document",
    tool_output="| Posten | Anschaffungswert |\n" * 60,   # no note in here
    output_truncated=True,
    output_chars=9000,
    output_notes=_AMBIGUITY_NOTE,
)


def test_the_engine_arms_the_truncation_ledger_from_the_field(agent_tree):
    engine = _engine(agent_tree, _client(["ok"], (_BIG_GRID_EVENT,)))
    engine.stream_response("lies die tabelle")
    assert engine._truncated_tools_turn == [
        "mcp__delfin-docs__read_document"]


def test_the_engine_arms_the_ambiguity_ledger_from_the_notes(agent_tree):
    engine = _engine(agent_tree, _client(["ok"], (_BIG_GRID_EVENT,)))
    engine.stream_response("lies die tabelle")
    assert engine._ambiguous_columns_turn == ["Anschaffungswert"]


def test_a_full_result_arms_nothing(agent_tree):
    plain = StreamEvent(type="tool_result", tool_name="list_files",
                        tool_output="a.pdf\nb.pdf\n", output_chars=12)
    engine = _engine(agent_tree, _client(["ok"], (plain,)))
    engine.stream_response("liste die dateien")
    assert engine._truncated_tools_turn == []
    assert engine._ambiguous_columns_turn == []


def test_a_backend_that_sends_the_whole_result_still_arms_it(agent_tree):
    """The CLI path passes tool output through unsliced and sets no field —
    the marker is genuinely in the text there."""
    whole = StreamEvent(
        type="tool_result", tool_name="Bash",
        tool_output="line\n... (stdout truncated, 4000 chars omitted)\n")
    engine = _engine(agent_tree, _client(["ok"], (whole,)))
    engine.stream_response("lauf")
    assert engine._truncated_tools_turn == ["Bash"]


def test_the_truncation_ledger_is_cleared_for_each_user_turn(agent_tree):
    """Documented as per-turn and never cleared: once anything truncated,
    every later answer with a two-digit count got the caveat."""
    engine = _engine(agent_tree, _client(["ok"], (_BIG_GRID_EVENT,)))
    engine.stream_response("lies die tabelle")
    assert engine._truncated_tools_turn
    engine.client.stream_message.side_effect = (
        lambda *a, **k: iter([StreamEvent(type="text_delta", text="ok")]))
    engine.stream_response("und jetzt etwas ganz anderes")
    assert engine._truncated_tools_turn == []
    assert engine._ambiguous_columns_turn == []


def test_the_count_caveat_reaches_a_real_answer(agent_tree):
    """End to end: a cut-short result, then an answer that counts from it."""
    answer = "Ich habe 31 Rechnungen geprüft."
    engine = _engine(agent_tree, _client([answer], (_BIG_GRID_EVENT,)))
    out = engine.stream_response("wie viele rechnungen?")
    assert "geschätzt und nicht gezählt" in out
    assert "31 Rechnungen" in out
