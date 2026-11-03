"""The evidence pool holds what the model saw, not what the UI previews.

Input: a tool result longer than the 2000-character UI preview but within
the model's context cap (DeepSeek: 5 KB), with a figure past the preview.
Before: the engine filled the observed-number pool from
``event.tool_output``, the 2000-character preview, while
``output_truncated`` reported only the model's own cut. The pool was
marked complete and lacked the figure, so a correct answer quoting it was
flagged as unsupported (measured: a 3229-character bash result, the figure
at about character 3000). A web_search result with eight hits is over
3.5 KB, so later snippets' numbers were lost the same way.
After: the event carries the text the model was given (``seen_text``) and
the pool is filled from it; the preview stays 2000 characters.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.api_client import StreamEvent

FIGURE = "E_gauche_minus_anti = 0.6234 kcal/mol"
LONG = "x" * 3000 + "\n" + FIGURE + "\n"


@pytest.fixture(autouse=True)
def _fresh_pool():
    vg.reset_observed_numbers()
    yield
    vg.reset_observed_numbers()


def _openai_client(ws: Path, monkeypatch):
    from delfin.agent import api_client as A
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc

    class _Delta:
        def __init__(self, content=None, tool_calls=None):
            self.content, self.tool_calls = content, tool_calls
            self.reasoning_content = None

    class _Choice:
        def __init__(self, delta, finish=None):
            self.index, self.delta, self.finish_reason = 0, delta, finish

    class _Chunk:
        def __init__(self, choices, usage=None):
            self.choices, self.usage = choices, usage

    class _Stream:
        def __init__(self, chunks):
            self._chunks = chunks

        def __iter__(self):
            return iter(self._chunks)

        def close(self):
            pass

    class _Call:
        def __init__(self, name, args):
            self.index, self.id, self.type = 0, "c0", "function"
            self.function = type("F", (), {"name": name, "arguments": args})()

    class _Registry:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    caps = mc.ModelCapabilities(model="kit.deepseek-v4-flash", provider="kit",
                                context_window=131072, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _Registry())
    client = A.create_client(backend="api", provider="kit",
                             model="kit.deepseek-v4-flash", cwd=str(ws),
                             api_key="dummy")
    client.set_permissions(A.KitToolPermissions(
        workspace=ws, mode="bypassPermissions", agent_role="solo_agent"))
    script = ws / "emit.py"
    script.write_text(f"print('x' * 3000)\nprint({FIGURE!r})\n")
    calls = {"n": 0}

    def _create(**_kw):
        calls["n"] += 1
        if calls["n"] == 1:
            args = json.dumps({"command": f"python3 {script.name}"})
            return _Stream([
                _Chunk([_Choice(_Delta(tool_calls=[_Call("bash", args)]))]),
                _Chunk([_Choice(_Delta(), finish="tool_calls")])])
        return _Stream([_Chunk([_Choice(_Delta(content="done"),
                                        finish="stop")])])

    client.client.chat.completions.create = _create
    return client


def test_the_event_carries_the_text_the_model_was_given(tmp_path, monkeypatch):
    client = _openai_client(tmp_path, monkeypatch)
    events = [ev for ev in client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=64)
        if ev.type == "tool_result"]
    assert events, "the scripted bash call produced no result"
    ev = events[0]
    assert len(ev.tool_output) <= 2000            # the preview is unchanged
    assert "0.6234" not in ev.tool_output
    assert "0.6234" in ev.seen_text and not ev.output_truncated


@pytest.fixture
def agent_tree(tmp_path):
    modes = tmp_path / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (tmp_path / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def test_a_figure_past_the_preview_grounds_the_answer(agent_tree):
    from delfin.agent.engine import AgentEngine
    result = StreamEvent(type="tool_result", tool_name="bash",
                         tool_output=LONG[:2000], output_chars=len(LONG))
    result.seen_text = LONG
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*_a, **_k):
        calls["n"] += 1
        if calls["n"] == 1:
            yield result
        yield StreamEvent(type="text_delta",
                          text="gauche lies 0.6234 kcal/mol above anti.")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                             pack_dir=agent_tree)
    out = engine.stream_response("How far apart are the conformers?")
    assert "[verify]" not in out, out
