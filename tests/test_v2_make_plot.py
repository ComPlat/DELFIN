"""make_plot: the tool result, the chat card, and what the model receives.

The result is ``DELFIN_CARD:`` + ``{"card": "plot", "path", "title",
"text"}`` (chat_plots.card_result). The dashboard draws the written file
with the chat's artifact preview (tab_agent._make_plot_card_html); the
model receives only ``text``, wrapped as untrusted.

The first version returned an ``<iframe srcdoc>`` card holding the base64
figure: a three-point line was a 39 kB result. The dashboard sees the first
2000 characters of a result and had no make_plot hook, so no card was ever
drawn, and the model was handed the base64 instead of the caption. The
tests below drive the real OpenAI loop and read the event the dashboard
receives and the message the model receives.
"""

from __future__ import annotations

import json

import pytest

import delfin.agent.chat_plots as cp
from delfin.agent import mcp_client as M
from delfin.agent.chat_plots import SpecError
from delfin.dashboard import tab_agent as T


def _spec() -> dict:
    return {
        "kind": "line",
        "title": "Energies per method",
        "units": "eV",
        "xlabel": "method",
        "ylabel": "energy",
        "x": ["a", "b", "c"],
        "y": [1.0, 2.0, 3.0],
    }


def _payload(out: str) -> dict:
    assert out.startswith("DELFIN_CARD:"), f"no card marker in {out[:60]!r}"
    return json.loads(out[len("DELFIN_CARD:"):])


def test_the_result_names_the_file_and_carries_no_figure(tmp_path):
    out = cp.make_plot(_spec(), out_dir=str(tmp_path))
    obj = _payload(out)
    assert set(obj) == {"card", "path", "title", "text"}
    assert obj["card"] == "plot"
    assert obj["path"] == str((tmp_path / "plot_line.svg").resolve())
    assert "<" not in out and "base64" not in out
    assert len(out) < 2000
    assert "Energies per method" in obj["text"]


def test_make_plot_bad_spec_is_refused_clearly(tmp_path):
    # Its own directory: tmp_path also holds the test HOME, where a
    # session-presence thread of an earlier test may write meanwhile.
    out = tmp_path / "plots"
    out.mkdir()
    with pytest.raises(SpecError):
        cp.make_plot({"kind": "pie"}, out_dir=str(out))
    assert list(out.rglob("*")) == [], "bad spec must not write files"


def test_make_plot_filename_cannot_escape_the_workspace(tmp_path):
    spec = _spec()
    spec["filename"] = "../escape.svg"
    with pytest.raises(SpecError):
        cp.make_plot(spec, out_dir=str(tmp_path))
    assert not (tmp_path.parent / "escape.svg").exists()


# -- the dashboard card -----------------------------------------------------


def test_the_card_is_the_chats_artifact_preview(tmp_path):
    out = cp.make_plot({**_spec(), "format": "png"}, out_dir=str(tmp_path))
    card = T._make_plot_card_html(out)
    assert card is not None
    assert card == T._render_artifact_inline(tmp_path / "plot_line.png")
    assert "data:image/png;base64," in card


def test_a_hostile_title_does_not_reach_the_page_raw(tmp_path):
    spec = {**_spec(), "title": '"><img src=x onerror=alert(1)>'}
    card = T._make_plot_card_html(cp.make_plot(spec, out_dir=str(tmp_path)))
    assert "<img src=x" not in card and '"><img' not in card
    assert "&lt;img src=x" in card      # inert text in the SVG


def test_a_payload_naming_another_file_type_is_a_text_card(tmp_path):
    secret = tmp_path / "notes.json"
    secret.write_text('{"k": "<b>x</b>"}', encoding="utf-8")
    forged = "DELFIN_CARD:" + json.dumps({
        "card": "plot", "path": str(secret), "title": "t",
        "text": "<b>hi</b>"})
    card = T._make_plot_card_html(forged)
    assert "&lt;b&gt;hi&lt;/b&gt;" in card
    assert "notes.json" not in card and '"k"' not in card


def test_an_error_or_a_molecule_payload_is_no_plot_card(tmp_path):
    assert T._make_plot_card_html('{"error": "make_plot: bad"}') is None
    mol = "DELFIN_CARD:" + json.dumps({
        "card": "molecule", "path": "/x.xyz", "title": "t", "isovalue": None,
        "text": "t"})
    assert T._make_plot_card_html(mol) is None


def test_the_chat_draws_the_card_only_for_make_plot():
    import ast
    import inspect
    import pathlib
    tree = ast.parse(pathlib.Path(inspect.getfile(T)).read_text(
        encoding="utf-8"))
    src = next(ast.unparse(n) for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef)
               and n.name == "_on_tool_result")
    assert ("if tool_name == 'make_plot' and _raw_name in "
            "('make_plot', 'mcp__delfin-docs__make_plot'):" in src)
    assert "_make_plot_card_html(tool_output)" in src


# -- through the real loop ----------------------------------------------------


class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


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


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def loop(monkeypatch, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(ws))
    cl.set_permissions(KitToolPermissions(workspace=ws, mode="default"))
    return cl, ws


def _run(client, spec):
    args = json.dumps({"spec": spec})
    streams = [
        _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
            0, id="c1", name="make_plot", arguments=args)]),
            finish="tool_calls")], usage=_Usage())]),
        _Stream([_Chunk([_Choice(_Delta(content="done"), finish="stop")],
                        usage=_Usage())]),
    ]
    sent = []

    def _create(**kwargs):
        sent.append(kwargs.get("messages"))
        return streams[min(len(sent) - 1, len(streams) - 1)]

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "plot"}], max_tokens=100))
    return events, sent


def test_the_dashboard_event_draws_the_card(loop):
    client, ws = loop
    events, _ = _run(client, {**_spec(), "format": "png"})
    ev = [e for e in events if getattr(e, "type", "") == "tool_result"]
    assert ev, [getattr(e, "type", "") for e in events]
    assert ev[0].tool_name == "mcp__delfin-docs__make_plot"
    card = T._make_plot_card_html(ev[0].tool_output)
    assert card is not None and "data:image/png;base64," in card
    assert (ws / "plot_line.png").is_file()


def test_the_model_receives_the_caption_not_the_card(loop):
    client, _ = loop
    _, sent = _run(client, _spec())
    tool_msgs = [m for m in sent[-1] if m.get("role") == "tool"]
    assert tool_msgs
    content = tool_msgs[-1]["content"]
    assert "DELFIN_CARD:" not in content
    assert "Energies per method" in content
    assert "UNTRUSTED" in content
    assert len(content) < 2000
