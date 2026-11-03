"""A file is shown once in a turn, whichever path would show it.

Three paths put a file in the dashboard chat in one turn: the workspace
diff after a file-writing tool, the card of an explicit ``show_molecule``,
``make_plot`` or ``view_image`` call, and the scan for files the answer
names. Each of them drew the file without asking the others, so one
``write_file`` of an .xyz and ONE ``show_molecule`` of it put the molecule
in the chat twice (field report, Raspberry Pi: [2.2]paracyclophane.xyz);
a shown file the answer then named, a plot followed by any later write,
and an image the turn wrote and then viewed were drawn twice as well.

Driven through the real tab: the scripted OpenAI-compatible endpoint
(tests/agent_terminal_probes/fake_model) answers one step per round, the
real tool loop runs the real tools in a temp workspace, and the views are
counted in what the turn appended to the chat. A file that a tool outside
the file-writing set produces (a geometry or MCP tool) is written by the
scripted endpoint while the round's request is pending, before the next
call runs.

Measured before the change: two views in every sequence below that
expects one, and four for a cube shown at the default level, at 0.05, at
0.05 again, and then named by the answer.
"""

from __future__ import annotations

import base64
import json
import re

import pytest

pytest.importorskip("ipywidgets")

from tests.agent_terminal_probes import fake_model  # noqa: E402

WATER = "3\nwater\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n"
H_CUBE = (
    "test cube\ngenerated\n"
    "1   0.0   0.0   0.0\n"
    "2   1.0   0.0   0.0\n"
    "2   0.0   1.0   0.0\n"
    "2   0.0   0.0   1.0\n"
    "1  0.0   0.0   0.0   0.0\n"
    "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n"
)
# A 3x2 PNG, byte for byte, so no imaging library is needed.
PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAMAAAACCAIAAAASFvFNAAAAEklEQVR4nGP4z8AARBAS"
    "TAAAzfAL9W1sCL0AAAAASUVORK5CYII=")
MARK = "run-scenario"


def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def _views(messages) -> list[tuple[str, str]]:
    """(kind, file name) of every card the messages draw, in order.

    A card is a ``delfin-artifact`` block; its name is the first
    ``<code>`` after it. kind: "3d" (the 3D viewer), "formula" (the
    structural formula card), "file" (any other preview)."""
    out = []
    for m in messages:
        text = str(m.get("content") or "")
        for at in [x.start() for x in re.finditer(r'class="delfin-artifact', text)]:
            tag = text[at:text.index(">", at)]
            name = re.search(r"<code>([^<]*)</code>", text[at:])
            kind = ("3d" if "delfin-chat-mol3d" in tag else
                    "formula" if "delfin-formula-card" in tag else "file")
            out.append((kind, name.group(1) if name else ""))
    return out


def _isovalues(messages) -> list[str]:
    return [iso for m in messages for iso in re.findall(
        r'class="delfin-artifact delfin-chat-mol3d"[^>]*data-mol3d-iso="([^"]*)"',
        str(m.get("content") or ""))]


@pytest.fixture
def run_turn(tmp_path, monkeypatch):
    """run_turn(steps, files={}) -> the chat messages one turn appended.

    steps: what the model answers, one per round -- {"tool", "args"},
    {"tools": [...]} for several calls in one response, {"text"}; a step's
    "before" callable runs while its request is pending (the file a tool
    outside the file-writing set leaves).

    Only the scripted endpoint is reachable: the suite strips every API
    key, and an agent CLI installed on the machine is hidden from the
    tab's provider list (found on PATH it is listed first, and a send
    reached the real CLI). Provider and model are asserted before the
    send, and the endpoint's log after it."""
    import shutil

    import ipywidgets as widgets
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    _which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda name, *a, **k: (
        None if name in ("claude", "codex") else _which(name, *a, **k)))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    ws = tmp_path / "ws"
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()
    built = []

    def _run(steps, files=None):
        for name, data in (files or {}).items():
            (ws / name).write_bytes(data if isinstance(data, bytes)
                                    else data.encode())

        def _plan(messages):
            starts = [i for i, m in enumerate(messages)
                      if m.get("role") == "user"
                      and MARK in fake_model._text_of(m)]
            done = sum(1 for m in messages[starts[-1] if starts else 0:]
                       if m.get("role") == "assistant")
            if done >= len(steps):
                return {"text": "Nothing more."}
            step = dict(steps[done])
            before = step.pop("before", None)
            if before is not None:
                before(ws)
            return step

        monkeypatch.setattr(fake_model, "plan", _plan)
        srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
        built.append(srv)
        monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
        monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
        path = user_settings.get_settings_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"agent": {
            "provider": "ollama", "model": fake_model.MODEL, "mode": "solo",
            "permission_profile": "all_free", "probe_endpoint": False}}))
        ctx = DashboardContext(calc_dir=tmp_path / "calc",
                               archive_dir=tmp_path / "archive",
                               office_dir=tmp_path / "office", agent_dir=ws)
        ctx.run_js = lambda script, **kw: None
        ctx.agent_workspace = str(ws)
        tab, refs = tab_agent.create_tab(ctx)
        built.append(refs)
        state = refs["state"]
        drop = {w.description: w for w in _walk(tab)
                if isinstance(w, widgets.Dropdown)}
        drop["Mode:"].value = "solo"
        assert drop["Provider:"].value == "ollama", drop["Provider:"].value
        assert drop["Model:"].value == fake_model.MODEL, drop["Model:"].value
        inp = next(w for w in _walk(tab) if isinstance(w, widgets.Textarea)
                   and "delfin-agent-input" in w._dom_classes)
        send = next(w for w in _walk(tab) if isinstance(w, widgets.Button)
                    and w.description == "Send")
        before_turn = len(state["chat_messages"])
        inp.value = f"{MARK} please"
        send.click()
        thread = state.get("_turn_thread")
        assert thread is not None, "the send started no turn"
        thread.join(90)
        assert not thread.is_alive(), "the turn did not end"
        log = (tmp_path / "model.log").read_text(encoding="utf-8")
        assert MARK in log, "the turn did not reach the scripted endpoint"
        return state["chat_messages"][before_turn:]

    yield _run
    for item in reversed(built):
        if isinstance(item, dict):
            item["shutdown"](save=False)
        else:
            item.shutdown()


def _write(name, content):
    return {"tool": "write_file", "args": {"path": name, "content": content}}


def _show(name, **extra):
    return {"tool": "show_molecule", "args": {"path": name, **extra}}


# -- molecules -----------------------------------------------------------------


def test_a_written_molecule_shown_by_the_next_round_is_one_view(run_turn):
    """The field sequence: the file is written, the next round shows it."""
    msgs = run_turn([_write("pcp.xyz", WATER), _show("pcp.xyz"),
                     {"text": "The structure pcp.xyz is shown above."}])
    assert _views(msgs) == [("3d", "pcp.xyz")]
    # The preview of an .xyz is the card show_molecule draws, so it is drawn
    # at once, under the write; the call then says where it is.
    card_at = next(i for i, m in enumerate(msgs)
                   if "delfin-artifact" in str(m.get("content")))
    call_at = next(i for i, m in enumerate(msgs)
                   if "show_molecule" in str(m.get("content")))
    assert card_at < call_at
    assert "shown above" in msgs[call_at]["content"]


def test_written_and_shown_in_one_round_is_one_view(run_turn):
    msgs = run_turn([{"tools": [_write("pcp.xyz", WATER), _show("pcp.xyz")]},
                     {"text": "Done: pcp.xyz."}])
    assert _views(msgs) == [("3d", "pcp.xyz")]


def test_text_between_the_write_and_the_call_is_still_one_view(run_turn):
    msgs = run_turn([_write("pcp.xyz", WATER),
                     dict(_show("pcp.xyz"), say="Now the 3D view."),
                     {"text": "Done: pcp.xyz."}])
    assert _views(msgs) == [("3d", "pcp.xyz")]


def test_a_shown_file_the_answer_names_is_one_view(run_turn):
    msgs = run_turn([_show("old.xyz"),
                     {"text": "This is old.xyz, the structure from before."}],
                    files={"old.xyz": WATER})
    assert _views(msgs) == [("3d", "old.xyz")]


def test_a_file_another_tool_wrote_is_one_view(run_turn):
    """A tool outside the file-writing set leaves geo.xyz, show_molecule
    draws it, a later write runs the workspace diff (which sees geo.xyz as
    new this turn), and the answer names it."""
    def _geometry_tool(ws):
        (ws / "geo.xyz").write_text(WATER)

    msgs = run_turn([dict(_show("geo.xyz"), before=_geometry_tool),
                     _write("notes.txt", "n\n"),
                     {"text": "geo.xyz is the built structure."}])
    assert _views(msgs) == [("3d", "geo.xyz")]


def test_a_formula_written_then_shown_is_the_formula_card(run_turn):
    """show_molecule draws .smi as the structural formula card; the plain
    preview of the written file is not drawn beside it."""
    pytest.importorskip("rdkit")
    msgs = run_turn([_write("ethanol.smi", "CCO ethanol\n"),
                     _show("ethanol.smi"), {"text": "Done."}])
    assert _views(msgs) == [("formula", "ethanol.smi")]
    card_at = next(i for i, m in enumerate(msgs)
                   if "delfin-artifact" in str(m.get("content")))
    call_at = next(i for i, m in enumerate(msgs)
                   if "show_molecule" in str(m.get("content")))
    assert call_at < card_at, "the card belongs under the call that drew it"


def test_a_written_formula_nobody_shows_is_still_drawn(run_turn):
    pytest.importorskip("rdkit")
    msgs = run_turn([_write("ethanol.smi", "CCO ethanol\n"),
                     {"text": "Wrote it."}])
    assert _views(msgs) == [("file", "ethanol.smi")]
    roles = [m["role"] for m in msgs]
    card_at = next(i for i, m in enumerate(msgs)
                   if "delfin-artifact" in str(m.get("content")))
    assert "assistant" in roles[card_at:], "the preview came after the answer"


def test_a_cube_written_then_shown_is_one_view_at_the_chosen_level(run_turn):
    msgs = run_turn([_write("h.cube", H_CUBE), _show("h.cube", isovalue=0.05),
                     {"text": "h.cube at 0.05."}])
    assert _views(msgs) == [("3d", "h.cube")]
    assert _isovalues(msgs) == ["0.05"]


def test_a_second_isovalue_is_a_second_view_and_a_repeat_is_not(run_turn):
    msgs = run_turn([_show("h.cube"), _show("h.cube", isovalue=0.05),
                     _show("h.cube", isovalue=0.05),
                     {"text": "h.cube at two isovalues."}],
                    files={"h.cube": H_CUBE})
    assert _views(msgs) == [("3d", "h.cube"), ("3d", "h.cube")]
    assert _isovalues(msgs) == ["0.02", "0.05"]
    assert sum("shown above" in str(m.get("content")) for m in msgs) == 1


def test_a_rewritten_file_is_new_content(run_turn):
    """One view per file CONTENT: the same name written again with other
    coordinates is shown again, as the workspace diff always did. (The
    write needs a read first: write_file refuses to overwrite a file the
    session has not read.)"""
    msgs = run_turn([_show("old.xyz"),
                     {"tool": "read_file", "args": {"path": "old.xyz"}},
                     _write("old.xyz", WATER.replace("0.96", "0.97")),
                     {"text": "Moved one atom in old.xyz."}],
                    files={"old.xyz": WATER})
    assert _views(msgs) == [("3d", "old.xyz"), ("3d", "old.xyz")]
    cards = [str(m.get("content")) for m in msgs
             if "delfin-chat-mol3d" in str(m.get("content"))]
    assert "0.97" not in cards[0] and "0.97" in cards[1]


# -- plots and images ------------------------------------------------------------


def test_a_plot_followed_by_a_write_is_one_view(run_turn):
    spec = {"kind": "line", "x": [1, 2, 3], "y": [1.0, 4.0, 9.0],
            "title": "t", "filename": "curve.svg"}
    msgs = run_turn([{"tool": "make_plot", "args": {"spec": spec}},
                     _write("notes.txt", "n\n"),
                     {"text": "The figure is curve.svg."}])
    assert _views(msgs) == [("file", "curve.svg")]


def test_an_image_the_turn_wrote_and_viewed_is_one_view(run_turn):
    """fig.png appears with the write (a script's figure), the workspace
    diff draws it, and view_image of it adds no second copy."""
    def _figure(ws):
        (ws / "fig.png").write_bytes(PNG)

    msgs = run_turn([dict(_write("notes.txt", "n\n"), before=_figure),
                     {"tool": "view_image", "args": {"path": "fig.png"}},
                     {"text": "fig.png looks right."}])
    assert _views(msgs) == [("file", "fig.png")]


def test_a_viewed_image_alone_is_still_drawn(run_turn):
    msgs = run_turn([{"tool": "view_image", "args": {"path": "fig.png"}},
                     {"text": "fig.png looks right."}],
                    files={"fig.png": PNG})
    assert _views(msgs) == [("file", "fig.png")]
