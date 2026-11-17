"""Both surfaces rerun a turn that cites a bare path that does not exist.

The dashboard used to own this check -- its inline Verify-Enforcement block
scanned a finished answer with ``scan_for_ungrounded_code_claims`` and forced
one correction turn when it found a ``kind == "nonexistent"`` flag. The engine
never did: its guard's ``scan_for_ungrounded_location_claims`` only flags
``file_line`` claims, so a bare (line-less) nonexistent path was invisible to
the terminal, which therefore answered once and never corrected.

The two tests characterise the rerun through the dashboard and through the
engine (the engine is what the terminal layer drives for a turn). The engine
test flipped from 1 stream / no feedback to 2 streams / one [Verify] row when
the decision moved into the engine's guard.

Two further pins guard the move's edges:

* A path the turn EXECUTED as a tool call is not an ungrounded citation --
  "python solver.py" grounds a claim about solver.py, and the guard
  must stay silent (one model call) even though the file exists nowhere.
* The [Verify] row for a nonexistent path must say the paths DO NOT EXIST
  (``code_claim_feedback``), not send the model to read a file that cannot
  be read (``location_claim_feedback`` wording).
* A cited path that EXISTS but was not read stays a dashboard-side SOFT
  hint ("cited but not read this turn") that forces nothing.
"""

from __future__ import annotations

import json
import textwrap

import pytest
from unittest.mock import MagicMock, patch

from tests.agent_terminal_probes import fake_model  # noqa: E402

# A path that does not exist in any repo.  ``verify_guard`` swallows paths
# that look like they might exist in a sibling root, so the name is
# deliberately inside a root that is not on disk.
NONEXISTENT = "delfin/src/nonexistent_file.py"
# A bare path citation with NO line number: the engine's location scanner
# only flags ``file_line`` claims, so this form slips past the terminal's
# guard while the dashboard's code scanner still flags it ``nonexistent``.
CITING = f"See `{NONEXISTENT}` for the loader."
#: A path that EXISTS but was not read this turn: the soft case.
SOFT_NAME = "exists_unread.py"
SOFT_CITING = f"See `{SOFT_NAME}` for the loader."


@pytest.fixture
def agent_tree(tmp_path):
    """A minimal agent pack, the shape the engine tests construct."""
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


# --------------------------------------------------------------------------
# terminal path -- today: no rerun
# --------------------------------------------------------------------------

def _client(replies):
    """A fake backend that returns ``replies`` in order and counts calls."""
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        yield from _text_delta(replies[min(i, len(replies) - 1)])

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


def _text_delta(text):
    from delfin.agent.api_client import StreamEvent
    yield StreamEvent(type="text_delta", text=text)
    yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)


def _term_engine(client, agent_tree):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def _tool_client(reply, tool_calls=()):
    """A fake backend that emits ``tool_calls`` ((name, input) pairs) with
    results before streaming ``reply`` -- the shape the engine's exec ledger
    is fed from (mirrors tests/test_grounding_and_budget.py::_tool_client)."""
    from delfin.agent.api_client import StreamEvent
    fake = MagicMock()
    fake._observed_files_session = set()

    def _stream(*a, **k):
        for name, payload in tool_calls:
            yield StreamEvent(type="tool_use", tool_name=name,
                              tool_input=payload)
            yield StreamEvent(type="tool_result", tool_name=name,
                              tool_output="ok")
        yield from _text_delta(reply)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


def _verify_rows(engine):
    return [m for m in engine.messages
            if m.get("_harness")
            and str(m.get("content", "")).startswith("[Verify]")]


def test_terminal_reruns_a_nonexistent_path(agent_tree):
    """The engine's grounding guard now treats a bare nonexistent citation
    as a correction trigger, so the terminal reruns the turn exactly as the
    dashboard does: two streams, one [Verify] feedback row appended."""
    client = _client([CITING, "Nothing more."])
    engine = _term_engine(client, agent_tree)
    engine.stream_response("What is the fix?")
    assert client.stream_message.call_count == 2, (
        "the terminal now forces one correction turn for a nonexistent path")
    feedback = _verify_rows(engine)
    assert len(feedback) == 1, f"terminal verification reruns: {len(feedback)}"
    assert NONEXISTENT in feedback[0]["content"]


def test_a_path_the_turn_executed_is_not_an_ungrounded_citation(agent_tree):
    """``solver.py`` exists nowhere, but the SAME turn ran it as a tool call
    and the result shows it ran. The claim is grounded by the turn's own
    action, so the guard stays silent: one model call, no [Verify] row."""
    client = _tool_client(
        "Das Skript solver.py läuft fehlerfrei.",
        tool_calls=(("Bash", '{"command": "python solver.py"}'),))
    engine = _term_engine(client, agent_tree)
    engine.stream_response("prüf das skript")
    assert client.stream_message.call_count == 1, (
        "an executed path must not cost a correction turn")
    assert _verify_rows(engine) == [], _verify_rows(engine)


def test_the_verify_row_says_the_paths_do_not_exist(agent_tree):
    """The [Verify] row for a nonexistent path must come from
    ``code_claim_feedback`` ('the cited paths ... do not exist'), not from
    ``location_claim_feedback`` ('not backed by a file read or grep') -- the
    latter orders the model to go and read a file that cannot be read."""
    client = _client([CITING, "Nothing more."])
    engine = _term_engine(client, agent_tree)
    engine.stream_response("What is the fix?")
    feedback = _verify_rows(engine)
    assert len(feedback) == 1, f"terminal verification reruns: {len(feedback)}"
    assert "do not exist" in feedback[0]["content"], feedback[0]["content"]
    assert "not backed by a file read" not in feedback[0]["content"]


# --------------------------------------------------------------------------
# dashboard path -- today: one rerun
# --------------------------------------------------------------------------

def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def _dashboard_tab(tmp_path, monkeypatch, plan, answer_file=None):
    """Drive one dashboard turn headlessly (the run_turn fixture shape);
    returns (state, refs) with the turn already finished."""
    import shutil
    import ipywidgets as widgets  # noqa: F401
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    MARK = "run-scenario"

    _which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda name, *a, **k: (
        None if name in ("claude", "codex") else _which(name, *a, **k)))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    ws = tmp_path / "ws"
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()
    if answer_file is not None:
        (ws / answer_file).write_text("x = 1\n")

    monkeypatch.setattr(fake_model, "plan", lambda m: plan(m))
    srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
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

    MARK = "run-scenario"
    drop = {w.description: w for w in _walk(tab)
            if isinstance(w, widgets.Dropdown)}
    drop["Mode:"].value = "solo"
    inp = next(w for w in _walk(tab)
               if isinstance(w, widgets.Textarea)
               and "delfin-agent-input" in w._dom_classes)
    send = next(w for w in _walk(tab)
                if isinstance(w, widgets.Button)
                and w.description == "Send")
    inp.value = f"{MARK} please"
    send.click()
    thread = refs["state"].get("_turn_thread")
    assert thread is not None, "the send started no turn"
    thread.join(120)
    assert not thread.is_alive(), "the turn did not end"
    return refs["state"], refs


def _dashboard_plan(*, first, later="Nothing more."):
    """Scripted plan: the first assistant round answers ``first``, any
    later round (the forced correction) answers ``later``."""
    MARK = "run-scenario"

    def _plan(messages):
        starts = [i for i, m in enumerate(messages)
                  if m.get("role") == "user" and MARK in fake_model._text_of(m)]
        done = sum(1 for m in messages[starts[-1] if starts else 0:]
                   if m.get("role") == "assistant")
        if done >= 1:
            return {"text": later}
        return {"text": first}

    return _plan


def test_dashboard_reruns_a_nonexistent_path(tmp_path, monkeypatch):
    """The dashboard still shows exactly one [Verify] citation correction
    and re-streams the turn (a second model call) when the answer cites a
    path that does not exist. The decision is the engine's; the tab renders
    the engine's verdict."""
    fake_plan = _dashboard_plan(first=CITING)
    state, refs = _dashboard_tab(tmp_path, monkeypatch, fake_plan)
    try:
        engine = state.get("engine")
        assert engine is not None, "the turn created no engine"
        verify = [m["content"] for m in engine.messages
                  if m.get("_harness") and m.get("role") == "user"
                  and str(m["content"]).startswith("[Verify]")]
        assert len(verify) == 1, f"dashboard verification reruns: {len(verify)}"
        assert NONEXISTENT in verify[0]
    finally:
        refs["shutdown"](save=False)


def test_the_dashboard_still_hints_at_a_cited_but_unread_path(
        tmp_path, monkeypatch):
    """A cited path that EXISTS but was not read this turn stays a SOFT
    dashboard hint ("cited but not read this turn"). It forces nothing and
    is not a rerun -- it went with the deleted dashboard block before this
    pin restored it."""
    fake_plan = _dashboard_plan(first=SOFT_CITING)
    state, refs = _dashboard_tab(tmp_path, monkeypatch, fake_plan,
                                 answer_file=SOFT_NAME)
    try:
        engine = state.get("engine")
        assert engine is not None, "the turn created no engine"
        said = [str(m.get("content", "")) for m in state["chat_messages"]
                if m.get("role") == "system"]
        hint = [s for s in said if "cited but not read this turn" in s]
        assert hint, said
        assert SOFT_NAME in hint[0], hint
        # A soft flag forces nothing: no correction turn.
        assert _verify_rows(engine) == [], _verify_rows(engine)
    finally:
        refs["shutdown"](save=False)
