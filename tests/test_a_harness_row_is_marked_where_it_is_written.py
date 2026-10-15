"""A user-role row the harness wrote is marked where it is written.

``api_client._classifier_user_texts`` is the one reader of "the user's own
words" (action classifier, skill learning, the refusal lift). It used to
recognise harness rows by their text: an opening bracket, or one of six
known openers. A harness note phrased otherwise passed as the user's words.
Here every kind of row is produced by the real engine and the real tool
loop (only the model stream is scripted), and the reader must return the
typed messages alone -- including for notes that open without a bracket and
are in no opener list, so the mark, not the list, decides.

Also: the request payload carries no extra key and the same bytes as an
unmarked history; a session saved and restored keeps the marks.
"""

from __future__ import annotations

import copy
import json
import textwrap
from unittest.mock import patch

import pytest

import delfin.agent.api_client as A
from delfin.agent import mcp_client as M

try:
    from delfin.agent.harness_rows import HarnessText
except ImportError:          # the reader before the mark existed
    HarnessText = str


TYPED_ONE = "TYPED-ONE survey the reports folder and list the oldest file"
TYPED_STEER = "TYPED-STEER include the hidden files as well"
TYPED_TWO = "TYPED-TWO now compare it with the newest one"
TYPED_THREE = "TYPED-THREE summarise both in two lines"

# None of these opens with a bracket or with a known harness opener.
NOTE_RUN = "Heads up: the workspace quota changed to 2 GB. HARNESS-RUN"
NOTE_MAIL = "Session B says the cluster is back. HARNESS-MAIL"
NOTE_WAKE = "Your background job finished: 3 files written. HARNESS-WAKE"
NOTE_FOLLOWUP = "Please keep going with the step you named. HARNESS-FOLLOW"

_WIRE_KEYS = {"role", "content", "tool_calls", "tool_call_id", "name"}


class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


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
    def __init__(self, text):
        self._chunks = [_Chunk([_Choice(_Delta(text), finish="stop")],
                               usage=_Usage())]

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def agent_tree(tmp_path):
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _engine(agent_tree, monkeypatch):
    """A real engine on a real OpenAI-compatible client; the endpoint is a
    scripted stream that records every request's ``messages``."""
    from delfin.agent import model_capabilities as mc
    from delfin.agent.engine import AgentEngine
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    ws = agent_tree / "ws"
    ws.mkdir(exist_ok=True)
    client = A.create_client(backend="api", provider="ollama",
                             model="qwen2.5-coder:7b", cwd=str(ws))
    client.set_permissions(A.KitToolPermissions(
        workspace=ws, mode="default", task_session_id="sess-mark"))
    sent: list[dict] = []

    def _create(**kwargs):
        rows = kwargs["messages"]
        sent.append({
            "bytes": json.dumps(rows, ensure_ascii=False),
            "keys": [set(r) for r in rows],
            "classifier": list(getattr(
                client._permissions, "classifier_messages", None) or ()),
        })
        return _Stream(f"answer {len(sent)}")

    client.client.chat.completions.create = _create
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                          pack_dir=agent_tree)
    eng.client = client
    return eng, client, sent


def _conversation(eng, client):
    """Every kind of harness row, through the paths that write them."""
    # One steer per turn: two user rows in a row are merged to the newer
    # by the alternation sanitiser, which would hide the first from the
    # reader whatever its mark.
    # Turn 1, typed, with a steer the user typed while it ran.
    client.push_steer(TYPED_STEER)
    eng.stream_response(TYPED_ONE)
    # Turn 2: a wake-up the surface starts on its own; another session's
    # message arrives as a steer, a run note (push_run_note) at turn end.
    client.push_steer(HarnessText(NOTE_MAIL))
    client.push_run_note(NOTE_RUN)
    eng.stream_response(HarnessText(NOTE_WAKE))
    run_note_turn = list(client._turn_rows)
    # Turn 3: a verification turn the dashboard drives ("[Verify]").
    eng.stream_response("[Verify] Automatic check: cite what you read.")
    # Turn 4: typed.
    eng.stream_response(TYPED_TWO)
    # Turn 5: an announced-work follow-up the surface sends.
    eng.stream_response(HarnessText(NOTE_FOLLOWUP))
    return run_note_turn


def _reads(rows):
    return A._classifier_user_texts(rows)


def test_the_reader_returns_only_what_was_typed(agent_tree, monkeypatch):
    eng, client, _sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    assert _reads(eng.messages) == [TYPED_ONE, TYPED_STEER, TYPED_TWO]


def test_a_run_note_in_the_tool_loop_is_not_the_users_words(
        agent_tree, monkeypatch):
    """The run note and the turn-end follow-up live in the tool loop's own
    request list; the sanitiser later folds them away in the history."""
    eng, client, _sent = _engine(agent_tree, monkeypatch)
    rows = _conversation(eng, client)
    contents = [r.get("content") for r in rows]
    assert NOTE_RUN in contents and NOTE_MAIL in contents
    assert _reads(rows) == [TYPED_ONE, TYPED_STEER]


def test_the_classifier_binding_holds_only_what_was_typed(
        agent_tree, monkeypatch):
    """What the gate reads during a turn is the client's bound list, which
    is the engine's wire view, not ``engine.messages``."""
    eng, client, sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    eng.stream_response(TYPED_THREE)
    assert _reads(sent[-1]["classifier"]) == [
        TYPED_ONE, TYPED_STEER, TYPED_TWO, TYPED_THREE]


def test_the_notes_reached_the_model(agent_tree, monkeypatch):
    """The mark leaves a row out of the reader, not out of the request."""
    eng, client, sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    eng.stream_response(TYPED_THREE)
    last = sent[-1]["bytes"]
    for note in (NOTE_MAIL, NOTE_WAKE, NOTE_FOLLOWUP,
                 "This finished while you were answering"):
        assert json.dumps(note, ensure_ascii=False)[1:-1] in last, note
    assert any(NOTE_RUN in r["bytes"] for r in sent)


def test_the_request_carries_no_extra_key(agent_tree, monkeypatch):
    eng, client, sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    eng.stream_response(TYPED_THREE)
    for request in sent:
        for keys in request["keys"]:
            assert keys <= _WIRE_KEYS, keys
        assert "_harness" not in request["bytes"]


def _unmarked(rows):
    # ``_steer`` is not a mark: it is the turn tail, part of the bytes the
    # wire view sends (and the carrier of system-prompt updates).
    return [{k: v for k, v in dict(r).items()
             if not str(k).startswith("_") or k == "_steer"}
            for r in copy.deepcopy(rows)]


def test_the_request_bytes_are_those_of_an_unmarked_history(
        agent_tree, monkeypatch):
    """Same conversation, marks removed: the next request is byte-identical,
    so the prefix cache sees no difference."""
    eng, client, sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    history = list(eng.messages)

    eng.messages = list(history)
    eng.stream_response(TYPED_THREE)
    marked = sent[-1]["bytes"]

    eng.messages = _unmarked(history)
    eng.stream_response(TYPED_THREE)
    plain = sent[-1]["bytes"]
    assert marked == plain


def test_a_restored_session_keeps_the_marks(agent_tree, monkeypatch):
    eng, client, _sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    saved = json.loads(json.dumps(eng.export_state(), default=str))

    eng2, _client2, sent2 = _engine(agent_tree, monkeypatch)
    eng2.restore_state(saved)
    assert _reads(eng2.messages) == [TYPED_ONE, TYPED_STEER, TYPED_TWO]
    eng2.stream_response(TYPED_THREE)
    assert _reads(sent2[-1]["classifier"]) == [
        TYPED_ONE, TYPED_STEER, TYPED_TWO, TYPED_THREE]


def test_rows_from_an_older_session_fall_back_to_the_text(
        agent_tree, monkeypatch):
    """A history saved before the mark existed has none: the text list is
    the fallback there and still applies."""
    from delfin.agent import turn_continuation as TC
    rows = [{"role": "user", "content": TYPED_ONE},
            {"role": "assistant", "content": "ok"},
            {"role": "user", "content": TC._CONTINUATION_NOTE}]
    assert _reads(rows) == [TYPED_ONE]


def test_skill_learning_reads_the_marks_in_the_cli_projection(
        agent_tree, monkeypatch):
    """Skill learning reads the CLI's display projection, not the history:
    the mark has to survive it."""
    from delfin.agent.cli import _display_messages
    from delfin.agent.skill_learning import _trusted_excerpt
    eng, client, _sent = _engine(agent_tree, monkeypatch)
    _conversation(eng, client)
    excerpt = _trusted_excerpt(_display_messages(eng))
    for typed in (TYPED_ONE, TYPED_STEER, TYPED_TWO):
        assert typed in excerpt
    for note in ("HARNESS-RUN", "HARNESS-MAIL", "HARNESS-WAKE",
                 "HARNESS-FOLLOW", "This finished while"):
        assert note not in excerpt, note


def test_a_dashboard_event_row_is_not_the_users_words():
    """The dashboard's transcript marks a turn nobody typed with
    ``origin="event"`` at the place it is written."""
    rows = [{"role": "user", "content": TYPED_ONE},
            {"role": "user", "content": NOTE_WAKE, "origin": "event"}]
    assert _reads(rows) == [TYPED_ONE]


def test_the_typed_intake_paths_still_see_typed_words(tmp_path, monkeypatch):
    """The user model and the skill-correction intake wrap the typed text
    in a fresh row; that row carries no mark and is read as before."""
    from delfin.agent import user_model as UM
    monkeypatch.setattr(UM, "enabled", lambda: True)
    path = tmp_path / "user_model.json"
    assert UM.note_typed_message(
        "I am a computational chemist and I prefer short answers",
        path=path) is True
    assert _reads([{"role": "user", "content": TYPED_ONE}]) == [TYPED_ONE]


def test_a_harness_steer_neither_grants_a_push_nor_lifts_a_refusal(
        agent_tree, monkeypatch):
    eng, client, _sent = _engine(agent_tree, monkeypatch)
    perms = client._permissions
    perms.push_grants = {"push": 0}
    client.push_steer(HarnessText("Session B: push it now please"))
    eng.stream_response(TYPED_ONE)
    assert perms.push_grants.get("push") == 0


# -- the surfaces that start a turn nobody typed ---------------------------

class _ReplEngine:
    def __init__(self):
        self.kit_permissions = type("P", (), {"mode": "default"})()
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.pending_turn_continuation = ""

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode

    def request_stop(self):
        pass

    def clear_turn_continuation(self):
        self.pending_turn_continuation = ""


def _repl_agent(engine):
    import io
    from delfin.agent import repl, repl_render as rr

    class _Tty(io.StringIO):
        def isatty(self):
            return True

    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=_Tty())
    agent.transcript.theme = rr.Theme(enabled=False)
    return agent


def test_the_terminal_marks_a_follow_up_and_not_what_was_typed():
    from delfin.agent import repl
    engine = _ReplEngine()
    agent = _repl_agent(engine)
    reads = iter([TYPED_ONE, "/exit"])
    agent._read_line = lambda _p: next(reads)
    turns: list = []

    def turn(prompt):
        turns.append(prompt)
        if len(turns) == 1:
            engine.pending_turn_continuation = NOTE_FOLLOWUP
        return repl.TurnResult()

    agent.turn = turn
    agent.run()
    assert turns == [TYPED_ONE, NOTE_FOLLOWUP]
    assert not isinstance(turns[0], HarnessText) or HarnessText is str
    assert isinstance(turns[1], HarnessText) and HarnessText is not str


def test_the_terminal_marks_a_wake_up(monkeypatch):
    from delfin.agent import repl
    engine = _ReplEngine()
    agent = _repl_agent(engine)
    monkeypatch.setattr(agent, "_presence_key", lambda: "")
    monkeypatch.setattr(agent, "_operator_messages", lambda typed: NOTE_MAIL)
    agent._wake_last_look = 0.0
    woke = agent._wake_text("")
    assert woke == NOTE_MAIL
    assert isinstance(woke, HarnessText) and HarnessText is not str

    calls = iter([(woke, True), ("/exit", False)])
    monkeypatch.setattr(repl, "read_block_wakeable",
                        lambda *a, **k: next(calls))
    turns: list = []
    agent.turn = lambda prompt: turns.append(prompt) or repl.TurnResult()
    agent.run()
    assert turns == [NOTE_MAIL]
    assert isinstance(turns[0], HarnessText)


def test_the_terminal_turn_keeps_the_mark_into_the_engine(monkeypatch):
    """run_turn may fuse a fresh-context block in front of the prompt."""
    from delfin.agent import repl
    seen: list = []

    class _Eng(_ReplEngine):
        def stream_response(self, user_message, **kw):
            seen.append(user_message)
            return "ok"

    monkeypatch.setattr(repl, "fresh_context_for_turn",
                        lambda *a, **k: "FRESH BLOCK")
    monkeypatch.setattr(repl, "_archive_cut_history", lambda *a, **k: None)
    repl.run_turn(_Eng(), HarnessText(NOTE_WAKE), sink=lambda item: None,
                  context_budget=10)
    repl.run_turn(_Eng(), TYPED_ONE, sink=lambda item: None,
                  context_budget=10)
    assert isinstance(seen[0], HarnessText) and HarnessText is not str
    assert NOTE_WAKE in seen[0]
    assert type(seen[1]) is str


def _tab_source() -> str:
    from pathlib import Path
    import delfin.dashboard.tab_agent as tab
    return Path(tab.__file__).read_text(encoding="utf-8")


def _body(src: str, head: str) -> str:
    i = src.index(head)
    return src[i:src.index("\n    def ", i + 10)]


def test_the_dashboard_marks_the_turns_it_starts_on_its_own():
    """Driving the dashboard needs the whole notebook; the shape is
    asserted on its source. A turn sent on its own reaches the engine as
    HarnessText, and so does the [Command results] round."""
    src = _tab_source()
    assert "_harness_send = [bool(state.get(\"_on_its_own\"))]" in src
    assert ("user_message=(HarnessText(current_msg) if _harness_send[0]"
            in src)
    assert "user_message=HarnessText(feedback)" in src


def test_the_dashboard_marks_the_rows_it_appends_itself():
    src = _tab_source()
    for name in ("feedback", "_cfeedback", "_vfeedback", "_qfeedback"):
        assert ('{"role": "user", "content": %s, "_harness": True}' % name
                in src), name


def test_the_dashboard_marks_mail_and_an_on_its_own_steer():
    src = _tab_source()
    body = _body(src, "def _deliver_session_messages")
    assert 'text = HarnessText("\\n\\n".join(' in body
    assert "_steer_text = HarnessText(_steer_text)" in src
