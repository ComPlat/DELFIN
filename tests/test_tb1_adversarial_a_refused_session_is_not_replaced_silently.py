"""Adversarial tests for TB1: a tab refused a session must not start another.

The defect (mine, verified on main ``7abca4483``): ``tab_agent.py:12088-12096``
appends "This conversation is held by a live session ..." and ``return``s
without recording anything. The tab then holds no session, and the next
``_on_send`` reaches ``_ensure_engine`` (``:21497`` -> ``:12516``), which mints
a brand-new session because ``state["engine"] is None``. The person's text
lands in a stranger session and the held conversation never sees it.

These tests are contract-based on purpose: they assert what the person
observes (a session is or is not created, the text is or is not kept, the
choice line is or is not shown) and never the patch's internal key names, so
they hold for any reasonable implementation.

Every test drives the real widgets of ``tab_agent.create_tab`` -- the real
Send button, the real input box, the real Load button -- never an inner
helper. Only ``AgentEngine`` is replaced; each minted engine carries a unique
id, which is what makes "a new session was created" observable at all.

Direction matters (builder checklist, red tests in both directions):
- tests marked MUST STAY GREEN pass on main today and must still pass after
  the patch -- they are the guard against a patch that blocks every send;
- ``test_typing_after_a_refusal...`` is RED on main and is the defect itself.

No host, user or account name appears here; HOME and the runtime dir are
redirected to ``tmp_path`` so no real ``~/.delfin`` is read.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

SID = "tb1a000000000001"
FREE = "tb1a0000000000f0"


def _walk(widget):
    yield widget
    children = getattr(widget, "children", ()) or ()
    for child in children:
        yield from _walk(child)


class _FakeEngine:
    """What the tab uses of an engine, and an id that can be counted."""

    def __init__(self, session_id):
        self.session_id = session_id
        self.messages = []
        self.route = ["solo_agent"]
        self.role_outputs = {}
        self.current_role = "solo_agent"
        self.current_role_index = 0
        self.is_cycle_complete = False
        self.cycle_phase = "idle"
        self.kit_permissions = None
        self.client = type("C", (), {"kill": lambda self: None})()

    def get_status(self):
        """The keys ``_render_status`` (:14070) reads positionally; a thin
        dict raises KeyError in the refresh thread and hides real failures."""
        return {"session_id": self.session_id, "mode": "solo",
                "backend": "api", "role": self.current_role,
                "role_index": 0, "role_total": 1,
                "input_tokens": 0, "output_tokens": 0, "cost_usd": 0.0,
                "provider": "", "model": "", "cached_tokens": 0}

    def export_state(self):
        return {"session_id": self.session_id, "engine_messages": []}

    def restore_state(self, data):
        """``_hand_state_to`` (:12468) is the one path a saved session takes
        into an engine; the tab aborts the restore when this raises, so a
        stub without it makes every load look refused."""
        self.restored = dict(data)
        self.session_id = str(data.get("session_id") or self.session_id)
        self.messages = list(data.get("engine_messages") or [])

    def build_handoff_message(self, task):
        return str(task)

    def compact_for_next_role(self):
        return None

    def advance_role(self):
        return False

    def stream_response(self, *a, **k):
        """No turn: these tests measure the session a send would start,
        never a model answer. Raising is deliberate -- a turn that ran here
        would write files behind the assertions."""
        raise AssertionError("a refused tab must not run a turn")


@pytest.fixture
def env(tmp_path, monkeypatch):
    """A private HOME, a private session store, and a fake engine factory."""
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.agent import session_store as SS
    from delfin.agent import engine as engine_mod

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "runtime"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "ws"))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setattr(SS, "_SESSIONS_DIR", tmp_path / "sessions")
    (tmp_path / "sessions").mkdir()
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "ollama", "model": "fake-model", "mode": "solo",
        "permission_profile": "all_free", "probe_endpoint": False}}))

    created = []

    def _factory(**kw):                      # noqa: ARG001 -- real ctor shape
        new = f"minted-{len(created) + 1:03d}"
        created.append(new)
        return _FakeEngine(new)

    monkeypatch.setattr(engine_mod, "AgentEngine", _factory)
    return {"tmp": tmp_path, "created": created, "SS": SS}


def _chat_text(state) -> str:
    """Everything the person can read in the chat, including the choice line."""
    return "\n".join(str(m.get("content") or "")
                     for m in state.get("chat_messages") or [])


def _session_files(env):
    """The sessions the store holds on disk -- what "a new session" means."""
    d = env["tmp"] / "sessions"
    return sorted(p.name for p in d.glob("*.json")
                  if not p.name.endswith(".turn.json"))


def _sessions_holding(env, text):
    """The session files whose content carries ``text``.

    The person's loss is not the file, it is the content: a session file for
    the new engine is written by the tab's own autosave on the first turn
    (``_auto_save_session``, ``:11943``), so "a session exists" is a weak
    signal while "the typed text is in a session" is the defect itself.
    """
    return sorted(p.name for p in (env["tmp"] / "sessions").glob("*.json")
                  if text in p.read_text(errors="replace"))


def _held_by_a_live_process(env):
    """A real second process holds SID's writer lock, as the field does.

    The pid must be alive and foreign: a lock held by this test's own pid
    reads as its own and would prove nothing (the lesson of
    ``test_a_reload_takes_over_from_its_dying_kernel``).
    """
    from tests.conftest import child_env
    holder_dir = env["tmp"] / "holder"
    holder_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                            env=child_env(holder_dir))
    lock = env["SS"]._SESSIONS_DIR / f"{SID}.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.write_text(json.dumps({"pid": proc.pid, "ts": time.time(),
                                "host": env["SS"]._this_host()}))
    assert env["SS"].session_lock_holder(SID) is not None, (
        "the fixture must produce a holder; without one this proves nothing")
    return proc


def _write_session(env, session_id, title, turn):
    (env["tmp"] / "sessions" / f"{session_id}.json").write_text(json.dumps({
        "session_id": session_id, "title": title, "mode": "solo",
        "workspace": str(env["tmp"] / "ws"),
        "created_at": 1.0, "updated_at": 1.0,
        "engine_messages": [
            {"role": "user", "content": turn},
            {"role": "assistant", "content": f"answer to {turn}"},
        ],
    }))


def _tab(env, initial_session_id=""):
    """The real agent tab, its real widgets, and its refs."""
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    ctx = DashboardContext(calc_dir=env["tmp"] / "calc",
                           agent_dir=env["tmp"] / "ws",
                           archive_dir=env["tmp"] / "archive",
                           office_dir=env["tmp"] / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(env["tmp"] / "ws")
    ctx.initial_session_id = initial_session_id
    ctx.presence_key = "tb1"
    tab, refs = tab_agent.create_tab(ctx)
    found = {"send": None, "input": None, "load": None, "session": None,
             "fork": None}
    for w in _walk(tab):
        if isinstance(w, widgets.Button):
            # Two buttons are described "Send": the chat's, and the hidden
            # bug-report one built earlier in the toolbar (``:9045``). The
            # class is what the page scripts use to find the chat's, so a
            # description match here would press the wrong one.
            classes = getattr(w, "_dom_classes", ()) or ()
            if "delfin-agent-send" in classes and found["send"] is None:
                found["send"] = w
            elif w.description == "Load" and found["load"] is None:
                found["load"] = w
            elif w.description == "Fork" and found["fork"] is None:
                found["fork"] = w
        elif isinstance(w, widgets.Textarea) and found["input"] is None:
            if "delfin-agent-input" in (getattr(w, "_dom_classes", ()) or ()):
                found["input"] = w
        elif isinstance(w, widgets.Dropdown) and w.description == "Session:":
            found["session"] = w
    assert all(found.values()), f"widgets not found: {sorted(found)}"
    found["refs"] = refs
    return found


def _open_held_conversation(tab, env):
    """Pick the held session in the list and press Load, as the person does."""
    session_ids = [v for _label, v in tab["session"].options]
    assert SID in session_ids, (
        "the held session must be offered in the list the person picks from")
    tab["session"].value = SID
    tab["load"].click()


def _type_and_send(tab, text):
    tab["input"].value = text
    tab["send"].click()


# ---------------------------------------------------------------------------
# The defect: red on main
# ---------------------------------------------------------------------------


def test_typing_after_a_refusal_starts_no_session_and_keeps_the_text(env):
    """The person's text after a refused open must not become a new session.

    Red on main: Load is refused for the held conversation, the tab keeps no
    session, and Send mints a fresh one -- the text lands in a stranger.
    """
    _write_session(env, SID, "held job", "held-turn")
    proc = _held_by_a_live_process(env)
    try:
        before = _session_files(env)
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            _open_held_conversation(tab, env)
            assert "held by a live session" in _chat_text(state), (
                "the fixture must reach the refusal; otherwise the rest "
                "measures the wrong path")
            assert env["SS"].session_lock_holder(SID) is not None, (
                "a live holder with its window keeps the lock")

            typed = "this belongs in the held conversation"
            _type_and_send(tab, typed)

            assert env["created"] == [], (
                "sending after a refusal minted a new session: the text went "
                "to a stranger instead of staying with the person")
            assert _sessions_holding(env, typed) == [], (
                "the typed text was written into a session after a refused "
                "open: it was sent instead of being kept for the person")
            assert tab["input"].value == typed, (
                "the typed text was taken out of the input box")
            said = _chat_text(state)
            assert "Start a new session from it" in said, (
                f"the person is not offered the fork choice; chat says:\n{said}")
            assert "Wait" in said, (
                f"the person is not offered the wait choice; chat says:\n{said}")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(timeout=30)


# ---------------------------------------------------------------------------
# The reuse point, pinned green on main
# ---------------------------------------------------------------------------


def test_the_fork_button_makes_exactly_one_session_from_the_held_history(env):
    """The "start a new session from it" path already exists and already works.

    Target 1 says reuse it, not build a second one. This pins what the builder
    must reuse: pressing Fork on the held conversation writes exactly one new
    session, its history is the held conversation's, and the tab then holds
    the fork. Green on main -- it is a contract, not a finding; if the patch
    breaks it, the patch changed the fork path.

    It also proves the fork works while the SOURCE is held: the new id has no
    lock of its own, which is why the refused state can safely offer it.
    """
    _write_session(env, SID, "held job", "held-turn")
    proc = _held_by_a_live_process(env)
    try:
        before = set(_session_files(env))
        tab = _tab(env)
        state = tab["refs"]["state"]
        try:
            _open_held_conversation(tab, env)
            assert "held by a live session" in _chat_text(state)

            tab["fork"].click()

            after = set(_session_files(env))
            forks = sorted(n for n in after - before
                           if n.startswith("fork-"))
            assert len(forks) == 1, (
                f"the fork must write exactly one session, wrote {forks}")
            forked = forks[0].removesuffix(".json")
            assert forked != SID
            forked_data = json.loads(
                (env["tmp"] / "sessions" / f"{forked}.json").read_text())
            assert forked_data.get("engine_messages") == [
                {"role": "user", "content": "held-turn"},
                {"role": "assistant", "content": "answer to held-turn"},
            ], "the fork must start from the held conversation's history"
            assert state.get("active_session_id") == forked, (
                "after forking, the tab must hold the fork, not stay empty")
        finally:
            tab["refs"]["shutdown"](save=False)
    finally:
        proc.kill()
        proc.wait(timeout=30)
