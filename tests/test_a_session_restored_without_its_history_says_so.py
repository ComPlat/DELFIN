"""A saved session whose model history is empty is opened as itself, and
says what it lacks.

Input: a session file with a chat transcript of real turns and
``engine_messages: []`` -- what the stopped fresh start wrote before it was
fixed, and what every such file on disk still holds. Output: the restored
tab says the model starts without the conversation the chat shows; a model
switch keeps the session's id, so the next save updates that file instead
of copying the transcript into a new one.

Measured before (field report): the reload showed 895 chat rows and the
model received one message, with nothing said. The model switch then
dropped the engine without a carry-over, because ``_drop_engine`` carried
only an engine that had messages -- the session id went with it, the next
engine minted a fresh one, and the first save wrote the transcript under
that new id. The report itself was filed from the copy.
"""

from __future__ import annotations

import json
import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent import scheduler as S                     # noqa: E402

_SID = "f0e1d2c3b4a5f6e7d8c9b0a1f2e3d4c5"


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    from delfin import user_settings
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))


def _save(tmp_path, engine_messages):
    from delfin.agent import session_store as ss
    chat = [{"role": "user", "content": "build the loop", "role_label": ""},
            {"role": "assistant", "content": "the plan is X", "role_label": ""}]
    ss.save_session(_SID, mode="solo", chat_messages=chat,
                    engine_messages=engine_messages, provider="kit",
                    model="kit.deepseek-v4-flash", perm_profile="all_free",
                    workspace=str(tmp_path / "ws"))
    ss.release_session_lock(_SID)


def _open(tmp_path):
    import ipywidgets as widgets
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.initial_session_id = _SID
    ctx.sessions_managed = True
    ctx.presence_key = _SID[:8]
    tab, refs = tab_agent.create_tab(ctx)
    found, stack = {}, [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Dropdown, widgets.Button,
                             widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        stack.extend(getattr(node, "children", ()) or ())
    return refs, found


def _system_lines(state):
    return [str(m.get("content") or "") for m in state["chat_messages"]
            if m.get("role") == "system"]


def test_a_restore_without_the_models_history_is_said(tmp_path):
    _save(tmp_path, engine_messages=[])
    refs, _w = _open(tmp_path)
    try:
        lines = _system_lines(refs["state"])
        assert any("starts without" in line for line in lines), lines[-4:]
    finally:
        refs["shutdown"](save=False)


def test_a_complete_restore_says_nothing_of_the_kind(tmp_path):
    _save(tmp_path, engine_messages=[
        {"role": "user", "content": "build the loop"},
        {"role": "assistant", "content": "the plan is X"}])
    refs, _w = _open(tmp_path)
    try:
        assert not any("starts without" in line
                       for line in _system_lines(refs["state"]))
    finally:
        refs["shutdown"](save=False)


def test_a_model_switch_keeps_the_id_of_a_session_without_history(
        tmp_path, monkeypatch):
    from types import SimpleNamespace

    from delfin.agent import api_client as A
    from delfin.agent import session_store as ss

    def _answer(self, system="", messages=None, **kw):
        yield SimpleNamespace(type="text_delta", text="an answer")

    monkeypatch.setattr(A.OpenAIClient, "stream_message", _answer)
    _save(tmp_path, engine_messages=[])
    refs, w = _open(tmp_path)
    state = refs["state"]
    try:
        w[("Dropdown", "Model:")].value = "kit.glm-5.3"
        w[("Textarea", "")].value = "Hi"
        w[("Button", "Send")].click()
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            if (not state.get("streaming") and any(
                    "an answer" in str(m.get("content"))
                    for m in state["chat_messages"])):
                break
            time.sleep(0.05)
        time.sleep(0.3)
        assert state["engine"].session_id == _SID
        assert state["active_session_id"] == _SID
        saved = sorted(p.stem for p in ss._SESSIONS_DIR.glob("*.json"))
        assert saved == [_SID], saved
    finally:
        refs["shutdown"](save=False)
