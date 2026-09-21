"""Sessions can write to each other.

Several agent sessions work side by side, sometimes in one repository. One
can tell another what it found, ask it to leave a file alone or say it has
pushed. A message waits until the receiving session takes it and always
reads as coming from another session, not from the user.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P
from delfin.agent.api_client import _DocToolExecutor


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def test_a_message_waits_until_it_is_taken_once():
    M.send("B", "I pushed the collector fix", from_key="A", from_title="collector")
    M.send("B", "leave tab_agent.py to me", from_key="A", from_title="collector")
    taken = M.take("B")
    assert [m["text"] for m in taken] == ["I pushed the collector fix",
                                         "leave tab_agent.py to me"]
    assert M.take("B") == []
    M.send("B", "one more", from_key="A")
    assert [m["text"] for m in M.take("B")] == ["one more"]


def test_a_message_never_reads_as_the_user():
    text = M.render({"from": "A", "from_title": "collector", "text": "done"})
    assert "not from the user" in text and 'to="A"' in text and text.endswith("done")


def _call(perms, **arguments):
    return json.loads(_DocToolExecutor()._execute_session_message(arguments, perms))


def test_the_tool_lists_the_other_sessions_and_writes_to_one(tmp_path):
    P.announce("A", title="collector", workspace=str(tmp_path))
    P.announce("B", title="docs", workspace=str(tmp_path))
    perms = SimpleNamespace(presence_key="A")

    listed = _call(perms)
    assert [s["key"] for s in listed["sessions"]] == ["B"]

    sent = _call(perms, to="B", message="the docs build is green")
    assert sent["status"] == "sent"
    (message,) = M.take("B")
    assert (message["from"], message["from_title"]) == ("A", "collector")

    assert "error" in _call(perms, to="nobody", message="x")
    assert "error" in _call(perms, to="B")
    assert "error" in _call(perms, to="A", message="to myself")


def test_an_idle_session_with_a_draft_keeps_its_messages_waiting(tmp_path):
    pytest.importorskip("ipywidgets")
    tab, refs = _tab(tmp_path)
    M.send("B", "hold off", from_key="A")
    _input(tab).value = "the user's own draft"
    refs["deliver_messages"]()
    assert [m["text"] for m in M.take("B")] == ["hold off"]
    refs["shutdown"]()


def test_a_running_turn_is_handed_the_message_between_rounds(tmp_path):
    pytest.importorskip("ipywidgets")
    tab, refs = _tab(tmp_path)
    steered = []
    refs["state"]["streaming"] = True
    refs["state"]["engine"] = SimpleNamespace(
        client=SimpleNamespace(push_steer=steered.append))
    M.send("B", "I pushed; rebase before you commit", from_key="A",
           from_title="collector")

    refs["deliver_messages"]()

    (text,) = steered
    assert "not from the user" in text and "rebase before you commit" in text
    assert "Message from collector" in str(refs["state"]["chat_messages"])
    refs["state"]["streaming"] = False
    refs["state"]["engine"] = None
    refs["shutdown"]()


def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.presence_key = "B"
    return tab_agent.create_tab(ctx)


def _input(tab):
    import ipywidgets as widgets
    stack = [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, widgets.Textarea) and "Message the agent" in (
                node.placeholder or ""):
            return node
        stack.extend(getattr(node, "children", ()) or ())
    raise AssertionError("no input box")
