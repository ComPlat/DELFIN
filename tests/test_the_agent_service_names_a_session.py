"""What a session is called and where it is reached are two things.

Its title (`/rename`, the agent view's ctrl+r) is what every window shows
and what the session store saves it under; its address -- the `-n` name it
started with, or its id -- is what mail, approvals and its presence record
reach it by. A rename changes the title only. The title is saved at every
turn's end (a crash keeps it), and a session resumed without a name comes
back under the one it was given -- never under the store's stand-in, the
start of its first message.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent.service import protocol
from tui_worker_client import WorkerClient


def _saved(w: WorkerClient, sid: str) -> dict | None:
    for root in (w.home, w.scratch):
        for path in root.rglob(f"{sid}.json"):
            try:
                return json.loads(path.read_text())
            except ValueError:
                continue
    return None


def _presence(w: WorkerClient) -> list[dict]:
    out = []
    for path in [*(w.home / ".delfin" / "session_presence").glob("*.json"),
                 *(w.scratch / ".delfin" / "session_presence").glob("*.json")]:
        try:
            out.append(json.loads(path.read_text()))
        except ValueError:
            continue
    return out


def _until(cond, timeout: float = 15.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


@pytest.mark.parametrize("text, kept", [
    ("  parser   fix ", "parser fix"),
    ("a\x1b[31mred\x07​name", "a red name"),
    ("\x1b]0;a title of its own\x07ok", "ok"),
    ("two\nlines here", "two lines here"),
    ("x" * 120, "x" * protocol.TITLE_MAX),
    ("\n\t \x1b", ""),
])
def test_a_title_is_one_plain_line_of_at_most_eighty_characters(text, kept):
    assert protocol.clean_title(text) == kept


def test_a_rename_changes_what_it_is_called_not_where_it_is_reached(tmp_path):
    with WorkerClient(tmp_path, "-n", "scan-a") as w:
        sid = w.init_result["session"]["sessionId"]
        assert w.init_result["session"]["title"] == "scan-a"
        mark = len(w.events)
        assert w.request("session.setTitle", {"title": "Parser\x1b[1m fix"}) == {"title": "Parser fix"}
        w.wait_for(lambda e: e.get("type") == "session.updated" and e.get("title") == "Parser fix", after=mark)
        mine = _until(lambda: [r for r in _presence(w) if r.get("session_id") == sid
                               and r.get("title") == "Parser fix"])
        # One record, under the address it started with, called what it is now.
        assert [r["key"] for r in mine] == ["scan-a"]
        assert [r["key"] for r in _presence(w) if r.get("session_id") == sid] == ["scan-a"]
        with pytest.raises(Exception):
            w.request("session.setTitle", {"title": " \x1b[0m "})


def test_the_title_is_saved_at_each_turn_and_a_resume_without_a_name_keeps_it(tmp_path):
    with WorkerClient(tmp_path, "--fake-save") as w:
        sid = w.init_result["session"]["sessionId"]
        assert w.init_result["session"].get("title") in (None, "")
        w.request("session.setTitle", {"title": "solvent scan"})
        w.turn("hello there")
        # Saved at the turn's end, before the session ends: a crash keeps it.
        saved = _until(lambda: (_saved(w, sid) or {}).get("title") == "solvent scan" and _saved(w, sid))
        assert saved and saved["title"] == "solvent scan"
    with WorkerClient(tmp_path, "--fake-save", "--resume", sid) as again:
        assert again.init_result["session"]["sessionId"] == sid
        assert again.init_result["session"]["title"] == "solvent scan"


def test_a_session_never_named_is_not_named_after_its_first_message_on_resume(tmp_path):
    with WorkerClient(tmp_path, "--fake-save") as w:
        sid = w.init_result["session"]["sessionId"]
        w.turn("hello there")
        assert _until(lambda: _saved(w, sid))
    # The store calls it "hello there" (its stand-in); that is no name.
    assert (_saved(w, sid) or {}).get("title") == "hello there"
    with WorkerClient(tmp_path, "--fake-save", "--resume", sid) as again:
        assert again.init_result["session"].get("title") in (None, "")


def test_a_title_given_at_the_start_is_kept_apart_from_the_address(tmp_path):
    with WorkerClient(tmp_path, "-n", "scan-b", "--title", "the B scan") as w:
        sid = w.init_result["session"]["sessionId"]
        assert w.init_result["session"]["title"] == "the B scan"
        w.turn("hello")
        mine = _until(lambda: [r for r in _presence(w) if r.get("session_id") == sid])
        assert [r["key"] for r in mine] == ["scan-b"]
        assert mine[0]["title"] == "the B scan"


def test_stored_title_tells_a_given_name_from_the_stand_in(tmp_path, monkeypatch):
    from delfin.agent import session_store
    from delfin.agent.service import worker

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", Path(tmp_path))
    long_prompt = "please read every log in this folder and " * 4

    def save(sid, title, first):
        (Path(tmp_path) / f"{sid}.json").write_text(json.dumps(
            {"session_id": sid, "title": title,
             "chat_messages": [{"role": "user", "content": first}, {"role": "assistant", "content": "ok"}]}))

    save("a", "hello there", "hello there")
    save("b", long_prompt[:80].strip() + "...", long_prompt)
    save("c", "solvent scan", "hello there")
    save("d", "", "hello there")
    assert worker.stored_title("a") == ""
    assert worker.stored_title("b") == ""
    assert worker.stored_title("c") == "solvent scan"
    assert worker.stored_title("d") == ""
    assert worker.stored_title("missing") == ""


def test_a_rename_of_an_idle_session_is_saved_at_once(tmp_path):
    # Not only at the next turn's end: a session killed before it keeps the name.
    with WorkerClient(tmp_path, "--fake-save") as w:
        sid = w.init_result["session"]["sessionId"]
        w.turn("hello there")
        w.request("session.setTitle", {"title": "renamed while idle"})
        saved = _until(lambda: (_saved(w, sid) or {}).get("title") == "renamed while idle")
        assert saved, _saved(w, sid)


def test_the_title_a_session_starts_with(monkeypatch):
    import argparse
    from delfin.agent.service import worker

    monkeypatch.setattr(worker, "stored_title", lambda sid: f"saved name of {sid}")
    engine = type("Engine", (), {"session_id": "s9"})()

    def args(**over):
        base = dict(title="", name="", resume="", continue_session=False, fork_session=False)
        return argparse.Namespace(**{**base, **over})

    assert worker.initial_title(args(title="given"), engine) == "given"
    assert worker.initial_title(args(name="scan-a"), engine) == "scan-a"
    assert worker.initial_title(args(), engine) == ""                       # new: no name
    for how in ({"resume": "s9"}, {"continue_session": True}, {"fork_session": True, "resume": "s9"}):
        assert worker.initial_title(args(**how), engine) == "saved name of s9", how


def test_rename_session_changes_the_title_and_nothing_else(tmp_path, monkeypatch):
    from delfin.agent import session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", Path(tmp_path))
    record = {"session_id": "s1", "title": "old", "chat_messages": [{"role": "user", "content": "hi"}],
              "cost_usd": 0.25, "workspace": "/w"}
    (Path(tmp_path) / "s1.json").write_text(json.dumps(record))
    assert session_store.rename_session("s1", "new name") is True
    after = json.loads((Path(tmp_path) / "s1.json").read_text())
    assert after == {**record, "title": "new name"}
    assert session_store.rename_session("missing", "x") is False
    assert session_store.rename_session("", "x") is False


def test_a_sleeping_session_renamed_is_saved_under_the_name_and_keeps_its_address(tmp_path):
    from tui_service import ServiceProcess

    svc = ServiceProcess(tmp_path, start=False)
    svc.env["DELFIN_AGENT_SLEEP_AFTER_S"] = "1"
    svc.start()
    try:
        a = svc.client()
        sid = svc.session(a, attach=False, name="first name")
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and not any(
                s["sessionId"] == sid and s["state"] == "sleeping" for s in a.request("sessions.list")["sessions"]):
            time.sleep(0.1)
        # A saved record of it, as a real session leaves one.
        store = svc.scratch / ".delfin" / "agent_sessions"
        store.mkdir(parents=True, exist_ok=True)
        (store / f"{sid}.json").write_text(json.dumps({"session_id": sid, "title": "first name",
                                                       "chat_messages": []}))
        a.request("session.setTitle", {"sessionId": sid, "title": "second name"})
        assert json.loads((store / f"{sid}.json").read_text())["title"] == "second name"
        roster = json.loads(next((svc.scratch / ".delfin" / "agent_service" / "run").glob("*.roster")).read_text())
        settings = next(r for r in roster["dormant"] if r["summary"]["sessionId"] == sid)["settings"]
        # Called by the new name; reached by the old address.
        assert settings["title"] == "second name" and settings["name"] == "first name"
        # Woken: still both.
        a.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
        mark = len(a.events)
        a.request("session.prompt", {"sessionId": sid, "text": "hello"})
        a.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark, timeout=60)
        presence = svc.scratch / ".delfin" / "session_presence"
        mine = _until(lambda: [json.loads(f.read_text()) for f in presence.glob("*.json")
                               if json.loads(f.read_text()).get("session_id") == sid])
        assert [(r["key"], r["title"]) for r in mine] == [("first name", "second name")]
        assert next(x for x in a.request("sessions.list")["sessions"] if x["sessionId"] == sid)["title"] == "second name"
    finally:
        svc.stop()
