"""A session opened in the window from the session store shows what was
said in it, not only its banner.

The window's replay comes from the agent service's event log, which holds
only what ran in the service. A conversation that ran in the dashboard or
the REPL, and every copy (a new id), has nothing there: the worker shows
the restored conversation as it opens -- a notice, then its last turns as
the window's own events (what the user wrote, the tools called and what
came back, the answers) -- once: they are durable, a later attach
replays them, and a session whose log has its turns shows them from there.
Rows the harness wrote are shown as the window shows them live (◇); the
brief a session about an attention item began with is left out of the
prompt it went in front of. A saved row it cannot read leaves the picture
incomplete, never the session unopenable.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.harness_rows import HARNESS_KEY
from delfin.agent.service.worker import BRIEF_END, BRIEF_HEADER, history_events
from tui_service import ServiceProcess, turn


def test_a_saved_conversation_becomes_the_windows_events():
    messages = [
        {"role": "system", "content": "you are DELFIN"},
        {"role": "user", "content": "optimise benzene"},
        {"role": "assistant", "content": "", "tool_calls": [
            {"id": "c1", "type": "function",
             "function": {"name": "bash", "arguments": json.dumps({"command": "orca opt.inp" + "x" * 900})}}]},
        {"role": "tool", "tool_call_id": "c1", "content": "TOTAL ENERGY -232.1\n" * 400},
        {"role": "assistant", "content": [{"type": "text", "text": "The geometry converged."}]},
        {"role": "user", "content": "[job 4711 finished]", HARNESS_KEY: True},
        {"role": "assistant", "content": "Noted."},
        {"role": "user", "content": f"{BRIEF_HEADER}Job 4711 failed.{BRIEF_END}What went wrong?"},
        {"role": "assistant", "content": "It ran out of memory."},
    ]
    events = history_events(messages)
    kinds = [k for k, _ in events]
    assert kinds == ["notice", "user.message", "tool.started", "tool.ended", "text.ended",
                     "user.message", "text.ended", "user.message", "text.ended"]
    assert events[0][1]["text"] == "earlier in this session, as it was saved"
    fields = [f for _, f in events]
    assert fields[1] == {"text": "optimise benzene", "source": "user"}
    assert fields[2]["name"] == "bash" and fields[2]["callId"] == "c1"
    assert len(fields[2]["input"]["command"]) == 400                  # a long input, cut
    assert len(fields[3]["output"]) == 2000                           # a long output, cut
    assert fields[4]["text"] == "The geometry converged."
    # The harness's row, as the window shows it live: not the user's words.
    assert fields[5] == {"text": "[job 4711 finished]", "source": "harness"}
    assert fields[6]["text"] == "Noted."
    assert fields[7] == {"text": "What went wrong?", "source": "user"}  # without the brief
    assert len({f["partId"] for f in fields if "partId" in f}) == 3   # each part its own


def test_a_brief_with_a_rule_in_it_and_rows_it_cannot_read():
    brief = "Job 4711 failed.\n\n---\n\nDetails: out of memory\n\n---\n\nmore"
    events = history_events([
        {"role": "user", "content": f"{BRIEF_HEADER}{brief}{BRIEF_END}What went wrong?"},
        {"role": "assistant", "content": "Memory.", "tool_calls": "not a list"},
        {"role": "assistant", "content": None, "tool_calls": [{"id": "c1", "function": "not a dict"},
                                                              "not a call", {"function": {"arguments": "{bad"}}]},
        {"role": "tool", "tool_call_id": "c1", "content": 42},
        {"role": "assistant", "content": [None, {"type": "text"}, {"type": "image_url"}]},
    ])
    fields = [f for _, f in events]
    assert fields[1] == {"text": "What went wrong?", "source": "user"}
    assert [k for k, _ in events].count("tool.started") == 2           # the two calls it could read
    assert any(k == "tool.ended" and f["output"] == "" for k, f in events)


def test_only_the_last_turns_and_how_many_are_left_out():
    messages = []
    for n in range(25):
        messages += [{"role": "user", "content": f"question {n}"}, {"role": "assistant", "content": f"answer {n}"}]
    events = history_events(messages)
    said = [f["text"] for k, f in events if k == "user.message"]
    assert said == [f"question {n}" for n in range(5, 25)]
    assert events[0][1]["text"] == "earlier in this session, as it was saved (5 earlier turns not shown)"
    assert history_events([]) == []
    assert history_events([{"role": "assistant", "content": "nobody asked"}]) == []
    assert history_events("not a list") == []


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _saved_session(service) -> tuple[str, str]:
    """A session that ran, was saved and ended: its id and workspace."""
    c = service.client()
    ws = str(service.workspace())
    sid = c.request("session.create", {"cwd": ws, "fake": True, "fakeDelay": 0.003, "fakeSave": True,
                                       "name": "saved-one"}, timeout=60)["sessionId"]
    turn(c, sid, "hello there")
    c.request("session.shutdown", {"sessionId": sid}, timeout=60)
    c.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid, timeout=60)
    return sid, ws


def test_a_copy_shows_the_conversation_it_was_copied_from(service):
    sid, ws = _saved_session(service)
    c = service.client()
    copy = c.request("session.create", {"cwd": ws, "resume": sid, "fork": True, "fake": True,
                                        "fakeDelay": 0.003, "attach": False}, timeout=60)["sessionId"]
    assert copy != sid
    attached = c.request("session.attach", {"sessionId": copy, "turns": 20})
    assert attached["session"]["title"] == "saved-one (fork)"           # the store's own copy
    replay = attached["events"]
    said = [(e["type"], e.get("text", "")) for e in replay
            if e.get("type") in ("notice", "user.message", "text.ended")]
    assert said[0] == ("notice", "earlier in this session, as it was saved")
    assert ("user.message", "hello there") in said
    assert any(t == "text.ended" and "(re: hello there)" in text for t, text in said)
    # Durable: attached again, the same rows once more -- not twice over.
    again = c.request("session.attach", {"sessionId": copy, "turns": 20})["events"]
    assert [e for e in again if e.get("type") == "user.message"] == \
        [e for e in replay if e.get("type") == "user.message"]


def test_a_resume_whose_log_has_the_turns_shows_them_once(service):
    sid, ws = _saved_session(service)
    c = service.client()
    resumed = c.request("session.create", {"cwd": ws, "resume": sid, "fake": True, "fakeDelay": 0.003,
                                           "attach": False}, timeout=60)["sessionId"]
    assert resumed == sid
    replay = c.request("session.attach", {"sessionId": sid, "turns": 20})["events"]
    said = [e["text"] for e in replay if e.get("type") == "user.message"]
    assert said == ["hello there"]
    assert not any(e.get("type") == "notice" and "as it was saved" in str(e.get("text"))
                   for e in replay)
