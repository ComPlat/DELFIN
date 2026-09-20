"""A session with nobody at its terminal can still ask somebody.

The terminal broker is bound only when ``stdin`` is a tty -- the comment
at the binding says so, and says why: without a callback a write inside
the workspace is allowed silently, "which is right for an unattended
run". But the one thing that always needs a human answer, the
Self-Modification Guard, then comes back as

    'write_file' targets the agent's own safety layer ('...') — this
    requires explicit user confirmation but no confirm_callback is
    configured, so it is refused.

Honest, and it means a headless session cannot touch the agent's own
safety layer at all. The most valuable work there is exactly what it
cannot do.

This broker gives such a session somebody to ask, over the filesystem:
the question is written to a directory, whoever is watching writes the
answer beside it. It is a security boundary, so every case below is about
REFUSING rather than guessing. What a supervisor may approve is what the
supervisor actually saw and actually answered -- nothing else.
"""

from __future__ import annotations

import json
import os
import threading
import time

import pytest

from delfin.agent import file_confirm as FC


@pytest.fixture()
def room(tmp_path, monkeypatch):
    place = tmp_path / "approvals"
    monkeypatch.setattr(FC, "requests_dir", lambda: place)
    return place


@pytest.fixture()
def broker(room):
    return FC.FileConfirmBroker(session_id="s1", timeout_s=3.0, poll_s=0.02)


def _answer_when_asked(room, decision, *, delay=0.05, mangle=None):
    """A supervisor in another thread."""
    def _run():
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            waiting = FC.pending(room)
            if waiting:
                rid = waiting[0]["id"]
                if mangle is None:
                    FC.answer(rid, decision, room=room, by="test")
                else:
                    mangle(room, rid)
                return
            time.sleep(0.02)
    t = threading.Thread(target=_run, daemon=True)
    time.sleep(delay) if False else None
    t.start()
    return t


# -- the contract the gate depends on ---------------------------------------

def test_the_callback_is_bound_so_the_gate_can_read_the_flag(broker):
    """The gate does perms.confirm_callback.__self__.last_timed_out."""
    assert getattr(broker.callback, "__self__", None) is broker


def test_the_expiry_flag_is_per_thread(broker):
    """Requests arrive on whichever thread runs tools and they overlap.
    One flag for all of them records an expiry as a refusal."""
    broker.last_timed_out = True
    seen = {}

    def _fresh():
        seen["value"] = broker.last_timed_out
    t = threading.Thread(target=_fresh, daemon=True)
    t.start()
    t.join(5)
    assert seen["value"] is False
    assert broker.last_timed_out is True


# -- asking and being answered ----------------------------------------------

def test_an_approval_comes_back_true(broker, room):
    t = _answer_when_asked(room, True)
    assert broker.callback("write_file", {"path": "api_client.py"}, "diff") is True
    t.join(5)


def test_a_refusal_comes_back_false_and_is_not_an_expiry(broker, room):
    t = _answer_when_asked(room, False)
    assert broker.callback("write_file", {"path": "x"}, "diff") is False
    assert broker.last_timed_out is False, (
        "a refusal is a decision; only absence sets the flag")
    t.join(5)


def test_the_question_carries_what_a_supervisor_needs(broker, room):
    seen = {}

    def _look(place, rid):
        seen.update(FC.pending(place)[0])
        FC.answer(rid, True, room=place)
    t = _answer_when_asked(room, True, mangle=_look)
    broker.callback("write_file", {"path": "delfin/agent/engine.py"},
                    "[SELF-MODIFICATION GUARD]\n--- a\n+++ b\n")
    t.join(5)
    assert seen["tool"] == "write_file"
    assert seen["path"] == "delfin/agent/engine.py"
    assert "SELF-MODIFICATION" in seen["preview"]
    assert seen["pid"] == os.getpid()
    assert seen["session_id"] == "s1"


# -- nobody there -----------------------------------------------------------

def test_silence_is_an_expiry_and_not_a_refusal(room):
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.3, poll_s=0.02)
    assert broker.callback("write_file", {"path": "x"}, "d") is False
    assert broker.last_timed_out is True, (
        "recorded as a refusal, a path the user never saw closes for the "
        "rest of the session")


def test_a_finished_exchange_leaves_the_pending_list(broker, room):
    t = _answer_when_asked(room, True)
    broker.callback("write_file", {"path": "x"}, "d")
    t.join(5)
    assert FC.pending(room) == []
    kept = list((room / "answered").glob("*.json"))
    assert kept, "what was approved, and by whom, is the record that matters"


# -- the refusing part ------------------------------------------------------

def _plant_answer(room, rid, payload, *, mode=0o600, mtime=None):
    """Put an answer in place in one step, the way the broker writes one.

    Built under a temporary name and moved into place. Written directly,
    the file existed for a moment with the CURRENT mtime and the default
    mode before chmod and utime reached it -- and a poller reading in
    that window saw a fresh, well-permissioned answer and took it. That
    is the whole of `test_an_answer_older_than_the_question_is_ignored`
    failing on a CI runner roughly one run in three while passing every
    time locally: the race was in this helper, not in the code it drives.

    Proved by widening the window: a 150 ms sleep between the write and
    the utime makes the old form return APPROVE three times out of three.
    """
    room.mkdir(parents=True, exist_ok=True)
    path = room / f"{rid}.answer.json"
    tmp = room / f".{rid}.answer.partial"
    tmp.write_text(json.dumps(payload), encoding="utf-8")
    os.chmod(tmp, mode)
    if mtime is not None:
        os.utime(tmp, (mtime, mtime))
    tmp.replace(path)
    return path


def test_an_answer_naming_another_request_is_not_an_answer(room):
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.4, poll_s=0.02)

    def _wrong(place, rid):
        _plant_answer(place, rid, {"id": "somebody-else",
                                   "decision": FC.APPROVE})
    t = _answer_when_asked(room, True, mangle=_wrong)
    assert broker.callback("write_file", {"path": "x"}, "d") is False
    assert broker.last_timed_out is True, "it expired; it was not refused"
    t.join(5)


def test_an_unrecognised_decision_is_not_an_approval(room):
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.4, poll_s=0.02)

    def _vague(place, rid):
        _plant_answer(place, rid, {"id": rid, "decision": True})
    t = _answer_when_asked(room, True, mangle=_vague)
    assert broker.callback("write_file", {"path": "x"}, "d") is False
    t.join(5)


def test_a_group_writable_answer_is_ignored(room):
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.4, poll_s=0.02)

    def _loose(place, rid):
        _plant_answer(place, rid, {"id": rid, "decision": FC.APPROVE},
                      mode=0o666)
    t = _answer_when_asked(room, True, mangle=_loose)
    assert broker.callback("write_file", {"path": "x"}, "d") is False, (
        "anyone on this machine could have written that file")
    t.join(5)


def test_an_answer_older_than_the_question_is_ignored(room):
    """A file placed in advance must not approve whatever comes next."""
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.4, poll_s=0.02)

    def _prewritten(place, rid):
        _plant_answer(place, rid, {"id": rid, "decision": FC.APPROVE},
                      mtime=time.time() - 3600)
    t = _answer_when_asked(room, True, mangle=_prewritten)
    assert broker.callback("write_file", {"path": "x"}, "d") is False
    t.join(5)


def test_a_symlinked_answer_is_ignored(room, tmp_path):
    broker = FC.FileConfirmBroker(session_id="s1", timeout_s=0.4, poll_s=0.02)
    target = tmp_path / "elsewhere.json"

    def _link(place, rid):
        target.write_text(json.dumps({"id": rid, "decision": FC.APPROVE}),
                          encoding="utf-8")
        (place / f"{rid}.answer.json").symlink_to(target)
    t = _answer_when_asked(room, True, mangle=_link)
    assert broker.callback("write_file", {"path": "x"}, "d") is False, (
        "the target of a link is not what the directory listing showed")
    t.join(5)


def test_the_directory_is_private(broker, room):
    t = _answer_when_asked(room, True)
    broker.callback("write_file", {"path": "x"}, "d")
    t.join(5)
    assert (os.stat(room).st_mode & 0o077) == 0, oct(os.stat(room).st_mode)


def test_answering_a_question_that_was_never_asked_fails(room):
    room.mkdir(parents=True, exist_ok=True)
    assert FC.answer("no-such-request", True, room=room) is False


def test_the_guard_actually_asks_through_it(room, tmp_path):
    """Driven through the executor, not through the broker alone.

    A broker nobody calls protects nothing, and the gate reaches the
    callback by a path of its own. Without one the same write comes back
    "no confirm_callback is configured, so it is refused" -- which is the
    state a headless run is in today.
    """
    from delfin.agent import api_client as A

    ws = tmp_path / "ws"
    (ws / "delfin" / "agent").mkdir(parents=True)
    (ws / "delfin" / "agent" / "api_client.py").write_text(
        "x = 1\n", encoding="utf-8")
    args = {"path": "delfin/agent/api_client.py", "content": "y = 2\n"}

    # 1. Nobody to ask.
    bare = A.KitToolPermissions(workspace=ws, mode="bypassPermissions")
    assert "no confirm_callback is configured" in str(
        A._doc_executor.execute("write_file", dict(args), bare))

    # 2. Somebody to ask, who says no.
    seen = {}

    def _supervisor():
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            waiting = FC.pending(room)
            if waiting:
                seen.update(waiting[0])
                FC.answer(waiting[0]["id"], False, room=room, by="test")
                return
            time.sleep(0.02)

    broker = FC.FileConfirmBroker(session_id="e2e", timeout_s=5, poll_s=0.02)
    perms = A.KitToolPermissions(workspace=ws, mode="bypassPermissions")
    perms.confirm_callback = broker.callback
    t = threading.Thread(target=_supervisor, daemon=True)
    t.start()
    out = str(A._doc_executor.execute("write_file", dict(args), perms))
    t.join(5)

    assert seen.get("path") == "delfin/agent/api_client.py"
    assert "SELF-MODIFICATION GUARD" in seen.get("preview", ""), (
        "the supervisor must see WHAT it is being asked to allow")
    assert "denied" in out
    assert (ws / "delfin" / "agent" / "api_client.py").read_text(
        encoding="utf-8") == "x = 1\n", "a refusal must not write"


def test_the_state_path_is_redirectable():
    """A test that approved something must not leave the answer in the
    user's real directory, where a later session would read it."""
    from delfin.agent import state_paths
    pairs = {(m, a) for m, a, _ in state_paths.USER_STATE_RESOLVERS}
    assert ("delfin.agent.file_confirm", "requests_dir") in pairs
