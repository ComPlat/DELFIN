"""An answer from outside the waiting process carries the code the person saw.

Input: a pending approval in each surface that publishes one -- the
headless ``FileConfirmBroker``, the terminal ``TerminalConfirmBroker``
and the dashboard ``KitConfirmBroker`` -- under a private HOME, each
driven on a worker thread through its real ``callback``. Answers are left
in the room the way any process of the same account can leave them: a
hand-written answer file, ``file_confirm.answer`` / ``answer_waiting`` /
``approval_answers.answer`` with a code, and ``delfin-agent approvals``
run as a child process.

Output: an answer file without the request's one-time code, or with a
wrong one, is not an answer -- it is moved aside, counted and recorded as
a security event; the third such answer expires the request. The right
code approves, and a refusal signed with it cannot be turned into an
approval. The code is shown on the human channel (the terminal frame,
the dashboard row, the headless session's terminal) and is in no file
under the temporary HOME or workspace after the exchange. The in-process
dialogs (a key in the terminal, a click in the dashboard) decide without
any code.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli_gate
from delfin.agent import file_confirm as fc
from delfin.agent import terminal_confirm as tc

from conftest import child_env


MARKER = "DELFIN_AGENT_SESSION"


@pytest.fixture()
def code_mod():
    from delfin.agent import approval_code
    return approval_code


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.delenv("DELFIN_SCRATCH_STATE", raising=False)
    monkeypatch.delenv(MARKER, raising=False)
    monkeypatch.setattr(fc, "requests_dir",
                        lambda: h / ".delfin" / "approvals")
    monkeypatch.setattr(tc, "_PENDING_DIR",
                        h / ".delfin" / "terminal_confirmations")
    monkeypatch.setattr(cli_gate, "state_dir", lambda: h / ".delfin" / "cli_gate")
    return h


@pytest.fixture()
def events(monkeypatch):
    seen = []
    monkeypatch.setattr(
        A, "_record_security_event",
        lambda kind, tool, detail, blocked=True: seen.append(
            (kind, tool, detail, blocked)))
    return seen


class Shown:
    """The headless session's human channel, captured."""

    def __init__(self):
        self.codes = {}

    def __call__(self, request_id, code):
        self.codes[request_id] = code


class Ask:
    """One request on a worker thread, as the engine makes it."""

    def __init__(self, call, *args):
        self.result = None
        self.done = threading.Event()

        def _run():
            try:
                self.result = call(*args)
            finally:
                self.done.set()
        threading.Thread(target=_run, daemon=True).start()


def _wait_for(predicate, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(0.01)
    raise AssertionError("condition not reached")


def _headless(home, shown, timeout_s=8.0):
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=timeout_s,
                                  poll_s=0.02, on_code=shown)
    ask = Ask(broker.callback, "write_file", {"path": "api_client.py"},
              "[SELF-MODIFICATION GUARD] diff")
    rows = _wait_for(lambda: fc.pending())
    return broker, ask, rows[0]["id"], Path(rows[0]["_file"])


def _forge(request_path: Path, body: dict) -> Path:
    """What any process of this account can write: owner, mode, id and age
    all pass the checks that existed before the code."""
    rid = request_path.name[:-len(".request.json")]
    answer = request_path.with_name(f"{rid}.answer.json")
    tmp = request_path.with_name(f".{rid}.forged")
    tmp.write_text(json.dumps(body), encoding="utf-8")
    os.chmod(tmp, 0o600)
    tmp.replace(answer)
    return answer


def _wrong(code: str) -> str:
    return "ZZZZZZ" if code.replace("-", "") != "ZZZZZZ" else "YYYYYY"


# -- the headless broker ----------------------------------------------------

def test_a_forged_answer_without_the_code_is_ignored(home, events):
    shown = Shown()
    broker, ask, rid, req = _headless(home, shown)
    answer = _forge(req, {"id": rid, "decision": "approve", "by": "x",
                          "at": time.time()})
    _wait_for(lambda: not answer.exists())
    assert not ask.done.wait(0.3), "a forged answer decided the request"
    assert any(e[0] == "approval_code_refused" and e[3] for e in events), events
    # The person, with the code, still gets through.
    assert fc.answer(rid, True, code=shown.codes[rid]) is True
    assert ask.done.wait(5) and ask.result is True
    assert broker.last_timed_out is False


def test_the_right_code_approves(home):
    shown = Shown()
    _, ask, rid, _ = _headless(home, shown)
    assert fc.answer(rid, True, code=shown.codes[rid].lower()) is True
    assert ask.done.wait(5) and ask.result is True


def test_a_wrong_code_is_refused_and_recorded(home, events):
    shown = Shown()
    _, ask, rid, req = _headless(home, shown)
    assert fc.answer(rid, True, code=_wrong(shown.codes[rid])) is True
    _wait_for(lambda: not req.with_name(f"{rid}.answer.json").exists())
    assert not ask.done.wait(0.3)
    assert [e for e in events if e[0] == "approval_code_refused"]
    assert all(shown.codes[rid] not in e[2] for e in events)
    assert fc.answer(rid, False, code=shown.codes[rid]) is True
    assert ask.done.wait(5) and ask.result is False


def test_three_wrong_codes_expire_the_request(home, events):
    shown = Shown()
    broker, ask, rid, req = _headless(home, shown, timeout_s=30.0)
    answer = req.with_name(f"{rid}.answer.json")
    started = time.monotonic()
    for _ in range(3):
        _wait_for(lambda: not answer.exists())
        if ask.done.is_set():
            break
        fc.answer(rid, True, code=_wrong(shown.codes[rid]))
    assert ask.done.wait(5), "three wrong codes did not end the request"
    assert ask.result is False
    assert time.monotonic() - started < 10
    refused = [e for e in events if e[0] == "approval_code_refused"]
    assert len(refused) == 3, refused


def test_a_signed_refusal_cannot_become_an_approval(home, code_mod):
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=30,
                                  on_code=Shown())
    rid, req = broker._write_request("write_file", {"path": "x"}, "diff")
    assert fc.answer(rid, False, code=code_mod.code_for(rid),
                     reason="not this one") is True
    answer = req.with_name(f"{rid}.answer.json")
    record = json.loads(answer.read_text())
    assert "proof" in record
    assert code_mod.code_for(rid) not in answer.read_text()
    record["decision"] = "approve"
    record.pop("reason", None)
    _forge(req, record)
    assert fc._read_answer(answer, rid, req.stat().st_mtime) is None


def test_an_answer_needs_a_code_the_answering_process_can_name(home, code_mod):
    """Another process has no registry: no code given, nothing written."""
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=30,
                                  on_code=Shown())
    rid, req = broker._write_request("write_file", {"path": "x"}, "diff")
    code = code_mod.code_for(rid)
    code_mod.forget(rid)
    assert fc.answer(rid, True) is False
    assert not req.with_name(f"{rid}.answer.json").exists()
    assert code


def test_the_request_file_says_a_code_is_needed_but_not_which(home, code_mod):
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=30,
                                  on_code=Shown())
    rid, req = broker._write_request("write_file", {"path": "x"}, "diff")
    record = json.loads(req.read_text())
    assert record.get("needs_code") is True
    assert code_mod.code_for(rid) not in req.read_text()


def test_the_cli_gate_relay_still_answers_its_own_hook(home):
    """The hook broker runs inside the external agent's process tree; it
    has no human channel of its own, so the relay's answer is accepted as
    before (a documented limit, not a code-less hole elsewhere)."""
    broker = fc.FileConfirmBroker(session_id="cli-gate-t", timeout_s=8,
                                  poll_s=0.02, require_code=False)
    ask = Ask(broker.callback, "write_file", {"path": "x"}, "diff")
    rid = _wait_for(lambda: fc.pending())[0]["id"]
    assert fc._write_answer(rid, True, by="delfin-cli-gate") is True
    assert ask.done.wait(5) and ask.result is True


# -- the terminal broker ----------------------------------------------------

def _terminal(home):
    broker = tc.TerminalConfirmBroker(session_id="t1", session_key="t1",
                                      poll_s=0.02)
    ask = Ask(broker.callback, "bash", {"command": "rm -rf build"},
              "rm -rf build")
    req = _wait_for(lambda: broker.take() if broker.take() and
                    broker.take().published else None)
    rid = Path(req.published).name[:-len(".request.json")]
    return broker, ask, req, rid


def test_the_terminal_key_still_decides_without_a_code(home):
    broker, ask, req, _ = _terminal(home)
    assert broker.resolve(req, True) is True
    assert ask.done.wait(5) and ask.result is True


def test_the_terminal_frame_shows_the_code(home, code_mod):
    _, _, req, rid = _terminal(home)
    frame = tc.render_request(req)
    assert code_mod.display(code_mod.code_for(rid)) in frame
    assert rid in frame


def test_a_terminal_question_ignores_a_forged_answer(home, events, code_mod):
    broker, ask, req, rid = _terminal(home)
    answer = _forge(Path(req.published), {"id": rid, "decision": "approve"})
    _wait_for(lambda: not answer.exists())
    assert not ask.done.wait(0.3)
    assert tc.answer_waiting(rid, True, code=code_mod.code_for(rid)) is True
    assert ask.done.wait(5) and ask.result is True


def test_a_choice_from_outside_needs_the_code(home, events, code_mod):
    from delfin.agent import approval_answers as ans
    broker = tc.TerminalConfirmBroker(session_id="t1", poll_s=0.02)
    payload = {"question": "which?", "options": [{"label": "red"},
                                                 {"label": "blue"}]}
    ask = Ask(broker.ask_user, payload)
    req = _wait_for(lambda: broker.take() if broker.take() and
                    broker.take().published else None)
    rid = Path(req.published).name[:-len(".request.json")]
    answer = _forge(Path(req.published), {"id": rid, "decision": "choose",
                                          "answers": ["blue"]})
    _wait_for(lambda: not answer.exists())
    assert not ask.done.wait(0.3)
    assert ans.answer(rid, "1", code=code_mod.code_for(rid)) == ["red"]
    assert ask.done.wait(5) and ask.result == {"answers": ["red"]}


# -- the dashboard broker ---------------------------------------------------

def _dashboard(home):
    from delfin.agent import kit_confirm as kc
    broker = kc.KitConfirmBroker(default_timeout_s=8.0)
    broker._ROOM_POLL_S = 0.02
    ask = Ask(broker.callback, "bash", {"command": "ls"}, "ls")
    _wait_for(lambda: tc.pending_for_supervisors())
    req = broker.pending_requests()[0]
    rid = Path(req.published).name[:-len(".request.json")]
    return broker, ask, req, rid


def test_the_dashboard_click_still_decides_without_a_code(home):
    broker, ask, req, _ = _dashboard(home)
    broker.resolve(req, True)
    assert ask.done.wait(5) and ask.result is True


def test_the_dashboard_room_needs_the_code(home, events, code_mod):
    broker, ask, req, rid = _dashboard(home)
    answer = _forge(Path(req.published), {"id": rid, "decision": "approve"})
    _wait_for(lambda: not answer.exists())
    assert not ask.done.wait(0.3)
    assert tc.answer_waiting(rid, True, code=code_mod.code_for(rid)) is True
    assert ask.done.wait(5) and ask.result is True


def test_the_dashboard_row_shows_the_code(home, code_mod):
    broker, _, req, rid = _dashboard(home)
    assert code_mod.display(code_mod.code_for(rid)) in \
        broker.outside_answer_hint(req)


# -- the command line, as another process ----------------------------------

def _cli(home, tmp_path, *argv):
    env = child_env(tmp_path)
    env["HOME"] = str(home)
    env.pop("DELFIN_SCRATCH_STATE", None)
    env.pop(MARKER, None)
    # The checkout under test, not whichever delfin the interpreter has
    # installed: the child runs from tmp_path.
    root = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in (root, env.get("PYTHONPATH", "")) if p)
    return subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", "approvals", *argv],
        env=env, capture_output=True, text=True, timeout=120, cwd=tmp_path,
        stdin=subprocess.DEVNULL)


def test_the_cli_without_a_code_is_refused(home, tmp_path):
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=30,
                                  on_code=Shown())
    rid, req = broker._write_request("write_file", {"path": "x"}, "diff")
    proc = _cli(home, tmp_path, "approve", rid)
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "code" in proc.stderr.lower()
    assert not req.with_name(f"{rid}.answer.json").exists()


def test_the_cli_with_the_code_writes_a_proof_not_the_code(home, tmp_path,
                                                         code_mod):
    broker = fc.FileConfirmBroker(session_id="s1", timeout_s=30,
                                  on_code=Shown())
    rid, req = broker._write_request("write_file", {"path": "x"}, "diff")
    code = code_mod.code_for(rid)
    proc = _cli(home, tmp_path, "approve", rid, "--code",
                code_mod.display(code))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    answer = req.with_name(f"{rid}.answer.json")
    assert code not in answer.read_text()
    assert fc._read_answer(answer, rid, req.stat().st_mtime) == (True, "")


# -- the code is on the human channel and in no file -----------------------

def test_the_code_reaches_no_file(home, tmp_path, code_mod):
    ws = tmp_path / "ws"
    ws.mkdir()
    shown = Shown()
    codes = []
    # headless: one wrong answer (recorded), then the right one
    broker, ask, rid, req = _headless(home, shown)
    fc.answer(rid, True, code=_wrong(shown.codes[rid]))
    _wait_for(lambda: not req.with_name(f"{rid}.answer.json").exists())
    fc.answer(rid, True, code=shown.codes[rid])
    assert ask.done.wait(5) and ask.result is True
    codes.append(shown.codes[rid])
    # terminal: answered from outside
    tbroker, task, treq, trid = _terminal(home)
    codes.append(code_mod.code_for(trid))
    tc.render_request(treq)
    assert tc.answer_waiting(trid, False, code=codes[-1], reason="no") is True
    assert task.done.wait(5)
    # dashboard: answered from outside
    dbroker, dask, dreq, drid = _dashboard(home)
    codes.append(code_mod.code_for(drid))
    assert tc.answer_waiting(drid, True, code=codes[-1]) is True
    assert dask.done.wait(5)

    forms = set()
    for c in codes:
        assert c
        # Upper case only: every id DELFIN writes is lower-case hex, and a
        # lower-case six-letter form could match one by chance.
        forms.update({c, code_mod.display(c)})
    hits = []
    for root in (home, ws, tmp_path):
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            try:
                text = p.read_bytes().decode("utf-8", "replace")
            except OSError:
                continue
            for f in forms:
                if f in text:
                    hits.append((str(p), f))
    assert not hits, hits
