"""What needs the user, from every DELFIN process, in the window's list.

The agent service reads the attention inbox -- questions a session parked,
confirmations, failed and finished runs, budget and security notices,
written by dashboards, the REPL, the job monitor and the scheduler -- and
the job monitor's findings, and sends the open items to every window:

* NEEDS YOU first (what blocks someone, then what failed), NOTICES after,
  the newest first in each; expired notices and closed entries left out;
* a failed job's entry carries its finding: the state, the diagnosis
  session, the prepared fix; a finding whose entry is closed is not shown;
* an entry of a session this service runs is left out (it says so itself);
* a session it belongs to is named, said to be open elsewhere or saved;
* a change to the records reaches every window within a few seconds;
* answering a parked question and dismissing an item go through the
  inbox's own helpers; a watch-only connection may do neither;
* a session started with a brief shows it as its first row and hands it
  to the model in front of the first prompt, not in what the window shows.

On the way: every inbox entry says which host wrote it, and the watchpost's
alerts arrive (their kind was refused, and the refusal swallowed).
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import state_paths
from delfin.agent.service import attention_feed, protocol
from delfin.agent.service.client import ServiceError
from tui_service import REPO, ServiceProcess, turn


def _entry(kind: str, title: str, *, detail: str = "", session_id: str = "", workspace: str = "",
           age: float = 60.0, status: str = "pending", host: str = "node1", options=None,
           now: float | None = None) -> dict:
    now = time.time() if now is None else now
    return {"id": f"att-{kind}-{abs(hash((title, age))) % 10**8:08d}", "kind": kind,
            "session_id": session_id, "title": title, "detail": detail,
            "options": options, "workspace": workspace, "host": host,
            "created_at": now - age, "status": status, "answer": None,
            "resolved_at": None, "acknowledged": False, "delivered": []}


def _write(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(r) + "\n" for r in records))


def _finding(job_id: str, folder: str, *, state: str = "OUT_OF_MEMORY", diagnosis: str = "",
             fix: dict | None = None, age: float = 120.0) -> dict:
    return {"job_id": job_id, "folder": folder, "state": state, "signatures": ["oom-kill"],
            "ts": time.time() - age, "diagnosis_session": diagnosis,
            "summary": "the SCF ran out of memory", "fix": fix or {}}


@pytest.fixture
def home(tmp_path):
    """This process's DELFIN state, redirected to a fresh directory (the
    suite's own redirect: nothing reaches the real ~/.delfin)."""
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


def _inbox(home: Path) -> Path:
    from delfin.agent import attention
    return attention._inbox_path()


def _findings(home: Path) -> Path:
    from delfin.agent import job_monitor
    return Path(job_monitor._FINDINGS_PATH)


# -- the items, from the records -------------------------------------------------------------

def test_the_items_come_in_the_agent_views_order(home):
    _write(_inbox(home), [
        _entry("run_finished", "Job 4702 finished", detail="watched job — COMPLETED", age=30),
        _entry("question_pending", "Question waiting: which functional?",
               detail="Which functional for the excited states?", options=["B3LYP", "CAM-B3LYP"],
               workspace="/calc/benzene", age=7200),
        _entry("run_failed", "Job 4711 failed (OUT_OF_MEMORY)", detail="out of memory",
               workspace="/calc/benzene", age=600),
        _entry("budget_warning", "Run budget exhausted — the agent is wrapping up", age=10),
        _entry("confirm_pending", "Confirmation required: bash", detail="$ rm -r build", age=5),
        _entry("question_pending", "Question waiting: answered already", status="resolved", age=1),
    ])
    items = attention_feed.snapshot()
    assert [(i["kind"], i["group"]) for i in items] == [
        ("confirm", "needs"), ("question", "needs"),         # what blocks, newest first
        ("failed", "needs"),                                 # then what failed
        ("budget", "notices"), ("finished", "notices")]      # then notices, newest first
    question = items[1]
    assert question["blocking"] is True
    assert question["text"] == "question: Which functional for the excited states?"
    assert question["options"] == ["B3LYP", "CAM-B3LYP"]
    assert question["workspace"] == "/calc/benzene" and question["host"] == "node1"
    assert items[0]["text"] == "confirm: bash"
    assert items[4]["name"] == "job 4702 finished"


def test_a_failed_jobs_entry_carries_its_finding(home):
    _write(_inbox(home), [_entry("run_failed", "Job 4711 failed (OUT_OF_MEMORY)",
                                 detail="out of memory", workspace="/calc/benzene")])
    _write(_findings(home), [
        # job_fix's assessment, as the monitor saves it: the prepared fix under "fix".
        _finding("4711", "/calc/benzene", diagnosis="diag-4711",
                 fix={"fix_class": "resource", "recommendation": "raise maxcore",
                      "fix": {"kind": "memory", "summary_line": "maxcore 3000 -> 4500",
                              "diff": "-%maxcore 3000\n+%maxcore 4500"}}),
        # A finding whose entry is closed (dismissed, expired): not shown.
        _finding("4690", "/calc/water"),
    ])
    items = attention_feed.snapshot()
    assert len(items) == 1
    job = items[0]
    assert job["name"] == "job 4711 failed"
    assert job["text"] == "OUT_OF_MEMORY · diagnosis ready · a fix is prepared"
    assert job["job"]["diagnosisSession"] == "diag-4711"
    assert job["job"]["fix"]["fix"]["summary_line"] == "maxcore 3000 -> 4500"
    assert job["job"]["folder"] == "/calc/benzene"


def test_an_assessment_without_a_prepared_fix_says_none_is_prepared(home):
    _write(_inbox(home), [_entry("run_failed", "Job 4712 failed (FAILED)", workspace="/calc/x")])
    _write(_findings(home), [_finding("4712", "/calc/x", state="FAILED",
                                      fix={"fix_class": "input", "recommendation": "check the input",
                                           "fix": None})])
    job = attention_feed.snapshot()[0]
    assert job["text"] == "FAILED"
    assert job["job"]["fix"]["recommendation"] == "check the input"


_EMITTER = """
import pathlib, sys
sys.path.insert(0, sys.argv[4])
from delfin.agent import attention
inbox = pathlib.Path(sys.argv[1])
attention._inbox_path = lambda: inbox          # the test's own inbox, whatever HOME says
attention._deliver = lambda event, **kw: ["desktop"]
for i in range(int(sys.argv[3])):
    attention.emit_attention("question_pending", title=f"Question waiting: {sys.argv[2]}-{i}",
                             detail=f"{sys.argv[2]}-{i}?", workspace="/calc/x")
"""


def test_the_list_holds_while_several_processes_write(home, tmp_path):
    """Dashboards, the REPL and the job monitor write the inbox at once:
    every read in between is a whole list (never an error, never a torn
    entry), and the last one has everything -- a torn last line aside."""
    from conftest import child_env
    procs, writes = 4, 25
    _inbox(home).parent.mkdir(parents=True, exist_ok=True)
    children = [subprocess.Popen([sys.executable, "-c", _EMITTER, str(_inbox(home)), f"p{k}", str(writes), str(REPO)],
                                 env=child_env(tmp_path)) for k in range(procs)]
    reads, most = 0, 0
    while any(c.poll() is None for c in children):
        items = attention_feed.snapshot()
        assert all(i["kind"] == "question" and i["text"].startswith("question: p") for i in items)
        most = max(most, len(items))
        reads += 1
    assert all(c.wait(timeout=120) == 0 for c in children)
    assert reads > 0
    with open(_inbox(home), "a", encoding="utf-8") as f:
        f.write('{"id": "att-torn", "kind": "question_pending", "title": "Quest')   # a writer died here
    items = attention_feed.snapshot()
    assert len(items) == procs * writes
    assert len({i["id"] for i in items}) == procs * writes and "att-torn" not in {i["id"] for i in items}


def test_expired_notices_go_and_nothing_that_blocks_ever_does(home):
    week = 8 * 86400
    _write(_inbox(home), [
        _entry("run_finished", "Job 1 finished", age=3 * 86400),        # 2 days to live
        _entry("question_pending", "Question waiting: an old one", age=week),
    ])
    items = attention_feed.snapshot()
    assert [i["kind"] for i in items] == ["question"]


def test_an_entry_of_a_session_this_service_runs_is_left_to_that_session(home):
    _write(_inbox(home), [
        _entry("budget_warning", "Run budget nearly spent", session_id="ours"),
        _entry("budget_warning", "Run budget nearly spent", session_id="theirs", age=61),
    ])
    items = attention_feed.snapshot(frozenset({"ours"}))
    assert [i["sessionId"] for i in items] == ["theirs"]


def test_the_session_an_item_belongs_to_is_named_and_where_it_is_said(home, monkeypatch):
    from delfin.agent import session_presence, session_store
    _write(_inbox(home), [
        _entry("question_pending", "Question waiting: open one", detail="open one?",
               session_id="s-open", host=""),
        _entry("question_pending", "Question waiting: saved one", detail="saved one?",
               session_id="s-saved", age=61),
        _entry("question_pending", "Question waiting: gone one", detail="gone one?",
               session_id="s-gone", age=62),
        _entry("question_pending", "Question waiting: nobody's", detail="nobody's?",
               workspace="/calc/water", age=63),
    ])
    monkeypatch.setattr(session_presence, "open_sessions", lambda **_: [
        {"session_id": "s-open", "title": "scan of benzene", "host": "login02"}])
    saved = {"s-saved": {"title": "relax Fe complex"}, "s-open": {"title": "scan of benzene"}}
    monkeypatch.setattr(session_store, "load_session", lambda sid: saved.get(sid))
    monkeypatch.setattr(session_store, "session_lock_holder", lambda sid: None)
    titles = attention_feed.Titles()
    by = {i["sessionId"]: i for i in attention_feed.snapshot(titles=titles)}
    assert by["s-open"]["name"] == "scan of benzene"
    # Open in the dashboard and saved: a window opens a copy of it.
    assert by["s-open"]["sessionOpen"] is True and by["s-open"]["sessionSaved"] is True
    assert by["s-open"]["host"] == "login02"
    assert by["s-saved"]["name"] == "relax Fe complex"
    assert by["s-saved"]["sessionSaved"] is True and by["s-saved"]["sessionOpen"] is False
    assert by["s-gone"]["sessionSaved"] is False
    # A question of no session is named for where it was asked.
    assert by[""]["name"] == "a question in water"
    # Saved later (its first turn ended): found at the next look.
    saved["s-gone"] = {"title": "late one"}
    by = {i["sessionId"]: i for i in attention_feed.snapshot(titles=titles)}
    assert by["s-gone"]["sessionSaved"] is True and by["s-gone"]["name"] == "late one"


def test_answer_and_dismiss_go_through_the_inbox(home):
    q = _entry("question_pending", "Question waiting: which?", detail="which?")
    n = _entry("run_finished", "Job 9 finished", age=61)
    _write(_inbox(home), [q, n])
    assert attention_feed.answer(q["id"], "CAM-B3LYP")["ok"] is True
    assert attention_feed.dismiss(n["id"])["ok"] is True
    records = {r["id"]: r for r in map(json.loads, _inbox(home).read_text().splitlines())}
    assert records[q["id"]]["status"] == "resolved" and records[q["id"]]["answer"] == "CAM-B3LYP"
    assert records[q["id"]]["acknowledged"] is False          # for the next session to take
    assert records[n["id"]]["acknowledged"] is True           # never handed to a model
    assert attention_feed.snapshot() == []


def test_every_entry_says_which_host_wrote_it(home, monkeypatch):
    import socket

    from delfin.agent import attention
    monkeypatch.setattr(attention, "_deliver", lambda *a, **k: [])
    attention.emit_attention("run_finished", title="Job 5 finished")
    record = json.loads(_inbox(home).read_text().splitlines()[-1])
    assert record["host"] == socket.gethostname()


def test_watchpost_alerts_reach_the_inbox(home, monkeypatch):
    from types import SimpleNamespace

    from delfin.agent import attention, scheduler_daemon
    from delfin.watchpost import scan
    monkeypatch.setattr(attention, "_deliver", lambda *a, **k: [])
    monkeypatch.setattr(scan, "run_scan", lambda *a, **k: [
        SimpleNamespace(severity="alert", check="ssh", what="a new authorized key")])
    scheduler_daemon.run_watchpost_entry(SimpleNamespace(workspace=str(home)))
    kinds = [json.loads(l)["kind"] for l in _inbox(home).read_text().splitlines()]
    assert kinds == ["security_alert"]


# -- the service ------------------------------------------------------------------------------

@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _service_inbox(svc: ServiceProcess) -> Path:
    # The service runs on scratch state (child_env): its ~/.delfin is there.
    return svc.scratch / ".delfin" / "attention_inbox.jsonl"


def _service_findings(svc: ServiceProcess) -> Path:
    return svc.scratch / ".delfin" / "monitor_findings.jsonl"


def test_the_service_lists_them_and_says_its_version(service):
    _write(_service_inbox(service), [
        _entry("question_pending", "Question waiting: which?", detail="which functional?",
               workspace=str(service.workspace())),
        _entry("run_failed", "Job 4711 failed (OUT_OF_MEMORY)", workspace="/calc", age=61)])
    _write(_service_findings(service), [_finding("4711", "/calc", diagnosis="diag-1")])
    c = service.client()
    assert c.hello["serviceVersion"] >= 3
    listed = c.request("attention.list")
    items = listed["items"]
    assert [i["kind"] for i in items] == ["question", "failed"]
    assert items[1]["job"]["diagnosisSession"] == "diag-1"
    # Versioned: a window keeps the newest of a list and an event that
    # cross on the socket.
    assert listed["version"] >= 1
    _write(_service_inbox(service), [_entry("question_pending", "Question waiting: another?")])
    # (The list's own look broadcast its version too: the next one is newer.)
    got, _ = c.wait_for(lambda e: e.get("type") == "attention.updated"
                        and e.get("version", 0) > listed["version"], timeout=10)
    assert got["version"] == listed["version"] + 1
    assert [i["text"] for i in got["items"]] == ["question: another?"]


def test_a_change_reaches_every_window_and_answers_and_dismissals_go_through(service):
    a, b = service.client(), service.client()
    assert a.request("attention.list")["items"] == []
    q = _entry("question_pending", "Question waiting: which?", detail="which functional?")
    n = _entry("run_finished", "Job 7 finished", age=61)
    _write(_service_inbox(service), [q, n])
    for c in (a, b):
        got, _ = c.wait_for(lambda e: e.get("type") == "attention.updated" and len(e["items"]) == 2,
                            timeout=10)
    start = len(b.events)
    assert a.request("attention.answer", {"id": q["id"], "text": "CAM-B3LYP"})["ok"] is True
    got, _ = b.wait_for(lambda e: e.get("type") == "attention.updated", after=start, timeout=10)
    assert [i["id"] for i in got["items"]] == [n["id"]]
    start = len(b.events)
    assert a.request("attention.dismiss", {"id": n["id"]})["ok"] is True
    got, _ = b.wait_for(lambda e: e.get("type") == "attention.updated", after=start, timeout=10)
    assert got["items"] == []
    # An item that is gone already says so, and nothing breaks.
    assert a.request("attention.dismiss", {"id": n["id"]})["ok"] is False
    with pytest.raises(ServiceError):
        a.request("attention.answer", {"id": q["id"], "text": "  "})


def test_a_session_of_this_service_says_its_own_and_its_items_come_back_when_it_ends(service):
    """Left out while the session runs here -- its row says it -- and
    listed once it has ended (nothing in the files changed: what the
    service runs is part of what it looks at)."""
    c = service.client()
    sid = c.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                       "fakeDelay": 0.003}, timeout=60)["sessionId"]
    _write(_service_inbox(service), [_entry("budget_warning", "Budget nearly spent", session_id=sid)])
    time.sleep(2.5)
    assert c.request("attention.list")["items"] == []
    start = len(c.events)
    c.request("session.shutdown", {"sessionId": sid}, timeout=60)
    got, _ = c.wait_for(lambda e: e.get("type") == "attention.updated" and e["items"], after=start, timeout=15)
    assert [i["sessionId"] for i in got["items"]] == [sid]


def test_what_goes_out_is_masked_as_a_sessions_events_are(home):
    token = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6Q7r8"
    _write(_inbox(home), [_entry("confirm_pending", "Confirmation required: bash",
                                 detail=f"curl -H 'Authorization: Bearer {token}' https://api.github.com/user")])
    _write(_findings(home), [])
    items = attention_feed.snapshot()
    assert items and token not in json.dumps(items)
    assert "redacted" in items[0]["detail"]
    assert items[0]["id"].startswith("att-")                      # identifiers as they are


def test_a_saved_session_renamed_since_is_named_anew(home):
    from delfin.agent import session_store
    session_store.save_session("s-renamed", workspace="/calc", title="first name",
                               chat_messages=[{"role": "user", "content": "go"}])
    _write(_inbox(home), [_entry("question_pending", "Question waiting: x?", detail="x?",
                                 session_id="s-renamed")])
    titles = attention_feed.Titles()
    assert attention_feed.snapshot(titles=titles)[0]["name"] == "first name"
    time.sleep(0.05)
    session_store.save_session("s-renamed", workspace="/calc", title="second name",
                               chat_messages=[{"role": "user", "content": "go"}, {"role": "assistant", "content": "ok"}])
    assert attention_feed.snapshot(titles=titles)[0]["name"] == "second name"


HEADER = "What this session was started to look into"


def _answer_of(c, sid, start, end) -> str:
    return " ".join(e.get("text", "") for e in c.events[start:end + 1]
                    if e.get("session") == sid and e.get("type") == "text.ended")


def test_a_prompt_cancelled_before_its_turn_leaves_the_brief_for_the_next(service):
    """Enter, then Esc at once: the turn never starts, and the brief goes
    with the next prompt (it was marked as given when the prompt was
    routed, and lost)."""
    c = service.client()
    for _ in range(6):
        # A fresh session each try: a cancel that comes too late lets the
        # turn take the brief.
        sid = c.request("session.create", {"cwd": str(service.workspace()), "fake": True, "fakeDelay": 0.003,
                                           "brief": "Job 4702 ended with OUT_OF_MEMORY.",
                                           "briefTitle": "job 4702"}, timeout=60)["sessionId"]
        start = len(c.events)
        c.request("session.prompt", {"sessionId": sid, "text": "Summarise it"})
        c.request("session.cancel", {"sessionId": sid})
        _, end = c.wait_for(lambda e, s=sid: e.get("type") == "turn.ended" and e.get("session") == s,
                            after=start, timeout=30)
        if c.events[end].get("error") == "cancelled before it started":
            break
    else:
        pytest.skip("the cancel never came before the line was taken")
    start, end = turn(c, sid, "Summarise it")
    assert HEADER in _answer_of(c, sid, start, end)


def test_a_session_about_an_item_that_sleeps_first_still_hands_its_brief_over(tmp_path):
    """It waits for the user by design, and may sleep (or move to a new
    service) before its first prompt: woken, it shows its context row once
    and the first prompt still carries the brief."""
    svc = ServiceProcess(tmp_path, start=False)
    svc.env["DELFIN_AGENT_SLEEP_AFTER_S"] = "1"
    svc.start()
    try:
        c = svc.client()
        sid = c.request("session.create", {"cwd": str(svc.workspace()), "fake": True, "fakeDelay": 0.003,
                                           "brief": "Job 4702 ended with OUT_OF_MEMORY.",
                                           "briefTitle": "job 4702"}, timeout=60)["sessionId"]
        # A first turn that fails before the model says anything (the fake
        # engine fails on "fail") has not taken it either.
        start, end = turn(c, sid, "please fail")
        assert c.events[end].get("error")
        c.request("session.detach", {"sessionId": sid})
        deadline = time.monotonic() + 40
        while time.monotonic() < deadline:
            row = next((r for r in c.request("sessions.list")["sessions"] if r["sessionId"] == sid), None)
            if row and row["state"] == "sleeping":
                break
            time.sleep(0.3)
        else:
            pytest.fail("the session never slept")
        attached = c.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
        assert len([e for e in attached["events"] if e.get("type") == "context"]) == 1
        start, end = turn(c, sid, "Summarise it")
        assert HEADER in _answer_of(c, sid, start, end)
        # Given once: the next prompt goes as typed.
        start, end = turn(c, sid, "and then?")
        assert HEADER not in _answer_of(c, sid, start, end)
    finally:
        svc.stop()


def test_a_watch_only_connection_may_not_answer_or_dismiss(service):
    q = _entry("question_pending", "Question waiting: which?", detail="which?")
    _write(_service_inbox(service), [q])
    probe = r"""
import json, sys
sys.path.insert(0, sys.argv[2])
from delfin.agent.service.client import ServiceClient, ServiceError
c = ServiceClient(sys.argv[1])
out = {"listed": len(c.request("attention.list")["items"])}
for method, params in (("attention.answer", {"id": sys.argv[3], "text": "yes"}),
                       ("attention.dismiss", {"id": sys.argv[3]})):
    try:
        c.request(method, params)
        out[method] = "done"
    except ServiceError as exc:
        out[method] = exc.code
print(json.dumps(out))
"""
    from conftest import child_env
    env = child_env(service.tmp)
    env.update(DELFIN_RUNTIME_DIR=service.env["DELFIN_RUNTIME_DIR"],
               PYTHONPATH=service.env["PYTHONPATH"], DELFIN_AGENT_SESSION="12345")
    run = subprocess.run([sys.executable, "-c", probe, str(service.socket), str(REPO), q["id"]],
                         env=env, capture_output=True, text=True, timeout=60)
    got = json.loads(run.stdout.strip().splitlines()[-1])
    assert got == {"listed": 1, "attention.answer": protocol.NOT_ALLOWED,
                   "attention.dismiss": protocol.NOT_ALLOWED}
    assert json.loads(_service_inbox(service).read_text().splitlines()[0])["status"] == "pending"


def test_a_session_started_with_a_brief_shows_it_and_hands_it_to_the_model(service):
    c = service.client()
    brief = "Job 4702 in /calc/water finished: COMPLETED after 2h14m."
    created = c.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                           "fakeDelay": 0.003, "brief": brief,
                                           "briefTitle": "job 4702 finished"}, timeout=60)
    sid = created["sessionId"]
    contexts = [e for e in created["events"] if e.get("type") == "context"]
    assert contexts and contexts[0]["text"] == brief and contexts[0]["title"] == "job 4702 finished"
    assert contexts[0]["seq"] > 0                                   # durable: replayed later
    start, end = turn(c, sid, "Summarise its results")
    said = [e for e in c.events[start:end + 1] if e.get("session") == sid]
    users = [e["text"] for e in said if e.get("type") == "user.message"]
    assert users == ["Summarise its results"]                      # as typed
    answer = " ".join(e.get("text", "") for e in said if e.get("type") == "text.ended")
    assert "What this session was started to look into" in answer  # the model had the brief first
    # Only the first prompt carries it.
    start, end = turn(c, sid, "and then?")
    answer = " ".join(e.get("text", "") for e in c.events[start:end + 1]
                      if e.get("session") == sid and e.get("type") == "text.ended")
    assert "started to look into" not in answer
    # A replay shows the brief as the session's first row again.
    replay = c.request("session.attach", {"sessionId": sid, "turns": 0})["events"]
    assert any(e.get("type") == "context" for e in replay)


def test_a_note_or_a_command_before_the_first_prompt_leaves_the_brief_for_it(service):
    c = service.client()
    brief = "Job 4711 in /calc/benzene ended with OUT_OF_MEMORY."   # (a fake turn fails on "fail")
    created = c.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                           "fakeDelay": 0.003, "brief": brief,
                                           "briefTitle": "job 4711 ended"}, timeout=60)
    sid = created["sessionId"]
    # A "#" note is remembered, not sent to the model: the brief waits.
    start = len(c.events)
    c.request("session.prompt", {"sessionId": sid, "text": "# the basis set is def2-TZVP"})
    c.wait_for(lambda e: e.get("session") == sid and e.get("type") == "notice"
               and "remember" in str(e.get("text", "")), after=start, timeout=30)
    start, end = turn(c, sid, "Summarise what happened")
    answer = " ".join(e.get("text", "") for e in c.events[start:end + 1]
                      if e.get("session") == sid and e.get("type") == "text.ended")
    assert "What this session was started to look into" in answer
