"""The agent service stops without leaving anything behind, and what it
tells a window agrees with what it shows.

* A session started while the service stops is refused or ended with the
  rest; no worker outlives the service. A stop noticed before the service
  serves is obeyed. Every session is asked to end at once, one still
  starting is ended at once, and one that stops reading holds nothing up.
* `service stop` returns once the service is gone, so the next window
  starts at once, and a window opened during a stop waits for it; a
  service that does not say hello (stuck, another protocol version) is
  named, and `service stop` still ends it -- only this shell's service.
* An attach answer and its replay stand at the same point, the replay
  comes from memory (no log read in the service's loop) and fits its
  budget; a resumed session replays how it is now, not how it ended; what
  a worker says right after it started reaches the window that created
  it. Granted folders are absolute, as the session's own folder is.
* A window that stops reading is dropped with its backlog; one that reads
  survives a burst.
"""

from __future__ import annotations

import asyncio
import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from conftest import child_env
from delfin.agent.service import eventlog, protocol, rundir, supervisor
from delfin.agent.service.client import ServiceError
from tui_service import ServiceProcess, of, turn


def _env(service: ServiceProcess) -> dict:
    env = child_env(service.tmp)
    env.update(DELFIN_RUNTIME_DIR=service.env["DELFIN_RUNTIME_DIR"],
               PYTHONPATH=service.env["PYTHONPATH"])
    return env


def _cli(service: ServiceProcess, *args: str, timeout: float = 120,
         prelude: str = "") -> subprocess.CompletedProcess:
    code = (f"{prelude}\nimport sys\nfrom delfin.agent import cli\n"
            f"sys.argv = ['delfin-agent', 'service', *{list(args)!r}]\n"
            "sys.exit(cli.main())")
    return subprocess.run([sys.executable, "-c", code], env=_env(service),
                          capture_output=True, text=True, timeout=timeout)


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


def _children_of(pid: int) -> list[int]:
    out = []
    for entry in Path("/proc").iterdir():
        if entry.name.isdigit():
            try:
                stat = (entry / "stat").read_text()
                if int(stat.rsplit(")", 1)[1].split()[1]) == pid:
                    out.append(int(entry.name))
            except (OSError, IndexError, ValueError):
                continue
    return out


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


# -- stopping -------------------------------------------------------------------------

def test_a_session_created_as_the_service_stops_does_not_outlive_it(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        a = svc.client()
        create = a.send("session.create", {"cwd": str(svc.workspace()), "fake": True,
                                           "fakeDelay": 0.003})
        a.send("service.shutdown")
        workers_seen: set[int] = set()
        deadline = time.monotonic() + 60
        while svc.proc.poll() is None and time.monotonic() < deadline:
            workers_seen.update(_children_of(svc.proc.pid))
            time.sleep(0.05)
        assert svc.proc.wait(timeout=60) == 0
        with pytest.raises(ServiceError) as refused:
            a.result(create, timeout=5)
        assert refused.value.code == protocol.NOT_ALLOWED
        assert not [p for p in workers_seen if _alive(p)]
        # Whatever worker started ended the orderly way, and saved.
        for log in (svc.scratch / ".delfin" / "agent_service" / "sessions").glob("*.jsonl"):
            reasons = [e.get("reason") for e in eventlog.read_events(log)
                       if e.get("type") == "session.ended"]
            assert reasons == ["shutdown"], reasons
    finally:
        svc.stop()


def test_nothing_acts_on_a_service_that_is_stopping(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    service_pid = service.proc.pid
    a.request("service.shutdown")            # answered once it is closing
    for method, params in (("session.create", {"cwd": str(service.workspace())}),
                           ("session.prompt", {"sessionId": sid, "text": "hi"})):
        with pytest.raises((ServiceError, OSError)) as refused:
            b.request(method, params, timeout=10)
        if isinstance(refused.value, ServiceError):
            assert refused.value.code in (protocol.NOT_ALLOWED, -1)
    assert service.proc.wait(timeout=60) == 0 and service_pid


def test_a_stop_noticed_before_the_service_serves_is_obeyed(tmp_path, monkeypatch):
    # In this process: its state goes to the test's own folders.
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "rt"))
    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    sup = supervisor.Supervisor(idle_exit_s=300)
    sup.stop()                           # e.g. the emergency stop, while starting
    path = Path(f"/tmp/dsvstop-{os.getpid()}.sock")

    async def run() -> float:
        started = time.monotonic()
        await asyncio.wait_for(sup.serve(path), timeout=10)
        return time.monotonic() - started

    assert asyncio.run(run()) < 5
    assert not path.exists()


def test_every_session_is_asked_to_end_at_once():
    asked: list[float] = []

    class Slow:
        id = "s"                              # started: it is asked, not signalled

        def __init__(self, loop):
            self.done = loop.create_future()
            self.proc = SimpleNamespace(returncode=None, send_signal=lambda s: None,
                                        kill=lambda: None)

        async def request(self, method, params=None, timeout=None):
            asked.append(time.monotonic())
            await asyncio.sleep(0.5)          # a session that takes its time
            self.proc.returncode = 0
            self.done.set_result(0)

    async def run() -> float:
        loop = asyncio.get_running_loop()
        sup = supervisor.Supervisor()
        sessions = [Slow(loop) for _ in range(4)]
        started = time.monotonic()
        await sup._end(sessions)
        return time.monotonic() - started

    took = asyncio.run(run())
    assert took < 1.5                         # one after another would be 2 s
    assert max(asked) - min(asked) < 0.2


def test_a_session_still_starting_is_ended_at_once():
    signals: list[int] = []

    class Starting:
        id = ""                               # its initialize not answered yet

        def __init__(self, loop):
            self.done = loop.create_future()
            self.proc = SimpleNamespace(returncode=None, send_signal=self._signal,
                                        kill=lambda: None)

        def _signal(self, sig):
            signals.append(sig)
            self.proc.returncode = -sig
            self.done.set_result(-sig)

        async def request(self, method, params=None, timeout=None):
            raise AssertionError("a worker still starting reads no request")

    async def run() -> float:
        sup = supervisor.Supervisor()
        started = time.monotonic()
        await sup._end([Starting(asyncio.get_running_loop())])
        return time.monotonic() - started

    assert asyncio.run(run()) < 1             # not the 5 s ask + 20 s grace
    assert signals == [signal.SIGTERM]


def test_a_worker_that_stops_reading_holds_no_request_and_no_stop():
    class FullPipe:
        def write(self, data):
            pass

        async def drain(self):
            await asyncio.Event().wait()      # never: nobody reads the pipe

    async def run() -> float:
        sup = supervisor.Supervisor()
        proc = SimpleNamespace(stdin=FullPipe(), returncode=None, pid=0)
        session = supervisor.Session(sup, proc, {})
        started = time.monotonic()
        with pytest.raises(TimeoutError):
            await session.request("session.shutdown", timeout=0.3)
        assert session._waiting == {}
        return time.monotonic() - started

    assert asyncio.run(run()) < 2


def test_service_stop_returns_once_the_service_is_gone(service):
    a = service.client()
    service.session(a)
    old = service.proc.pid
    stopped = _cli(service, "stop")
    assert stopped.returncode == 0 and "stopped" in stopped.stdout, stopped.stderr
    assert not _alive(old)
    # The next window starts a service at once, never "exited before it answered".
    started = subprocess.run(
        [sys.executable, "-c", "from delfin.agent.service import launcher; "
                               "print(launcher.ensure_service()['pid'])"],
        env=_env(service), capture_output=True, text=True, timeout=60)
    assert started.returncode == 0, started.stderr
    new = int(started.stdout.strip())
    assert new != old
    os.kill(new, signal.SIGTERM)


def test_a_service_that_does_not_answer_is_named_and_stopped(service):
    pid = service.proc.pid
    os.kill(pid, signal.SIGSTOP)
    try:
        status = _cli(service, "status")
        assert status.returncode == 1
        assert f"(pid {pid}), but it did not answer" in status.stdout
        stopped = _cli(service, "stop", prelude=(
            "from delfin.agent.service import cli_service\n"
            "cli_service.UNANSWERING_TERM_WAIT_S = 1.0"))
        assert stopped.returncode == 0, stopped.stdout + stopped.stderr
        assert "SIGKILL" in stopped.stdout
    finally:
        try:
            os.kill(pid, signal.SIGCONT)
        except ProcessLookupError:
            pass
    assert service.proc.wait(timeout=10) != 0


def test_a_service_of_another_protocol_version_is_named(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    code = textwrap.dedent("""
        import importlib, sys
        from delfin.agent.service import protocol
        protocol.PROTOCOL_VERSION = 99           # a release from the future
        sys.exit(importlib.import_module("delfin.agent.service.supervisor").main())
    """)
    proc = subprocess.Popen([sys.executable, "-c", code], env=svc.env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            start_new_session=True)
    try:
        deadline = time.monotonic() + 30
        while not svc.socket.exists() and time.monotonic() < deadline:
            time.sleep(0.1)
        # Busy, by its roster: a launcher of this version must not cut it.
        roster = svc.scratch / ".delfin" / "agent_service" / "run" / f"{rundir.hostname()}.roster"
        deadline = time.monotonic() + 30
        while not roster.exists() and time.monotonic() < deadline:
            time.sleep(0.1)
        roster.write_text(json.dumps({"pid": proc.pid, "sessions": [
            {"sessionId": "busy1", "state": "running", "pending": 0, "workspace": str(svc.tmp)}]}))
        probe = subprocess.run(
            [sys.executable, "-c",
             "from delfin.agent.service import launcher\n"
             "try:\n    launcher.ensure_service(timeout=5)\n"
             "except launcher.ServiceDidNotStart as e:\n    print(e)\n"],
            env=svc.env, capture_output=True, text=True, timeout=60)
        assert "speaks another protocol version" in probe.stdout, probe.stdout + probe.stderr
        assert f"(pid {proc.pid})" in probe.stdout
        status = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", "status"],
                                env=svc.env, capture_output=True, text=True, timeout=60)
        assert f"(pid {proc.pid}), but it speaks another protocol version" in status.stdout
    finally:
        proc.terminate()
        proc.wait(timeout=60)
        svc.stop()


def test_a_service_of_another_runtime_dir_is_not_this_shells(service):
    # The same HOME, another DELFIN_RUNTIME_DIR: the run file names the
    # service, but not this shell's socket -- nothing to name or stop here.
    other = _env(service)
    other["DELFIN_RUNTIME_DIR"] = str(service.tmp / "rt-other")
    stopped = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", "stop"],
                             env=other, capture_output=True, text=True, timeout=60)
    assert "no agent service runs" in stopped.stdout, stopped.stdout + stopped.stderr
    time.sleep(0.5)
    assert service.proc.poll() is None
    assert service.client().request("service.status")


def test_a_window_waits_for_a_service_that_is_stopping(tmp_path):
    # A session that takes its time to save (its worker stopped, SIGSTOP)
    # holds the service's stop up; a window opened meanwhile waits for the
    # stop and then starts the next service -- it never fails half-way
    # with "run `service stop`". The grace is cut short for the test.
    svc = ServiceProcess(tmp_path, start=False)
    code = textwrap.dedent("""
        import sys
        import delfin.agent.service.supervisor as supervisor   # (its command line says so)
        supervisor.SHUTDOWN_GRACE_S = 4.0
        sys.exit(supervisor.main(sys.argv[1:]))
    """)
    svc.proc = subprocess.Popen([sys.executable, "-c", code, "--idle-exit", "120"],
                                env=svc.env, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL, start_new_session=True)
    try:
        svc.wait_until_up()
        a = svc.client()
        sid = svc.session(a)
        worker = next(s["pid"] for s in a.request("sessions.list")["sessions"]
                      if s["sessionId"] == sid)
        assert worker in _children_of(svc.proc.pid)
        os.kill(worker, signal.SIGSTOP)
        old = svc.proc.pid
        a.request("service.shutdown")
        status = _cli(svc, "status")
        assert "starting or stopping" in status.stdout, status.stdout
        opened = subprocess.run(
            [sys.executable, "-c", "from delfin.agent.service import launcher; "
                                   "print(launcher.ensure_service(timeout=5)['pid'])"],
            env=_env(svc), capture_output=True, text=True, timeout=120)
        assert opened.returncode == 0, opened.stdout + opened.stderr
        new = int(opened.stdout.strip())
        assert new != old and not _alive(old) and not _alive(worker)
        os.kill(new, signal.SIGTERM)
    finally:
        try:
            os.kill(worker, signal.SIGCONT)
        except (ProcessLookupError, UnboundLocalError):
            pass
        svc.stop()


def test_the_cli_says_when_this_shell_may_not_act(service):
    service.session(service.client())
    env = _env(service)
    env["DELFIN_AGENT_SESSION"] = "1"            # as in an agent session's shell
    out = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "service", "stop"],
                         env=env, capture_output=True, text=True, timeout=60)
    assert out.returncode == 1
    assert "Traceback" not in out.stderr
    assert "may watch, not act" in out.stderr
    assert service.proc.poll() is None


def test_status_counts_windows_not_itself(service):
    a = service.client()
    service.session(a)
    status = _cli(service, "status")
    assert "1 session, 0 windows" in status.stdout


# -- what a window is told -----------------------------------------------------------------

def test_a_session_folder_must_be_an_absolute_existing_path(service):
    a = service.client()
    for cwd in ("", "project", str(service.tmp / "missing")):
        with pytest.raises(ServiceError) as refused:
            a.request("session.create", {"cwd": cwd, "fake": True})
        assert refused.value.code == protocol.INVALID_PARAMS


def test_a_folder_granted_to_a_session_must_be_absolute(service):
    a = service.client()
    for key in ("addDirs", "readDirs"):
        # Taken relative to the service's own folder ($HOME), it would grant that.
        with pytest.raises(ServiceError) as refused:
            a.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                         key: ["data"]})
        assert refused.value.code == protocol.INVALID_PARAMS
        assert key in refused.value.message


def test_a_resumed_session_replays_how_it_is_now_not_how_it_ended(service):
    a = service.client()
    sid = service.session(a, permissionMode="default")
    a.request("session.setMode", {"sessionId": sid, "mode": "plan"})
    turn(a, sid, "hello")
    mark = len(a.events)
    a.request("session.shutdown", {"sessionId": sid})
    a.wait_for(lambda e: e.get("type") == "session.removed" and e.get("session") == sid,
               after=mark, timeout=60)
    created = a.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                           "fakeDelay": 0.003, "resume": sid,
                                           "permissionMode": "acceptEdits",
                                           "attach": False}, timeout=120)
    assert created["sessionId"] == sid
    got = a.request("session.attach", {"sessionId": sid, "turns": 4})
    kinds = [e["type"] for e in got["events"]]
    # The earlier run's end and posture are history: a window told
    # "ended" closes the session, one told "plan" shows the wrong posture.
    assert "session.ended" not in kinds and "session.updated" not in kinds
    assert kinds.count("session.info") == 1 and kinds[0] == "session.info"
    assert got["events"][0]["approval"] == "acceptEdits" == got["session"]["approval"]
    assert "user.message" in kinds and got["state"] == "idle"
    # One session, one worker: a second resume while it runs is refused.
    with pytest.raises(ServiceError) as refused:
        a.request("session.create", {"cwd": str(service.workspace()), "fake": True,
                                     "resume": sid}, timeout=60)
    assert refused.value.code == protocol.NOT_ALLOWED


def test_the_creator_hears_of_its_session_in_the_answer_first(service):
    a = service.client()
    mark = len(a.events)
    msg = a.send("session.create", {"cwd": str(service.workspace()), "fake": True,
                                    "fakeDelay": 0.003})
    created = a.result(msg, timeout=60)
    # What the worker said since it started comes with the answer.
    assert [e["type"] for e in created["events"]][:1] == ["session.info"]
    assert all(isinstance(e.get("seq"), int) for e in created["events"])
    a.wait_for("session.summary", after=mark)
    # The summary is sent after the answer (it would name an unknown session).
    with a._cond:
        assert msg not in a._responses
    sid = created["sessionId"]
    start, end = turn(a, sid, "hello")
    seqs = [e["seq"] for e in of(a, sid, start) if "seq" in e]
    assert seqs == sorted(seqs) and seqs[0] > max(e["seq"] for e in created["events"])


def test_an_attach_and_its_replay_stand_at_the_same_point(service):
    a, b = service.client(), service.client()
    sid = service.session(a, permissionMode="default")
    problems = []
    for _ in range(6):
        mark = len(a.events)
        a.request("session.prompt", {"sessionId": sid, "text": "please edit the parser"})
        asked, _ = a.wait_for(lambda e: e.get("type") == "permission.asked"
                              and e.get("session") == sid, after=mark)
        for _ in range(8):
            got = b.request("session.attach", {"sessionId": sid, "turns": 2})
            open_in_replay = {}
            for e in got["events"]:
                if e["type"] == "permission.asked":
                    open_in_replay[e["requestId"]] = e
                elif e["type"] == "permission.resolved":
                    open_in_replay.pop(e["requestId"], None)
            if set(open_in_replay) != {p["requestId"] for p in got["pending"]}:
                problems.append((sorted(open_in_replay), [p["requestId"] for p in got["pending"]]))
        a.request("permission.reply", {"sessionId": sid, "requestId": asked["requestId"],
                                       "choice": "y"})
        for _ in range(8):
            got = b.request("session.attach", {"sessionId": sid, "turns": 2})
            open_in_replay = set()
            for e in got["events"]:
                if e["type"] == "permission.asked":
                    open_in_replay.add(e["requestId"])
                elif e["type"] == "permission.resolved":
                    open_in_replay.discard(e["requestId"])
            if open_in_replay != {p["requestId"] for p in got["pending"]}:
                problems.append((sorted(open_in_replay), [p["requestId"] for p in got["pending"]]))
        a.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid,
                   after=mark)
    assert problems == []


def test_an_attach_replays_from_memory_not_from_the_log(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    turn(a, sid, "hello")
    log = service.scratch / ".delfin" / "agent_service" / "sessions" / f"{sid}.jsonl"
    log.rename(log.with_suffix(".moved"))
    got = b.request("session.attach", {"sessionId": sid, "turns": 1})
    assert [e["type"] for e in got["events"]][:2] == ["session.info", "user.message"]
    assert any(e["type"] == "turn.ended" for e in got["events"])


def test_a_peek_is_answered_from_memory_too(service):
    # The agent view peeks (session.replay) on every key: the service
    # answers from what it holds, as for an attach, and reads no log.
    a, b = service.client(), service.client()
    sid = service.session(a)
    turn(a, sid, "hello")
    log = service.scratch / ".delfin" / "agent_service" / "sessions" / f"{sid}.jsonl"
    log.rename(log.with_suffix(".moved"))
    got = b.request("session.replay", {"sessionId": sid, "turns": 1})
    kinds = [e["type"] for e in got["events"]]
    assert kinds[:2] == ["session.info", "user.message"] and "turn.ended" in kinds
    with pytest.raises(ServiceError) as refused:
        b.request("session.replay", {"sessionId": sid, "turns": "two"})
    assert refused.value.code == protocol.INVALID_PARAMS


def _turn_events(n: int, size: int) -> list[dict]:
    events, seq = [], 0
    for t in range(1, n + 1):
        for kind, extra in (("user.message", {"text": f"turn {t}", "source": "user"}),
                            ("turn.started", {"turnId": f"t{t}", "source": "user"}),
                            ("text.ended", {"partId": f"p{t}", "text": "x" * size}),
                            ("turn.ended", {"turnId": f"t{t}"})):
            seq += 1
            events.append({"type": kind, "session": "s", "seq": seq, **extra})
    return events


def test_a_replay_keeps_the_newest_whole_turns_that_fit_its_budget():
    events = _turn_events(10, 1000)
    picked = eventlog.tail_turns(events, 0, max_bytes=3500)
    assert picked[0]["type"] == "notice" and "not shown" in picked[0]["text"]
    kept = [e for e in picked if e["type"] == "user.message"]
    assert [e["text"] for e in kept] == ["turn 9", "turn 10"]
    assert eventlog.tail_turns(events, 0, max_bytes=10**9) == events
    window = eventlog.RecentEvents(3500)
    for e in events:
        window.add(e)
    assert window.dropped and window.events()[0]["type"] == "user.message"
    assert window.turns() == 2
    # One turn larger than the budget is kept from where it fits.
    huge = eventlog.RecentEvents(500)
    for e in _turn_events(1, 2000):
        huge.add(e)
    assert huge.events() and huge.events()[-1]["type"] == "turn.ended"


def _log(tmp_path: Path, events: list[dict], tail: bytes = b"") -> Path:
    path = tmp_path / "s.jsonl"
    path.write_bytes(b"".join(protocol.encode(e) for e in events) + tail)
    return path


def test_older_turns_merge_under_the_ones_seen_live(tmp_path):
    window = eventlog.RecentEvents(10**6)
    events = _turn_events(3, 10)
    for e in events[4:]:
        window.add(e)                        # seen live: turns 2 and 3
    older, cut = eventlog.read_newest(_log(tmp_path, events), 10**6)
    window.merge_older_sized(older, cut)     # the log: all of it
    assert [e["seq"] for e in window.events()] == [e["seq"] for e in events]
    assert not window.dropped


def test_a_resumed_log_is_read_from_its_end_within_the_budget(tmp_path):
    events = _turn_events(10, 1000)
    # A line still being written at the end, and one damaged on the way.
    path = _log(tmp_path, events, tail=b'{"type": "text.ended", "te')
    lines = path.read_bytes().split(b"\n")
    lines[3] = b"garbage"
    path.write_bytes(b"\n".join(lines))
    whole = eventlog.read_events(path)
    for chunk in (7, 1000, 1 << 20):         # lines across the read blocks
        got, cut = eventlog.read_newest(path, 10**9, chunk=chunk)
        assert [e for e, _ in got] == whole and not cut
        assert all(size == len(protocol.encode(e)) for e, size in got)
        newest, cut = eventlog.read_newest(path, 3500, chunk=chunk)
        assert cut and [e for e, _ in newest] == whole[-len(newest):]
        assert sum(size for _, size in newest) <= 3500
    assert eventlog.read_newest(tmp_path / "missing.jsonl", 100) == ([], False)


def test_a_replay_under_the_current_header_drops_what_it_supersedes():
    header = {"type": "session.info", "session": "s", "seq": 0, "approval": "acceptEdits"}
    run1 = [{"type": "session.info", "session": "s", "seq": 1, "approval": "default"},
            {"type": "session.updated", "session": "s", "seq": 2, "approval": "plan"},
            *[{**e, "seq": e["seq"] + 2} for e in _turn_events(1, 10)],
            {"type": "session.ended", "session": "s", "seq": 7, "reason": "exit"}]
    run2 = [{"type": "session.info", "session": "s", "seq": 8, "approval": "acceptEdits"}]
    picked = eventlog.tail_turns(run1 + run2, 0, header=header)
    assert picked[0] is header
    assert [e["type"] for e in picked[1:]] == ["user.message", "turn.started",
                                                "text.ended", "turn.ended"]
    # The end of the run that ends now is the window's to hear.
    ending = run1 + run2 + [{"type": "session.ended", "session": "s", "seq": 9,
                             "reason": "shutdown"}]
    assert eventlog.tail_turns(ending, 0, header=header)[-1]["seq"] == 9


def test_a_window_that_stops_reading_is_dropped_with_its_backlog(monkeypatch):
    monkeypatch.setattr(supervisor, "CLIENT_STUCK_S", 0.1)
    monkeypatch.setattr(supervisor, "CLIENT_BACKLOG_BYTES", 100)
    gone = []

    class Transport:
        def __init__(self):
            self.size = 0
            self.aborted = False

        def get_write_buffer_size(self):
            return self.size

        def abort(self):
            self.aborted = True

    class Writer:
        def __init__(self):
            self.transport = Transport()

        def write(self, data):
            self.transport.size += len(data)

        def close(self):
            pass

    async def run():
        sup = SimpleNamespace(client_gone=gone.append)
        stuck = supervisor.Client(sup, None, Writer(), True, 1)
        reading = supervisor.Client(sup, None, Writer(), True, 2)
        stuck.send_bytes(b"x" * 500)
        reading.send_bytes(b"x" * 500)
        reading.writer.transport.size = 0          # it reads: caught up
        await asyncio.sleep(0.3)
        return stuck, reading

    stuck, reading = asyncio.run(run())
    assert stuck.closed and stuck.writer.transport.aborted
    assert not reading.closed
    assert gone == [stuck]


def test_what_a_worker_says_right_after_it_starts_reaches_its_creator(tmp_path):
    # A worker that is busy from its first moment and asks a question that
    # is settled at once: the service's view must be the worker's.
    stub = tmp_path / "busy_worker.py"
    stub.write_text(textwrap.dedent("""
        import json, sys, time
        from delfin.agent.service import protocol
        out = sys.stdout.buffer
        def send(m):
            out.write(protocol.encode(m)); out.flush()
        seq = 0
        def event(type_, **f):
            global seq
            if type_ in protocol.DURABLE_EVENTS:
                seq += 1
                f["seq"] = seq
            send(protocol.event({"type": type_, "session": "stub1", "ts": 0, **f}))
        info = {"sessionId": "stub1", "workspace": "/w", "provider": "fake", "model": "m",
                "mode": "solo", "approval": "default", "banner": [], "notices": []}
        event("session.info", **info)
        event("session.state", state="idle")
        for line in sys.stdin.buffer:
            msg = json.loads(line)
            if msg.get("method") == "initialize":
                send(protocol.result(msg["id"], {"protocolVersion": 1, "pid": 1,
                     "session": info, "state": "idle", "pending": []}))
                event("session.state", state="running")
                event("permission.asked", requestId="p1", kind="confirm", tool="bash",
                      preview="$ ls", options=[], flags={})
                event("permission.resolved", requestId="p1", decision="expired", by="expired")
                event("notice", text="started at once")
            elif msg.get("method") == "session.shutdown":
                send(protocol.result(msg["id"], {"closing": True}))
                event("session.ended", reason="shutdown", exitCode=0)
                break
    """))
    main = tmp_path / "stub_main.py"
    main.write_text(textwrap.dedent(f"""
        import sys
        from delfin.agent.service import supervisor as s
        _init = s.Supervisor.__init__
        def init(self, *a, **k):
            k["worker_cmd"] = [sys.executable, {str(stub)!r}]
            _init(self, *a, **k)
        s.Supervisor.__init__ = init
        sys.exit(s.main(["--idle-exit", "120"]))
    """))
    svc = ServiceProcess(tmp_path, start=False)
    proc = subprocess.Popen([sys.executable, str(main)], env=svc.env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            start_new_session=True)
    svc.proc = proc
    try:
        svc.wait_until_up()
        a = svc.client()
        created = a.request("session.create", {"cwd": str(svc.workspace())}, timeout=60)
        listed = a.request("sessions.list")["sessions"][0]
        # Either the answer already stands after the worker's first events,
        # or the summary that follows it does; never the older state.
        a.wait_for(lambda e: e.get("type") == "session.summary"
                   and e["summary"]["state"] == "running", timeout=10)
        assert listed["pending"] == 0
        assert created["pending"] == [] or listed["pending"] == 0
        said = [e["type"] for e in created["events"]] + \
            [e["type"] for e in a.events if e.get("session") == "stub1"]
        assert "notice" in said
    finally:
        svc.stop()
