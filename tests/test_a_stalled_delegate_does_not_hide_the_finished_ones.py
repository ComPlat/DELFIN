"""Collecting delegates was head-of-line blocked, and background ones
could not be found at all.

**Fan-out collected in emission order.** The pool submits every
``subagent`` call concurrently, but the results were taken one by one in
the order the model emitted them, each with its OWN timeout of
``max_wall_s + 120`` (1020 s at the default). A delegate that finished in
5 seconds went unnoticed for up to that whole budget behind a stalled
predecessor, and N stalled delegates cost N times it. Collection now runs
under ONE deadline for the whole fan-out, taking results as they land.

**A thread that died before the registry entry existed reported
``unknown``.** The live entry was written inside the run, so anything
that failed before it -- the runner raising on the way in, a foreign
runner that stores nothing -- left the reserved id resolving to "no such
id", which this module's own docstring says the model reads as "I made
the id up". The entry is now reserved when the id is, before the thread
starts, and the background runner records a ``died`` entry instead of
swallowing the exception. "Started and crashed" and "never existed" are
different answers, and only one of them tells the parent to start the
work again.

**The note pointed at the wrong id.** For ``background=true`` with a
``resume_id`` the collectable id is the resumed one, and the returned
``sa_id`` said so -- while the human-readable note next to it still named
the reserved id the parent must not poll.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import re
import threading
import time

import pytest

from delfin.agent import api_client as ac
from delfin.agent import subagents as sa


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")


# ---------------------------------------------------------------------------
# One deadline for the fan-out, results taken as they land
# ---------------------------------------------------------------------------

@pytest.fixture
def pool():
    ex = cf.ThreadPoolExecutor(max_workers=8)
    try:
        yield ex
    finally:
        ex.shutdown(wait=False)


def _stall(seconds: float = 30.0):
    def _fn():
        time.sleep(seconds)
        return "late"
    return _fn


def test_a_finished_delegate_is_collected_while_a_sibling_stalls(pool):
    futures = {
        "stalled": pool.submit(_stall()),
        "quick": pool.submit(lambda: "REPORT"),
    }
    out = ac._collect_subagent_futures(futures, 0.5)
    assert out["quick"] == "REPORT"
    assert "stalled" not in out


def test_the_deadline_covers_the_whole_fan_out_not_each_item(pool):
    futures = {f"s{i}": pool.submit(_stall()) for i in range(4)}
    t0 = time.monotonic()
    out = ac._collect_subagent_futures(futures, 0.4)
    elapsed = time.monotonic() - t0
    assert out == {}
    assert elapsed < 4 * 0.4, elapsed


def test_a_delegate_that_raised_reports_its_failure(pool):
    def _boom():
        raise RuntimeError("child exploded")

    out = ac._collect_subagent_futures({"x": pool.submit(_boom)}, 5.0)
    assert "child exploded" in json.loads(out["x"])["error"]


def test_no_futures_is_not_a_wait(pool):
    assert ac._collect_subagent_futures({}, 5.0) == {}


def test_the_tool_loop_collects_them_this_way():
    """The loop must not fall back to per-item waiting."""
    import inspect
    src = inspect.getsource(ac.OpenAIClient.stream_message)
    assert "_collect_subagent_futures" in src


# ---------------------------------------------------------------------------
# A background id resolves from the moment it is handed out
# ---------------------------------------------------------------------------

class _Perms:
    def __init__(self, runner):
        self.subagent_runner = runner
        self.subagent_depth = 0


def _spawn(runner, **extra) -> dict:
    args = {"subagent_type": "explore", "description": "probe the thing",
            "prompt": "look around carefully for the bug in module x",
            "background": True}
    args.update(extra)
    return json.loads(ac._doc_executor._execute_subagent(args, _Perms(runner)))


def test_a_background_id_resolves_while_its_run_is_still_starting():
    release = threading.Event()
    started = threading.Event()

    def _runner(**kw):
        started.set()
        release.wait(timeout=5)
        return {}

    out = _spawn(_runner)
    try:
        assert started.wait(timeout=5)
        assert sa.get_subagent_result(out["sa_id"])["status"] == "running"
    finally:
        release.set()


def test_a_thread_that_dies_early_is_not_reported_as_a_made_up_id():
    done = threading.Event()

    def _runner(**kw):
        done.set()
        raise RuntimeError("no workspace")

    out = _spawn(_runner)
    assert done.wait(timeout=5)
    for _ in range(50):
        status = sa.get_subagent_result(out["sa_id"])["status"]
        if status != "running":
            break
        time.sleep(0.05)
    assert status == "died"


def test_the_died_record_says_what_happened():
    def _runner(**kw):
        raise RuntimeError("no workspace")

    out = _spawn(_runner)
    for _ in range(50):
        res = sa.get_subagent_result(out["sa_id"])
        if res["status"] != "running":
            break
        time.sleep(0.05)
    assert "no workspace" in res.get("error", "")


def test_a_run_that_stored_no_report_does_not_stay_running_forever():
    out = _spawn(lambda **kw: {"ok": True})
    for _ in range(50):
        status = sa.get_subagent_result(out["sa_id"])["status"]
        if status != "running":
            break
        time.sleep(0.05)
    assert status == "died"


def test_a_died_entry_is_not_shown_as_running(tmp_path):
    sa.reserve_running("abc123", subagent_type="explore", description="d")
    assert "abc123" in sa.read_running()
    sa.mark_running_died("abc123", "thread ended")
    assert "abc123" not in sa.read_running()
    assert "abc123" in sa.read_running(include_dead=True)


def test_a_finished_run_still_reads_as_finished(monkeypatch):
    sa.reserve_running("abc123", subagent_type="explore", description="d")
    monkeypatch.setattr(sa, "load_subagent_session", lambda i: {
        "subagent_type": "explore", "description": "d",
        "final_report": "REPORT", "error": "", "interactions": []})
    monkeypatch.setattr(sa, "read_running", lambda **kw: {})
    assert sa.get_subagent_result("abc123")["status"] == "finished"


# ---------------------------------------------------------------------------
# The note names the id the parent has to poll
# ---------------------------------------------------------------------------

def test_the_note_names_the_collectable_id_on_a_resume():
    out = _spawn(lambda **kw: {}, resume_id="deadbeef")
    assert out["sa_id"] == "deadbeef"
    assert "deadbeef" in out["note"]


def test_the_note_does_not_name_an_id_nothing_will_answer_for():
    out = _spawn(lambda **kw: {}, resume_id="deadbeef")
    named = re.findall(r"sa_id='([^']*)'", out["note"])
    assert named == ["deadbeef"]


def test_a_fresh_background_run_still_names_its_own_id():
    out = _spawn(lambda **kw: {})
    assert out["sa_id"] in out["note"]
