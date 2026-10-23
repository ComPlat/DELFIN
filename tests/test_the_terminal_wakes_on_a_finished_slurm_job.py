"""The terminal wakes when a SLURM calculation IT submitted ends.

Reviewer adversarial tests (package C, welle 11). The terminal idle wake
(``repl._wake_text``) today reports finished *shells* only; a SLURM job that
finished while the session idled left it silent for up to 90 min. These tests
exercise the promised fix: ``job_wake.finished_watched_jobs(workspace, seen,
session_id, run_fn=None)`` — one ``check_agent_jobs(consume=False,
marker="wake_notified")`` call, session-scoped, rendered through
``wake_prompt`` — and its wiring into ``repl._wake_text``.

Fakes only: an injected ``run_fn`` stands in for squeue/sacct, and the watch
file lives under a tmp-path workspace. No real scheduler, no real sleep.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import job_wake
from delfin.agent import job_monitor as jm


# -- a fake scheduler backend -------------------------------------------------

_SCHEDULE = {}        # job_id -> terminal state | "RUNNING" | "" (aged out)
_MODE = "squeue"      # "squeue" answers; "down" makes every call return None
_RUNNING = "RUNNING"  # the one state squeue still lists


def _scheduler(cmd):
    """Fake scheduler answering exactly the commands the monitor issues
    (``squeue -j … -o '%i %T'`` lists jobs still in the queue; ``sacct``
    journals everything that has ended). A job whose squeue no longer lists
    it is asked of sacct; one neither knows is gone from accounting.

    Three answers, matching the monitor's tri-state contract: ``""`` when the
    scheduler answered but does not list the job (aged out), None when the
    scheduler could not be asked at all (down). Only the latter degrades.
    """
    if _MODE == "down":
        return None
    if cmd[0] == "squeue":
        lines = [f"{jid} {_RUNNING}" for jid, state in _SCHEDULE.items()
                 if state == _RUNNING]
        return "\n".join(lines)  # "" = answered, nothing running
    if cmd[0] == "sacct":
        lines = [f"{jid}   {state}" for jid, state in _SCHEDULE.items()
                 if state and state != _RUNNING]
        return "\n".join(lines)  # "" = answered, no finished jobs
    raise AssertionError(f"unexpected command: {cmd!r}")


@pytest.fixture(autouse=True)
def _reset(monkeypatch, tmp_path):
    global _SCHEDULE, _MODE
    _SCHEDULE = {}
    _MODE = "squeue"
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH",
                        tmp_path / "index.json")


@pytest.fixture
def ws(tmp_path):
    """The per-test workspace: a fresh tmp_path with the agent watch index
    redirected into it by the autouse _reset fixture."""
    return tmp_path


def _register(ws, jid, *, session_id, state, kind="slurm"):
    """Put a job on the agent watch file exactly the way a submission does."""
    jm.register_agent_job(ws, jid, description=f"job {jid}",
                          extra={"kind": kind, "session_id": session_id})
    _SCHEDULE[jid] = state


# -- the core gap: SLURM completion wakes the terminal that submitted it -----

def test_a_finished_slurm_job_wakes_the_terminal_session(ws):
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    out = job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler)
    assert len(out) == 1
    assert out[0]["kind"] == "slurm"
    assert out[0]["ok"] is True
    assert out[0]["job_id"] == "12345"


def test_the_result_is_rendered_in_delfins_own_words(ws):
    _register(ws, "12345", session_id="session-a", state="FAILED")
    out = job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler)
    text = job_wake.wake_prompt(out)
    assert "12345" in text and "FAILED" in text


# -- session scoping: only the owning session is woken -----------------------

def test_another_sessions_job_does_not_wake_this_session(ws):
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    assert job_wake.finished_watched_jobs(
        ws, set(), session_id="session-b", run_fn=_scheduler) == []


def test_a_job_with_no_session_still_reaches_its_session(ws):
    # A job from before watches had an owner is nobody's yet.
    _register(ws, "12345", session_id="", state="COMPLETED")
    assert job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler) == []


def test_without_a_session_nothing_is_filtered(ws):
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    _register(ws, "67890", session_id="session-b", state="COMPLETED")
    out = job_wake.finished_watched_jobs(ws, set(), run_fn=_scheduler)
    assert {d["job_id"] for d in out} == {"12345", "67890"}


# -- exactly once ------------------------------------------------------------

def test_a_finished_job_is_reported_exactly_once(ws):
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    seen: set = set()
    assert job_wake.finished_watched_jobs(
        ws, seen, session_id="session-a", run_fn=_scheduler) != []
    assert job_wake.finished_watched_jobs(
        ws, seen, session_id="session-a", run_fn=_scheduler) == []


def test_the_daemons_peek_is_not_taken_away(ws):
    # A consume=False + marker peek by the DAEMON (session_id=None, its own
    # marker) does NOT consume: it still leaves the completion for the owning
    # session's turn to report once. (The old body peeked as session-a and
    # then expected session-a's OWN turn to hear it again -- that same-session
    # double delivery was the V3-1 defect removed by the fix.)
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    assert jm.check_agent_jobs(ws, run_fn=_scheduler,
                               consume=False, marker="daemon_notified") != []
    assert jm.check_agent_jobs(ws, run_fn=_scheduler, session_id="session-a",
                               consume=True) != []


# -- negatives and degraded state -------------------------------------------

def test_a_running_job_is_not_reported(ws):
    _register(ws, "12345", session_id="session-a", state="RUNNING")
    assert job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler) == []


def test_nothing_watched_wakes_nobody(ws):
    assert job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler) == []


def test_an_unanswerable_scheduler_is_reported_as_degraded(ws):
    _register(ws, "12345", session_id="session-a", state="COMPLETED")
    global _MODE
    _MODE = "down"
    out = job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler)
    assert out and out[0].get("degraded"), "a down scheduler must degrade, not wait silently"


def test_a_job_that_left_squeue_still_reports_its_outcome(ws):
    # A finished job is no longer in squeue — sacct journals it. Asking only
    # squeue would miss it, so the outcome must still reach the terminal.
    _register(ws, "12345", session_id="session-a", state="FAILED")
    out = job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler)
    assert len(out) == 1 and out[0]["job_id"] == "12345" and not out[0]["ok"]


def test_a_job_gone_from_accounting_is_not_reported(ws):
    # Both squeue and sacct answer but do not list the id (""): aged out of
    # accounting. That is not a completion — nothing wakes; the entry is
    # left for the pruner. (Distinct from a down scheduler, which degrades.)
    _register(ws, "12345", session_id="session-a", state="")
    assert job_wake.finished_watched_jobs(
        ws, set(), session_id="session-a", run_fn=_scheduler) == []


# -- it does not open its own scheduler loop --------------------------------

def test_no_short_period_poll_loop_is_added():
    import ast, inspect
    src = pathlib.Path(inspect.getfile(job_wake)).read_text(encoding="utf-8")
    tree = ast.parse(src)
    sleeps = [n for n in ast.walk(tree)
              if isinstance(n, ast.Call)
              and isinstance(n.func, ast.Attribute) and n.func.attr == "sleep"]
    assert not sleeps, "the terminal wake must not add its own polling loop"


def test_status_is_sourced_from_the_watch_file_and_scheduler_only():
    import ast, inspect
    src = pathlib.Path(inspect.getfile(job_wake)).read_text(encoding="utf-8")
    # The one change the fix is allowed: a single check_agent_jobs read with
    # the wake marker, inside finished_watched_jobs. No bare 'squeue' string.
    assert "check_agent_jobs" in src
    assert "wake_notified" in src
    # No direct scheduler call of its own: no subprocess/shell/Popen suffice
    # to run squeue itself, only the injected run_fn handed to check_agent_jobs.
    tree = ast.parse(src)
    direct = [n for n in ast.walk(tree)
              if isinstance(n, (ast.Subscript, ast.Attribute))
              and getattr(n, "attr", "") in ("Popen", "system", "run", "squeue")]
    assert not direct, "finished_watched_jobs must not open its own scheduler call"


# -- the terminal actually calls it, session-scoped --------------------------

def test_the_terminal_wake_text_includes_slurm_watches():
    import pathlib
    root = pathlib.Path(job_wake.__file__).resolve().parents[2]
    src = (root / "delfin/agent/repl.py").read_text(encoding="utf-8")
    assert "finished_watched_jobs" in src
    at = src.index("finished_watched_jobs")
    assert "session_id=" in src[at:at + 200], (
        "repl._wake_text must scope finished_watched_jobs by session")
