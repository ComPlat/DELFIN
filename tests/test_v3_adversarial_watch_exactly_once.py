"""Adversarial tests: a watched job's completion is delivered exactly once.

Package V3 promises (phase 2): one completion record per finished
job/watch, delivered exactly once to the session that started it.

A watched job lives in ``<ws>/.delfin/agent_watched_jobs.json``. It is
read by two delivery surfaces that use *different* read semantics:

* the idle peek  -- ``check_agent_jobs(consume=False, marker="wake_notified")``
  (repl.py idle wake, tab_agent.py idle tick) -- KEEPS the entry and only
  sets a flag;
* the busy-turn  -- ``check_agent_jobs(consume=True)`` (engine.py
  ``_build_finished_jobs_block``) -- POPS the entry no matter what flags
  are set.

If a completion is first seen by the idle peek and then by the busy turn
of the SAME session, the current code delivers it twice. These tests pin
the "delivered exactly once" contract across that sequence.

No real scheduler: a fake ``run_fn`` returns a terminal state, like the
existing ``tests/test_job_monitor.py::_fake_run`` idiom.
"""

from __future__ import annotations

from delfin.agent import job_monitor as JM


def _fake_run(squeue: str = "", sacct: str = "") -> callable:
    def run(cmd):
        if cmd[0] == "squeue":
            return squeue
        if cmd[0] == "sacct":
            return sacct
        return ""
    return run


def _deliver_once_peek_then_drain(ws, *, session_id="session-a"):
    """The realistic sequence: idle peek first, then a busy turn drains.

    Returns (peek_deliveries, drain_deliveries). The contract is that the
    sum is 1 -- one completion reaches the session, exactly once."""
    run = _fake_run(sacct="111  COMPLETED\n")
    peek = JM.check_agent_jobs(ws, run_fn=run, consume=False,
                               marker="wake_notified",
                               session_id=session_id)
    drain = JM.check_agent_jobs(ws, run_fn=run, consume=True,
                                session_id=session_id)
    return peek, drain


def test_idle_peek_then_busy_drain_delivers_a_completion_once(tmp_path):
    JM.register_agent_job(tmp_path, "111", "prod run",
                          extra={"session_id": "session-a"})
    peek, drain = _deliver_once_peek_then_drain(tmp_path)

    # The completion reached the session...
    assert (len(peek) + len(drain)) == 1, (
        f"one completion delivered {len(peek) + len(drain)} times "
        f"(peek={len(peek)}, drain={len(drain)}) -- the idle peek and the "
        "busy turn both announced the same finished job")
    assert (peek + drain)[0]["job_id"] == "111"
    assert (peek + drain)[0]["state"] == "COMPLETED"


def test_idle_peek_alone_delivers_once(tmp_path):
    """Peek only (session goes idle and stays idle): one delivery, and the
    entry stays for the daemon -- no premature retirement."""
    JM.register_agent_job(tmp_path, "222", "opt freq run")
    run = _fake_run(sacct="222  COMPLETED\n")
    first = JM.check_agent_jobs(tmp_path, run_fn=run, consume=False,
                                marker="wake_notified")
    second = JM.check_agent_jobs(tmp_path, run_fn=run, consume=False,
                                 marker="wake_notified")
    assert len(first) == 1
    assert second == []          # same peek type is told once
    assert "222" in JM.load_watched(
        JM._agent_watch_path(tmp_path))["jobs"]   # still watched


def test_busy_drain_alone_delivers_once_and_retires(tmp_path):
    """Busy turn only: one delivery and the entry is retired (popped), so a
    restart or a later peek sees nothing."""
    JM.register_agent_job(tmp_path, "333", "bench")
    run = _fake_run(sacct="333  COMPLETED\n")
    first = JM.check_agent_jobs(tmp_path, run_fn=run, consume=True)
    second = JM.check_agent_jobs(tmp_path, run_fn=run, consume=True)
    assert len(first) == 1
    assert second == []          # popped, nothing left to report
