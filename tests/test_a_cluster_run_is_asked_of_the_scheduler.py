"""A cluster run's state was never asked of SLURM, and cancel never checked.

``Runtime.submit_application(backend="slurm")`` records the job id and sets
the run PENDING. The only writer that could move it on was the job itself,
from the compute node. So every way a batch job can end WITHOUT running its
payload -- OUT_OF_MEMORY, TIMEOUT, NODE_FAIL, an admin ``scancel``, a held
job that never starts -- left the record PENDING for ever: no timeout, no
expiry, no sweep, and nothing anywhere that asked ``squeue`` or ``sacct``.
On a shared node the run store is the thing that says whether an allocation
is still being spent, and a record that can only ever say "pending" says
nothing at all.

``cancel()`` was worse than silent. A missing ``scancel`` was skipped and
the record marked CANCELLED anyway; a present one had its return code
discarded and the record marked CANCELLED anyway. Because CANCELLED is
terminal, nothing would ever look at that record again -- while the job ran
its full allocation on a node other people were queued for, and the caller
had been told True.

And ``RunHandle.wait()`` defaulted to no timeout with a 0.05 s poll: a 20 Hz
busy-wait on the login node that, for exactly the records above, could never
return.

No real queue is touched here -- the scheduler probe and ``scancel`` are
injected.
"""

from __future__ import annotations

import time

import pytest

from delfin.tools import _runtime as rt


@pytest.fixture
def store(tmp_path) -> rt.RunStore:
    return rt.RunStore(tmp_path / "runs")


def _slurm_record(store: rt.RunStore, job_id: str = "4242",
                  status: str = rt.RunStatus.PENDING.value) -> rt.RunRecord:
    rec = rt.RunRecord(
        id="run1", kind="application", name="opt", created_at=rt._now(),
        status=status, work_dir="/tmp/wd",
        metrics={"cores": 32, "backend": "slurm", "slurm_job_id": job_id},
    )
    store.save(rec)
    return rec


def _answer(state, monkeypatch, calls: list | None = None):
    def _query(job_id):
        if calls is not None:
            calls.append(job_id)
        return state
    monkeypatch.setattr(rt, "_slurm_query", _query)


# ---------------------------------------------------------------------------
# Reading a run asks the scheduler
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("state,expected_status", [
    ("OUT_OF_MEMORY", rt.RunStatus.FAILED.value),
    ("TIMEOUT", rt.RunStatus.FAILED.value),
    ("NODE_FAIL", rt.RunStatus.FAILED.value),
    ("FAILED", rt.RunStatus.FAILED.value),
    ("PREEMPTED", rt.RunStatus.FAILED.value),
    ("CANCELLED", rt.RunStatus.CANCELLED.value),
])
def test_a_job_that_died_in_the_queue_stops_reading_as_pending(
        store, monkeypatch, state, expected_status):
    """Every one of these used to be indistinguishable from "still queued"."""
    _slurm_record(store)
    _answer(state, monkeypatch)
    runtime = rt.Runtime(store)

    rec = runtime.get("run1")

    assert rec.status == expected_status
    assert rec.done is True
    assert rec.finished_at, "a terminal record must say when it ended"


def test_the_failure_names_the_slurm_state(store, monkeypatch):
    """'It failed' is not actionable; OUT_OF_MEMORY tells you what to change."""
    _slurm_record(store)
    _answer("OUT_OF_MEMORY", monkeypatch)

    rec = rt.Runtime(store).get("run1")

    assert "OUT_OF_MEMORY" in rec.error
    assert "4242" in rec.error


def test_a_started_job_advances_from_pending_to_running(store, monkeypatch):
    _slurm_record(store)
    _answer("RUNNING", monkeypatch)

    rec = rt.Runtime(store).get("run1")

    assert rec.status == rt.RunStatus.RUNNING.value
    assert rec.done is False


def test_a_queued_job_stays_pending(store, monkeypatch):
    _slurm_record(store)
    _answer("PENDING", monkeypatch)

    rec = rt.Runtime(store).get("run1")

    assert rec.status == rt.RunStatus.PENDING.value


def test_completed_without_a_result_is_a_failure_not_a_success(
        store, monkeypatch):
    """The compute-node executor writes the outcome. If SLURM says the job is
    over and no outcome was written, the payload never ran or died before
    saving -- reporting success would invent outputs that do not exist."""
    _slurm_record(store)
    _answer("COMPLETED", monkeypatch)

    rec = rt.Runtime(store).get("run1")

    assert rec.status == rt.RunStatus.FAILED.value
    assert rec.outputs == {}
    assert "COMPLETED" in rec.error
    assert "without the run recording a result" in rec.error


def test_a_result_the_compute_node_wrote_is_never_overwritten(
        store, monkeypatch):
    """The job's own record of what happened always wins over our inference."""
    rec = _slurm_record(store, status=rt.RunStatus.SUCCESS.value)
    rec.outputs = {"energy_Eh": -76.4}
    store.save(rec)
    _answer("TIMEOUT", monkeypatch)

    got = rt.Runtime(store).get("run1")

    assert got.status == rt.RunStatus.SUCCESS.value
    assert got.outputs == {"energy_Eh": -76.4}


def test_a_local_run_never_asks_the_scheduler(store, monkeypatch):
    """No job id, no query -- reading a local run must stay free."""
    store.save(rt.RunRecord(id="local1", kind="application", name="opt",
                            created_at=rt._now(), metrics={"backend": "local"}))
    calls: list = []
    _answer("FAILED", monkeypatch, calls)

    rec = rt.Runtime(store).get("local1")

    assert calls == []
    assert rec.status == rt.RunStatus.PENDING.value


def test_listing_runs_reconciles_too(store, monkeypatch):
    _slurm_record(store)
    _answer("NODE_FAIL", monkeypatch)

    runs = rt.Runtime(store).list_runs()

    assert [r.status for r in runs] == [rt.RunStatus.FAILED.value]


# ---------------------------------------------------------------------------
# "Could not be asked" is not "nothing changed"
# ---------------------------------------------------------------------------

def test_an_unreachable_scheduler_is_reported_as_unknown_not_running(
        store, monkeypatch):
    """The whole point of the tri-state. A record whose scheduler could not
    be reached is unverified, and a caller that treats that as "still going"
    waits for ever on a job that was killed hours ago."""
    _slurm_record(store)
    _answer(None, monkeypatch)
    runtime = rt.Runtime(store)
    handle = rt.RunHandle("run1", runtime)

    state = handle.state()

    assert state["known"] is False
    assert state["status"] == rt.RunStatus.PENDING.value
    assert state["slurm_state"] == "unavailable"
    assert "could not be asked" in state["detail"]


def test_a_reachable_scheduler_reports_a_known_state(store, monkeypatch):
    _slurm_record(store)
    _answer("RUNNING", monkeypatch)

    state = rt.RunHandle("run1", rt.Runtime(store)).state()

    assert state["known"] is True
    assert state["status"] == rt.RunStatus.RUNNING.value


def test_an_id_the_scheduler_does_not_know_closes_only_after_the_grace(
        store, monkeypatch):
    """Accounting can lag a fresh submission by seconds. A job still missing
    a quarter of an hour later is not going to appear."""
    _slurm_record(store)
    _answer("", monkeypatch)
    runtime = rt.Runtime(store)

    rec = runtime.get("run1")
    assert rec.status == rt.RunStatus.PENDING.value      # inside the grace
    assert rec.metrics["slurm_state"] == "absent"

    rec.metrics["slurm_absent_since"] = time.time() - rt._SLURM_ABSENT_GRACE_S - 1
    store.save(rec)

    rec = runtime.get("run1")
    assert rec.status == rt.RunStatus.FAILED.value
    assert "unknown to the scheduler" in rec.error


# ---------------------------------------------------------------------------
# cancel() reports what actually happened
# ---------------------------------------------------------------------------

def test_cancel_without_scancel_returns_false_and_touches_nothing(
        store, monkeypatch):
    """A 32-core job kept running while its record said CANCELLED, so nothing
    would ever look at it again."""
    _slurm_record(store)
    monkeypatch.setattr("shutil.which", lambda name: None)

    ok = rt.Runtime(store).cancel("run1")

    assert ok is False
    assert store.get("run1").status == rt.RunStatus.PENDING.value


def test_a_failing_scancel_returns_false_and_touches_nothing(
        store, monkeypatch):
    _slurm_record(store)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/scancel")

    class _Proc:
        returncode = 1
        stdout = ""
        stderr = "Invalid job id specified"

    monkeypatch.setattr("subprocess.run", lambda *a, **kw: _Proc())

    ok = rt.Runtime(store).cancel("run1")

    assert ok is False
    assert store.get("run1").status == rt.RunStatus.PENDING.value


def test_a_successful_scancel_is_confirmed_by_a_state_read(
        store, monkeypatch):
    """scancel exiting 0 means the request was accepted, not that the job is
    over. The terminal status is written from what the scheduler then says."""
    _slurm_record(store)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/scancel")

    class _Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    monkeypatch.setattr("subprocess.run", lambda *a, **kw: _Proc())
    _answer("CANCELLED", monkeypatch)

    ok = rt.Runtime(store).cancel("run1")

    assert ok is True
    assert store.get("run1").status == rt.RunStatus.CANCELLED.value


def test_a_cancel_the_scheduler_has_not_applied_yet_stays_non_terminal(
        store, monkeypatch):
    """Better a record that says "still going" and gets reconciled next read
    than a terminal one that closes the file on a live job."""
    _slurm_record(store)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/scancel")

    class _Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    monkeypatch.setattr("subprocess.run", lambda *a, **kw: _Proc())
    _answer("RUNNING", monkeypatch)

    ok = rt.Runtime(store).cancel("run1")

    assert ok is True
    rec = store.get("run1")
    assert rec.done is False
    assert any(e.get("event") == "slurm_cancel_requested" for e in rec.events)


# ---------------------------------------------------------------------------
# wait() is bounded and does not spin
# ---------------------------------------------------------------------------

def test_wait_has_a_default_bound_instead_of_looping_for_ever(
        store, monkeypatch):
    _slurm_record(store)
    _answer("PENDING", monkeypatch)
    monkeypatch.setattr(rt, "_DEFAULT_WAIT_TIMEOUT_S", 0.2)
    handle = rt.RunHandle("run1", rt.Runtime(store))

    t0 = time.monotonic()
    rec = handle.wait()                       # no timeout argument at all
    elapsed = time.monotonic() - t0

    assert elapsed < 5
    assert rec is not None and rec.done is False


def test_an_expired_wait_can_be_told_from_a_finished_one(store, monkeypatch):
    _slurm_record(store)
    _answer(None, monkeypatch)
    monkeypatch.setattr(rt, "_DEFAULT_WAIT_TIMEOUT_S", 0.2)
    handle = rt.RunHandle("run1", rt.Runtime(store))

    rec = handle.wait()

    assert rec.done is False
    assert handle.state()["known"] is False


def test_waiting_backs_off_instead_of_polling_at_twenty_hertz(
        store, monkeypatch):
    """0.05 s with no bound is a busy core on the login node, held beside the
    queue it is waiting on."""
    _slurm_record(store)
    _answer("PENDING", monkeypatch)
    monkeypatch.setattr(rt, "_DEFAULT_WAIT_TIMEOUT_S", 1.0)
    monkeypatch.setattr(rt, "_WAIT_POLL_MAX_S", 0.08)
    sleeps: list[float] = []
    real_sleep = time.sleep

    def _sleep(seconds):
        sleeps.append(seconds)
        real_sleep(0)

    monkeypatch.setattr(time, "sleep", _sleep)
    rt.RunHandle("run1", rt.Runtime(store)).wait(poll=0.01)

    assert len(sleeps) > 2, "wait must actually sleep between polls"
    assert sleeps[1] > sleeps[0], "the poll interval must grow, not stay flat"
    assert max(sleeps) <= 0.08, "and it must stop growing at the ceiling"
    # The default first poll is a second-scale wait, not the 0.05 s that
    # made an unbounded wait a 20 Hz spin.
    assert rt._WAIT_POLL_START_S >= 0.5
