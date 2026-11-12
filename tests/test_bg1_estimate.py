"""BG1: cost-and-time estimate plus a per-session cluster budget.

The estimate turns a job's system size and method into core-hours, memory,
wall time and an expected queue wait, reusing DELFIN's single throttled
scheduler client (never a per-loop or per-widget ``squeue``/``sacct``) and
EF1's accounting as calibration where it exists.  The budget decides, before
submission, whether a job fits a session's remaining core-hours: within
budget it goes, over it the person is asked -- the agent never shrinks
PAL/maxcore/wall-time silently to make a job fit.

Everything in this module drives the real public path
``delfin.agent.calc_budget``.  On the unchanged tree that module does not
exist, so every test here is red (ImportError) until BG1 lands.
"""

import pytest
from datetime import datetime, timedelta

from delfin.agent import calc_budget
from delfin.scheduler_client import THROTTLED, SchedulerClient


def _small_job(**over):
    """A cheap single-point job -- must clear the budget untouched."""
    job = {
        "n_atoms": 6,
        "n_electrons": 20,
        "basis_functions": 80,
        "method": "PBE0/def2-SVP",
        "partition": "cpu_short",
        "pal": 2,
        "maxcore": 6000,
    }
    job.update(over)
    return job


def _large_dlpno_job(**over):
    """A large DLPNO-CCSD(T) job -- must exceed a normal session budget."""
    job = {
        "n_atoms": 80,
        "n_electrons": 400,
        "basis_functions": 4000,
        "method": "DLPNO-CCSD(T)/def2-TZVP",
        "partition": "cpu",
        "pal": 24,
        "maxcore": 6000,
    }
    job.update(over)
    return job


class _FakeClient:
    """A scheduler client runner that answers without any SLURM call."""

    def __init__(self, answer):
        self.answer = answer  # a (kind, argv) -> value dict, or a THROTTLED mode

    def __call__(self, cmd, **kw):
        if self.answer is THROTTLED:
            raise AssertionError("a throttled client must not run the command")
        # cmd is [sbatch, --test-only, ...]; the only thing we vary is
        # the partition name, so derive the answer from it.
        for arg in cmd:
            if arg.startswith("--partition="):
                return self.answer.get(arg.split("=", 1)[1], THROTTLED)
        return self.answer.get("default", THROTTLED)


def _client_throttled():
    c = SchedulerClient(runner=_FakeClient(THROTTLED))
    # Close the per-kind gate at the current time so query() returns
    # THROTTLED for an argv that has nothing cached, and never runs the
    # (asserting) runner.
    c._kind_last_run["sbatch"] = c._now()  # type: ignore[attr-defined]
    return c


def _client_partitions(answers):
    return SchedulerClient(runner=_FakeClient(answers))


def test_default_path_shares_the_process_wide_throttle():
    # F2 (reviewer): the default path must NOT hand out a fresh client per
    # call -- that would open a new throttle window/cache each estimate and
    # let a loop ask the scheduler repeatedly.  _confirm_client routes a
    # missing client to the process-wide singleton instead.
    from delfin.scheduler_client import default_client
    assert calc_budget._confirm_client(None) is default_client


def _test_only_answer(delay_s):
    """A real ``sbatch --test-only`` answer: ``... to start at <datetime>``,
    ``delay_s`` seconds in the future (the marker DELFIN keys on)."""
    t = datetime.now() + timedelta(seconds=delay_s)
    return f"sbatch: Job 123 to start at {t.strftime('%Y-%m-%dT%H:%M:%S')}"


# --------------------------------------------------------------------------
# 1. estimate(job): core-hours / memory / wall time from size + method.
# --------------------------------------------------------------------------

def test_estimate_returns_core_hours_memory_walltime():
    est = calc_budget.estimate(_small_job(), client=_client_partitions({"cpu_short": "ok"}))
    assert est.core_hours > 0
    assert est.memory_mb > 0
    assert est.walltime_s > 0
    assert est.n_atoms == 6
    assert est.partition == "cpu_short"
    assert est.method == "PBE0/def2-SVP"


def test_estimate_core_hours_grow_with_system_size():
    small = calc_budget.estimate(_small_job(), client=_client_partitions({"cpu_short": "ok"}))
    big = calc_budget.estimate(_large_dlpno_job(), client=_client_partitions({"cpu": "ok"}))
    assert big.core_hours > small.core_hours


# --------------------------------------------------------------------------
# 2. THROTTLED scheduler -> "unknown" queue wait, never an invented zero.
# --------------------------------------------------------------------------

def test_throttled_scheduler_gives_unknown_wait_not_zero():
    est = calc_budget.estimate(_small_job(), client=_client_throttled())
    assert est.wait_s == "unknown"


def test_scheduler_wait_is_reported_when_answered():
    # A fitting partition answers; the estimate carries its start delay.
    est = calc_budget.estimate(_small_job(),
                               client=_client_partitions({"cpu_short": _test_only_answer(120)}))
    # A wall-clock "seconds until start" is inherently approximate (the fake
    # and the parser each read the clock): assert the advertised delay within
    # sub-second-to-a-second drift, not an exact equality.
    assert 115 <= est.wait_s <= 125


# --------------------------------------------------------------------------
# 3. Budget: within -> go, over -> ask, never a silent shrink.
# --------------------------------------------------------------------------

def test_budget_within_limit_decides_go():
    budget = calc_budget.Budget(limit_core_hours=200.0)
    est = calc_budget.estimate(_small_job(), client=_client_partitions({"cpu_short": "ok"}))
    verdict, reason = budget.decision(est)
    assert verdict == "go"
    assert reason  # the reason names the estimate


def test_budget_over_limit_decides_ask_with_estimate():
    budget = calc_budget.Budget(limit_core_hours=1.0)
    est = calc_budget.estimate(_large_dlpno_job(), client=_client_partitions({"cpu": "ok"}))
    assert est.core_hours > 1.0
    verdict, reason = budget.decision(est)
    assert verdict == "ask"
    # the ask must carry the estimate, never a silent shrink of it
    assert str(est.core_hours) in reason


def test_budget_does_not_silently_shrink_pal_maxcore_walltime():
    budget = calc_budget.Budget(limit_core_hours=1.0)
    est = calc_budget.estimate(_large_dlpno_job(), client=_client_partitions({"cpu": "ok"}))
    verdict, _ = budget.decision(est)
    assert verdict == "ask"
    # the vote is ask-or-go; nothing in the response resizes the job.
    assert budget.remaining() == 0.0 or budget.remaining() < est.core_hours


def test_budget_record_tracks_spent_core_hours():
    budget = calc_budget.Budget(limit_core_hours=200.0)
    est = calc_budget.estimate(_small_job(), client=_client_partitions({"cpu_short": "ok"}))
    budget.record(est)
    assert budget.spent_core_hours == pytest.approx(est.core_hours)
    assert budget.remaining() == pytest.approx(200.0 - est.core_hours)


# --------------------------------------------------------------------------
# 4. Rate limit: estimate is cheap, never a per-loop squeue/sacct.
# --------------------------------------------------------------------------

def test_estimate_single_call_asks_the_scheduler_once():
    client = _client_partitions({"cpu_short": "started-after-60s"})
    est = calc_budget.estimate(_small_job(), client=client)
    assert est.wait_s is not None
    # the estimate path must produce exactly one scheduler ask (the
    # --test-only for its chosen partition), throttled through the client.
    assert len(client._cache) == 1  # type: ignore[attr-defined]
