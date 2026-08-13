"""Near the end of a walltime, splitting the pool can finish nothing at all.

The archived TADF job had two jobs ready — an excited-state optimisation each
for S1 and T1 — and both were on the critical path for the nine jobs behind
them. The scheduler gave each 20 of 40 cores, and 60 h later neither had
finished, so nothing downstream ever started.

Running them one after another at full width would have finished one. These
tests pin when that trade is taken and, just as importantly, when it is not:
splitting is the better default, so serialising may only happen when both the
remaining walltime and a runtime estimate are actually known.

Also covered is the walltime parsing itself, because a wrong number here would
make the scheduler refuse jobs that do fit.
"""

from types import SimpleNamespace

import pytest

from delfin.cluster_utils import _parse_walltime
from delfin.workflows.engine.classic import _WorkflowManager


def _job(job_id: str, *, cores_min=10, cores_optimal=10, cores_max=40):
    return SimpleNamespace(
        job_id=job_id,
        cores_min=cores_min,
        cores_optimal=cores_optimal,
        cores_max=cores_max,
    )


def _manager(hints: dict, monkeypatch, remaining):
    """A stand-in exposing only what the decision actually reads."""
    monkeypatch.setattr(
        "delfin.workflows.engine.classic.get_walltime_limit", lambda: remaining
    )
    return SimpleNamespace(
        label="classic",
        _get_duration_hint=lambda job: hints.get(job.job_id),
        _scaled_duration=_WorkflowManager._scaled_duration,
    )


def _decide(manager, jobs, total=40, available=40):
    return _WorkflowManager._serialize_for_walltime(manager, jobs, total, available)


# --------------------------------------------------------------------------
# when to serialise
# --------------------------------------------------------------------------

def test_two_jobs_that_only_fit_alone_are_run_one_at_a_time(monkeypatch):
    """Both measured at 30 h on 10 cores. Split over 20 cores each that is
    15 h and neither lands inside the 12 h left; alone on 40 it is 7.5 h and
    one does."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    hints = {"esd_S1": 30 * 3600, "esd_T1": 30 * 3600}  # measured at 10 cores
    manager = _manager(hints, monkeypatch, remaining=12 * 3600)

    decision = _decide(manager, jobs)

    assert decision is not None, "one of them fits alone — that beats neither"
    chosen, cores = decision
    assert chosen.job_id in {"esd_S1", "esd_T1"}
    assert cores == 40, "the chosen job gets the whole pool"


def test_a_split_that_still_lands_is_left_alone(monkeypatch):
    """The boundary case: 20 h on 10 cores is 10 h split over 20, which fits
    in 12 h. Serialising here would be a pure loss."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    hints = {"esd_S1": 20 * 3600, "esd_T1": 20 * 3600}
    manager = _manager(hints, monkeypatch, remaining=12 * 3600)

    assert _decide(manager, jobs) is None


def test_the_job_that_lands_soonest_is_the_one_chosen(monkeypatch):
    """Neither fits split (20 h and 16 h against 12 h left); alone they need
    10 h and 8 h, so the 8 h one is the one with room to spare."""
    jobs = [_job("slow"), _job("quick")]
    hints = {"slow": 40 * 3600, "quick": 32 * 3600}
    manager = _manager(hints, monkeypatch, remaining=12 * 3600)

    decision = _decide(manager, jobs)

    assert decision is not None
    assert decision[0].job_id == "quick"


def test_plenty_of_walltime_keeps_the_jobs_parallel(monkeypatch):
    """Splitting is the better default and must stay untouched."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    hints = {"esd_S1": 2 * 3600, "esd_T1": 2 * 3600}
    manager = _manager(hints, monkeypatch, remaining=48 * 3600)

    assert _decide(manager, jobs) is None


def test_nothing_fitting_either_way_keeps_the_jobs_parallel(monkeypatch):
    """Partial progress is resumable now, so spreading it across both jobs is
    worth more than forcing one that will not finish either."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    hints = {"esd_S1": 500 * 3600, "esd_T1": 500 * 3600}
    manager = _manager(hints, monkeypatch, remaining=12 * 3600)

    assert _decide(manager, jobs) is None


# --------------------------------------------------------------------------
# when the scheduler must not interfere
# --------------------------------------------------------------------------

@pytest.mark.parametrize("remaining", [None, 0])
def test_an_unknown_walltime_changes_nothing(monkeypatch, remaining):
    """Off a cluster, or when squeue cannot answer, the previous behaviour has
    to survive exactly."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    hints = {"esd_S1": 10 * 3600, "esd_T1": 10 * 3600}
    manager = _manager(hints, monkeypatch, remaining=remaining)

    assert _decide(manager, jobs) is None


def test_a_missing_runtime_estimate_changes_nothing(monkeypatch):
    """A first-time job has no history. Guessing would be worse than splitting."""
    jobs = [_job("esd_S1"), _job("esd_T1")]
    manager = _manager({"esd_S1": 10 * 3600}, monkeypatch, remaining=12 * 3600)

    assert _decide(manager, jobs) is None


def test_a_single_ready_job_is_not_this_decision(monkeypatch):
    manager = _manager({"esd_S1": 10 * 3600}, monkeypatch, remaining=1 * 3600)

    assert _decide(manager, [_job("esd_S1")]) is None


def test_a_walltime_lookup_that_raises_changes_nothing(monkeypatch):
    """Scheduling must never fail because a cluster query did."""
    def boom():
        raise OSError("squeue unavailable")

    monkeypatch.setattr("delfin.workflows.engine.classic.get_walltime_limit", boom)
    manager = SimpleNamespace(
        label="classic",
        _get_duration_hint=lambda job: 3600.0,
        _scaled_duration=_WorkflowManager._scaled_duration,
    )

    assert _decide(manager, [_job("a"), _job("b")]) is None


def test_the_core_projection_is_optimistic_by_design(monkeypatch):
    """Doubling the cores is assumed to halve the runtime. Real speed-up is
    worse, so the projection serialises less often than reality warrants —
    the safe direction for a change to scheduling."""
    assert _WorkflowManager._scaled_duration(3600.0, 10, 20) == 1800.0
    assert _WorkflowManager._scaled_duration(3600.0, 10, 10) == 3600.0
    assert _WorkflowManager._scaled_duration(3600.0, 10, 0) == 3600.0


# --------------------------------------------------------------------------
# walltime parsing
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "text,seconds",
    [
        ("2-12:00:00", 216000),   # SLURM TimeLimit for the archived job
        ("1-06:30:00", 109800),
        ("11:54:52", 42892),      # RunTime without a day part
        ("29:30", 1770),          # squeue %L drops to MM:SS near the end
        ("3-04:15", 274500),      # D-HH:MM
        ("48", 172800),           # PBS bare hours
    ],
)
def test_every_slurm_time_format_parses(text, seconds):
    assert _parse_walltime(text) == seconds


@pytest.mark.parametrize(
    "text", ["UNLIMITED", "NOT_SET", "INVALID", "N/A", "", "   ", "garbage", "1-xx:00:00"]
)
def test_no_limit_never_reads_as_no_time_left(text):
    """None means unknown. Returning 0 here would make the scheduler treat an
    unlimited job as one that is out of time."""
    assert _parse_walltime(text) is None
