"""BG1 phase 2: the estimate learns from EF1's real sacct records.

After a job finishes, the agent records estimate vs actual (EF1 core-hours)
so that the next estimate for the same method family is corrected.  The
MEASURED value must come from the real writer -- ``cluster_efficiency.account``
with a fake sacct runner producing a genuine ``sacct -P`` line -- never a
dict typed by hand.  On the unchanged tree ``record_outcome`` does not exist,
so every test here is red (AttributeError) until BG1 phase 2 lands.
"""

import pytest

from delfin.agent import calc_budget
from delfin.agent.cluster_efficiency import account


@pytest.fixture(autouse=True)
def _clear_learned():
    """Keep the learned-factor registry per-test: the learning tests are
    order-independent (one test's recorded EF1 factor must not leak into
    another's fresh-family assertion)."""
    calc_budget._LEARNED.clear()
    yield
    calc_budget._LEARNED.clear()


# A real ``sacct -P -X`` row for a completed 2-core job that ran one wall
# hour:  elapsed 01:00:00, AllocCPUS 2  =>  core_hours_used 2.0.
_SACCT_OK = (
    "12345|calc_small|COMPLETED|0:0|01:00:00|02:00:00|cpu=1.00|2|1048576K|"
    "2026-11-07T12:00:00|2026-11-07T12:10:00|2026-11-07T13:10:00|02:00:00\n"
)
# A row whose run state is UNKNOWN / incomplete: it must NOT steer the model.
_SACCT_UNKNOWN = (
    "12346|calc_unknown|UNKNOWN|||02:00:00|cpu=1.00|2|1048576K|"
    "2026-11-07T12:00:00|2026-11-07T12:00:00|2026-11-07T12:30:00|00:30:00\n"
)


def _run_sacct(text: str):
    """A fake ``run_fn`` so ``account`` parses real sacct text, not shell."""
    def runner(argv):
        return text
    return runner


def _real_record(text: str):
    """The real-shaped EF1 record, produced by the real writer."""
    # account() keys its output by the IDS PASSED IN (the sacct text's own
    # JobID is the match key), so request the id that the fixture line carries.
    jid = text.split("|", 1)[0].strip()
    report = account([jid], run_fn=_run_sacct(text))
    assert report is not None and report["jobs"], (
        "sacct fixture must parse into at least one job record")
    return next(iter(report["jobs"].values()))


def _fake_sbatch(delay_s=5):
    """A hermetic sbatch --test-only runner (real format) -- never a real call."""
    from datetime import datetime, timedelta
    from delfin.scheduler_client import SchedulerClient
    t = (datetime.now() + timedelta(seconds=delay_s)).strftime("%Y-%m-%dT%H:%M:%S")
    return SchedulerClient(runner=lambda argv: f"sbatch: Job to start at {t}")


def _pbe0_estimate(client=None):
    """A small PBE0 estimate (family ``dft``); core-hours > 0 and hermetic."""
    return calc_budget.estimate({
        "n_atoms": 6, "n_electrons": 20, "basis_functions": 80,
        "method": "PBE0", "partition": "cpu_short", "pal": 2, "maxcore": 6000,
    }, client=client if client is not None else _fake_sbatch())


def test_learn_scales_next_estimate_for_the_same_family():
    # A real EF1 record says the job really used 2.0 core-hours.
    rep = _real_record(_SACCT_OK)
    est_before = _pbe0_estimate()
    calc_budget.record_outcome(rep, est_before)
    est_after = _pbe0_estimate()
    # 2.0 actual vs est_before.core_hours -- the learned factor must apply.
    assert est_after.core_hours > 0
    assert est_after.core_hours != est_before.core_hours
    assert "learned" not in est_before.source
    assert "learned" in est_after.source


def test_unmeasured_job_does_not_steer_the_model():
    # A record without a measured core-hours value (state UNKNOWN) is not a
    # zero -- it must be skipped, leaving the model unchanged.
    rep = _real_record(_SACCT_UNKNOWN)
    est = _pbe0_estimate()
    calc_budget.record_outcome({"job_id": rep["job_id"], "core_hours_used": None,
                                "state": rep["state"]}, est)
    assert est.core_hours == calc_budget.estimate({
        "n_atoms": 6, "n_electrons": 20, "basis_functions": 80,
        "method": "PBE0", "partition": "cpu_short", "pal": 2, "maxcore": 6000,
    }).core_hours


def test_a_fresh_family_starts_unscaled():
    # Before any record for a family, the estimate is the plain model.
    assert calc_budget.learned_scale("never_seen") is None
    est = _pbe0_estimate()
    assert "learned" not in est.source
