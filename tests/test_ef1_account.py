"""EF1 phase-1 red tests: cluster_efficiency.account().

Accounting from RECORDED sacct text only — no live cluster in tests.
The module routes its single sacct call through DELFIN's shared
scheduler client (delfin.scheduler_client), so a repeated call within
the interval re-serves the cache instead of hitting the runner.
"""
import pytest

from delfin.scheduler_client import SchedulerClient, THROTTLED


SACCT_COLUMNS = ("JobID,JobName,State,ExitCode,Elapsed,Timelimit,"
                 "ReqTRES,AllocCPUS,MaxRSS,Submit,Start,End")

# One completed 16-core, 16 GB job: 10 min queue wait, 2 h run.
SACCT_OK = (
    "1001|et1_opt|COMPLETED|0:0|02:00:00|24:00:00|"
    "billing=16,cpu=16,mem=16000M,node=1|16|9200M|"
    "2026-10-17T08:00:00|2026-10-17T08:10:00|2026-10-17T10:10:00\n"
)


def _fake_client(runner, now, interval=25.0):
    return SchedulerClient(min_interval_s=interval, clock=(lambda: now[0]),
                           runner=runner)


def _sacct_runner(text, calls):
    def run(argv):
        calls.append(list(argv))
        return text
    return run


def test_account_reports_requested_vs_used():
    from delfin.agent import cluster_efficiency as ce
    calls = []
    run = _sacct_runner(SACCT_OK, calls)
    report = ce.account(["1001"], run_fn=run)
    job = report["jobs"]["1001"]
    assert job["state"] == "COMPLETED"
    assert job["requested"]["cores"] == 16
    assert job["requested"]["mem_mb"] == 16000
    assert job["used"]["cores"] == 16
    assert job["used"]["mem_mb"] == 9200
    assert job["used"]["wall_s"] == 7200
    assert job["queue_wait_s"] == 600
    assert job["core_hours_used"] == pytest.approx(32.0)
    assert job["core_hours_wasted"] == pytest.approx(0.0)


def test_failed_job_counts_core_hours_as_waste():
    from delfin.agent import cluster_efficiency as ce
    text = (
        "2001|et2_opt|FAILED|1:0|01:00:00|24:00:00|"
        "billing=8,cpu=8,mem=8000M,node=1|8|1200M|"
        "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T09:00:00\n"
    )
    calls = []
    run = _sacct_runner(text, calls)
    job = ce.account(["2001"], run_fn=run)["jobs"]["2001"]
    assert job["state"] == "FAILED"
    assert job["core_hours_wasted"] == pytest.approx(8.0)
    assert job["core_hours_used"] == pytest.approx(0.0)


def test_repeated_job_failed_then_completed_counts_the_failed_run_as_waste():
    from delfin.agent import cluster_efficiency as ce
    text = (
        "3001|et3_opt|FAILED|1:0|01:00:00|24:00:00|"
        "billing=8,cpu=8,mem=8000M,node=1|8|300M|"
        "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T09:00:00\n"
        "3002|et3_opt|COMPLETED|0:0|03:00:00|24:00:00|"
        "billing=8,cpu=8,mem=8000M,node=1|8|7000M|"
        "2026-10-17T09:05:00|2026-10-17T09:05:00|2026-10-17T12:05:00\n"
    )
    report = ce.account(["3001", "3002"], run_fn=_sacct_runner(text, []))
    failed = report["jobs"]["3001"]
    ok = report["jobs"]["3002"]
    assert failed["is_repeated_waste"] is True
    assert failed["core_hours_wasted"] == pytest.approx(8.0)
    assert ok["is_repeated_waste"] is False
    # Summary-level waste counts the failed run, not the success.
    assert report["total"]["core_hours_wasted"] == pytest.approx(8.0)
    assert report["total"]["core_hours_used"] == pytest.approx(24.0)


def test_throttled_sacct_returns_none_never_zeros():
    from delfin.agent import cluster_efficiency as ce
    calls = []

    def run(argv):
        calls.append(list(argv))
        return THROTTLED

    assert ce.account(["4001"], run_fn=run) is None
    assert len(calls) == 1  # asked once through the gate, no retry


def test_repeated_call_within_interval_serves_cache_one_sacct_call():
    from delfin.agent import cluster_efficiency as ce
    calls = []
    run = _sacct_runner(SACCT_OK, calls)
    now = [1000.0]
    client = _fake_client(run, now)
    first = ce.account(["1001"], run_fn=run, client=client)
    now[0] += 5.0
    second = ce.account(["1001"], run_fn=run, client=client)
    assert len(calls) == 1
    assert second["jobs"]["1001"]["core_hours_used"] == (
        first["jobs"]["1001"]["core_hours_used"])


def test_scheduler_unavailable_is_reported_not_fabricated():
    from delfin.agent import cluster_efficiency as ce
    calls = []

    def run(argv):
        calls.append(list(argv))
        return None

    assert ce.account(["5001"], run_fn=run) is None


def test_one_sacct_call_for_many_ids_not_a_loop():
    from delfin.agent import cluster_efficiency as ce
    calls = []
    text = SACCT_OK.replace("et1_opt", "et4_opt")
    report = ce.account(["1001", "1002", "1003"],
                        run_fn=_sacct_runner(text, calls))
    assert len(calls) == 1
    assert set(report["jobs"]) == {"1001", "1002", "1003"}
