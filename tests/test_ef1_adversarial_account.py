"""EF1 phase-1 reviewer adversarial tests for cluster_efficiency.account().

Written by nacht-s2 (reviewer) BEFORE any of these behaviours is fixed;
every test here encodes a defect the reviewer found in 49eae57c.
Recorded sacct text only — no live cluster.
"""
import pytest

from delfin.agent import cluster_efficiency as ce
from delfin.scheduler_client import SchedulerClient, THROTTLED

SACCT_COLUMNS = ("JobID,JobName,State,ExitCode,Elapsed,Timelimit,"
                 "ReqTRES,AllocCPUS,MaxRSS,Submit,Start,End,TotalCPU")


def _sacct_runner(text, calls):
    def run(argv):
        calls.append(list(argv))
        return text
    return run


# ---------------------------------------------------------------- D1
# elided line:
# 2026-10-18T14:00:00  (submit field with a date separator instead of T)
# The brief says queue_wait from Start-Submit; a wrong parse must stay
# None, never a fabricated number.
def test_unparseable_submit_yields_queue_wait_none():
    text = (
        "6001|et6_opt|COMPLETED|0:0|00:30:00|24:00:00|"
        "billing=4,cpu=4,mem=4000M,node=1|4|100M|"
        "2026-10-18 14:00:00|2026-10-18T14:00:00|2026-10-18T14:30:00\n"
    )
    job = ce.account(["6001"], run_fn=_sacct_runner(text, []))["jobs"]["6001"]
    assert job["queue_wait_s"] is None


# ---------------------------------------------------------------- D2
# elided line:
# 2026-10-17T25:61:00
def test_impossible_iso_timestamp_is_none():
    text = (
        "6002|et7_opt|COMPLETED|0:0|00:30:00|24:00:00|"
        "billing=4,cpu=4,mem=4000M,node=1|4|100M|"
        "2026-10-17T25:61:00|2026-10-17T08:10:00|2026-10-17T08:40:00\n"
    )
    job = ce.account(["6002"], run_fn=_sacct_runner(text, []))["jobs"]["6002"]
    assert job["queue_wait_s"] is None


# ---------------------------------------------------------------- D3
# A sacct field can be empty. None-coalescing: no zeros invented.
def test_empty_fields_stay_none_not_zero():
    text = ("6003|et8_opt|RUNNING|0:0||||16||"
            "2026-10-17T08:00:00|2026-10-17T08:05:00|Unknown\n")
    job = ce.account(["6003"], run_fn=_sacct_runner(text, []))["jobs"]["6003"]
    assert job["requested"]["cores"] == 16  # from AllocCPUS fallback
    assert job["requested"]["mem_mb"] is None
    assert job["requested"]["wall_s"] is None
    assert job["used"]["mem_mb"] is None
    assert job["core_hours_used"] == pytest.approx(0.0)  # elapsed unknown -> 0


# ---------------------------------------------------------------- D4
# A name that appears in rows for BOTH ids: account(["1001","1002"])
# with the same job_name must not let the repeat check read rows of the
# other id as a later run under the FIRST id's own rows (cross-id scan
# is over rows.items(), so it must key correctly).
def test_repeated_check_scans_other_ids_not_own_rows():
    text = (
        "7001|same_name|FAILED|1:0|01:00:00|24:00:00|"
        "billing=8,cpu=8,mem=8000M,node=1|8|300M|"
        "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T09:00:00\n"
        "7002|same_name|COMPLETED|0:0|03:00:00|24:00:00|"
        "billing=8,cpu=8,mem=8000M,node=1|8|7000M|"
        "2026-10-17T09:05:00|2026-10-17T09:05:00|2026-10-17T12:05:00\n"
    )
    report = ce.account(["7001", "7002"],
                        run_fn=_sacct_runner(text, []))
    # the failed run of the SAME id can never be its own "later run"
    assert report["jobs"]["7001"]["is_repeated_waste"] is True
    assert report["jobs"]["7002"]["is_repeated_waste"] is False


# ---------------------------------------------------------------- D5
# Step rows ("1001.batch") are filtered by the .split(".")[0] key, but
# an -X output for an array job id "1001_3" keys by array element. The
# split is on "." only — array sub-jobs stay separate ids, which is
# correct. What must hold: an id with ONLY step rows still reports.
def test_id_keyed_by_array_element_not_lost():
    text = (
        "8001_3|et9_opt|COMPLETED|0:0|01:00:00|24:00:00|"
        "billing=4,cpu=4,mem=4000M,node=1|4|500M|"
        "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T09:00:00\n"
    )
    job = ce.account(["8001_3"],
                     run_fn=_sacct_runner(text, []))["jobs"]["8001_3"]
    assert job["state"] == "COMPLETED"
    assert job["core_hours_used"] == pytest.approx(4.0)


# ---------------------------------------------------------------- D6
# allow-list check: the argv an account() call builds is exactly the
# documented one — no option abbreviation, no extra flags.
def test_sacct_argv_is_the_documented_full_form():
    calls = []
    ce.account(["9001", "9002"], run_fn=_sacct_runner("", calls))
    assert calls == [["sacct", "-j", "9001,9002", "-n", "-X", "-P",
                      "-o", SACCT_COLUMNS]]


# ---------------------------------------------------------------- D7
# THROTTLED surfacing when the gate is a real SchedulerClient with a
# fresh argv (different ids) inside the interval: second distinct call
# is blocked -> None, and the runner is hit exactly once across BOTH.
def test_second_distinct_call_within_interval_blocked_once():
    calls = []
    now = [1000.0]
    client = SchedulerClient(min_interval_s=25.0,
                             clock=(lambda: now[0]),
                             runner=_sacct_runner(SACCT_OK_TEXT, calls))
    first = ce.account(["1101"], client=client)
    now[0] += 1.0
    second = ce.account(["1102"], client=client)
    assert first is not None
    assert second is None
    assert len(calls) == 1


SACCT_OK_TEXT = (
    "1101|etA_opt|COMPLETED|0:0|01:00:00|24:00:00|"
    "billing=8,cpu=8,mem=8000M,node=1|8|800M|"
    "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T09:00:00\n"
)


# ---------------------------------------------------------------- D8
# A CANCELLED+ or TIMEOUT+ state keeps its base state after split.
def test_state_with_suffix_maps_to_base_wasted_state():
    for state in ("CANCELLED by 848213", "TIMEOUT+", "OUT_OF_MEMORY+",
                  "NODE_FAIL+", "PREEMPTED+"):
        text = (
            "9500|etB_opt|COMPLETED|0:0|00:30:00|24:00:00|"
            "billing=4,cpu=4,mem=4000M,node=1|4|100M|"
            "2026-10-17T08:00:00|2026-10-17T08:00:00|2026-10-17T08:30:00\n"
        )
        text = text.replace("COMPLETED", state)
        job = ce.account(["9500"], run_fn=_sacct_runner(text, []))["jobs"]["9500"]
        assert job["state"] == state.split("+")[0].split(" ")[0]
        assert job["core_hours_wasted"] == pytest.approx(2.0)

