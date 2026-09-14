"""A waiting job says why it waits and when SLURM expects it to start.

The table showed squeue's raw "(Priority)" and nothing else; the estimate sat
in a separate box fed by a second squeue on every refresh.
"""

from datetime import datetime

from delfin.dashboard.backend_base import JobInfo
from delfin.dashboard.tab_job_status import _relative_start, build_slurm_job_table


def _pending(**extra):
    base = {
        "partition": "cpu,cpu_il", "time_used": "0:00", "nodes": "1",
        "reason": "(Priority)", "reason_code": "Priority", "time_limit": "2-00:00:00",
        "cpus": "40", "memory": "240000M", "start_estimate": "2026-09-16T21:20:09",
    }
    base.update(extra)
    return JobInfo(job_id="900", name="TADF 1", status="PD", extra=base)


def test_a_waiting_job_shows_its_reason_and_expected_start():
    options = []
    html = build_slurm_job_table([_pending()], options, now=datetime(2026, 9, 14, 14, 0))

    assert "<b>Priority</b>" in html
    assert "higher priority" in html
    assert "2026-09-16 21:20" in html
    assert "in 2 d 7 h" in html
    assert "cpu, cpu_il" in html and "first to free wins" in html
    assert options == [("900 - TADF 1", "900")]


def test_a_job_without_an_estimate_says_so():
    html = build_slurm_job_table([_pending(start_estimate="")], [])
    assert "not estimated yet" in html


def test_a_running_job_shows_its_node_and_a_name_cannot_inject_markup():
    job = JobInfo(job_id="901", name="<b>x</b>", status="R",
                  extra={"reason": "uc3n024", "start_estimate": "2026-09-14T08:00:00"})
    html = build_slurm_job_table([job], [])
    assert "uc3n024" in html
    assert "&lt;b&gt;x&lt;/b&gt;" in html and "<b>x</b>" not in html
    assert "not estimated yet" not in html


def test_the_relative_start_reads_like_a_person_would_say_it():
    now = datetime(2026, 9, 14, 12, 0)
    assert _relative_start("2026-09-14T12:00:30", now) == "any moment now"
    assert _relative_start("2026-09-14T12:45:00", now) == "in 45 min"
    assert _relative_start("2026-09-14T15:10:00", now) == "in 3 h 10 min"
    assert _relative_start("", now) == ""
