"""A waiting job says why it waits and when SLURM expects it to start.

The table showed squeue's raw "(Priority)" and nothing else; the estimate sat
in a separate box fed by a second squeue on every refresh.
"""

from datetime import datetime

from delfin.dashboard.backend_base import JobInfo
from delfin.dashboard.tab_job_status import (
    _fmt_duration,
    _fmt_memory,
    _fmt_when,
    _relative_start,
    build_slurm_job_table,
    empty_queue_html,
)

NOW = datetime(2026, 9, 14, 14, 0)


def _pending(**extra):
    base = {
        "partition": "cpu,cpu_il", "time_used": "0:00", "nodes": "1",
        "reason": "(Priority)", "reason_code": "Priority", "time_limit": "2-00:00:00",
        "cpus": "40", "memory": "240000M", "start_estimate": "2026-09-16T21:20:09",
    }
    base.update(extra)
    return JobInfo(job_id="900", name="TADF 1", status="PD", extra=base)


def _running(**extra):
    base = {
        "partition": "cpu_il", "time_used": "1:02:03", "nodes": "1", "reason": "uc2n872",
        "reason_code": "None", "time_limit": "1-00:00:00", "cpus": "40",
        "memory": "240000M", "start_estimate": "2026-09-14T12:57:57",
    }
    base.update(extra)
    return JobInfo(job_id="901", name="Opt Fe", status="R", extra=base)


def test_a_waiting_job_shows_its_reason_and_expected_start():
    options = []
    html = build_slurm_job_table([_pending()], options, now=NOW)

    assert "Waiting" in html
    assert "Priority" in html and "higher priority" in html
    assert "16 Sep 21:20" in html and "in 2 d 7 h" in html
    assert "Expected start Wed 16 Sep, 21:20, in 2 d 7 h" in html
    assert ">cpu<" in html and ">cpu_il<" in html
    assert "whichever partition frees first" in html
    assert "40 CPUs · 240 GB" in html and "limit 2 d" in html
    assert 'title="Priority: Waiting for jobs with higher priority to start first."' in html
    assert "About estimated starts" in html
    assert options == [("900 - TADF 1", "900")]


def test_a_raw_reason_code_gets_a_short_label_and_keeps_the_code_in_the_tooltip():
    html = build_slurm_job_table(
        [_pending(reason="(QOSMaxCpuPerUserLimit)", reason_code="QOSMaxCpuPerUserLimit")], [], now=NOW)
    assert ">Queue limit<" in html
    assert "(SLURM reason QOSMaxCpuPerUserLimit)" in html

    held = build_slurm_job_table([_pending(reason="(JobHeldUser)", reason_code="JobHeldUser")], [], now=NOW)
    assert ">Held by you<" in held and "djs-stuck" in held, "a job that will not start by waiting stands out"


def test_a_job_without_an_estimate_says_so():
    html = build_slurm_job_table([_pending(start_estimate="")], [], now=NOW)
    assert "not estimated yet" in html


def test_a_running_job_shows_its_node_and_how_far_it_is():
    html = build_slurm_job_table([_running()], [], now=NOW)

    assert "Running" in html
    assert "uc2n872" in html
    assert "1 h 2 min / 1 d" in html
    assert "since</span> 14 Sep 12:57" in html
    assert 'title="Started Mon 14 Sep, 12:57"' in html
    assert "width:4%" in html, "3723 s of 86400 s"
    assert "About estimated starts" not in html, "nothing is waiting"


def test_a_job_close_to_its_limit_is_marked():
    html = build_slurm_job_table([_running(time_used="22:00:00")], [], now=NOW)
    assert "djs-warn" in html


def test_the_summary_counts_the_whole_queue_even_when_the_table_is_filtered():
    jobs = [_pending(), _running(), _running(cpus="16")]
    html = build_slurm_job_table([jobs[0]], [], now=NOW, summary_jobs=jobs)

    assert '>2</div><div class="djs-stat-label">Running' in html
    assert '>1</div><div class="djs-stat-label">Waiting' in html
    assert '>56</div><div class="djs-stat-label">CPUs in use' in html
    assert ">in 2 d 7 h</div><div class=\"djs-stat-label\">Next expected start" in html


def test_a_name_cannot_inject_markup():
    job = JobInfo(job_id="902", name='<b>x</b>" onmouseover="y', status="R",
                  extra={"reason": "uc3n024"})
    html = build_slurm_job_table([job], [], now=NOW)
    assert "&lt;b&gt;x&lt;/b&gt;" in html and "<b>x</b>" not in html
    assert '" onmouseover="' not in html, "a quote in a name must not leave its title attribute"


def test_a_job_stays_on_one_line_however_long_its_text():
    """Long names and explanations are cut with an ellipsis, not wrapped."""
    long_name = "TADF_" + "very_long_ligand_name_" * 8
    html = build_slurm_job_table([_pending(), _running()], [], now=NOW)
    wide = build_slurm_job_table([_running()], [], now=NOW)

    assert "table-layout: fixed" in html
    assert "white-space: nowrap" in html and "text-overflow: ellipsis" in html
    assert "<br" not in html and "<div class=\"djs-sub\"" not in html, "no second line inside a row"
    assert "@container" in html, "columns give way on narrow panes"
    named = build_slurm_job_table([_running()], [], now=NOW).replace("Opt Fe", long_name)
    assert long_name in named and "<br" not in named
    assert wide.count("<tr>") == 2, "header and one row"


def test_an_empty_queue_says_what_will_appear_there():
    assert "No jobs in the queue" in empty_queue_html()
    assert "No waiting jobs" in empty_queue_html("pending")


def test_durations_memory_and_times_read_like_a_person_would_say_them():
    assert _fmt_duration(172800) == "2 d"
    assert _fmt_duration(129600) == "1 d 12 h"
    assert _fmt_duration(3723) == "1 h 2 min"
    assert _fmt_duration(7200) == "2 h"
    assert _fmt_duration(30) == "< 1 min"
    assert _fmt_memory("240000M") == "240 GB"
    assert _fmt_memory("240G") == "240 GB"
    assert _fmt_memory("1500M") == "1.5 GB"
    assert _fmt_memory("500M") == "500 MB"
    assert _fmt_when("2026-09-16T21:20:09") == "Wed 16 Sep, 21:20"
    assert _fmt_when("N/A") == ""


def test_the_relative_start_reads_like_a_person_would_say_it():
    now = datetime(2026, 9, 14, 12, 0)
    assert _relative_start("2026-09-14T12:00:30", now) == "any moment now"
    assert _relative_start("2026-09-14T12:45:00", now) == "in 45 min"
    assert _relative_start("2026-09-14T15:10:00", now) == "in 3 h 10 min"
    assert _relative_start("", now) == ""
