"""RS1 phase-1 red test: CPU efficiency per job from real sacct data.

The wave-21 lesson (SOLVATOR at pal 2 ran at ~1.5 % CPU) is a number
EF1's :func:`account` does not report today: it measures requested vs
used cores, wall time and core-hours, but never the fraction of the
allocated cores that actually computed (TotalCPU vs AllocCPUS *
Elapsed). Package RS1 needs exactly that fraction to give pre-submit
advice from measurements instead of guesses.

Red on the previous commit: ``account`` has no CPU-efficiency field.
"""
from delfin.agent import cluster_efficiency as ce


# A real-shaped sacct -P -X -n output: one SOLVATOR-like job asking 2
# cores for 01:00:00, running 00:59:40 with 120 CPU-seconds total CPU
# (1.5 cores per core-second would be 0; here ~1.7 % efficiency).
SACCT_SOLVATOR_LIKE = (
    "101|xtb_solvator|COMPLETED|0:0|00:59:40|01:00:00|"
    "cpu=2,mem=4G|2|1.5G|2026-10-20T22:00:00|2026-10-20T22:01:00|"
    "2026-10-20T23:00:40|00:02:00\n"
    "102|dft_opt|COMPLETED|0:0|00:30:00|01:00:00|"
    "cpu=8,mem=8G|8|3.2G|2026-10-20T22:00:00|2026-10-20T22:05:00|"
    "2026-10-20T22:35:00|07:12:00\n"
)


def _sacct_runner(text):
    def run(argv):
        return text
    return run


def _account() -> dict:
    return ce.account(
        ["101", "102"],
        run_fn=_sacct_runner(SACCT_SOLVATOR_LIKE),
    )


def test_cpu_efficiency_is_reported_per_job():
    rows = _account()
    assert "cpu_efficiency" in rows["jobs"]["101"], (
        "account() must report the fraction of allocated cores that "
        "actually computed, from TotalCPU / (AllocCPUS * Elapsed)")
    eff = rows["jobs"]["101"]["cpu_efficiency"]
    assert eff is not None
    # Elapsed 00:59:40 = 3580 s; TotalCPU 00:02:00 = 120 CPU-s;
    # efficiency = 120 / (2 * 3580) ~ 0.017
    assert abs(eff - 120.0 / (2 * 3580.0)) < 1e-3, eff


def test_well_parallelised_dft_step_shows_high_efficiency():
    rows = _account()
    eff = rows["jobs"]["102"]["cpu_efficiency"]
    # 00:30:00 = 1800 s at 8 cores; TotalCPU 07:12:00 = 25920 CPU-s
    # => efficiency = 25920 / (8 * 1800) = 0.9 — a near-fully-used run.
    assert eff is not None
    assert abs(eff - 25920.0 / (8 * 1800.0)) < 1e-3, eff
    assert eff > 0.5, eff


def test_missing_total_cpu_stays_none_not_zero():
    """A scheduler row without TotalCPU must not invent 0.0 %."""
    out = (
        "201|some_job|COMPLETED|0:0|00:10:00|01:00:00|cpu=4,mem=2G|4|1G|"
        "2026-10-20T22:00:00|2026-10-20T22:00:10|2026-10-20T22:10:00\n"
    )
    rows = ce.account(["201"], run_fn=_sacct_runner(out))
    assert rows["jobs"]["201"]["cpu_efficiency"] is None


# --------------------------------------------------------------------
# Operator cautions (wave-22, RS1): two artefacts of the scheduler's
# TotalCPU number must never become a fabricated efficiency.
# --------------------------------------------------------------------

def test_timeout_job_reports_unknown_efficiency_never_zero():
    """A job killed at its time limit computed, but this cluster's
    TotalCPU for TIMEOUT rows is a known artefact (0 although the job
    computed). The efficiency must be UNKNOWN (None), never 0.0, and
    the row must carry a note saying why."""
    out = (
        "301|shell_run|TIMEOUT|0:9|00:30:00|00:30:00|cpu=12,mem=8G|12|2G|"
        "2026-10-20T22:00:00|2026-10-20T22:01:00|2026-10-20T22:31:00|"
        "00:00:00\n"
    )
    rows = ce.account(["301"], run_fn=_sacct_runner(out))
    job = rows["jobs"]["301"]
    assert job["cpu_efficiency"] is None, job
    assert any("TIMEOUT" in note for note in job.get("notes", [])), job


def test_cancelled_job_reports_unknown_efficiency_never_zero():
    """CANCELLED (also by the backup) is the same artefact family."""
    out = (
        "302|killed_run|CANCELLED|0:0|00:15:00|01:00:00|cpu=4,mem=2G|4|1G|"
        "2026-10-20T22:00:00|2026-10-20T22:01:00|2026-10-20T22:16:00|"
        "00:00:00\n"
    )
    rows = ce.account(["302"], run_fn=_sacct_runner(out))
    job = rows["jobs"]["302"]
    assert job["cpu_efficiency"] is None, job
    assert any("TotalCPU" in note for note in job.get("notes", [])), job


def test_zero_artefact_never_reported_as_low_efficiency_on_wasted_state():
    """The red direction: an arbitrary tiny-but-plausible TotalCPU on a
    TIMEOUT row is still the artefact, not a measurement — efficiency
    stays None even when the number would have computed to 1 %."""
    out = (
        "303|shell_run|TIMEOUT|0:9|00:30:00|00:30:00|cpu=2,mem=4G|2|1G|"
        "2026-10-20T22:00:00|2026-10-20T22:01:00|2026-10-20T22:31:00|"
        "00:00:30\n"
    )
    rows = ce.account(["303"], run_fn=_sacct_runner(out))
    assert rows["jobs"]["303"]["cpu_efficiency"] is None
