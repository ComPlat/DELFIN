"""RS1 adversarial tests, phase 1b (reviewer nacht-s15).

Target: the TotalCPU-artefact handling on TIMEOUT/CANCELLED rows
(cluster_efficiency.py:49, :160-169, merged @92549fad).

Attack angles:
1. Compound sacct state strings: real scheduler rows say
   "CANCELLED by 12345", "TIMEOUT+" (kill-on-invalid-dep), lower-case
   or suffixed variants. The normalised state must still be recognised
   as the artefact family — a look-alike spelling that slips through
   would report a fabricated ~0 % as advice-grade efficiency.
2. A genuinely low measured efficiency on a healthy COMPLETED row must
   stay MEASURED (not swept into the artefact rule) — the artefact
   rule is state-gated, never a general excuse.
3. Healthy rows carry an empty notes list (shape unchanged for other
   consumers).
4. PREEMPTED / NODE_FAIL / OUT_OF_MEMORY keep measured efficiency
   (builder said FAILED etc. keep theirs — pinned).
"""
from delfin.agent import cluster_efficiency as ce


def _run(text):
    def run(argv):
        return text
    return run


def _row(job_id, state, total_cpu, elapsed="00:30:00", alloc="2"):
    return (f"{job_id}|job_{job_id}|{state}|0:0|{elapsed}|01:00:00|"
            f"cpu={alloc},mem=4G|{alloc}|1G|2026-10-20T22:00:00|"
            f"2026-10-20T22:00:00|2026-10-20T22:30:00|{total_cpu}\n")


def test_compound_cancelled_by_n_is_artefact_family():
    out = _row("401", "CANCELLED by 12345", "00:00:00")
    job = ce.account(["401"], run_fn=_run(out))["jobs"]["401"]
    assert job["cpu_efficiency"] is None, job
    assert any("TotalCPU" in n for n in job.get("notes", [])), job


def test_compound_timeout_plus_suffix_is_artefact_family():
    out = _row("402", "TIMEOUT+", "00:00:15")
    job = ce.account(["402"], run_fn=_run(out))["jobs"]["402"]
    assert job["cpu_efficiency"] is None, job
    assert any("TIMEOUT" in n for n in job.get("notes", [])), job


def test_low_efficiency_on_completed_row_stays_measured():
    """The artefact rule is STATE-gATED. A COMPLETED row with a real
    1.5 % efficiency is a measurement, never reported UNKNOWN."""
    out = _row("403", "COMPLETED", "00:00:30")
    job = ce.account(["403"], run_fn=_run(out))["jobs"]["403"]
    eff = job["cpu_efficiency"]
    assert eff is not None, job
    assert abs(eff - 30.0 / (2 * 1800.0)) < 1e-6, eff
    assert job.get("notes") == [], job


def test_other_wasted_states_keep_measured_efficiency():
    out = _row("404", "PREEMPTED", "00:10:00")
    job = ce.account(["404"], run_fn=_run(out))["jobs"]["404"]
    assert job["cpu_efficiency"] is not None, job
    # 600 CPU-s / (2 cores * 1800 s) = 0.1667 — measured, whatever the value.
    assert abs(job["cpu_efficiency"] - 600.0 / 3600.0) < 1e-3, job
    out = _row("405", "NODE_FAIL", "00:20:00")
    job = ce.account(["405"], run_fn=_run(out))["jobs"]["405"]
    assert job["cpu_efficiency"] is not None, job
    assert job.get("notes") == [], job


def test_failed_state_keeps_measured_efficiency():
    out = _row("406", "FAILED", "01:00:00", elapsed="00:10:00", alloc="8")
    job = ce.account(["406"], run_fn=_run(out))["jobs"]["406"]
    eff = job["cpu_efficiency"]
    assert eff is not None, job
    assert abs(eff - 3600.0 / (8 * 600.0)) < 1e-3, eff
