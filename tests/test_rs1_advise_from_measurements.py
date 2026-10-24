"""RS1 phase-2 red tests: pre-submit advice from measured efficiency.

advise() sizes PAL/maxcore through DELFIN's own logic but knows
nothing about how a step kind actually used its cores on this
cluster — the wave-21 lesson (SOLVATOR shell at pal 2 measured
1.53 % CPU, .gate/practice/sv1_tiny_solvator/delfin_7922181.out,
'CPU Utilized' 00:01:31) is invisible before a submit kills the job
at its time limit. Phase 2 lets the caller hand in the account()
records EF1 already produces and turns a measured low-parallel
step into advice. Advice only: the user's PAL/maxcore/time are
never changed.

Red on the previous commit: advise() has no efficiency_records
parameter and no measured-efficiency advice.
"""
import pytest

from delfin.agent import cluster_efficiency as ce


def _rec(jid, name, state, eff, wall_s, cores=2):
    """One account() job row, the shape EF1's account() returns."""
    return {"job_id": jid, "job_name": name, "state": state,
            "cpu_efficiency": eff, "notes": [],
            "used": {"cores": cores, "wall_s": wall_s}}


def test_low_efficiency_step_gets_measured_advice():
    low = _rec("101", "xtb_solvator", "COMPLETED", 0.015, 3580)
    adv = ce.advise("xtb_solvator", 12, efficiency_records=[low])
    base = ce.advise("xtb_solvator", 12)
    # Advice only: nothing about the user's sizing changes.
    assert adv["PAL"] == base["PAL"]
    assert adv["maxcore_mb"] == base["maxcore_mb"]
    eff = adv["efficiency"]
    assert eff["low_parallel"] is True
    assert eff["measured_runs"] == 1
    assert eff["mean_cpu_efficiency"] == pytest.approx(0.015)
    assert eff["max_measured_wall_s"] == 3580
    # Every number in the advice names where it comes from.
    assert "account()" in eff["source"]
    joined = "\n".join(adv["notes"])
    assert "1.5" in joined          # the measured efficiency
    assert "xtb_solvator" in joined  # the step it was measured on


def test_well_parallelised_dft_step_gets_no_advice():
    good = _rec("102", "dft_opt", "COMPLETED", 0.9, 1800, cores=8)
    adv = ce.advise("opt", 20, efficiency_records=[good])
    assert adv["efficiency"]["low_parallel"] is False
    assert adv["efficiency"]["advice"] is None


def test_artefact_rows_are_excluded_never_zero():
    """TIMEOUT/CANCELLED rows carry the scheduler's TotalCPU artefact
    (operator caution, wave 22): they must be excluded from the mean,
    never read as 0 %."""
    killed = _rec("301", "xtb_solvator", "TIMEOUT", None, 1800)
    adv = ce.advise("xtb_solvator", 12, efficiency_records=[killed])
    eff = adv["efficiency"]
    assert eff["measured_runs"] == 0
    assert eff["low_parallel"] is False
    assert eff["advice"] is None
    assert eff["excluded_runs"] == 1


def test_step_name_collision_is_exact_match():
    """'crest_preopt' is not a measurement of 'crest': substring
    matching would mis-attribute one step's efficiency to another."""
    other = _rec("103", "crest_preopt", "COMPLETED", 0.02, 600)
    adv = ce.advise("crest", 30, efficiency_records=[other])
    assert adv["efficiency"]["measured_runs"] == 0
    assert adv["efficiency"]["advice"] is None


def test_no_records_behaves_like_today():
    adv = ce.advise("opt", 20)
    base = ce.advise("opt", 20)
    assert adv["efficiency"] is None
    assert adv["PAL"] == base["PAL"]


def test_boundary_threshold_from_measured_populations():
    """The boundary sits between the two recorded populations
    (~1.5 % SOLVATOR vs ~90 % DFT); exactly 10 % is still fine."""
    at = _rec("104", "opt", "COMPLETED", ce.LOW_PARALLEL_EFFICIENCY, 600)
    adv = ce.advise("opt", 20, efficiency_records=[at])
    assert adv["efficiency"]["low_parallel"] is False
    below = _rec("105", "opt", "COMPLETED",
                 ce.LOW_PARALLEL_EFFICIENCY - 0.001, 600)
    adv = ce.advise("opt", 20, efficiency_records=[below])
    assert adv["efficiency"]["low_parallel"] is True


def test_hostile_step_names_do_not_crash():
    other = _rec("106", "opt", "COMPLETED", 0.9, 600)
    for step in ("OPT", "x" * 5000, "schiff_basesüß", "", None):
        adv = ce.advise(step, 12, efficiency_records=[other])
        assert "efficiency" in adv


# --------------------------------------------------------------------
# Reviewer findings F1/F2 (s15, on c0cfb0f7). Red on c0cfb0f7.
# --------------------------------------------------------------------

def test_nan_efficiency_is_not_silent_low_parallel_false():
    """F1: a NaN record passes isinstance(float) and poisons the mean;
    NaN < threshold is False, so a genuinely low-parallel step would
    silently get low_parallel=False and no advice. NaN/inf must be
    excluded like the TIMEOUT/CANCELLED rows, with a count."""
    low = _rec("101", "xtb_solvator", "COMPLETED", 0.015, 3580)
    poisoned = _rec("102", "xtb_solvator", "COMPLETED", float("nan"),
                    3500)
    adv = ce.advise("xtb_solvator", 12,
                    efficiency_records=[low, poisoned])
    eff = adv["efficiency"]
    assert eff["low_parallel"] is True, eff
    assert eff["measured_runs"] == 1, eff
    assert eff["excluded_runs"] == 1, eff
    assert eff["mean_cpu_efficiency"] == pytest.approx(0.015), eff


def test_inf_efficiency_is_excluded_like_nan():
    """F1: +inf passes the same isinstance filter."""
    low = _rec("101", "xtb_solvator", "COMPLETED", 0.015, 3580)
    poisoned = _rec("102", "xtb_solvator", "COMPLETED", float("inf"),
                    3500)
    adv = ce.advise("xtb_solvator", 12,
                    efficiency_records=[low, poisoned])
    eff = adv["efficiency"]
    assert eff["low_parallel"] is True, eff
    assert eff["excluded_runs"] == 1, eff


def test_summary_always_carries_the_advice_key():
    """F2: healthy low/normal summaries and zero-run summaries all
    carry "advice" (None when not low) — one shape for consumers."""
    low = _rec("101", "xtb_solvator", "COMPLETED", 0.015, 3580)
    good = _rec("102", "dft_opt", "COMPLETED", 0.9, 1800, cores=8)
    low_adv = ce.advise("xtb_solvator", 12, efficiency_records=[low])
    assert "advice" in low_adv["efficiency"], low_adv
    good_adv = ce.advise("opt", 20, efficiency_records=[good])
    assert good_adv["efficiency"]["advice"] is None, good_adv
    zero = ce.advise("crest", 4, efficiency_records=[good])
    assert zero["efficiency"]["advice"] is None, zero
