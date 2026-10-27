"""RS1 adversarial tests, phase 2 (reviewer nacht-s15).

Target: _efficiency_summary / advise(efficiency_records=...)
(cluster_efficiency.py:342-396, 398-411, merged @c0cfb0f7).

Attack angles:
1. NaN cpu_efficiency sneaks past isinstance(float); the mean becomes
   NaN and NaN < 0.10 is False — advice silently swallowed for a
   poisoned record. The summary must not report low_parallel=False as
   if measured.
2. bool cpu_efficiency (True/False) passes isinstance(int).
3. efficiency_records as a non-list iterable (a dict) — .get() on the
   string keys must not crash silently into nonsense.
4. Empty records list: zero-run summary, never None, no advice.
5. A well-parallelised step keeps low_parallel False and advice None
   (both directions pinned).
6. advise() must never change PAL/maxcore/partitions when records are
   present (advice only).
"""
import math

from delfin.agent import cluster_efficiency as ce


def _rec(name, state="COMPLETED", eff=0.9, wall=600, **extra):
    rec = {"job_id": "1", "job_name": name, "state": state,
           "requested": {}, "used": {"wall_s": wall},
           "cpu_efficiency": eff, "notes": [], "queue_wait_s": None,
           "core_hours_used": 0.0, "core_hours_wasted": 0.0,
           "is_repeated_waste": False}
    rec.update(extra)
    return rec


def test_nan_efficiency_is_not_silent_low_parallel_false():
    s = ce._efficiency_summary("sv", [_rec("sv", eff=float("nan"))])
    assert not (s["measured_runs"] == 1 and s["mean_cpu_efficiency"] is not None
                and math.isnan(s["mean_cpu_efficiency"])
                and s["low_parallel"] is False), s


def test_bool_efficiency_is_rejected_or_counted_honestly():
    s = ce._efficiency_summary("sv", [_rec("sv", eff=True)])
    # A bool is not a measurement: either excluded, or the mean is
    # exactly 1.0 — never an opaque 1 in a float mean presented as data.
    assert s["mean_cpu_efficiency"] in (None, 1.0), s


def test_dict_records_do_not_produce_silent_nonsense():
    s = ce._efficiency_summary("sv", {"a": _rec("sv")})
    # A dict iterates its keys ('a'), which match no step: the honest
    # outcome is the zero-run summary, not a crash, not advice.
    assert s["low_parallel"] is False and s["advice"] is None, s
    assert s["measured_runs"] == 0, s


def test_empty_records_zero_run_summary():
    s = ce._efficiency_summary("sv", [])
    assert s["measured_runs"] == 0 and s["low_parallel"] is False
    assert s["advice"] is None and s["mean_cpu_efficiency"] is None


def test_well_parallelised_step_gets_no_advice():
    s = ce._efficiency_summary("opt", [_rec("opt", eff=0.92),
                                       _rec("opt", eff=0.88)])
    assert s["low_parallel"] is False
    assert s.get("advice") is None  # FINDING for builder: healthy rows
    # lack the "advice" key entirely while zero-run rows carry None —
    # an inconsistent shape for s["advice"] consumers. Pinned via .get
    # here; the key-presence inconsistency is reported to nacht-s13.
    assert abs(s["mean_cpu_efficiency"] - 0.9) < 1e-9


def test_low_parallel_step_gets_advice_with_source():
    s = ce._efficiency_summary("sv", [_rec("sv", eff=0.015, wall=3580)])
    assert s["low_parallel"] is True
    assert "0.015" in s["advice"].replace("1.5", "0.015") or "%" in s["advice"]
    assert s["max_measured_wall_s"] == 3580
    assert "account()" in s["source"]


def test_advise_leaves_sizing_untouched_with_records():
    base = ce.advise("sv", 12)
    recs = [_rec("sv", eff=0.015, wall=3580)]
    adv = ce.advise("sv", 12, efficiency_records=recs)
    assert adv["PAL"] == base["PAL"], (adv["PAL"], base["PAL"])
    assert adv.get("maxcore") == base.get("maxcore")
    assert adv.get("partitions") == base.get("partitions")
