"""RS1 adversarial tests, phase 1 (reviewer nacht-s15).

Target of the attack: cpu_efficiency in EF1's account()
(delfin/agent/cluster_efficiency.py:144-157, merged @029c7f10).

Angles:
1. Red control: account() at the wave base a4c9ebaf had no
   cpu_efficiency field at all (KeyError) — this file's field checks
   are red there and green on the fix.
2. Shape attacks on the computation: empty TotalCPU, unparseable
   TotalCPU ("UNK", "UNLIMITED"), zero AllocCPUS, zero elapsed,
   negative-looking values, and a row where efficiency exceeds 1.0
   (scheduler accounting artifacts) must all stay None or the honest
   number — never crash, never fabricate 0.0 %.
3. Real sacct columns: TotalCPU can carry a per-node format or end in
   different formats; the parser must survive any string.
4. A 12-field row (no TotalCPU column) keeps every other field intact.
"""
from delfin.agent import cluster_efficiency as ce


def _run(text):
    def run(argv):
        return text
    return run


def _row(job_id, total_cpu, alloc="2", elapsed="00:10:00", state="COMPLETED"):
    return (f"{job_id}|job_{job_id}|{state}|0:0|{elapsed}|01:00:00|"
            f"cpu={alloc},mem=4G|{alloc}|1G|2026-10-20T22:00:00|"
            f"2026-10-20T22:00:00|2026-10-20T22:10:00|{total_cpu}\n")


def test_efficiency_present_and_correct_on_merged_tree():
    """Red control at a4c9ebaf: no cpu_efficiency field (KeyError)."""
    rows = ce.account(["1"], run_fn=_run(_row("1", "00:05:00")))
    eff = rows["jobs"]["1"]["cpu_efficiency"]
    # 300 CPU-s / (2 cores * 600 s) = 0.25
    assert eff is not None
    assert abs(eff - 0.25) < 1e-6, eff


def test_empty_and_unparseable_total_cpu_stay_none():
    for junk in ("", "UNK", "UNLIMITED", "1-02:03:04:05", "not-a-time"):
        out = _row("7", junk)
        rows = ce.account(["7"], run_fn=_run(out))
        assert rows["jobs"]["7"]["cpu_efficiency"] is None, junk
        # The rest of the accounting must survive intact.
        assert rows["jobs"]["7"]["used"]["cores"] == 2, junk
        assert rows["jobs"]["7"]["used"]["wall_s"] == 600, junk


def test_zero_alloc_cpus_or_zero_elapsed_stay_none():
    rows = ce.account(["8"], run_fn=_run(
        _row("8", "00:05:00", alloc="0")))
    assert rows["jobs"]["8"]["cpu_efficiency"] is None
    rows = ce.account(["9"], run_fn=_run(
        _row("9", "00:00:10", elapsed="00:00:00")))
    assert rows["jobs"]["9"]["cpu_efficiency"] is None


def test_efficiency_over_one_is_reported_not_clamped():
    """Scheduler artifacts can exceed 1.0 (TimeLimitTotalCPU) — report
    the honest number, never silently clamp."""
    rows = ce.account(["10"], run_fn=_run(
        _row("10", "01:00:00", elapsed="00:10:00", alloc="2")))
    eff = rows["jobs"]["10"]["cpu_efficiency"]
    assert eff is not None
    assert eff > 1.0, eff


def test_twelve_field_row_keeps_fields_and_none_efficiency():
    out = ("202|old_shape|COMPLETED|0:0|00:10:00|01:00:00|"
           "cpu=4,mem=2G|4|1G|2026-10-20T22:00:00|"
           "2026-10-20T22:00:10|2026-10-20T22:10:00\n")
    rows = ce.account(["202"], run_fn=_run(out))
    job = rows["jobs"]["202"]
    assert job["cpu_efficiency"] is None
    assert job["used"]["cores"] == 4
    assert job["used"]["wall_s"] == 600
    assert job["state"] == "COMPLETED"
    assert job["requested"]["cores"] == 4


def test_parse_sacct_survives_per_node_total_cpu_shape():
    """sacct -P prints per-node TotalCPU as a pipe-free string; any
    string content must parse to a row, never raise."""
    out = _row("11", "arbitrary text with spaces")
    rows = ce.account(["11"], run_fn=_run(out))
    assert rows["jobs"]["11"]["cpu_efficiency"] is None
    assert rows["jobs"]["11"]["job_name"] == "job_11"
