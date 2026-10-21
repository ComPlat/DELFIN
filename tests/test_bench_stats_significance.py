"""A compare must label each per-task pass-rate difference as within or
beyond noise, and the `bench compare` markdown must surface the statistical
statement -- not just a direction.

The data model already computes a pooled two-sided Fisher p
(``compare_runs`` → ``summary.significant``/``significance_p``), but:

1. per-task rows carry only raw rate deltas -- nothing says whether a single
   task's 3/3 → 0/3 is anything but a sampling accident at N=3;
2. ``format_compare_markdown`` -- what ``bench compare`` prints -- never shows
   the p-value, the alpha, or the Wilson intervals, so a reader cannot tell a
   real change from run-to-run noise.

Known values (two-sided Fisher, confirmed against ``_fisher_exact_2x2_pvalue``):
  3/3 vs 0/3 -> p = 0.100 (not resolvable at N=3, alpha .05)
  4/4 vs 0/4 -> p = 2/70 ~ 0.0286 (resolvable)
  9/10 vs 1/10 -> p = 0.00109 (resolvable)
  2/3 vs 1/3 -> p = 1.0 (not resolvable)
"""

import pytest

from delfin.agent.benchmark import classify_pass_delta, compare_runs, format_compare_markdown


# ── pure helper: within / beyond noise for one task's pass delta ──────────

def test_classify_pass_delta_known_values():
    # 3/3 vs 0/3: the classic N=3 trap -- not resolvable at alpha .05.
    d = classify_pass_delta(3, 3, 0, 3)
    assert d["label"] == "within"
    assert d["p"] == pytest.approx(0.1, abs=1e-9)

    # 4/4 vs 0/4: one more repeat per arm and the same swing resolves.
    d = classify_pass_delta(4, 4, 0, 4)
    assert d["label"] == "beyond"
    assert d["p"] == pytest.approx(2 / 70, abs=1e-9)

    # A large, clean effect resolves comfortably.
    d = classify_pass_delta(9, 10, 1, 10)
    assert d["label"] == "beyond"
    assert d["p"] < 0.05

    # 2/3 vs 1/3 is indistinguishable: p = 1.0.
    d = classify_pass_delta(2, 3, 1, 3)
    assert d["label"] == "within"
    assert d["p"] == pytest.approx(1.0)


def test_classify_pass_delta_defaults_to_alpha_05():
    assert classify_pass_delta(4, 4, 0, 4)["label"] == "beyond"
    # At a stricter alpha the same data is not resolvable.
    assert classify_pass_delta(4, 4, 0, 4, alpha=0.01)["label"] == "within"


def test_classify_identical_rates_are_within_noise():
    d = classify_pass_delta(5, 10, 5, 10)
    assert d["label"] == "within"
    assert d["p"] == pytest.approx(1.0)


# ── compare_runs exposes the label on each per-task row ───────────────────

def _row(task_id: str, n_pass: int, n: int) -> dict:
    flags = [True] * n_pass + [False] * (n - n_pass)
    return {
        "task_id": task_id,
        "success": n_pass * 2 >= n,
        "success_rate": n_pass / n,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 50,
        "cost_usd": 0.1,
        "duration_s": 10.0,
        "tool_calls": 5,
        "input_tokens": 400,
        "output_tokens": 200,
    }


def test_compare_runs_per_task_carries_noise_label_and_p():
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    row = out["per_task"][0]
    assert row["delta_noise"] == "within"
    assert row["delta_p"] == pytest.approx(0.1, abs=1e-9)
    # The observed delta is still reported -- direction is not suppressed.
    assert row["success_rate_old"] == 1.0
    assert row["success_rate_new"] == 0.0


def test_compare_runs_per_task_beyond_when_resolvable():
    out = compare_runs([_row("t", 4, 4)], [_row("t", 0, 4)])
    assert out["per_task"][0]["delta_noise"] == "beyond"


# ── the printed markdown surfaces the statistical statement ───────────────

def test_markdown_prints_noise_label_and_p_for_insufficient_n():
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    md = format_compare_markdown(out)
    assert "within noise" in md
    assert "0.1" in md            # the p that decides it
    assert "pooled" in md.lower()


def test_markdown_prints_beyond_noise_when_resolvable():
    out = compare_runs([_row("t", 9, 10)], [_row("t", 1, 10)])
    md = format_compare_markdown(out)
    assert "beyond noise" in md
    assert "0.0011" in md


def test_markdown_prints_wilson_intervals_from_repeats():
    # Wilson interval is the honest spread from N repeats; the compare
    # output must carry it, not just the point estimate.
    out = compare_runs([_row("t", 3, 3)], [_row("t", 0, 3)])
    md = format_compare_markdown(out)
    assert "wilson" in md.lower()
    # 3/3 is NOT 100%: the interval must show an upper bound below 100%.
    assert "%" in md


def test_markdown_identical_runs_call_it_within_noise():
    out = compare_runs([_row("t", 2, 3)], [_row("t", 2, 3)])
    assert "within noise" in format_compare_markdown(out)
