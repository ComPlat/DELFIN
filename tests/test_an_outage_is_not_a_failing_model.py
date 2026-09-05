"""A run that never reached the model was scored as a model that failed.

First live office run, 2026-08-12, against the KIT-hosted Qwen with
three repeats. Two of eleven tasks came back at quality 35, pass rate 0,
zero tool calls, about 1.4 seconds. Their entire recorded model output
was the harness's own retry banners::

    ⏳ Transient API error (BadRequestError); retrying 1/3 in 2s…
    ⏳ Transient API error (BadRequestError); retrying 2/3 in 3s…
    ⏳ Transient API error (BadRequestError); retrying 3/3 in 6s…

and the error field said "No deployments available for selected model".
The endpoint had no capacity. Nothing about the model was observed — and
the suite booked it as a model that failed, taking 9/11 down to 7/11 and
writing that into the file baselines are compared against.

A baseline poisoned this way is worse than no baseline: it lowers the
bar permanently and silently, so the next real regression passes under
it. It is the shape this project keeps removing, sitting inside the
instrument that is supposed to detect the shape.

The distinction that matters is narrow on purpose. A model that answered
badly is a failure. A tool that broke inside a live turn is a failure. A
request that never reached the model is neither — it is a measurement
that did not happen, and it is excluded from the rate and the average
and named in the summary.
"""

from __future__ import annotations

import pytest

from delfin.agent import benchmark as B


_RETRY_BANNERS = (
    "⏳ Transient API error (BadRequestError); retrying 1/3 in 2s…\n\n"
    "⏳ Transient API error (BadRequestError); retrying 2/3 in 3s…\n\n"
    "⏳ Transient API error (BadRequestError); retrying 3/3 in 6s…\n"
)
_NO_DEPLOYMENT = (
    "Error code: 400 - {'detail': \"No deployments available for selected "
    "model, Try again in 5 seconds. Passed model=kit.qwen3.5-397b-A17b\"}"
)


def _task(tid="t"):
    return B.Task(id=tid, task_class="office", mode="office", prompt="x")


# ---------------------------------------------------------------------------
# The incident, as it happened
# ---------------------------------------------------------------------------

def test_the_run_that_produced_only_retry_banners_is_unmeasured():
    traj = B.Trajectory(text=_RETRY_BANNERS, error=_NO_DEPLOYMENT,
                        duration_s=1.5)
    assert B.is_unmeasured(traj) is True


def test_it_is_marked_on_the_result():
    traj = B.Trajectory(text=_RETRY_BANNERS, error=_NO_DEPLOYMENT)
    assert B.score_outcome(_task(), traj, model="m").unmeasured is True


@pytest.mark.parametrize("err", [
    "Error code: 503 - Service Unavailable",
    "Open WebUI: Server Connection Error",
    "connection refused",
    "Error code: 429 - rate limit exceeded",
    "Error code: 401 - invalid api key",
    "engine init failed: no credential",
])
def test_every_shape_of_not_reaching_the_model_counts(err):
    assert B.is_unmeasured(B.Trajectory(text="", error=err)) is True


# ---------------------------------------------------------------------------
# ...and a real failure is still a real failure
# ---------------------------------------------------------------------------

def test_a_model_that_answered_badly_is_a_failure_not_an_outage():
    traj = B.Trajectory(text="Die Summe beträgt 1.986,40 EUR.", error="")
    assert B.is_unmeasured(traj) is False


def test_a_turn_that_did_work_and_then_broke_is_still_scored():
    """The narrow half of the rule: a transport error near the END of a
    turn that ran tools and produced an answer says plenty about the
    model, and must not be excused."""
    traj = B.Trajectory(
        text="Ich habe die Datei gelesen und komme auf 1.986,40 EUR.",
        tool_calls=[{"name": "read_document", "input": "{}"}],
        error="Error code: 503 - Service Unavailable")
    assert B.is_unmeasured(traj) is False


def test_a_tool_error_inside_a_live_turn_is_not_an_outage():
    traj = B.Trajectory(text="Der Lesevorgang schlug fehl.",
                        error="OfficeError: could not open workbook")
    assert B.is_unmeasured(traj) is False


def test_an_error_free_run_is_never_unmeasured():
    assert B.is_unmeasured(B.Trajectory(text="", error="")) is False


# ---------------------------------------------------------------------------
# The summary excludes it, and says so
# ---------------------------------------------------------------------------

def _rows():
    good = B.score_outcome(
        _task("ok"), B.Trajectory(text="fine", duration_s=1.0), model="m")
    bad = B.score_outcome(
        _task("bad"), B.Trajectory(text="wrong", duration_s=1.0), model="m")
    out = B.score_outcome(
        _task("outage"),
        B.Trajectory(text=_RETRY_BANNERS, error=_NO_DEPLOYMENT), model="m")
    return [good, bad, out]


def test_the_rate_counts_only_what_was_measured():
    s = B.summarise_run(_rows())
    assert s["n_tasks"] == 2, "the outage is not one of the tasks measured"
    assert s["n_unmeasured"] == 1


def test_the_average_does_not_include_the_outage():
    rows = _rows()
    s = B.summarise_run(rows)
    scored = [r.quality_0_100 for r in rows if not r.unmeasured]
    assert s["avg_quality"] == pytest.approx(sum(scored) / len(scored))


def test_the_summary_names_which_task_was_not_measured():
    assert B.summarise_run(_rows())["unmeasured_tasks"] == ["outage"]


def test_a_run_of_nothing_but_outages_reports_no_rate_rather_than_zero():
    out = B.score_outcome(
        _task("a"), B.Trajectory(text=_RETRY_BANNERS, error=_NO_DEPLOYMENT),
        model="m")
    s = B.summarise_run([out])
    assert s["n_tasks"] == 0 and s["n_unmeasured"] == 1
    assert s["pass_rate"] == 0.0 and s["avg_quality"] == 0.0


def test_an_empty_run_still_has_the_keys():
    s = B.summarise_run([])
    assert s["n_unmeasured"] == 0 and s["unmeasured_tasks"] == []


def test_the_cost_and_duration_still_count_the_outage():
    """It cost wall-clock and possibly tokens; hiding that would be the
    opposite error."""
    rows = _rows()
    s = B.summarise_run(rows)
    assert s["total_duration_s"] == pytest.approx(
        sum(r.duration_s for r in rows))


def test_the_summary_line_says_it_out_loud():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "cli.py").read_text(encoding="utf-8")
    i = src.index("passed \"")
    block = src[i:i + 1500]
    assert "NOT MEASURED" in block
    assert "n_unmeasured" in block
