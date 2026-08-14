"""A run conclusion is only useful while one thing decides it.

CI is the line anybody glances at to see whether main is healthy, and for six
commits it could not answer. The heavy suite ran as a job inside CI and
cancelled its own predecessor on every push — right for a ninety-minute
informational run, since the verdict on a commit five pushes back is not
worth a runner. But GitHub marks the whole RUN cancelled when any job in it
is cancelled, `continue-on-error` included. So CI on main read "cancelled"
every time while `tests` and `lint` were passing every time, and the two
states were indistinguishable from the outside.

The fix was to separate the audiences: CI answers "may this merge", and
slow-tests.yml answers "does the heavy suite still hold". Neither can
overwrite the other's answer any more.

That is a structural property, not a line of code, and the natural way to
undo it is to add one more job to CI because it belongs there thematically.
These tests hold the shape: whatever gates stays gating, and nothing that
cancels or is allowed to fail sits in the file whose conclusion people read.
"""

from pathlib import Path

import pytest
import yaml

_WORKFLOWS = Path(__file__).resolve().parents[1] / ".github" / "workflows"

# The contexts GitHub branch protection requires on main. Renaming a job
# renames its check and silently drops the requirement.
_REQUIRED_CHECKS = {"tests (py3.11)", "lint (ruff)"}


def _load(name: str) -> dict:
    path = _WORKFLOWS / name
    assert path.exists(), f"{name} is missing"
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _jobs(name: str) -> dict:
    return _load(name).get("jobs", {})


# --------------------------------------------------------------------------
# nothing in the gating workflow may colour its conclusion
# --------------------------------------------------------------------------

def test_no_job_in_ci_cancels_a_sibling_run():
    """A job-level concurrency group with cancel-in-progress is what did it.
    The workflow-level group on CI is a different thing and stays: it is
    keyed per commit on main and never cancels."""
    offenders = {
        name: job["concurrency"]
        for name, job in _jobs("ci.yml").items()
        if isinstance(job.get("concurrency"), dict)
        and job["concurrency"].get("cancel-in-progress") in (True, "true")
    }

    assert not offenders, (
        f"{sorted(offenders)} would make CI's run conclusion read 'cancelled' "
        "on a push that supersedes it, whatever the required checks did. A job "
        "that needs to cancel its predecessor belongs in its own workflow."
    )


def test_no_job_in_ci_is_allowed_to_fail():
    """continue-on-error says "this result does not count". A result that
    does not count does not belong beside the ones that do — it is reporting,
    and reporting has its own file."""
    allowed = [n for n, j in _jobs("ci.yml").items()
               if j.get("continue-on-error") in (True, "true")]

    assert not allowed, (
        f"{allowed} is reporting, not gating, and sits in the workflow whose "
        "conclusion is read as the gate's verdict."
    )


def test_ci_has_no_job_long_enough_to_outlive_the_gate():
    """The gate is waited on. A four-hour job in the same file makes the run
    take four hours to reach a conclusion, whatever the fast suite said."""
    too_long = {n: j["timeout-minutes"] for n, j in _jobs("ci.yml").items()
                if int(j.get("timeout-minutes", 0)) > 60}

    assert not too_long, f"{too_long} do not belong in the gating workflow"


# --------------------------------------------------------------------------
# the checks branch protection names must keep existing
# --------------------------------------------------------------------------

def test_the_required_checks_are_still_produced_by_ci():
    """Branch protection matches on the rendered job name. Renaming one drops
    the requirement without any error anywhere — main simply stops being
    protected by it."""
    rendered = set()
    for name, job in _jobs("ci.yml").items():
        display = job.get("name", name)
        versions = (job.get("strategy", {})
                       .get("matrix", {})
                       .get("python-version", []))
        if "${{ matrix.python-version }}" in display:
            for v in versions:
                rendered.add(display.replace("${{ matrix.python-version }}", str(v)))
        else:
            rendered.add(display)

    missing = _REQUIRED_CHECKS - rendered
    assert not missing, (
        f"{sorted(missing)} is required on main but no job in ci.yml produces "
        f"it any more. Produced: {sorted(rendered)}"
    )


# --------------------------------------------------------------------------
# the heavy suite still runs somewhere
# --------------------------------------------------------------------------

def test_the_slow_suite_did_not_just_disappear():
    """Splitting it out and forgetting to run it would look identical from
    CI's side: green, quiet, and testing half of what it used to."""
    jobs = _jobs("slow-tests.yml")
    assert jobs, "slow-tests.yml defines no job"

    steps = [s for job in jobs.values() for s in job.get("steps", [])]
    runs = " ".join(s.get("run", "") for s in steps)

    assert '-m "slow"' in runs or "-m 'slow'" in runs, (
        "no step runs the slow-marked tests"
    )


def test_the_slow_suite_still_runs_without_anybody_asking():
    on = _load("slow-tests.yml")[True]  # YAML parses the key `on:` as True

    assert "schedule" in on, "nightly run is gone"
    assert "workflow_dispatch" in on, "cannot be started by hand"


def test_the_slow_suite_reports_what_it_found():
    """It carried continue-on-error while it lived in CI, for one reason: to
    keep a failure from colouring the gating run. Nothing gates on it now, so
    the flag would only guarantee a green tick over a red result — which is
    the state it spent months in, when it was schedule-only and its silence
    read like health."""
    for name, job in _jobs("slow-tests.yml").items():
        assert job.get("continue-on-error") not in (True, "true"), (
            f"{name} would report success whatever the suite did"
        )


def test_the_slow_suite_keeps_one_run_at_a_time():
    """Without this it stacks: main takes a push every few minutes and a
    ninety-minute job piles up five deep within the hour."""
    wf = _load("slow-tests.yml")
    conc = wf.get("concurrency")

    assert isinstance(conc, dict), "no concurrency group"
    assert conc.get("cancel-in-progress") in (True, "true")


def test_a_slow_test_has_a_deadline_of_its_own():
    """One call joins a thread with timeout=None. Without a per-test deadline
    it holds the job until the runner's limit and pytest's failure list —
    printed only at the end — never appears: the run dies with the names
    still unknown."""
    steps = [s for job in _jobs("slow-tests.yml").values()
             for s in job.get("steps", [])]
    runs = " ".join(s.get("run", "") for s in steps)

    assert "--timeout=" in runs, "a hang would be a stuck run, not a named failure"
    assert "--timeout-method=signal" in runs, (
        "the thread method takes the process down with it, losing every "
        "remaining test as well"
    )


# --------------------------------------------------------------------------
# both files stay loadable
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(p.name for p in _WORKFLOWS.glob("*.yml")))
def test_every_workflow_parses_and_declares_its_permissions(name):
    wf = _load(name)

    assert wf.get("permissions") is not None, (
        f"{name} inherits whatever the repository default grants; say it here"
    )
    assert wf.get("jobs"), f"{name} defines no job"
