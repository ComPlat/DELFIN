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

# The contexts the `main-required-checks` ruleset requires on main, by the
# workflow that produces them (docs/MAIN_PROTECTION.md). Renaming a job
# renames its check and silently drops the requirement.
_REQUIRED_CHECKS = {
    "ci.yml": {
        "tests (py3.11)",
        "lint (ruff)",
        "container isolation (docker, podman)",
        "generalist toolchains (node, cargo, go)",
        "reference environment (py3.12, numpy 2.2)",
    },
    "codeql.yml": {"analyze (python)", "analyze (actions)"},
    "dependency-review.yml": {"dependency review"},
    "tui.yml": {"window and service"},
}


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
# the checks the ruleset names must keep existing, and keep reporting
# --------------------------------------------------------------------------

def _rendered_check_names(workflow: str) -> set[str]:
    """Job names as GitHub renders them: one per matrix value of every
    `${{ matrix.<key> }}` the name uses."""
    rendered = set()
    for key, job in _jobs(workflow).items():
        names = [job.get("name", key)]
        matrix = job.get("strategy", {}).get("matrix", {})
        for var, values in matrix.items():
            token = f"${{{{ matrix.{var} }}}}"
            if isinstance(values, list):
                names = [n.replace(token, str(v)) if token in n else n
                         for n in names for v in values]
        rendered.update(names)
    return rendered


@pytest.mark.parametrize("workflow", sorted(_REQUIRED_CHECKS))
def test_the_required_checks_are_still_produced(workflow):
    """The ruleset matches on the rendered job name. Renaming one drops the
    requirement without any error anywhere -- main simply stops being
    protected by it."""
    rendered = _rendered_check_names(workflow)
    missing = _REQUIRED_CHECKS[workflow] - rendered
    assert not missing, (
        f"{sorted(missing)} is required on main but no job in {workflow} "
        f"produces it any more. Produced: {sorted(rendered)}"
    )


@pytest.mark.parametrize("workflow", sorted(_REQUIRED_CHECKS))
def test_a_required_check_reports_on_every_pull_request(workflow):
    """A required check that never starts is not a failed check, it is a
    pending one, and GitHub waits for it forever: a `paths:` filter on
    pull_request keeps every pull request outside those paths from merging.
    A workflow that applies to part of the tree decides inside a job and
    skips with `if:` -- a skipped job counts as passed."""
    on = _load(workflow)[True]  # YAML parses the key `on:` as True
    assert "pull_request" in on, f"{workflow} does not run on pull requests"
    pr = on["pull_request"] or {}
    filtered = sorted(k for k in ("paths", "paths-ignore", "branches",
                                  "branches-ignore") if k in pr)
    assert not filtered, (
        f"{workflow} filters pull_request by {filtered}; its required checks "
        "would never report on a pull request outside the filter"
    )


_DOC = _WORKFLOWS.parents[1] / "docs" / "MAIN_PROTECTION.md"


def _documented_checks() -> list[tuple[str, str]]:
    """(check, workflow) for every table row ``| `check` | wf.yml |`` and
    every bullet ``- `check` (wf.yml...`` in docs/MAIN_PROTECTION.md."""
    import re
    out = []
    for line in _DOC.read_text(encoding="utf-8").splitlines():
        row = re.match(r"\|\s*(`[^|]+`)\s*\|\s*([\w-]+\.yml)\s*\|", line)
        bullet = re.match(r"-\s*`([^`]+)`\s*\(([\w-]+\.yml)\b", line)
        if row:
            for check in re.findall(r"`([^`]+)`", row.group(1)):
                out.append((check, row.group(2)))
        elif bullet:
            out.append((bullet.group(1), bullet.group(2)))
    return out


def test_the_protection_document_names_the_workflow_of_each_check():
    """The document is the reference a maintainer reads when a check is
    missing; a check filed under the wrong workflow (the dependency audit
    moved from ci.yml to audit.yml) sends them to a file that does not
    produce it. A name ending in ``(...)`` stands for its matrix."""
    documented = _documented_checks()
    assert ("dependency audit", "audit.yml") in documented, documented
    wrong = []
    for check, workflow in documented:
        rendered = _rendered_check_names(workflow)
        if check.endswith("(...)"):
            prefix = check[:-len("(...)")]
            ok = any(r.startswith(prefix) for r in rendered)
        else:
            ok = check in rendered
        if not ok:
            wrong.append((check, workflow, sorted(rendered)))
    assert not wrong, wrong
    for workflow, checks in _REQUIRED_CHECKS.items():
        for check in checks:
            assert (check, workflow) in documented, (
                f"{check} is required ({workflow}) and not in the document")


def test_the_window_tests_are_skipped_only_on_a_definite_no():
    """A job whose `needs` failed is skipped, and a skipped required check
    passes. So the window tests must run when the job deciding whether they
    apply did not succeed, and be skipped only when it said "not covered"."""
    cond = _jobs("tui.yml")["test"].get("if", "")
    assert "!cancelled()" in cond, cond
    assert "needs.changes.result != 'success'" in cond, cond


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
