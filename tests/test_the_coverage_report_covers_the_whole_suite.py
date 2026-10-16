"""The nightly coverage figure is over the whole fast suite, or it is red.

Input: the `coverage-shard` and `coverage` jobs in .github/workflows/ci.yml,
and the `test-shard` / `test` jobs they are built after.

Output: the properties that make the combined figure mean "the fast suite":
- the coverage shards split the suite the way test-shard does (same matrix,
  same plugin options, the same pytest command plus --cov), so
  ci_shards.py's plan assigns every file to exactly one of them;
- each shard keeps its data file AND its shard report, also after a
  failure;
- `coverage` runs whatever the shards concluded (otherwise a failed shard
  skips it and the night has no report), combines the data, and ends with
  `ci_shards.py verify` over the same shard count, so a missing shard, an
  unfinished test or a failing one turns it red;
- nothing in it runs on a push or a pull request, and no job outlives the
  gate's 60-minute ceiling (test_the_gate_can_still_say_it_passed.py).

Before: one job ran the fast suite under tracing in 60 minutes while the
suite alone measured 47.7 minutes of tests (2026-10-15), with
`|| true` after pytest: at best a figure, at worst a timeout and none.
"""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _jobs() -> dict:
    return yaml.safe_load((ROOT / ".github" / "workflows" / "ci.yml")
                          .read_text(encoding="utf-8"))["jobs"]


def _step(job: dict, needle: str) -> dict:
    return next(s for s in job["steps"] if needle in s.get("name", ""))


_NIGHTLY = "github.event_name == 'schedule' || github.event_name == 'workflow_dispatch'"


def test_the_coverage_shards_split_the_suite_like_the_gate():
    jobs = _jobs()
    gate, cov = jobs["test-shard"], jobs["coverage-shard"]
    assert cov["strategy"]["matrix"] == gate["strategy"]["matrix"]
    assert cov["strategy"]["fail-fast"] is False

    gate_run = _step(gate, "Run fast test suite")
    cov_run = _step(cov, "Measure coverage")
    assert cov_run["env"] == gate_run["env"]
    assert " ".join(cov_run["run"].split()) == (
        " ".join(gate_run["run"].split()) + " --cov=delfin --cov-report=")
    assert '"numpy<2"' in _step(cov, "Install")["run"]
    assert int(cov["timeout-minutes"]) <= 60


def test_each_shard_keeps_its_data_and_its_report_after_a_failure():
    cov = _jobs()["coverage-shard"]
    keep = _step(cov, "Keep the shard's data")
    upload = _step(cov, "Upload the shard's data")
    for step in (keep, upload):
        assert step["if"].replace(" ", "") == "${{always()}}"
    assert "coverage-data/" in upload["with"]["path"]
    assert "shard-report/" in upload["with"]["path"]
    assert upload["with"]["name"] == "coverage-shard-${{ matrix.shard }}"


def test_the_combined_report_is_red_unless_the_whole_suite_was_measured():
    jobs = _jobs()
    combine = jobs["coverage"]
    assert combine["needs"] == "coverage-shard"
    assert "!cancelled()" in combine["if"]
    assert "coverage combine" in _step(combine, "Combine")["run"]

    verify = combine["steps"][-1]
    n = len(jobs["coverage-shard"]["strategy"]["matrix"]["shard"])
    assert "ci_shards.py verify" in verify["run"]
    assert f"--shards {n} " in verify["run"] + " "
    assert '--sha "$SHA"' in verify["run"]
    assert verify["if"].replace(" ", "") == "${{!cancelled()}}"
    for name in ("coverage-shard", "coverage"):
        assert "|| true" not in " ".join(s.get("run", "") for s in jobs[name]["steps"])


def test_coverage_runs_only_at_night_or_on_demand():
    jobs = _jobs()
    for name in ("coverage-shard", "coverage"):
        assert _NIGHTLY in jobs[name]["if"], jobs[name]["if"]
