"""The fast suite runs as shards, and the required check proves they ran it all.

Input: .github/ci/ci_shards.py (plan, verify), .github/ci/shard_plugin.py
(the pytest side) and the `test-shard` / `test` jobs in ci.yml.
Output: the properties the required check `tests (py3.11)` rests on:

- every collected test file is assigned to exactly one shard, whatever the
  durations file knows, and the plan does not depend on input order;
- the shards' reports pass verification only when together they ran every
  collected test exactly once, for this commit, with exit status 0; every
  other case, a missing report included, fails;
- the plugin selects after ``-m`` deselection and leaves no shard variable
  behind for a test's subprocess;
- the workflow runs the aggregator whatever the shards concluded, fails it
  unless they all succeeded, and verifies the same shard count the matrix
  defines;
- the shard step's command is still the plain fast-suite command, which
  delfin.agent.verify_recipe reads as this project's test command.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CI_DIR = ROOT / ".github" / "ci"
FAST = 'pytest -q -m "not slow" --timeout=300 --timeout-method=signal'


def _load_ci_shards():
    spec = importlib.util.spec_from_file_location(
        "ci_shards_under_test", CI_DIR / "ci_shards.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cs = _load_ci_shards()


def _ci_jobs() -> dict:
    return yaml.safe_load(
        (ROOT / ".github" / "workflows" / "ci.yml").read_text(
            encoding="utf-8"))["jobs"]


# --------------------------------------------------------------------------
# the plan
# --------------------------------------------------------------------------

def _repo_test_files() -> list[str]:
    return sorted(p.relative_to(ROOT).as_posix()
                  for p in (ROOT / "tests").rglob("test_*.py"))


@pytest.mark.parametrize("shards", [1, 2, 3, 4, 5, 8])
def test_every_file_lands_in_exactly_one_shard(shards):
    durations = cs.load_durations()
    files = _repo_test_files() + ["tests/test_not_measured_yet.py"]
    plan = cs.plan(files, durations, shards)

    assert set(plan) == set(files)
    assert set(plan.values()) <= set(range(1, shards + 1))
    assert len(set(plan.values())) == min(shards, len(files))


def test_a_file_without_a_measured_duration_is_still_assigned():
    plan = cs.plan(["tests/test_a.py", "tests/test_new.py"],
                   {"tests/test_a.py": 10.0}, 4)
    assert set(plan) == {"tests/test_a.py", "tests/test_new.py"}


def test_the_plan_is_a_function_of_its_inputs_only():
    durations = cs.load_durations()
    files = _repo_test_files()
    a = cs.plan(files, durations, 4)
    b = cs.plan(list(reversed(files)), dict(reversed(durations.items())), 4)
    assert a == b
    assert cs.plan_digest(a) == cs.plan_digest(b)


def test_the_measured_durations_balance_four_shards():
    """702 s per shard on the 2026-10-15 measurement; a plan that leaves one
    shard 5 % above the mean wastes what the split was for."""
    durations = cs.load_durations()
    files = sorted(durations)
    plan = cs.plan(files, durations, 4)
    load = [sum(durations[f] for f in files if plan[f] == k)
            for k in range(1, 5)]
    assert max(load) <= 1.05 * (sum(load) / 4), load


# --------------------------------------------------------------------------
# verification of the shard reports
# --------------------------------------------------------------------------

_UNIVERSE = [f"tests/test_{c}.py::test_{i}" for c in "abcdef" for i in range(3)]


def _reports(root: Path, shards: int = 3, sha: str = "abc") -> Path:
    files = sorted({cs.file_of(n) for n in _UNIVERSE})
    plan = cs.plan(files, {}, shards)
    for k in range(1, shards + 1):
        d = root / f"tests-shard-{k}"
        d.mkdir(parents=True)
        sel = [n for n in _UNIVERSE if plan[cs.file_of(n)] == k]
        (d / "universe.txt").write_text("\n".join(sorted(_UNIVERSE)) + "\n")
        (d / "selected.txt").write_text("".join(n + "\n" for n in sel))
        (d / "ran.txt").write_text("".join(n + "\n" for n in sel))
        (d / "meta.json").write_text(json.dumps({
            "shard": k, "shards": shards, "sha": sha, "exitstatus": 0,
            "plan_digest": cs.plan_digest(plan)}))
    return root


def test_complete_reports_pass(tmp_path):
    problems, summary = cs.verify(_reports(tmp_path / "r"), 3, "abc")
    assert problems == []
    assert any("union of shards: 18" in line for line in summary)


def _meta(d: Path, **change):
    meta = json.loads((d / "meta.json").read_text())
    meta.update(change)
    (d / "meta.json").write_text(json.dumps(meta))


def _drop_line(path: Path, index: int = 0):
    lines = path.read_text().splitlines()
    del lines[index]
    path.write_text("".join(l + "\n" for l in lines))


def _append(path: Path, line: str):
    path.write_text(path.read_text() + line + "\n")


_TAMPER = {
    "a shard's artefact is missing":
        lambda r: shutil.rmtree(r / "tests-shard-2"),
    "every artefact is missing":
        lambda r: [shutil.rmtree(d) for d in list(r.iterdir())],
    "a report lacks its meta.json (the run died)":
        lambda r: (r / "tests-shard-1" / "meta.json").unlink(),
    "a report lacks ran.txt":
        lambda r: (r / "tests-shard-3" / "ran.txt").unlink(),
    "a report lacks universe.txt":
        lambda r: (r / "tests-shard-3" / "universe.txt").unlink(),
    "a report lacks selected.txt":
        lambda r: (r / "tests-shard-3" / "selected.txt").unlink(),
    "a shard index twice, one missing":
        lambda r: _meta(r / "tests-shard-2", shard=1),
    "a shard index outside 1..N":
        lambda r: _meta(r / "tests-shard-2", shard=7),
    "a report from another split":
        lambda r: _meta(r / "tests-shard-1", shards=4),
    "a report from another commit":
        lambda r: _meta(r / "tests-shard-1", sha="def"),
    "a shard whose pytest failed":
        lambda r: _meta(r / "tests-shard-3", exitstatus=1),
    "the shards computed different plans":
        lambda r: _meta(r / "tests-shard-3", plan_digest="0" * 64),
    "an unreadable meta.json":
        lambda r: (r / "tests-shard-1" / "meta.json").write_text("{"),
    "a selected test that never finished":
        lambda r: _drop_line(r / "tests-shard-1" / "ran.txt"),
    "a test that ran twice":
        lambda r: _append(r / "tests-shard-1" / "ran.txt",
                          (r / "tests-shard-1" / "ran.txt")
                          .read_text().splitlines()[0]),
    "a test nobody selected":
        lambda r: _drop_line(r / "tests-shard-2" / "selected.txt")
        or _drop_line(r / "tests-shard-2" / "ran.txt"),
    "a test selected by two shards":
        lambda r: [_append(r / "tests-shard-2" / n,
                           (r / "tests-shard-1" / "selected.txt")
                           .read_text().splitlines()[0])
                   for n in ("selected.txt", "ran.txt")],
    "a shard collected a different suite":
        lambda r: _append(r / "tests-shard-2" / "universe.txt",
                          "tests/test_z.py::test_0"),
    "an extra report directory":
        lambda r: shutil.copytree(r / "tests-shard-1", r / "stray"),
}


@pytest.mark.parametrize("case", sorted(_TAMPER))
def test_an_incomplete_split_fails(tmp_path, case):
    root = _reports(tmp_path / "r")
    _TAMPER[case](root)
    problems, _ = cs.verify(root, 3, "abc")
    assert problems, f"{case}: verify passed"


def test_no_download_at_all_fails(tmp_path):
    problems, _ = cs.verify(tmp_path / "never-created", 4, "abc")
    assert problems


def test_a_test_split_from_its_file_fails(tmp_path):
    """No file may run in two shards, even if no node id repeats."""
    root = _reports(tmp_path / "r")
    one = root / "tests-shard-1"
    moved = one.joinpath("selected.txt").read_text().splitlines()[0]
    _drop_line(one / "selected.txt")
    _drop_line(one / "ran.txt")
    two = root / "tests-shard-2"
    _append(two / "selected.txt", moved)
    _append(two / "ran.txt", moved)
    problems, _ = cs.verify(root, 3, "abc")
    assert any("split between" in p for p in problems), problems


def test_the_command_line_exits_1_on_an_incomplete_split(tmp_path):
    root = _reports(tmp_path / "r")
    ok = subprocess.run([sys.executable, str(CI_DIR / "ci_shards.py"),
                         "verify", str(root), "--shards", "3", "--sha", "abc"],
                        capture_output=True, text=True,
                        env={**os.environ, "GITHUB_STEP_SUMMARY": ""})
    assert ok.returncode == 0, ok.stdout + ok.stderr
    shutil.rmtree(root / "tests-shard-3")
    bad = subprocess.run([sys.executable, str(CI_DIR / "ci_shards.py"),
                          "verify", str(root), "--shards", "3", "--sha", "abc"],
                         capture_output=True, text=True,
                         env={**os.environ, "GITHUB_STEP_SUMMARY": ""})
    assert bad.returncode == 1
    assert "shard 3 has no report" in bad.stdout


# --------------------------------------------------------------------------
# the plugin, driven through pytest itself
# --------------------------------------------------------------------------

_CONFTEST = '''
import pytest
def pytest_configure(config):
    config.addinivalue_line("markers", "slow: heavy")
def pytest_collection_modifyitems(config, items):
    for item in items:
        if "heavy" in item.nodeid:
            item.add_marker(pytest.mark.slow)
'''

_ENV_TEST = '''
import os
def test_the_shard_left_nothing_behind():
    assert "PYTEST_ADDOPTS" not in os.environ
    assert ".github" not in os.environ.get("PYTHONPATH", "")
'''


def _tree(root: Path) -> Path:
    root.mkdir()
    (root / "conftest.py").write_text(_CONFTEST)
    for i in range(6):
        (root / f"test_f{i}.py").write_text(
            "import pytest\n"
            "@pytest.mark.parametrize('x', [1, 2])\n"
            "def test_a(x):\n    pass\n")
    (root / "test_heavy.py").write_text("def test_heavy():\n    pass\n")
    (root / "test_env.py").write_text(_ENV_TEST)
    return root


def _run_shard(tree: Path, k: int, n: int, report: Path, *extra: str):
    env = {**os.environ, "GITHUB_SHA": "abc",
           "PYTHONPATH": str(CI_DIR),
           "PYTEST_ADDOPTS": f"-p shard_plugin --shard={k}/{n} "
                             f"--shard-report={report}"}
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         *extra],
        cwd=tree, env=env, capture_output=True, text=True, timeout=120)


def test_three_shards_run_the_suite_once(tmp_path):
    tree = _tree(tmp_path / "suite")
    reports = tmp_path / "reports"
    for k in (1, 2, 3):
        proc = _run_shard(tree, k, 3, reports / f"tests-shard-{k}",
                          "-m", "not slow")
        assert proc.returncode == 0, proc.stdout + proc.stderr

    problems, _ = cs.verify(reports, 3, "abc")
    assert problems == []
    universe = (reports / "tests-shard-1" / "universe.txt").read_text()
    assert "test_heavy" not in universe, "selection ran before -m"
    assert universe.count("\n") == 13
    durations = json.loads(
        (reports / "tests-shard-1" / "durations.json").read_text())
    assert durations and all(f.startswith("test_") for f in durations)


def test_a_shard_that_saw_another_suite_fails_verification(tmp_path):
    """Shard 3 without -m: it collects the slow test too."""
    tree = _tree(tmp_path / "suite")
    reports = tmp_path / "reports"
    for k in (1, 2):
        _run_shard(tree, k, 3, reports / f"tests-shard-{k}", "-m", "not slow")
    _run_shard(tree, 3, 3, reports / "tests-shard-3")

    problems, _ = cs.verify(reports, 3, "abc")
    assert any("different suite" in p for p in problems), problems


def test_a_failing_test_reaches_the_report(tmp_path):
    tree = _tree(tmp_path / "suite")
    (tree / "test_f0.py").write_text("def test_a():\n    assert False\n")
    reports = tmp_path / "reports"
    for k in (1, 2):
        _run_shard(tree, k, 2, reports / f"tests-shard-{k}", "-m", "not slow")

    problems, _ = cs.verify(reports, 2, "abc")
    assert any("exit status 1" in p for p in problems), problems


def test_without_a_shard_the_plugin_does_nothing(tmp_path):
    tree = _tree(tmp_path / "suite")
    env = {**os.environ, "PYTHONPATH": str(CI_DIR),
           "PYTEST_ADDOPTS": "-p shard_plugin"}
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "-m", "not slow", "--deselect", "test_env.py"],
        cwd=tree, env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "12 passed" in proc.stdout, proc.stdout


# --------------------------------------------------------------------------
# the workflow
# --------------------------------------------------------------------------

def test_the_required_check_is_the_aggregator():
    jobs = _ci_jobs()
    named = [k for k, j in jobs.items() if j.get("name") == "tests (py3.11)"]
    assert named == ["test"], named
    agg = jobs["test"]
    assert agg["needs"] in ("test-shard", ["test-shard"])
    # Without !cancelled() a failed shard would skip the aggregator, and a
    # skipped required check counts as passed.
    assert agg["if"].replace(" ", "") == "${{!cancelled()}}"


def _step(job: dict, name_part: str) -> dict:
    hits = [s for s in job["steps"] if name_part in s.get("name", "")]
    assert len(hits) == 1, [s.get("name") for s in job["steps"]]
    return hits[0]


@pytest.mark.parametrize("result,passes", [
    ("success", True), ("failure", False), ("cancelled", False),
    ("skipped", False), ("", False)])
def test_the_aggregator_fails_unless_every_shard_succeeded(result, passes):
    step = _step(_ci_jobs()["test"], "Every shard succeeded")
    assert step["env"]["SHARDS_RESULT"].replace(" ", "") == \
        "${{needs.test-shard.result}}"
    proc = subprocess.run(["bash", "-e", "-c", step["run"]],
                          env={"PATH": os.environ["PATH"],
                               "SHARDS_RESULT": result},
                          capture_output=True, text=True)
    assert (proc.returncode == 0) is passes


def test_the_aggregator_verifies_the_split_the_matrix_defines():
    jobs = _ci_jobs()
    shard_job = jobs["test-shard"]
    matrix = shard_job["strategy"]["matrix"]
    # strategy.job-total is the shard count only while shard is the sole
    # matrix dimension.
    assert list(matrix) == ["shard"], matrix
    n = len(matrix["shard"])
    assert matrix["shard"] == list(range(1, n + 1))
    assert shard_job["strategy"]["fail-fast"] is False

    verify = _step(jobs["test"], "exactly once")
    assert f"--shards {n} " in verify["run"] + " "
    assert '--sha "$SHA"' in verify["run"]
    assert verify["env"]["SHA"].replace(" ", "") == "${{github.sha}}"
    assert verify["if"].replace(" ", "") == "${{!cancelled()}}"

    upload = _step(shard_job, "Keep the shard's report")
    assert upload["if"].replace(" ", "") == "${{always()}}"
    assert upload["uses"].startswith("actions/upload-artifact@")
    assert upload["with"]["name"] == "tests-shard-${{ matrix.shard }}"
    assert upload["with"]["if-no-files-found"] == "error"
    download = _step(jobs["test"], "Fetch the shard reports")
    assert download["uses"].startswith("actions/download-artifact@")
    assert download["with"]["pattern"] == "tests-shard-*"
    assert download["with"]["path"] in verify["run"]


def test_a_shard_runs_what_the_unsharded_job_ran():
    shard_job = _ci_jobs()["test-shard"]
    run = _step(shard_job, "Run fast test suite")
    assert run["run"].strip() == FAST
    addopts = run["env"]["PYTEST_ADDOPTS"]
    assert "-p shard_plugin" in addopts
    assert "--shard=${{ matrix.shard }}/${{ strategy.job-total }}" in addopts
    assert run["env"]["PYTHONPATH"].endswith("/.github/ci")

    install = _step(shard_job, "Install package + test deps")["run"]
    assert 'pip install -e ".[agent,docs,dev]" "numpy<2"' in install
    checkout = _step(shard_job, "Checkout")
    assert checkout["with"]["fetch-depth"] == 0
    assert "3.11" in str(_step(shard_job, "Set up Python")["with"])
    assert int(shard_job["timeout-minutes"]) <= 60


def test_the_project_test_command_is_still_the_fast_suite():
    """verify_recipe renders the first CI test step to the agent as the
    command to run on changed files."""
    from delfin.agent import verify_recipe

    tests = [s for s in verify_recipe.discover(ROOT).steps if s.kind == "test"]
    assert tests and tests[0].command == FAST, tests[:1]
