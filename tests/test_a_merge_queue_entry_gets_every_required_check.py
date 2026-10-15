"""A merge-queue entry receives every check the ruleset requires.

Input: .github/workflows/*.yml and the required-check list held in
tests/test_the_gate_can_still_say_it_passed.py.
Output: the properties a merge queue rests on. The queue builds a commit on
a gh-readonly-queue/main/* ref, fires `merge_group`, and waits for every
required check on that commit; a check that never reports holds the entry
until check_response_timeout_minutes and then removes it. So:

- every required check is produced by a workflow that has a `merge_group`
  trigger, by a job not switched off for that event;
- `merge_group` carries no filter (GitHub knows only `types` there);
- no concurrency group cancels a queue run (ci.yml groups it by commit);
- the window tests run for a queue entry: tui.yml either runs them
  unconditionally or, with its `changes` job, decides "covered" for any
  event other than pull_request;
- the `binary` matrix (not required, macOS runners) stays out of the queue.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"


def _required_checks() -> set[str]:
    """The names the gate test holds, whatever shape it keeps them in (a
    set of names, or names per workflow)."""
    path = Path(__file__).with_name("test_the_gate_can_still_say_it_passed.py")
    spec = importlib.util.spec_from_file_location("gate_checks", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    req = mod._REQUIRED_CHECKS
    if isinstance(req, dict):
        return set().union(*req.values())
    return set(req)


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _rendered(job_key: str, job: dict) -> set[str]:
    names = [job.get("name", job_key)]
    matrix = job.get("strategy", {}).get("matrix", {}) or {}
    for var, values in matrix.items():
        token = f"${{{{ matrix.{var} }}}}"
        if isinstance(values, list):
            names = [n.replace(token, str(v)) if token in n else n
                     for n in names for v in values]
    return set(names)


def _producers() -> dict[str, tuple[str, str, dict]]:
    out = {}
    for path in sorted(WORKFLOWS.glob("*.yml")):
        wf = _load(path)
        for key, job in (wf.get("jobs") or {}).items():
            for name in _rendered(key, job):
                out[name] = (path.name, key, job)
    return out


@pytest.mark.parametrize("check", sorted(_required_checks()))
def test_a_required_check_reports_on_a_queue_entry(check):
    producers = _producers()
    assert check in producers, f"no workflow produces {check!r}"
    workflow, key, job = producers[check]
    on = _load(WORKFLOWS / workflow)[True]  # YAML reads `on:` as True
    assert "merge_group" in on, (
        f"{check!r} comes from {workflow}, which does not run on "
        "merge_group: a queue entry would wait for it until the timeout")
    assert "merge_group" not in str(job.get("if", "")), (
        f"{workflow}:{key} is conditioned on merge_group")


def _queue_workflows() -> list[str]:
    producers = _producers()
    return sorted({producers[c][0] for c in _required_checks()
                   if c in producers})


@pytest.mark.parametrize("workflow", _queue_workflows())
def test_the_queue_trigger_is_not_filtered(workflow):
    trigger = _load(WORKFLOWS / workflow)[True].get("merge_group") or {}
    assert set(trigger) <= {"types"}, trigger
    assert trigger.get("types", ["checks_requested"]) == ["checks_requested"]


@pytest.mark.parametrize("workflow", _queue_workflows())
def test_no_concurrency_group_cancels_a_queue_run(workflow):
    wf = _load(WORKFLOWS / workflow)
    conc = wf.get("concurrency")
    if isinstance(conc, dict):
        cancel = str(conc.get("cancel-in-progress", False))
        assert cancel in ("False", "false") or \
            "github.event_name != 'merge_group'" in cancel, cancel
    for key, job in wf["jobs"].items():
        jc = job.get("concurrency")
        assert not (isinstance(jc, dict)
                    and jc.get("cancel-in-progress") not in (None, False)), (
            f"{workflow}:{key} has a job-level cancelling group")


def test_a_queue_run_of_ci_is_grouped_by_its_commit():
    group = _load(WORKFLOWS / "ci.yml")["concurrency"]["group"]
    assert "github.event_name == 'merge_group' && format('queue-{0}', " \
           "github.sha)" in group, group


def test_a_queue_entry_runs_the_window_tests(tmp_path):
    jobs = _load(WORKFLOWS / "tui.yml")["jobs"]
    if "changes" not in jobs:
        assert "if" not in jobs["test"], jobs["test"].get("if")
        return
    # With the `changes` job: drive its script for a merge_group event on a
    # commit that touches nothing the window tests cover.
    step = next(s for s in jobs["changes"]["steps"] if s.get("id") == "diff")
    repo = tmp_path / "repo"
    repo.mkdir()
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c",
           "user.email=t@example.invalid"]
    subprocess.run(git[:3] + ["init", "-q"], check=True)
    for text in ("a\n", "b\n"):
        (repo / "README.md").write_text(text)
        subprocess.run(git + ["add", "README.md"], check=True)
        subprocess.run(git + ["commit", "-q", "-m", "x"], check=True)
    out = tmp_path / "out"
    env = {"PATH": os.environ["PATH"], "GITHUB_OUTPUT": str(out),
           "HOME": str(tmp_path)}
    for event, expected in (("merge_group", "true"),
                            ("pull_request", "false")):
        out.write_text("")
        subprocess.run(["bash", "-e", "-c", step["run"]], cwd=repo,
                       env={**env, "EVENT": event}, check=True,
                       capture_output=True)
        assert f"covered={expected}" in out.read_text(), (event,
                                                          out.read_text())


def test_the_binary_builds_stay_out_of_the_queue():
    binary = _load(WORKFLOWS / "tui.yml")["jobs"]["binary"]
    assert binary["if"].replace(" ", "") == \
        "${{github.event_name!='merge_group'}}"
