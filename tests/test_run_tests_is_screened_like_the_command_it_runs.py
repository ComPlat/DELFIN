"""Under Bypass, run_tests is screened like the caged pytest it runs.

run_tests starts pytest in the cage the shell runs in (``_test_run_argv``).
The same run spelled as bash passes unasked under Bypass where the cage
holds it (``_bypass_screened_action``); spelled as run_tests it went to the
action watcher, and one refusal ended a test run the task asked for.
Measured 2026-11-15 in two unattended swarm runs: "Running existing tests
does not match the task" and "Testing unrelated modules", both runs their
brief named.

The whole-suite run that run_tests hands to a background job ran with the
argv as built, outside every cage; it now gets the shell's cage too.

Driven through the real gate with a watcher that refuses everything, so a
pass means the model was not asked.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ
from delfin.agent.api_client import KitToolPermissions, _doc_executor

REFUSE = "DECISION: refuse\nREASON: nobody asked for this."
TASK = [{"role": "user", "content": "move the helpers and run the tests"}]


class Runner:
    def __init__(self):
        self.prompts = []

    def __call__(self, prompt, model):
        self.prompts.append(prompt)
        return REFUSE


def _settings(monkeypatch, broad=False):
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": "m", "rules": [], "broad": broad})


def _session(tmp_path, monkeypatch, *, cage="bwrap", mode="bypassPermissions"):
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (cage, "auto"))
    ws = tmp_path / "wt"
    ws.mkdir(exist_ok=True)
    perms = KitToolPermissions(workspace=str(ws), mode=mode)
    perms.classifier_runner = Runner()
    perms.classifier_messages = TASK
    perms.confirm_callback = None
    return perms


ARGS = [{"target": "tests/test_shell_words_leaf_identity.py"},
        {"target": "tests/test_a.py", "pytest_args": ["tests/test_b.py", "-q"]}]


@pytest.mark.parametrize("args", ARGS)
def test_a_caged_test_run_needs_no_model_under_bypass(tmp_path, monkeypatch, args):
    _settings(monkeypatch)
    perms = _session(tmp_path, monkeypatch)
    assert _doc_executor._run_permission_gate("run_tests", dict(args), perms) is None
    assert perms.classifier_runner.prompts == []


def test_the_same_run_as_bash_was_already_screened_without_the_model(tmp_path,
                                                                    monkeypatch):
    """The control: the shell spelling of the run, unchanged."""
    _settings(monkeypatch)
    perms = _session(tmp_path, monkeypatch)
    cmd = "python -m pytest -q tests/test_shell_words_leaf_identity.py"
    assert _doc_executor._run_permission_gate("bash", {"command": cmd}, perms) is None
    assert perms.classifier_runner.prompts == []


def test_without_a_cage_the_model_still_judges(tmp_path, monkeypatch):
    _settings(monkeypatch)
    perms = _session(tmp_path, monkeypatch, cage="none")
    assert _doc_executor._run_permission_gate("run_tests", dict(ARGS[0]), perms)
    assert perms.classifier_runner.prompts


def test_the_broad_scope_keeps_the_model_for_a_test_run(tmp_path, monkeypatch):
    """The user asking for a second opinion on everything."""
    _settings(monkeypatch, broad=True)
    perms = _session(tmp_path, monkeypatch)
    assert _doc_executor._run_permission_gate("run_tests", dict(ARGS[0]), perms)
    assert perms.classifier_runner.prompts


def test_the_whole_suite_job_runs_in_the_shells_cage(tmp_path, monkeypatch):
    perms = _session(tmp_path, monkeypatch)
    perms.task_session_id = "sess-suite"
    started = {}

    class _Job:
        job_id = "j1"

    class _Registry:
        def start(self, command, **kw):
            started.update(kw, command=command)
            return _Job()

    monkeypatch.setattr(BJ, "get_registry", lambda: _Registry())
    monkeypatch.setattr(A, "_bash_isolation_argv",
                        lambda cmd, cwd, p, **kw: ["CAGED", cmd])
    out = json.loads(_doc_executor._run_tests_in_background(
        "", [], 300, 600, {}, perms))
    assert out["status"] == "started", out
    assert started["argv"] == ["CAGED", started["command"]], started
    assert "pytest" in started["command"]
    assert started["session_id"] == "sess-suite"
