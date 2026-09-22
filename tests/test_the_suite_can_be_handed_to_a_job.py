"""`run_tests` can hand the suite over, instead of only saying to.

The previous pass made it REFUSE a run longer than the shell's ceiling
and name `bash_background`. That was right and incomplete: the caller
then has to build the pytest invocation itself, which is the one thing
this tool exists to know.

It now starts the run as a job and hands back the id — on request, and
always above the ceiling. The command is the SAME one the foreground
path runs, because test_runner owns it (`build_command`). A second
spelling is how a producer and a renderer came to disagree elsewhere in
this codebase, and the cost was six live wake-ups that named nothing.

The other half of this file is the idle prompt's tick. It was 0.1 s —
ten wake-ups a second where readline blocked. `read_ready` is a select:
a key returns at once and a signal interrupts the wait, so the timeout
governs nothing but how often an idle prompt wakes to do nothing.
"""

from __future__ import annotations

import ast
import inspect
import json

import pytest

from delfin.agent import test_runner as TR
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


# -- one command, two ways to run it ----------------------------------------

def test_the_command_is_built_in_one_place(tmp_path):
    cmd = TR.build_command(tmp_path, target="tests/test_x.py",
                           pytest_args=["-q", "--tb=line"])
    assert cmd[1:3] == ["-m", "pytest"]
    assert "tests/test_x.py" in " ".join(cmd)
    assert "--tb=line" in cmd
    assert cmd[-1] == "-q"


def test_the_foreground_path_uses_it():
    """Not a copy: run_tests calls the same builder."""
    src = inspect.getsource(TR.run_tests)
    assert "build_command(" in src


def test_the_background_path_uses_it_too():
    from delfin.agent import api_client as A
    src = inspect.getsource(A._DocToolExecutor._run_tests_in_background)
    assert "build_command(" in src
    assert "get_registry().start(" in src


# -- handing it over --------------------------------------------------------

@pytest.fixture()
def run_tests(tmp_path, monkeypatch):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    started: dict = {}

    class _Job:
        job_id = "01e5b151"

    class _Registry:
        def start(self, command, **kw):
            started["command"] = command
            started["kw"] = kw
            return _Job()

    from delfin.agent import bash_jobs as bj
    monkeypatch.setattr(bj, "get_registry", lambda: _Registry())

    def _call(**args):
        return json.loads(eng._execute_run_tests(dict(args), perms)), started
    return _call


def test_asking_for_background_starts_a_job(run_tests):
    out, started = run_tests(pytest_args=["-q"], background=True)
    assert out.get("status") == "started"
    assert out.get("job_id") == "01e5b151"
    assert "pytest" in started["command"]
    assert "do not wait on it" in out.get("note", "").lower()
    assert "bash_output" in out.get("note", "")


def test_a_run_over_the_ceiling_goes_there_without_asking(run_tests):
    out, started = run_tests(pytest_args=["-q"], timeout_s=1800)
    assert out.get("status") == "started", out
    assert out.get("requested_timeout_s") == 1800
    assert started["command"], "the job really was started"


def test_a_short_run_still_runs_here(run_tests):
    out, started = run_tests(target="nothing_here.py", pytest_args=["-q"],
                             timeout_s=60)
    assert out.get("status") != "started"
    assert not started, "a short run must not be handed away"


def test_a_registry_that_cannot_start_says_so(tmp_path, monkeypatch):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms
    from delfin.agent import bash_jobs as bj

    class _Bad:
        def start(self, *a, **k):
            raise RuntimeError("no room")
    monkeypatch.setattr(bj, "get_registry", lambda: _Bad())
    out = json.loads(eng._execute_run_tests(
        {"pytest_args": ["-q"], "background": True}, perms))
    assert "could not start" in out.get("error", "")


# -- the idle tick ----------------------------------------------------------

def test_the_idle_prompt_does_not_wake_ten_times_a_second():
    from delfin.agent import repl as R

    tree = ast.parse(inspect.getsource(R))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "read_boxed")
    src = ast.unparse(fn)
    assert "_IDLE_TICK_S" in src, "the tick is a named number, not a literal"
    ticks = [n.value for n in ast.walk(tree)
             if isinstance(n, ast.Assign)
             and getattr(n.targets[0], "id", "") == "_IDLE_TICK_S"]
    assert ticks and ticks[0].value >= 0.5, (
        "0.1 s is ten wake-ups a second where readline blocked")
