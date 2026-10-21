"""A whole test suite is not waited out in the foreground.

Measured on the shared bug archive, 2026-09-16 .. 10-14: 38 whole-suite
invocations in 8 sessions (bash 17, bash_background 16, run_tests with
no file 5), at least 144 minutes of blocking tool time, 7 cut at the
600 s ceiling. One session ran the same suite five times in a row with
``timeout`` 540, 580, 560, 550 in front of it; the repeat check compared
the exact text and saw five different commands. A user wrote three
messages while a 30-minute ``run_tests`` held the turn, and none of them
was read until it ended.

  a whole-suite bash run          started as a watched background job
  a whole-suite run_tests run     the same, with a note naming the
                                  foreground alternative (a target)
  a small suite, a file, -k       untouched: runs in the foreground
  a sub-agent                     untouched
  timeout N / nice / env / X=1    stripped before a repeat is compared
  a command that timed out        not re-run in the foreground in the
                                  same turn under another prefix
  a user message during a run     ends the run after a grace period and
                                  says why
"""

from __future__ import annotations

import json
import subprocess
import time
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as bj
from delfin.agent import contained_run
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _suite(root, n):
    tests = root / "tests"
    tests.mkdir(exist_ok=True)
    for i in range(n):
        (tests / f"test_m{i}.py").write_text(
            "def test_x():\n    assert True\n", encoding="utf-8")
    return root


class _Job:
    job_id = "0b5e7a11"

    class proc:  # noqa: N801 - mirrors the registry job's attribute
        pid = 4242


@pytest.fixture
def started(monkeypatch, tmp_path):
    """The registry records what would have been started; nothing runs."""
    calls: list[dict] = []

    class _Registry:
        def start(self, command=None, **kw):
            calls.append({"command": command, **kw})
            return _Job()

        def get(self, job_id, ws=None):
            return _Job()

    monkeypatch.setattr(bj, "get_registry", lambda: _Registry())
    from delfin.agent import job_monitor as jm
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH",
                        tmp_path / "watch-index.json")
    return calls


@pytest.fixture
def executor(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms
    return eng, perms


# -- recognising the run -----------------------------------------------------

FIELD_SPELLINGS = [
    "python -m pytest tests/ -q",
    "timeout 580 python -m pytest tests/ -q --tb=line --deselect "
    "tests/test_a.py::test_b 2>&1 | tail",
    "python3 -m pytest tests/ -q -x --ignore=tests/test_slurm 2>/dev/null "
    "| tail -3; echo \"EXIT=${PIPESTATUS[0]}\"",
    "set -o pipefail; python -m pytest tests/ -q --tb=line -x 2>&1 | tail -4",
    "pytest -q",
]

NOT_THE_SUITE = [
    "python -m pytest tests/test_m1.py -q",
    "python -m pytest tests/test_m1.py::test_x",
    "python -m pytest tests/ -q -k co2",
    "python -m pytest tests/ -q --collect-only",
    "grep -rn pytest tests/",
]


@pytest.mark.parametrize("cmd", FIELD_SPELLINGS)
def test_each_field_spelling_is_a_whole_suite_run(tmp_path, cmd):
    from delfin.agent import long_runs as L
    _suite(tmp_path, 45)
    assert L.suite_files_in_command(cmd, tmp_path) >= L.LARGE_SUITE_MIN_FILES


@pytest.mark.parametrize("cmd", NOT_THE_SUITE)
def test_a_chosen_part_is_not(tmp_path, cmd):
    from delfin.agent import long_runs as L
    _suite(tmp_path, 45)
    assert L.suite_files_in_command(cmd, tmp_path) == 0


def test_a_small_suite_keeps_the_foreground(tmp_path):
    from delfin.agent import long_runs as L
    _suite(tmp_path, 3)
    assert L.suite_files_in_command("python -m pytest tests/ -q", tmp_path) == 0


@pytest.mark.parametrize("a,b", [
    ("timeout 540 python -m pytest tests/ -q",
     "timeout 580 python -m pytest tests/ -q 2>&1 | tail -5"),
    ("nice -n 10 make check", "make check"),
    ("FOO=1 BAR=2 make check", "make check"),
    ("env LC_ALL=C make check", "make check"),
    ("set -o pipefail; make check 2>&1 | tail -3", "make check"),
    ("cd /w && timeout 30 make check", "make check"),
])
def test_waiting_prefixes_do_not_make_another_command(a, b):
    from delfin.agent import long_runs as L
    assert L.strip_command_prefixes(a) == L.strip_command_prefixes(b)


def test_a_different_command_stays_different():
    from delfin.agent import long_runs as L
    assert (L.strip_command_prefixes("timeout 30 make check")
            != L.strip_command_prefixes("timeout 30 make lint"))


def test_the_repeat_guard_sees_the_respelled_command():
    """Three failures under three timeouts are one failing command."""
    from delfin.agent.action_grounding import FailureBudget
    budget = FailureBudget()
    for t in (540, 580, 560):
        budget.record("bash", {"command": f"timeout {t} python -m pytest "
                                          "tests/ -q 2>&1 | tail"}, ok=False)
    assert budget.repeat_hint(
        "bash", {"command": "timeout 550 python -m pytest tests/ -q"}) is not None


# -- run_tests ----------------------------------------------------------------

def test_run_tests_over_a_large_suite_is_handed_to_a_watched_job(
        executor, started, tmp_path, monkeypatch):
    from delfin.agent import test_runner
    monkeypatch.setattr(test_runner, "run_tests", MagicMock(
        side_effect=AssertionError("the whole suite ran in the foreground")))
    _suite(tmp_path, 45)
    eng, perms = executor
    out = json.loads(eng._execute_run_tests({"pytest_args": ["-q"]}, perms))
    assert out.get("status") == "started", out
    assert out.get("job_id") == _Job.job_id
    assert "target" in out.get("note", "")
    from delfin.agent import job_monitor as jm
    watched = jm.load_watched(jm._agent_watch_path(tmp_path))["jobs"]
    assert _Job.job_id in watched


@pytest.mark.parametrize("args", [
    {"target": "tests/test_m1.py"},
    {"pytest_args": ["-q", "tests/", "-k", "co2"]},
])
def test_run_tests_over_a_chosen_part_runs_in_the_foreground(
        executor, started, tmp_path, monkeypatch, args):
    from delfin.agent import test_runner
    ran = MagicMock(return_value={"status": "ok"})
    monkeypatch.setattr(test_runner, "run_tests", ran)
    _suite(tmp_path, 45)
    eng, perms = executor
    out = json.loads(eng._execute_run_tests(dict(args), perms))
    assert out == {"status": "ok"}
    assert not started


def test_a_sub_agent_still_runs_it_itself(executor, started, tmp_path,
                                          monkeypatch):
    from delfin.agent import test_runner
    ran = MagicMock(return_value={"status": "ok"})
    monkeypatch.setattr(test_runner, "run_tests", ran)
    _suite(tmp_path, 45)
    eng, perms = executor
    perms.subagent_depth = 1
    out = json.loads(eng._execute_run_tests({"pytest_args": ["-q"]}, perms))
    assert out == {"status": "ok"}


# -- bash ---------------------------------------------------------------------

def test_a_whole_suite_bash_run_is_started_not_waited_for(
        executor, started, tmp_path, monkeypatch):
    monkeypatch.setattr(contained_run, "run", MagicMock(
        side_effect=AssertionError("the whole suite ran in the foreground")))
    _suite(tmp_path, 45)
    eng, perms = executor
    cmd = ("timeout 580 python -m pytest tests/ -q --tb=line 2>&1 | tail -5")
    out = json.loads(eng._execute_bash({"command": cmd}, perms))
    assert out.get("status") == "started", out
    assert out.get("moved_to_background") is True
    assert started and started[0]["command"] == cmd
    assert "run_tests" in out.get("note", "")


def test_a_small_suite_bash_run_is_untouched(executor, started, tmp_path,
                                             monkeypatch):
    fake = subprocess.CompletedProcess(["x"], 0, "3 passed\n", "")
    monkeypatch.setattr(contained_run, "run", MagicMock(return_value=fake))
    _suite(tmp_path, 3)
    eng, perms = executor
    out = json.loads(eng._execute_bash(
        {"command": "python -m pytest tests/ -q"}, perms))
    assert out.get("exit_code") == 0
    assert not started


def test_a_command_that_timed_out_is_not_rerun_under_another_prefix(
        executor, monkeypatch):
    runs = []

    def _run(*a, **kw):
        runs.append(a)
        raise subprocess.TimeoutExpired(cmd="x", timeout=kw.get("timeout", 1))

    monkeypatch.setattr(contained_run, "run", _run)
    eng, perms = executor
    first = json.loads(eng._execute_bash(
        {"command": "timeout 580 make check 2>&1 | tail"}, perms))
    assert "timed out" in first["error"]
    t0 = time.monotonic()
    second = json.loads(eng._execute_bash(
        {"command": "timeout 560 make check"}, perms))
    assert time.monotonic() - t0 < 2
    assert len(runs) == 1, "the respelled command ran again"
    assert "bash_background" in second["error"]
    # A new turn is a new chance: the user may have changed something.
    perms.turn_token = A._next_turn_token()
    json.loads(eng._execute_bash({"command": "make check"}, perms))
    assert len(runs) == 2


# -- a user message during a long run ----------------------------------------

def test_a_user_message_ends_a_long_command(executor, monkeypatch):
    monkeypatch.setattr(A, "USER_INTERRUPT_AFTER_S", 0.5)
    eng, perms = executor
    perms.user_message_waiting = lambda: True
    t0 = time.monotonic()
    out = json.loads(eng._execute_bash({"command": "sleep 30"}, perms))
    assert time.monotonic() - t0 < 10
    assert out.get("exit_code") == contained_run.STOPPED_RETURNCODE
    assert out.get("interrupted_by_user_message") is True
    assert "message" in out.get("note", "")


def test_a_short_command_finishes_although_a_message_waits(executor,
                                                           monkeypatch):
    monkeypatch.setattr(A, "USER_INTERRUPT_AFTER_S", 30)
    eng, perms = executor
    perms.user_message_waiting = lambda: True
    out = json.loads(eng._execute_bash({"command": "echo done"}, perms))
    assert out.get("exit_code") == 0
    assert "interrupted_by_user_message" not in out


def test_the_client_knows_a_user_message_is_waiting():
    import threading
    c = A.OpenAIClient.__new__(A.OpenAIClient)
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    assert c._user_steer_waiting() is False
    c.push_steer(A._HarnessText("[watch] a job finished"))
    assert c._user_steer_waiting() is False
    c.push_steer("push it now")
    assert c._user_steer_waiting() is True


def test_the_engine_installs_the_probe(tmp_path):
    from delfin.agent import engine as E
    from delfin.agent.api_client import StreamEvent
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")
    seen = {}

    def _one_round(**kw):
        probe = getattr(eng.client.kit_permissions, "user_message_waiting",
                        None)
        seen["probe"] = probe() if callable(probe) else None
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=1,
                          output_tokens=1, cost_usd=0.0)

    eng.client._user_steer_waiting = lambda: True
    eng.client.stream_message = MagicMock(side_effect=_one_round)
    eng.stream_response("hello")
    assert seen.get("probe") is True
