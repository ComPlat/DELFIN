"""Hooks and the test runner run under a process budget (phase 4).

Both go through contained_run.run today without budget=True, so a hook
or a pytest run could fan out without bound. The public call paths —
run_hooks (fires each hook command) and test_runner.run_tests (starts
pytest) — must arm a guard with their own profile: "hook" and "tests".
"""
import time

import pytest

from delfin.agent import contained_run, hooks as hooks_mod
from delfin.agent import process_budget, test_runner


def _wait_for(condition, timeout_s=20.0, what="condition"):
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        value = condition()
        if value:
            return value
        time.sleep(0.2)
    raise AssertionError(f"timed out waiting for {what}")


def _fanout(n: int) -> str:
    return f"for i in $(seq 1 {n}); do sleep 15 & done; sleep 15"


def _shrink_profiles(monkeypatch, name, max_processes):
    monkeypatch.setattr(
        process_budget, "limits_for_profile",
        lambda n: process_budget.BudgetLimits(
            max_processes=max_processes if n == name else 4,
            max_rss_mb=1e9, max_cpu_seconds=1e9))


# --- contained_run.run: the profile plumbing itself ------------------------

def test_contained_run_accepts_a_budget_profile(monkeypatch):
    """A breach through contained_run.run(budget_profile=...) surfaces in
    stderr — the call path hooks and tests actually use."""
    _shrink_profiles(monkeypatch, "hook", 4)
    proc = contained_run.run(
        _fanout(12), shell=True, timeout=30,
        budget_profile="hook",
    )
    assert proc.returncode != 0
    assert "process budget" in (proc.stderr or "").lower()


def test_budget_false_still_means_no_guard():
    """The old callers (budget=False, no profile) get exactly what they
    had: no guard, so a fan-out command runs unmolested until timeout."""
    proc = contained_run.run(
        "for i in $(seq 1 8); do sleep 3 & done; sleep 3",
        shell=True, timeout=20)
    # No budget ended it; every sleeper finished, so the shell exits 0.
    assert proc.returncode == 0


# --- hooks: the public run_hooks path --------------------------------------

def _hook_config(tmp_path, command, event="PreToolUse"):
    hk = hooks_mod.HookCommand(
        type="command", command=command, matcher="", timeout_s=30)
    cfg = hooks_mod.HooksConfig(
        by_event={event: [hk]})
    return cfg


def test_run_hooks_passes_the_hook_profile(monkeypatch):
    seen = {}
    real_run = contained_run.run

    def spy_run(*args, **kw):
        seen.update(kw)
        return real_run(*args, **kw)

    import delfin.agent.hooks as H
    monkeypatch.setattr("delfin.agent.contained_run.run", spy_run)
    cfg = _hook_config(None, "echo hi")
    H.run_hooks("PreToolUse", cfg, tool_name="bash",
                arguments={"command": "true"})
    assert seen.get("budget_profile") == "hook"


def test_a_hook_fanout_is_killed_by_the_budget(monkeypatch, tmp_path):
    """A hook whose command fans out past the hook profile is killed and
    reported through the public run_hooks result."""
    _shrink_profiles(monkeypatch, "hook", 4)
    cfg = _hook_config(tmp_path, _fanout(12))
    t0 = time.monotonic()
    results = hooks_mod.run_hooks("PreToolUse", cfg, tool_name="bash",
                                  arguments={"command": "true"})
    dur = time.monotonic() - t0
    assert len(results) == 1
    r = results[0]
    assert r.exit_code != 0
    # It was the budget, not the 30 s timeout: a timeout had cost ~30 s.
    assert dur < 20, f"hook took {dur:.1f}s — looks like a timeout, not a breach"
    assert "process budget" in (r.stderr or "").lower()


# --- test runner: the public run_tests path --------------------------------

def test_run_tests_passes_the_tests_profile(monkeypatch, tmp_path):
    seen = {}
    real_run = contained_run.run

    def spy_run(*args, **kw):
        seen.update(kw)
        return real_run(*args, **kw)

    import delfin.agent.test_runner as TR
    monkeypatch.setattr("delfin.agent.contained_run.run", spy_run)
    result = test_runner.run_tests(str(tmp_path), pytest_args=["-q",
                           "--collect-only"])
    assert seen.get("budget_profile") == "tests"


def test_run_tests_kills_a_fanout_suite(monkeypatch, tmp_path):
    """A test file that spawns 12 sleepers past the (shrunken) tests
    profile is killed, and the result says the run was ended by the
    budget, not merely timed out."""
    _shrink_profiles(monkeypatch, "tests", 4)
    test_file = tmp_path / "test_fanout.py"
    test_file.write_text(
        "import subprocess\n"
        "def test_fans_out():\n"
        "    subprocess.run('for i in $(seq 1 12); do sleep 15 & done; "
        "sleep 15', shell=True)\n")
    t0 = time.monotonic()
    result = test_runner.run_tests(
        str(tmp_path), target=str(test_file), pytest_args=["-q", "-x"],
        timeout_s=60)
    dur = time.monotonic() - t0
    # The pytest process was killed by the guard: not a clean pass.
    assert result.get("status") != "passed"
    assert dur < 45, f"run took {dur:.1f}s — looks like a timeout, not a breach"
    tail = (result.get("raw_stderr_tail") or "") + (result.get(
        "raw_stdout_tail") or "") + str(result.get("error") or "")
    assert "process budget" in tail.lower()
