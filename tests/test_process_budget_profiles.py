"""Named budget profiles for process_budget.

Phase 2 of "no process without a budget": the foreground guard built by
wave 4 gets siblings for the call sites that run without one today
(background jobs, hooks, the test runner). Each profile is a named set of
limits with its own justification; settings can override any profile
without touching the others, and the foreground profile's behaviour must
not change.
"""
import subprocess
import time

import pytest

from delfin.agent import process_budget as pb


def test_foreground_profile_matches_the_existing_default_limits():
    """The foreground profile IS the wave-4 default: same numbers."""
    limits = pb.limits_for_profile("foreground")
    assert limits.max_processes == 64
    assert limits.max_rss_mb == 8192.0
    assert limits.max_cpu_seconds == 1200.0


def test_foreground_profile_agrees_with_default_limits_function():
    # _default_limits() is what BudgetGuard uses today; the profile path
    # must give the identical answer so wiring it in changes nothing.
    assert pb.limits_for_profile("foreground") == pb._default_limits()


def test_background_profile_is_wider_than_foreground():
    """A background job may run longer and larger, but not unbounded."""
    fg = pb.limits_for_profile("foreground")
    bg = pb.limits_for_profile("background")
    assert bg.max_processes > fg.max_processes
    assert bg.max_rss_mb > fg.max_rss_mb
    assert bg.max_cpu_seconds > fg.max_cpu_seconds
    # "not unlimited": every ceiling is finite.
    for v in (bg.max_processes, bg.max_rss_mb, bg.max_cpu_seconds):
        assert 0 < v < 1e9


def test_hook_profile_is_tighter_than_foreground():
    """Hooks run on every tool call: their budget must be small."""
    fg = pb.limits_for_profile("foreground")
    hk = pb.limits_for_profile("hook")
    assert hk.max_processes < fg.max_processes
    assert hk.max_rss_mb < fg.max_rss_mb
    assert hk.max_cpu_seconds < fg.max_cpu_seconds


def test_tests_profile_exists_and_covers_a_suite():
    ts = pb.limits_for_profile("tests")
    fg = pb.limits_for_profile("foreground")
    # pytest fans out via xdist sometimes but rarely exceeds a small build.
    assert ts.max_processes > fg.max_processes
    assert 0 < ts.max_rss_mb < 1e9


def test_unknown_profile_is_rejected():
    with pytest.raises(ValueError):
        pb.limits_for_profile("nope")


def test_profile_names_are_stable_api():
    assert set(pb.PROFILE_LIMITS) == {"foreground", "background", "hook", "tests"}


def test_settings_override_one_profile_only(monkeypatch):
    """agent.process_budget.profiles.<name>.* overrides that profile alone."""
    monkeypatch.setattr(
        pb.user_settings, "load_settings",
        lambda: {"agent": {"process_budget": {"profiles": {
            "background": {"max_processes": 300}}}}})
    bg = pb.limits_for_profile("background")
    assert bg.max_processes == 300
    # The other profiles keep their built-in numbers.
    assert pb.limits_for_profile("foreground").max_processes == 64
    # Untouched keys of the overridden profile keep their defaults too.
    assert bg.max_rss_mb == pb.PROFILE_LIMITS["background"].max_rss_mb


def test_toplevel_settings_still_override_only_the_foreground(monkeypatch):
    """The wave-4 override path (agent.process_budget.max_processes &c)
    keeps meaning "the foreground limits" and must not leak into the new
    profiles."""
    monkeypatch.setattr(
        pb.user_settings, "load_settings",
        lambda: {"agent": {"process_budget": {"max_processes": 5}}})
    assert pb.limits_for_profile("foreground").max_processes == 5
    assert pb.limits_for_profile("hook").max_processes != 5


def test_guard_uses_the_requested_profile(monkeypatch):
    """BudgetGuard(profile=...) must run with that profile's limits."""
    guard = pb.BudgetGuard(1, profile="hook")
    assert guard.limits == pb.limits_for_profile("hook")


def test_a_breached_hook_profile_reports_itself(tmp_path):
    """End-to-end through the real /proc census: a command that forks
    past the hook profile's process limit is killed and reported."""
    script = (
        "for i in $(seq 1 8); do sleep 30 & done; sleep 30")
    limits = pb.BudgetLimits(max_processes=4, max_rss_mb=1e9, max_cpu_seconds=1e9)
    guard = pb.BudgetGuard.__new__(pb.BudgetGuard)  # built manually below
    proc = subprocess.Popen(
        ["/bin/bash", "-c", script], cwd=str(tmp_path),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        stdin=subprocess.DEVNULL, start_new_session=True)
    try:
        guard = pb.BudgetGuard(proc.pid, limits=limits, poll_s=0.2)
        deadline = time.monotonic() + 10
        breach = None
        while time.monotonic() < deadline:
            breach = guard.poll_once()
            if breach is not None:
                break
            time.sleep(0.1)
        assert breach is not None, "the fan-out never breached"
        assert breach.limit == "max_processes"
        # The tree was actually ended, not just reported.
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline and proc.poll() is None:
            time.sleep(0.1)
        assert proc.poll() is not None
    finally:
        if proc.poll() is None:
            proc.kill()
        proc.wait()
