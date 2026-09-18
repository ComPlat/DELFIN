"""A suite that takes half an hour is started, not sat on.

Measured across the four sessions of 2026-09-18: eighteen tool calls ran
longer than five minutes, 154 minutes between them, and every session hit
the same wall in the same way.

    bash                 capped at 600 s, killed, ten minutes gone — all
                         four sessions, once or twice each
    run_tests            accepts timeout_s up to 1800 and runs in the
                         FOREGROUND: one call held a turn for 1803 s
    bash_background      the one that works, and the one they reached for
                         only after losing the ten minutes

`bash` learned this in the morning ("a long run is started, not waited
out"); `run_tests` did not. It does now — and it refuses rather than
clamps, because clamping spends the ceiling before saying the same
thing, which is exactly what cost each session its ten minutes.

The way out is real as of today: a finished background job wakes the
session by itself, in the dashboard and at the terminal prompt.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture()
def run_tests(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    def _call(**args):
        return json.loads(eng._execute_run_tests(dict(args), perms)), perms
    return _call


def test_a_half_hour_suite_is_refused_at_once(run_tests):
    out, perms = run_tests(pytest_args=["-q", "tests/"], timeout_s=1800)
    err = out.get("error", "")
    assert "started, not waited out" in err
    assert "bash_background" in err
    assert "wakes this session" in err, "the way out has to be named"
    assert out.get("requested_timeout_s") == 1800
    assert out.get("max_timeout_s") == perms.bash_max_timeout_s


def test_it_refuses_before_spending_the_time(run_tests):
    """The whole point: the answer costs nothing. A clamp would spend the
    ceiling first and then say the same sentence."""
    import time
    t0 = time.monotonic()
    run_tests(pytest_args=["-q"], timeout_s=1800)
    assert time.monotonic() - t0 < 5.0


def test_the_ceiling_is_the_shell_s_own(run_tests):
    """One number for both, so the agent does not have to learn two."""
    out, perms = run_tests(pytest_args=["-q"], timeout_s=100000)
    assert out.get("max_timeout_s") == perms.bash_max_timeout_s


@pytest.mark.parametrize("timeout", [5, 60, 300, 600])
def test_a_run_inside_the_ceiling_is_untouched(run_tests, timeout):
    """Below it nothing changes — the tool still runs the tests."""
    out, _ = run_tests(target="nothing_here.py", pytest_args=["-q"],
                       timeout_s=timeout)
    assert "started, not waited out" not in str(out.get("error", ""))


def test_the_default_is_untouched(run_tests):
    out, _ = run_tests(target="nothing_here.py", pytest_args=["-q"])
    assert "started, not waited out" not in str(out.get("error", ""))
