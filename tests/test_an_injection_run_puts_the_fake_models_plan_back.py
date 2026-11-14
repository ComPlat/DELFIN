"""An injection-suite run puts the scripted model's plan back.

Input: ``tests.injection_suite.runner.run_case`` with the scripted provider,
which swaps ``fake_model.plan`` for the case's obedient planner.

Before: the swap stayed after the run. Every later test in the same process
that served fake_model got the planner that answers "ok" -- in CI,
test_real_run_webfetch_plan_mode failed ("no tool result reached the
model") whenever a shard ran it after
test_a_tainted_turn_needs_the_users_words_to_send_out (2026-11-14, the shard
split of two unrelated PRs put them together).
Now the module's own plan is restored after the run, also when the child
times out. The child process is a stand-in here; the order the CI hit is
pinned too, with the real CLI.
"""
from __future__ import annotations

import subprocess

import pytest

from tests.agent_terminal_probes import fake_model
from tests.injection_suite import cases as C
from tests.injection_suite import runner


def _case():
    return next(c for c in C.build() if c.id == "web_fetch_exfil")


@pytest.mark.parametrize("child", ["answers", "times_out"])
def test_the_plan_is_the_modules_own_after_a_run(tmp_path, monkeypatch, child):
    own = fake_model.plan

    def fake_run(argv, **kw):
        if child == "times_out":
            raise subprocess.TimeoutExpired(argv, 1)
        return subprocess.CompletedProcess(argv, 0, "{}", "")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    try:
        runner.run_case(_case(), tree=tmp_path, base=tmp_path, timeout=5)
    except Exception:
        pass  # what collect() makes of a stand-in output is not the point
    assert fake_model.plan is own


def test_a_scripted_run_first_does_not_change_the_next_real_session(tmp_path):
    """The CI order, through the real entry points: one injection case,
    then the plan-mode web_fetch session that used to fail after it."""
    from tests import test_real_run_webfetch_plan_mode as later
    tree = later.REPO
    r = runner.run_case(_case(), tree=tree, base=tmp_path, timeout=240)
    assert not r.attack_success
    mp = pytest.MonkeyPatch()
    try:
        (tmp_path / "next").mkdir()
        later.test_real_session_refuses_webfetch_in_plan_mode(tmp_path / "next", mp)
    finally:
        mp.undo()
