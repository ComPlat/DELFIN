"""Control tests for the per-call process and memory budget.

Red on the commit before ``delfin/agent/process_budget.py`` existed:
every case imports the module, so the whole file fails at collection.
"""
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent.process_budget import (  # noqa: E402
    BudgetLimits,
    BudgetGuard,
    count_descendants,
    sample_descendants,
    verdict,
)


def _spawn(cmd):
    return subprocess.Popen(
        cmd, shell=True, start_new_session=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def _wait_for(predicate, timeout=10.0, interval=0.05):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(interval)
    return None


class TestCounting:
    def test_counts_a_tree_of_sleeps(self):
        # A shell that starts 5 sleeps; the root itself does not count,
        # so the census must see exactly 5 descendants.
        root = _spawn("sleep 30 & sleep 30 & sleep 30 & sleep 30 & sleep 30 & wait")
        try:
            n = _wait_for(lambda: (count_descendants(root.pid) == 5
                                   and count_descendants(root.pid)) or None)
            assert n == 5, f"expected 5 descendants, saw {count_descendants(root.pid)}"
        finally:
            os.killpg(root.pid, signal.SIGKILL)
            root.wait()

    def test_root_pid_missing_counts_zero(self):
        # A pid that has certainly exited: the count must be 0, not raise.
        gone = _spawn("true")
        gone.wait()
        assert count_descendants(gone.pid) == 0

    def test_zombies_do_not_count_as_running(self):
        # Children that exit while the parent is SIGSTOPped cannot be
        # reaped and turn zombie; they must not count as running. Only
        # the one live sleep counts: 1 descendant, not 3.
        root = _spawn("true & true & sleep 30 & kill -STOP $self"
                      .replace("$self", "$$"))
        try:
            n = _wait_for(lambda: count_descendants(root.pid) == 1)
            assert n == 1, f"zombies leaked into the count: {n}"
        finally:
            os.killpg(root.pid, signal.SIGKILL)
            root.wait()


class TestVerdict:
    def test_within_limits_is_none(self):
        sample = {"processes": 4, "rss_mb": 100.0, "cpu_s": 10.0}
        limits = BudgetLimits(max_processes=8, max_rss_mb=4096.0,
                              max_cpu_seconds=600.0)
        assert verdict(sample, limits) is None

    def test_process_breach_names_the_limit(self):
        sample = {"processes": 9, "rss_mb": 100.0, "cpu_s": 10.0}
        limits = BudgetLimits(max_processes=4, max_rss_mb=4096.0,
                              max_cpu_seconds=600.0)
        breach = verdict(sample, limits)
        assert breach is not None
        assert "process" in breach.message
        assert "9" in breach.message and "4" in breach.message

    def test_rss_and_cpu_breaches(self):
        limits = BudgetLimits(max_processes=4, max_rss_mb=50.0,
                              max_cpu_seconds=5.0)
        assert verdict({"processes": 1, "rss_mb": 51.0, "cpu_s": 1.0},
                       limits) is not None
        assert verdict({"processes": 1, "rss_mb": 1.0, "cpu_s": 6.0},
                       limits) is not None


class TestGuard:
    def test_breach_kills_the_whole_group_but_not_outsiders(self):
        # Inside the group: a shell that keeps 8 sleeps alive (limit 4).
        inside = _spawn(
            "while true; do sleep 30 & sleep 30 & sleep 30 & sleep 30 & "
            "sleep 30 & sleep 30 & sleep 30 & sleep 30 & wait; done"
        )
        # Outside the group: an unrelated sleep in its own session.
        outside = _spawn("sleep 30")
        guard = BudgetGuard(
            inside.pid, limits=BudgetLimits(
                max_processes=4, max_rss_mb=4096.0, max_cpu_seconds=600.0),
            poll_s=0.05, term_grace_s=2.0)
        try:
            guard.start()
            finished = _wait_for(
                lambda: inside.poll() is not None, timeout=15.0)
            assert finished is not None, "the group survived the breach"
            # The outside process must still be alive.
            assert outside.poll() is None, "the guard killed outside its group"
        finally:
            guard.stop()
            for p in (inside, outside):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                p.wait()

    def test_no_breach_leaves_the_group_running(self):
        inside = _spawn("sleep 30")
        guard = BudgetGuard(
            inside.pid, limits=BudgetLimits(
                max_processes=8, max_rss_mb=4096.0, max_cpu_seconds=600.0),
            poll_s=0.05, term_grace_s=2.0)
        try:
            guard.start()
            time.sleep(0.5)
            assert inside.poll() is None, "guard killed a compliant call"
        finally:
            guard.stop()
            os.killpg(inside.pid, signal.SIGKILL)
            inside.wait()

    def test_sample_sees_the_tree(self):
        root = _spawn("sleep 30 & sleep 30 & wait")
        try:
            sample = _wait_for(
                lambda: (sample_descendants(root.pid)["processes"] == 2
                         and sample_descendants(root.pid)) or None)
            assert sample is not None
            assert set(sample) >= {"processes", "rss_mb", "cpu_s"}
        finally:
            os.killpg(root.pid, signal.SIGKILL)
            root.wait()
