"""Controls for session_load: the load an open session puts on this host.

Every case here was written against the previous commit, where
``delfin/agent/session_load.py`` did not exist: the whole file is red
there (ImportError), which is the control line of the commit that adds
the module.
"""

from __future__ import annotations

import os
import time

import pytest

from delfin.agent import session_load


def _write_proc(root, pid, ppid, uid, rss_kb, utime_ticks=0, stime_ticks=0,
                state="S", comm="python3"):
    """One process directory under an artificial /proc root."""
    d = root / str(pid)
    d.mkdir(parents=True, exist_ok=True)
    # comm may contain spaces and parens; the parser must find the fields
    # after the LAST ')', exactly like /proc does. Real layout: state,
    # ppid, then nine fields before utime/stime (fields 14/15).
    stat = (f"{pid} ({comm}) {state} {ppid} 0 0 0 0 0 0 0 0 0 "
            f"{utime_ticks} {stime_ticks} 0 0 0 0 0 0 0 0 "
            f"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")
    (d / "stat").write_text(stat, encoding="utf-8")
    (d / "status").write_text(
        f"Name:\t{comm}\nUid:\t{uid}\nVmRSS:\t{rss_kb} kB\n",
        encoding="utf-8")
    return d


@pytest.fixture
def proc_root(tmp_path):
    root = tmp_path / "proc"
    root.mkdir()
    return root


class TestProcessTree:
    def test_descendants_only_root_included(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 100)
        _write_proc(proc_root, 11, 10, os.getuid(), 100)
        _write_proc(proc_root, 12, 11, os.getuid(), 100)
        _write_proc(proc_root, 99, 1, os.getuid(), 100)   # unrelated
        assert session_load.process_tree(10, proc_root=proc_root) == [10, 11, 12]

    def test_other_uid_excluded(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 100)
        _write_proc(proc_root, 11, 10, os.getuid() + 1, 100)  # not ours
        assert session_load.process_tree(10, proc_root=proc_root) == [10]

    def test_missing_pid_is_empty(self, proc_root):
        assert session_load.process_tree(424242, proc_root=proc_root) == []

    def test_real_tree_of_this_test_process(self):
        # The pytest process itself is a tree of at least one process.
        pids = session_load.process_tree(os.getpid())
        assert os.getpid() in pids
        assert len(pids) >= 1

    def test_bad_stat_line_is_skipped_not_fatal(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 100)
        d = proc_root / "11"
        d.mkdir()
        (d / "stat").write_text("garbage without parens", encoding="utf-8")
        assert session_load.process_tree(10, proc_root=proc_root) == [10]


class TestLoadOf:
    def test_counts_procs_rss_and_d_state(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 100)
        _write_proc(proc_root, 11, 10, os.getuid(), 300, state="D")
        _write_proc(proc_root, 12, 10, os.getuid(), 600, state="R")
        load = session_load.load_of(10, proc_root=proc_root)
        assert load["procs"] == 3
        assert load["rss_bytes"] == (100 + 300 + 600) * 1024
        assert load["d_state"] == 1

    def test_cpu_cores_from_difference(self, proc_root):
        # utime+stime = 3.0 s of CPU; 1.0 s of wall clock between two
        # measurements means one core busy on average.
        _write_proc(proc_root, 10, 1, os.getuid(), 1,
                    utime_ticks=200, stime_ticks=100)
        load = session_load.load_of(
            10, previous={"cpu_seconds": 2.0,
                          "monotonic": time.monotonic() - 1.0},
            proc_root=proc_root)
        assert load["cpu_cores"] == pytest.approx(1.0, abs=0.05)

    def test_no_previous_means_no_cores(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 1)
        assert session_load.load_of(10, proc_root=proc_root)["cpu_cores"] is None

    def test_other_uid_not_counted(self, proc_root):
        _write_proc(proc_root, 10, 1, os.getuid(), 100)
        _write_proc(proc_root, 11, 10, os.getuid() + 1, 99999)
        load = session_load.load_of(10, proc_root=proc_root)
        assert load["procs"] == 1
        assert load["rss_bytes"] == 100 * 1024


class TestAlarms:
    def test_defaults_quiet_under_limits(self):
        load = {"procs": 10, "rss_bytes": 1 << 30, "cpu_cores": 1.0,
                "d_state": 0}
        assert session_load.alarms(load) == []

    def test_default_thresholds(self):
        load = {"procs": 61, "rss_bytes": (12 << 30) + 1,
                "cpu_cores": 4.5, "d_state": 3}
        found = session_load.alarms(load)
        assert any("61" in a for a in found)
        assert any("12" in a for a in found)
        assert any("4" in a for a in found)
        assert any("D" in a for a in found)

    def test_custom_limits(self):
        load = {"procs": 5, "rss_bytes": 1 << 30, "cpu_cores": 0.5,
                "d_state": 0}
        limits = {"max_procs": 4, "max_rss_gb": 0.5, "max_cores": 0.1,
                  "max_d_state": 0}
        assert len(session_load.alarms(load, limits)) == 3
