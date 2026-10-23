"""Controls for `delfin-agent sessions` showing each session's load.

Written against the commit before the CLI wiring: there, `sessions`
has no --load/--watch flags, and the open-sessions block carries no
procs/rss/cpu columns — every case below is red on that commit.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import cli, session_load


@pytest.fixture
def fake_proc(tmp_path, monkeypatch):
    """An artificial /proc with one session tree of 3 processes."""
    root = tmp_path / "proc"
    root.mkdir()

    def write(pid, ppid, rss_kb, state="S", utime=0, stime=0, uid=None):
        uid = os.getuid() if uid is None else uid
        d = root / str(pid)
        d.mkdir(parents=True, exist_ok=True)
        stat = (f"{pid} (python3) {state} {ppid} 0 0 0 0 0 0 0 0 0 "
                f"{utime} {stime} 0 0 0 0 0 0 0 0 "
                f"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")
        (d / "stat").write_text(stat, encoding="utf-8")
        (d / "status").write_text(
            f"Name:\tpython3\nUid:\t{uid}\nVmRSS:\t{rss_kb} kB\n",
            encoding="utf-8")

    write(10, 1, 100)
    write(11, 10, 300, state="D")
    write(12, 10, 600)
    write(99, 1, 99999)  # another user's, same name
    monkeypatch.setattr(session_load, "_MY_UID", os.getuid())
    return root


def _rows_with_load(proc_root, records, previous=None, limits=None):
    """What cmd_sessions will call: one load row per open session."""
    out = []
    for rec in records:
        pid = int(rec.get("pid") or 0)
        load = session_load.load_of(pid, previous, proc_root=proc_root)
        out.append((rec, load, session_load.alarms(load, limits)))
    return out


class TestLoadRows:
    def test_one_row_per_session_with_pid(self, fake_proc):
        records = [{"key": "s1", "pid": 10, "host": "this-host"}]
        rows = _rows_with_load(fake_proc, records)
        assert len(rows) == 1
        rec, load, found = rows[0]
        assert rec["key"] == "s1"
        assert load["procs"] == 3
        assert load["rss_bytes"] == 1000 * 1024
        assert load["d_state"] == 1
        assert found == []  # 3 procs / ~1 MB is quiet

    def test_other_host_pid_zero_skipped(self, fake_proc):
        records = [{"key": "s2", "pid": 0, "host": "other-host"},
                   {"key": "s3", "pid": 10, "host": "this-host"}]
        rows = [r for r in _rows_with_load(fake_proc, records)
                if r[1]["procs"] > 0]
        assert [r[0]["key"] for r in rows] == ["s3"]

    def test_pid_of_dead_session_is_quiet(self, fake_proc):
        records = [{"key": "gone", "pid": 424242, "host": "this-host"}]
        rows = _rows_with_load(fake_proc, records)
        assert rows[0][1]["procs"] == 0


class TestRendering:
    def test_format_columns(self, fake_proc):
        from delfin.agent import cli as _h
        line = _h.format_load_row("s1", {"procs": 3, "rss_bytes": 3 << 30,
                                         "cpu_cores": 0.5, "d_state": 0}, [])
        assert "s1" in line
        assert "3" in line
        assert "3.0G" in line
        assert "0.5" in line

    def test_alarms_highlighted(self, fake_proc):
        from delfin.agent import cli as _h
        line = _h.format_load_row(
            "s1", {"procs": 80, "rss_bytes": 1 << 30, "cpu_cores": None,
                   "d_state": 0},
            ["80 processes (limit 60)"])
        assert "!" in line
        assert "80 processes" in line

    def test_cpu_none_renders_as_dash(self, fake_proc):
        from delfin.agent import cli as _h
        line = _h.format_load_row("s1", {"procs": 1, "rss_bytes": 0,
                                         "cpu_cores": None, "d_state": 0}, [])
        assert "-" in line


class TestParser:
    def test_sessions_has_load_and_watch_flags(self):
        parser = cli.build_parser()
        args = parser.parse_args(["sessions", "--load", "--watch", "5"])
        assert getattr(args, "load", False) is True
        assert int(getattr(args, "watch", 0)) == 5

    def test_sessions_plain_has_no_load(self):
        parser = cli.build_parser()
        args = parser.parse_args(["sessions"])
        assert getattr(args, "load", False) is False
        assert getattr(args, "watch", None) is None

    def test_cmd_sessions_load_prints_one_line_per_session(self, fake_proc,
                                                           capsys, monkeypatch):
        # The real session_presence dir belongs to the running sessions;
        # the record here is the same shape announce writes, pointed at
        # THIS test process so load_of reads a real, live tree.
        import socket

        records = [{"key": "selftest", "pid": os.getpid(),
                    "host": socket.gethostname()}]
        from delfin.agent import session_presence as _pres
        monkeypatch.setattr(_pres, "open_sessions", lambda **kw: records)
        args = cli.build_parser().parse_args(["sessions", "--load"])
        assert cli.cmd_sessions(args) == 0
        out = capsys.readouterr().out
        assert "selftest" in out
        assert "procs" in out
        assert "rss" in out

    def test_cmd_sessions_load_without_load_still_lists_history(
            self, capsys):
        # Plain `sessions` keeps its old behaviour: history table, no
        # load block. A regression guard, green before and after.
        args = cli.build_parser().parse_args(["sessions", "--limit", "1"])
        assert cli.cmd_sessions(args) == 0

    def test_plain_sessions_shows_load_columns_for_open_sessions(
            self, fake_proc, capsys, monkeypatch):
        # The default listing names what each open session on this host
        # costs: procs, rss, cpu next to the open block.
        import socket

        records = [{"key": "selftest", "pid": os.getpid(),
                    "host": socket.gethostname()}]
        from delfin.agent import session_presence as _pres
        monkeypatch.setattr(_pres, "open_sessions", lambda **kw: records)
        args = cli.build_parser().parse_args(["sessions"])
        assert cli.cmd_sessions(args) == 0
        out = capsys.readouterr().out
        assert "selftest" in out
        assert "procs" in out and "rss" in out and "cpu" in out
