"""What is pending is counted at the ``bug ls`` door and by nothing else.

Input: an archive directory holding freshly written reports. Output: the
lines ``cli.cmd_bug`` prints for the ``ls`` action -- one row per un-triaged
report and the count -- driven through the real entry point
(``delfin-agent bug ls``), not by calling a helper.

Semantics: the pending count arrives at that door; the daemon's liveness
already has its own answer in the ops report (``stop_all.status`` and its
``daemon_pid_files``, printed by ``cli.py``). ``bug_watcher.watcher_status``
-- a second daemon-status reader that no caller reached -- is deleted.
Rejected alternative: wiring ``watcher_status`` into the dashboard: it would
be a third view beside the ops report and ``bug ls``, a second answer that
can disagree with the first.

The invariant test fails at the baseline (the second reader is still
defined); the counting tests pass at both arms -- they are the never-worse
guard that the one door still answers.
"""
from __future__ import annotations

import argparse
import io
from contextlib import redirect_stdout

from delfin.agent import bug_report as br
from delfin.agent import bug_watcher
from delfin.agent import cli


def _ls(archive):
    out = io.StringIO()
    with redirect_stdout(out):
        code = cli.cmd_bug(argparse.Namespace(bug_action="ls",
                                              archive=str(archive)))
    return code, out.getvalue()


def test_the_ls_door_counts_what_is_pending(tmp_path):
    archive = tmp_path / "archive"
    report = br.write_bug_report(chat_messages=[], settings={},
                                 archive_dir=archive, session_id="pending1")
    code, text = _ls(archive)
    assert code == 0
    assert report.name in text
    assert "1 pending in" in text
    empty = tmp_path / "empty"
    empty.mkdir()
    code, text = _ls(empty)
    assert code == 0
    assert "No un-triaged reports in" in text


def test_no_second_daemon_status_reader_is_defined():
    assert not hasattr(bug_watcher, "watcher_status")