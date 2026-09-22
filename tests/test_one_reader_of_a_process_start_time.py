"""A pid is a number the system hands out again — and four places asked.

``where.py`` asked first, about the dashboard note. Then ``lifeline``,
then ``bash_jobs``, each with its own copy of the same eleven lines. Then
the job status tab, which wrote its own and wrote it wrongly:

    stat = f.read().split()
    start_ticks = int(stat[21])          # "field 22"

Field 2 of ``/proc/<pid>/stat`` is the command name, in brackets, and it
may contain spaces. Every other reader in this codebase splits after the
LAST closing bracket for exactly that reason. Measured against a binary
named ``a b``: the whole-line split answered 0 — the boot moment — where
the true value was 26055418, and the tab turned that into an elapsed time
of three days for a process one second old.

So there is one reader now, and this file pins it. The duplication was
not the bug; it was what let the bug live in only one of the copies.
"""

from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import time

import pytest

from delfin.agent import proc_identity


def _start_time_readers(read):
    """Files that take field 22 out of a stat line, however they spell it.

    Not ``split()[19]``: two of the four readers parsed in two steps and
    a pattern tied to one spelling found only the other two. The field
    index is what they all have in common, so the index is what is looked
    for — inside the passage that has just read a stat line.
    """
    found = []
    for name, text in read():
        for match in re.finditer(r"/proc/.{0,24}stat", text):
            window = text[match.end():match.end() + 400]
            if re.search(r"\[\s*(19|21)\s*\]", window):
                found.append(name)
                break
    return sorted(set(found))


def test_there_is_one_reader_of_a_process_start_time():
    root = pathlib.Path(proc_identity.__file__).parent.parent      # delfin/

    def _read():
        for path in sorted(root.rglob("*.py")):
            yield path.name, path.read_text(encoding="utf-8",
                                            errors="replace")

    readers = _start_time_readers(_read)
    assert readers == ["proc_identity.py"], (
        f"more than one place reads a process start time: {readers}")


def test_where_still_answers_through_it():
    from delfin.agent import where
    pid = os.getpid()
    assert where._process_start(pid) == proc_identity.process_start(pid)


def test_the_job_ledger_and_the_lifeline_answer_through_it():
    from delfin.agent import bash_jobs, lifeline
    assert bash_jobs._proc_start_ticks is proc_identity.start_ticks
    assert lifeline._start_ticks is proc_identity.start_ticks


def test_a_missing_process_has_no_start_time(gone_pid):
    assert proc_identity.start_ticks(gone_pid) is None
    assert proc_identity.start_ticks(0) is None
    assert proc_identity.start_ticks("not a pid") is None


def test_a_live_process_is_alive_and_a_gone_one_is_not(gone_pid):
    assert proc_identity.alive(os.getpid(),
                               proc_identity.process_start(os.getpid()))
    assert proc_identity.alive(gone_pid, "1") is False


def test_a_pid_from_another_machine_is_not_answered():
    assert proc_identity.alive(os.getpid(), "", "somewhere-else") is None


def test_a_command_name_with_a_space_does_not_shift_the_field(tmp_path):
    """The case the job status tab got wrong, driven for real."""
    if not os.path.isdir("/proc/self"):
        pytest.skip("no /proc to read")
    sleeper = shutil.which("sleep")
    if not sleeper:
        pytest.skip("no sleep to run")
    spaced = tmp_path / "a b"
    shutil.copy(sleeper, spaced)

    proc = subprocess.Popen([str(spaced), "30"])
    try:
        for _ in range(50):
            if os.path.exists(f"/proc/{proc.pid}/stat"):
                break
            time.sleep(0.02)
        line = open(f"/proc/{proc.pid}/stat", encoding="utf-8").read()
        assert " " in line[line.index("(") + 1:line.rindex(")")], (
            f"this process was supposed to have a space in its name: "
            f"{line[:60]!r}")
        truth = int(line[line.rindex(")") + 1:].split()[19])
        assert proc_identity.start_ticks(proc.pid) == truth
        assert str(truth) != line.split()[21], (
            "the whole-line split was supposed to read a different field "
            "here — without that the case proves nothing")
    finally:
        proc.kill()
        proc.wait(timeout=5)
