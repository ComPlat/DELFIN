"""IO2: a script piped INTO sbatch must be judged like one passed by file.

`echo '<walking script>' | sbatch` is a submit of that text: sbatch reads
its batch script from stdin, and a real bash runs what a pipe delivers
(measured: `echo 'echo PIPED-PROGRAM' | bash` prints PIPED-PROGRAM). The
gate split the line into `echo '<text>'` and `sbatch` segments, saw no
script argument on the second, and allowed it.

The producers the gate can actually read are resolved:
  * ``echo``/``printf`` with a literal argument -> judge that text;
  * ``cat <file>``                              -> judge that file;
  * anything else                              -> refuse, because the text
    cannot be known (fail closed, the same principle as a path that is not
    a file).

A ``cat script.sh | sbatch`` that is clean must still run, and a name
given as an ARGUMENT is unaffected: ``echo x | sbatch batch.sh`` runs
batch.sh, so the pipe is that command's stdin and the FILE is judged.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "du -sh /pfs/data6/workspace/"
CLEAN = "sbatch --version"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def test_an_echoed_walk_into_sbatch_is_refused(lustre):
    verdict, text, may_run = _v(f"echo '{WALK}' | sbatch")
    assert verdict == "refuse", (
        f"FAIL-OPEN: a piped walking script was not judged "
        f"(verdict={verdict!r})")
    assert "tree_walk_on_lustre" in text
    assert may_run is False


def test_a_printf_walk_into_sbatch_is_refused(lustre):
    verdict, text, may_run = _v(f"printf '%s' '{WALK}' | sbatch")
    assert verdict == "refuse", (
        f"FAIL-OPEN: a printf'd walking script was not judged "
        f"(verdict={verdict!r})")
    assert "tree_walk_on_lustre" in text


def test_a_catted_walking_file_into_sbatch_is_refused(lustre, tmp_path):
    script = tmp_path / "walk.sh"
    script.write_text("#!/bin/bash\n" + WALK + "\n", encoding="utf-8")
    verdict, text, may_run = _v(f"cat {script} | sbatch")
    assert verdict == "refuse", (
        f"FAIL-OPEN: a catted walking script was not judged "
        f"(verdict={verdict!r})")
    assert "tree_walk_on_lustre" in text


def test_a_catted_clean_file_into_sbatch_runs(lustre, tmp_path):
    script = tmp_path / "clean.sh"
    script.write_text("#!/bin/bash\n" + CLEAN + "\n", encoding="utf-8")
    verdict, text, may_run = _v(f"cat {script} | sbatch")
    assert verdict != "refuse", (
        f"OVER-REFUSAL: a clean catted script was refused: {text[:150]}")


def test_a_clean_echo_into_sbatch_runs(lustre):
    verdict, text, may_run = _v(f"echo '{CLEAN}' | sbatch")
    assert verdict != "refuse", (
        f"OVER-REFUSAL: a clean echoed script was refused: {text[:150]}")


def test_an_unreadable_producer_fails_closed(lustre):
    """`generate_it | sbatch`: the text cannot be known -> refuse."""
    verdict, text, may_run = _v("generate_it | sbatch")
    assert verdict == "refuse", (
        f"FAIL-OPEN: an unreadable producer was allowed "
        f"(verdict={verdict!r})")
    assert may_run is False


def test_a_script_argument_beats_the_pipe(lustre, tmp_path):
    """`echo x | sbatch batch.sh` runs batch.sh -- the file is judged."""
    script = tmp_path / "batch.sh"
    script.write_text("#!/bin/bash\n" + CLEAN + "\n", encoding="utf-8")
    verdict, text, may_run = _v(f"echo hello | sbatch {script}")
    assert verdict != "refuse", (
        f"OVER-REFUSAL: the named script should have been judged and is "
        f"clean: {text[:150]}")


def test_a_pipe_not_into_sbatch_is_untouched(lustre):
    """A pipe somewhere else in the line must not be refused."""
    verdict, text, may_run = _v("du -sh /pfs/data6/workspace/ | head -1")
    assert verdict == "allow", (
        f"OVER-REFUSAL: an ordinary pipe was refused: {text[:150]}")
