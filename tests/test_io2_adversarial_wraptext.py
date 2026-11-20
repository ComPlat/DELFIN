"""Adversarial round 8: --wrap TEXT is fed through check_submit's PATH path.

The phase-3 hook checks `sbatch --wrap=<text>` as text -- that is the
package's own deliverable 2. But check_submit decides path-vs-text by
SHAPE (hpc_io_check.py:540-546): a newline-free str ending in a script
suffix is treated as a PATH and read from disk. A --wrap text is
frequently exactly that shape ("bash run.sh", "python3 analysis.py"),
so the hook raises FileNotFoundError instead of returning a verdict.

Consequence on the bash route: an ordinary `sbatch --wrap='bash run.sh'`
does not get a check at all -- it crashes the gate. That is neither
"refuse" nor "run with a warning"; it is an exception in the one place
that decides whether a submit is allowed.

Attribution (measured, not assumed): the same probe on a scratch branch
at 44e0ab4ef (before the every-segment fix) and at 8c97b4085 gives the
IDENTICAL result for all 11 shapes, so this predates the segment fix; it
is a defect of the phase-3 hook as delivered, not a regression.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check


ENV = {"HOME": "/home/u"}
MOUNTS = (
    "proc /proc proc rw 0 0\n"
    "/dev/sda1 / ext4 rw 0 0\n"
    "10.0.0.1@o2ib:/lfs /pfs lustre rw 0 0\n"
)
WALK = "du -sh /pfs/data6/workspace"


@pytest.fixture()
def mounts():
    return hpc_io_check.MountsSource(text=MOUNTS)


# Each of these is a legal --wrap text whose shape ends in a script suffix
# or contains one, and none of them names a file that exists.
SUFFIX_SHAPED = [
    "run.sh",
    "./run.sh",
    "analysis.py",
    "submit_delfin.sh",
    "bash run.sh",
    "sh ./run.sh",
    "python3 analysis.py",
]


@pytest.mark.parametrize("text", SUFFIX_SHAPED)
def test_a_wrap_text_is_not_read_from_disk(mounts, tmp_path, monkeypatch, text):
    """The hook must return a verdict for wrap TEXT, never raise."""
    monkeypatch.chdir(tmp_path)
    verdict, finding, may_run = hpc_io_check.sbatch_gate(
        f"sbatch --wrap={text!r}", env=ENV, mounts_path=mounts)
    assert verdict == "allow", (
        f"a clean wrap text {text!r} should be allowed; got {verdict!r} "
        f"with {finding!r}")


@pytest.mark.parametrize("text", SUFFIX_SHAPED)
def test_a_wrap_text_ending_in_a_suffix_does_not_crash(
        mounts, tmp_path, monkeypatch, text):
    """FileNotFoundError from check_submit's shape heuristic (H1) leaks."""
    monkeypatch.chdir(tmp_path)
    try:
        hpc_io_check.sbatch_gate(
            f"sbatch --wrap={text!r}", env=ENV, mounts_path=mounts)
    except OSError as exc:
        pytest.fail(
            "sbatch_gate raised on a --wrap TEXT shaped like a script "
            f"name ({text!r}): {type(exc).__name__}: {exc}. check_submit "
            "reads it as a PATH (hpc_io_check.py:540-546); the --wrap "
            "route must pass it as TEXT.")


def test_the_wrapped_walk_is_still_caught(mounts, tmp_path, monkeypatch):
    """Control: the defect must not be 'fixed' by refusing every wrap."""
    monkeypatch.chdir(tmp_path)
    verdict, finding, may_run = hpc_io_check.sbatch_gate(
        f"sbatch --wrap={WALK!r}", env=ENV, mounts_path=mounts)
    assert verdict == "refuse"
    assert may_run is False
    assert "tree_walk_on_lustre" in finding


def test_a_clean_wrap_without_a_suffix_still_runs(mounts, tmp_path, monkeypatch):
    """Control: 'make -j8' has no script suffix and was never affected."""
    monkeypatch.chdir(tmp_path)
    verdict, finding, may_run = hpc_io_check.sbatch_gate(
        "sbatch --wrap='make -j8'", env=ENV, mounts_path=mounts)
    assert verdict == "allow"
    assert may_run is True


def test_a_wrap_naming_a_walking_script_is_still_refused(
        mounts, tmp_path, monkeypatch):
    """Control: when the named file DOES exist its content is judged.

    This is the behaviour that makes the shape heuristic tempting -- it
    must keep working once the text path stops crashing.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "run.sh").write_text(
        "#!/bin/bash\n" + WALK + "\n", encoding="utf-8")
    verdict, finding, may_run = hpc_io_check.sbatch_gate(
        "sbatch --wrap='run.sh'", env=ENV, mounts_path=mounts)
    assert verdict == "refuse"
    assert may_run is False
