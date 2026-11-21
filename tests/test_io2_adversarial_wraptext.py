"""Guards on the `sbatch --wrap=` route: the green half of rounds 8/9.

The red half of this investigation lives in the reviewer's uncommitted
probe, because a test that is red against the code cannot be committed:

- a wrap TEXT shaped like a script name (`hpc_io_check.py:540-546` reads
  it as a PATH and FileNotFoundError escapes at :651);
- a walk that sits behind a separator inside the wrap text (the split at
  `hpc_io_check.py:632` is not quote-aware).

Both findings are with the builder and the QS with file:line and the fix
shape.

What IS committed here is the half a fix must not destroy: a wrap text
that carries a Lustre walk must keep being refused, and a harmless wrap
must keep running. Both pass on 8c97b4085, so a fix that breaks one is
caught rather than celebrated.

A third case was dropped after the QS ruled requirement (a) for
`sbatch --wrap='bash run.sh'`. It had asserted that a wrap text NAMING an
existing walking script must be refused, which is requirement (b) --
resolving the named file -- and the brief never asks for it. The ruling
is that the wrap TEXT is the unit: a wrap whose text carries no walk is
allowed even when a file of that name walks. That assertion is red
today, so it sits in the probe with the rest of the acceptance criteria.
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


def test_the_wrapped_walk_is_still_caught(mounts, tmp_path, monkeypatch):
    """The crash must not be 'fixed' by refusing or ignoring every wrap."""
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
