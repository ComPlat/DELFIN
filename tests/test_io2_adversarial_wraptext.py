"""Guards on the `sbatch --wrap=` route: the green half of round 8/9.

The red half of this investigation -- a wrap TEXT shaped like a script
name (hpc_io_check.py:540-546 reads it as a PATH and FileNotFoundError
escapes at :651) -- lives in the reviewer's uncommitted probe, because a
test that is red against the code cannot be committed. The finding is
with the builder and the QS with file:line and the fix shape.

What IS committed here is the other half: the three behaviours a fix for
that crash must NOT destroy. Each passes on 8c97b4085, so a fix that
breaks one is caught rather than celebrated.
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


@pytest.mark.parametrize("spelling", ["run.sh", "./run.sh"])
def test_a_wrap_naming_a_walking_script_is_still_refused(
        mounts, tmp_path, monkeypatch, spelling):
    """The fix must not turn every wrap text into a pure TEXT check.

    Today the shape heuristic gives this for free: the wrap text names a
    file, the file is read, its content walks Lustre, the submit is
    refused. A fix that simply passes is_text=True from
    hpc_io_check.py:651 keeps the crash away but LOSES this refusal --
    the check then sees the literal characters "run.sh" and allows.

    So this is the guard that separates "checks the text" from "checks
    the text AND, when it names a real script, that script". A correct
    fix refuses both spellings; a naive is_text fix allows both.

    Deliberately NOT parametrised over "bash run.sh": that wrap text
    MENTIONS a file rather than naming one, so resolving it would need
    shell-word parsing the brief never asked for. I raised it with the
    QS as an open design question instead of asserting a requirement
    the package does not state.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "run.sh").write_text(
        "#!/bin/bash\n" + WALK + "\n", encoding="utf-8")

    # What the two readings do with the same string, measured: the pure
    # TEXT reading finds nothing (that is the trap), the FILE reading
    # finds the walk. If check() ever starts refusing a bare "run.sh" as
    # text, this test no longer discriminates and must be revisited.
    assert hpc_io_check.check(spelling, env=ENV, mounts_path=mounts) == []

    verdict, finding, may_run = hpc_io_check.sbatch_gate(
        f"sbatch --wrap={spelling!r}", env=ENV, mounts_path=mounts)
    assert verdict == "refuse", (
        f"wrap text {spelling!r} names an existing walking script; the "
        "submit must be judged on that script's content")
    assert may_run is False
