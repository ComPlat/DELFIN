"""Guards on the `sbatch` gate's segment handling: the green half of round 7.

The red half -- a shell separator INSIDE a quoted argument (`--wrap='…; …'`
or a script path containing `;`) still hides the submit, because
`sbatch_gate` splits on a plain regex at hpc_io_check.py:632 -- lives in
the reviewer's uncommitted probe. A test that is red against the code
cannot be committed.

What IS committed here is the pair of behaviours a fix for that must not
break: the control shapes without a separator must keep refusing, and a
wrapped text that walks nothing must keep running. Both pass on
8c97b4085.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
WRAP_WALK = "du -sh /pfs/data6/workspace"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


@pytest.fixture
def env():
    return {"HOME": "/home/u"}


def _verdict(cmd: str, env):
    verdict, _text, may_run = hpc_io_check.sbatch_gate(cmd, env=env)
    return verdict, may_run


def test_the_wrapped_text_alone_is_refused(lustre, env):
    """Control: without a separator the wrapped walk is already caught."""
    verdict, may_run = _verdict(f"sbatch --wrap='{WRAP_WALK}'", env)
    assert verdict == "refuse", "the control shape must refuse"
    assert may_run is False


def test_the_script_alone_is_refused(lustre, tmp_path, env):
    """Control: the script form must keep refusing a walking script."""
    script = tmp_path / "run.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, may_run = _verdict(f"sbatch {script}", env)
    assert verdict == "refuse"
    assert may_run is False


def test_a_clean_wrapped_text_with_a_separator_still_runs(lustre, env):
    """The fix must not start refusing wrapped text that walks nothing."""
    verdict, may_run = _verdict("sbatch --wrap='cd /tmp; echo hi'", env)
    assert verdict == "allow" and may_run is True, (
        "the fix must not start refusing wrapped text that walks nothing")
