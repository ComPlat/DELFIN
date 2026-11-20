"""IO2 phase 3 (px) — adversarial round 7: the segment split is not
quote-aware, so a separator INSIDE a quoted argument still hides a submit.

Round 6 closed the suffix bypass (`sbatch run.sh && echo ok`) and the
quoted-path hole. Both fixes were measured green. This file pins the two
shapes that share one root cause and are still open: ``sbatch_gate``
splits the command with a plain regex on ``;`` ``&&`` ``||`` ``|``
(``hpc_io_check.py:632``, ``re.split``), which cannot see that the
separator sits inside a quoted string. What follows is therefore handled
as two shell segments when the shell would treat it as one argument.

Measured on 8c97b4085 and confirmed identical on 44e0ab4ef (a scratch
branch at the pre-fix commit), so these are PRE-EXISTING holes that the
round-6 fix neither introduced nor closed -- not regressions:

    sbatch --wrap='du -sh /pfs/.../workspace'              -> refuse
    sbatch --wrap='du -sh /pfs/.../workspace; echo ok'     -> allow   (Q1)
    sbatch "a;b/run.sh"  (walking script)                  -> allow   (Q2)

Q1 is in the package's own scope: deliverable 2 says ``sbatch
--wrap=...`` checks the WRAPPED TEXT. A wrapped text containing an
ordinary ``;`` is submitted unchecked.
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


# --- Controls: without a separator, both shapes are refused ---------------

def test_the_wrapped_text_alone_is_refused(lustre, env):
    verdict, may_run = _verdict(f"sbatch --wrap='{WRAP_WALK}'", env)
    assert verdict == "refuse", "the control shape must refuse"
    assert may_run is False


def test_the_script_alone_is_refused(lustre, tmp_path, env):
    script = tmp_path / "run.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, may_run = _verdict(f"sbatch {script}", env)
    assert verdict == "refuse"
    assert may_run is False


# --- Q1: a separator INSIDE the wrapped text -----------------------------
# Each of these wraps text the check refuses when it stands alone. The
# separator is part of the shell command sbatch will run, not a second
# command in the caller's line, and splitting on it turns the wrapped
# text into a fragment with no absolute path left to judge.

@pytest.mark.parametrize("wrapped", [
    f"{WRAP_WALK}; echo ok",
    f"cd /pfs/work; {WRAP_WALK}",
    f"echo hi && {WRAP_WALK}",
    f"{WRAP_WALK} | tee log",
    f"true || {WRAP_WALK}",
])
def test_a_separator_inside_the_wrapped_text_is_checked(lustre, env, wrapped):
    cmd = f"sbatch --wrap='{wrapped}'"
    verdict, may_run = _verdict(cmd, env)
    assert verdict == "refuse", (
        "the segment split is not quote-aware, so the separator inside "
        f"the --wrap text ended the segment and the walk was never "
        f"checked: {cmd}")
    assert may_run is False


def test_a_clean_wrapped_text_with_a_separator_still_runs(lustre, env):
    verdict, may_run = _verdict("sbatch --wrap='cd /tmp; echo hi'", env)
    assert verdict == "allow" and may_run is True, (
        "the fix must not start refusing wrapped text that walks nothing")


# --- Q2: a separator inside a quoted script path -------------------------
# A job folder may legitimately be named with a semicolon. The split cuts
# the path mid-string, the last token names nothing, and the missing-file
# rule returns allow -- the same fail-open as round 6's B2, one quoting
# level deeper.

def test_a_separator_inside_a_quoted_script_path_is_checked(
        lustre, tmp_path, env):
    folder = tmp_path / "a;b"
    folder.mkdir()
    script = folder / "run.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, may_run = _verdict(f'sbatch "{script}"', env)
    assert verdict == "refuse", (
        "the semicolon inside the quoted path was treated as a shell "
        "separator, so the walking script was never checked")
    assert may_run is False
