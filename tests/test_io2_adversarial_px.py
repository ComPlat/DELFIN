"""IO2 phase 3 (px) — adversarial: the agent's bash route must not be
bypassable by ordinary shell spelling.

The hook ``hpc_io_check.sbatch_gate`` decides whether an ``sbatch`` a
model typed into its bash tool may run. ``tests/test_io2_px_sbatch_gate.py``
pins the shapes the builder chose; this file pins the shapes a model
actually types, which are the same command with something appended.

Every case here is measured on the merged builder tip. A case is red
when the hook allows a script it refuses when the command stands alone,
so the bypass is proven by comparison, not by assertion:

    sbatch bad.sh                     -> refuse   (control, green)
    sbatch bad.sh && echo submitted   -> allow    (B1, red)
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
# A script the check refuses when it stands alone: a recursive stat of a
# whole workspace tree on Lustre, which is what the two operator mails to
# the cluster team described.
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


@pytest.fixture
def env():
    return {"HOME": "/home/u"}


@pytest.fixture
def walk_script(tmp_path):
    path = tmp_path / "run.sh"
    path.write_text(WALK, encoding="utf-8")
    return path


def _verdict(cmd: str, env):
    verdict, _text, may_run = hpc_io_check.sbatch_gate(cmd, env=env)
    return verdict, may_run


# --- Control: standing alone, the same script is refused ------------------

def test_the_script_alone_is_refused(lustre, walk_script, env):
    verdict, may_run = _verdict(f"sbatch {walk_script}", env)
    assert verdict == "refuse", "the control shape must refuse"
    assert may_run is False


def test_a_clean_script_alone_is_allowed(lustre, tmp_path, env):
    clean = tmp_path / "clean.sh"
    clean.write_text("#!/bin/bash\necho hello\n", encoding="utf-8")
    verdict, may_run = _verdict(f"sbatch {clean}", env)
    assert verdict == "allow" and may_run is True


# --- B1: appending a second command bypasses the check --------------------
# ``segments[-1]`` is the only segment inspected (hpc_io_check.py:600), so
# a separator AFTER the submit hands the gate a non-sbatch segment. These
# are the spellings a model writes to see whether its job was queued.

@pytest.mark.parametrize("suffix", [
    " && echo submitted",
    " && squeue -u $USER",
    "; echo done",
    " | tee submit.log",
    " > submit.log",
    " 2>&1",
])
def test_a_command_after_the_submit_does_not_hide_it(
        lustre, walk_script, env, suffix):
    cmd = f"sbatch {walk_script}{suffix}"
    verdict, may_run = _verdict(cmd, env)
    assert verdict == "refuse", (
        f"the check saw only the last segment, so `{suffix.strip()}` after "
        f"the submit let a walking script through: {cmd}")
    assert may_run is False, "a refused submit must not run"


# --- B2: a quoted script path with a space --------------------------------
# ``target.split()`` breaks the quoted path into two tokens, the last of
# which names nothing, and a FileNotFoundError is read as "let sbatch
# report it" -- so a walking script in a directory with a space is never
# checked.

def test_a_quoted_script_path_with_a_space_is_checked(lustre, tmp_path, env):
    spaced = tmp_path / "my dir"
    spaced.mkdir()
    script = spaced / "run.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, may_run = _verdict(f'sbatch "{script}"', env)
    assert verdict == "refuse", (
        "a quoted path is split on whitespace, so the walking script was "
        "not found and the gate fell through to allow")
    assert may_run is False
