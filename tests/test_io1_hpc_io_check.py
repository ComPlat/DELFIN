"""IO1: delfin.hpc_io_check — Lustre metadata-load check for job scripts.

The cluster operators measured 43 million open/close calls from one Python
start and a test suite run straight from HOME on the shared Lustre file
systems, both from job scripts that looked harmless (operator notes,
14.11. and 16.11.2026). The fix rule: unpack the environment to the
node-local $TMPDIR, keep scratch on the node, copy results back once.
Nothing in DELFIN checked a script for this.

A check that fires without Lustre mounted would reject clean scripts, so
every rule test here runs twice: with a Lustre /proc/mounts fixture and
with a mounts file that has no lustre entry — the finding must exist only
in the first. check() is pure: no file walk, no subprocess.
"""
from __future__ import annotations

import textwrap

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = """\
/dev/sda1 / ext4 rw 0 0
mgmt@o2ib:fs /pfs lustre rw,lazystatfs 0 0
"""


def _check(script: str, *, lustre: bool = True):
    env = {"HOME": "/pfs/home/ka_user"} if lustre else {"HOME": "/home/ka_user"}
    return hpc_io_check.check(script, env=env,
                              mounts_path=hpc_io_check.MountsSource(
                                  LUSTRE_MOUNTS if lustre else "/dev/sda1 / ext4 rw 0 0"))


def _rules(findings):
    return {f.rule for f in findings}


# --- Lustre detection: exactly one read of /proc/mounts, type lustre -----

def test_unreadable_mounts_file_does_not_silently_disable_the_check():
    # QS finding (nacht-s8, IO1): read() returning "" on OSError makes
    # lustre_paths() -> [] and check() exits early — the check disables
    # itself exactly when /proc/mounts is unreadable. A default mounts
    # file that cannot be read must be a hard failure, not a silent no-op.
    script = "#!/bin/bash\ncd /pfs/home/ka_user/delfin\npytest -q tests/\n"
    env = {"HOME": "/pfs/home/ka_user"}
    src = hpc_io_check.MountsSource(None, path="/nonexistent")
    with pytest.raises(RuntimeError):
        hpc_io_check.check(script, env=env, mounts_path=src)


def test_lustre_detection_from_mounts_text():
    src = hpc_io_check.MountsSource(LUSTRE_MOUNTS)
    assert hpc_io_check.lustre_paths(src) == ["/pfs"]


def test_no_lustre_mount_means_no_lustre():
    src = hpc_io_check.MountsSource("/dev/sda1 / ext4 rw 0 0")
    assert hpc_io_check.lustre_paths(src) == []


# --- Rule: interpreter / pip from Lustre inside a loop -> warn ----------

PY_LOOP = textwrap.dedent("""\
    #!/bin/bash
    for f in /pfs/data/*.txt; do
        python /pfs/env/bin/analyse.py "$f"
    done
    """)


def test_interpreter_from_lustre_in_loop_warns():
    findings = _check(PY_LOOP)
    assert "lustre_interpreter_in_loop" in _rules(findings)
    f = [x for x in findings if x.rule == "lustre_interpreter_in_loop"][0]
    assert f.severity == "warn"
    assert f.line == 3
    assert "$TMPDIR" in f.instead
    assert f.why


def test_interpreter_from_lustre_in_loop_is_silent_without_lustre():
    assert _rules(_check(PY_LOOP, lustre=False)) == set()


# --- Rule: per-item interpreter start (the 43M open/close shape) ---------

PIP_LOOP = "#!/bin/bash\nfor pkg in a b c; do pip install \"$pkg\"; done\n"


def test_pip_loop_warns():
    findings = _check(PIP_LOOP)
    f = [x for x in findings if x.rule == "lustre_interpreter_in_loop"][0]
    assert f.severity == "warn"


# --- Rule: tree walks over workspace / HOME / mount -> refuse -----------

WALK = "#!/bin/bash\ndu -sh /pfs/data/workspace/\n"


def test_tree_walk_over_workspace_refuses():
    findings = _check(WALK)
    f = [x for x in findings if x.rule == "tree_walk_on_lustre"][0]
    assert f.severity == "refuse"


def test_tree_walk_over_job_own_dir_is_fine():
    script = "#!/bin/bash\ndu -sh \"$TMPDIR/out\"\n"
    assert _check(script) == []


def test_find_without_recursion_is_fine():
    script = "#!/bin/bash\nfind /pfs/data -maxdepth 1 -name '*.inp'\n"
    assert _check(script) == []


# --- Rule: scratch on Lustre -> warn ------------------------------------

SCRATCH = "#!/usr/bin/env python\n# %scratch /pfs/data/scratch\nimport os\nos.environ['XTB_TMPDIR'] = '/pfs/data/scratch'\n"


def test_scratch_on_lustre_warns():
    findings = _check(SCRATCH)
    f = [x for x in findings if x.rule == "scratch_on_lustre"][0]
    assert f.severity == "warn"
    assert "$TMPDIR" in f.instead


# --- Rule: pytest straight from a HOME/Lustre checkout -> refuse --------

PYTEST = "#!/bin/bash\ncd /pfs/home/ka_user/delfin\npytest -q tests/\n"


def test_pytest_from_home_checkout_refuses():
    findings = _check(PYTEST)
    f = [x for x in findings if x.rule == "pytest_from_lustre_checkout"][0]
    assert f.severity == "refuse"
    assert "git archive" in f.instead or "tar" in f.instead


def test_pytest_from_tmpdir_checkout_is_fine():
    script = "#!/bin/bash\ncd \"$TMPDIR/delfin-copy\"\npytest -q tests/\n"
    assert _check(script) == []


# --- Shape of the API ----------------------------------------------------

def test_check_returns_findings_with_all_fields():
    findings = _check(PY_LOOP)
    for f in findings:
        assert f.rule and f.why and f.instead
        assert f.severity in ("refuse", "warn")
        assert f.line >= 1


def test_empty_script_has_no_findings():
    assert _check("# nothing here\n") == []
