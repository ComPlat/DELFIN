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
                                  text=LUSTRE_MOUNTS if lustre
                                  else "/dev/sda1 / ext4 rw 0 0"))


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
    src = hpc_io_check.MountsSource(path="/nonexistent")
    with pytest.raises(RuntimeError):
        hpc_io_check.check(script, env=env, mounts_path=src)


def test_lustre_detection_from_mounts_text():
    src = hpc_io_check.MountsSource(text=LUSTRE_MOUNTS)
    assert hpc_io_check.lustre_paths(src) == ["/pfs"]


def test_no_lustre_mount_means_no_lustre():
    src = hpc_io_check.MountsSource(text="/dev/sda1 / ext4 rw 0 0")
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


# --- Operator live-run round (compute node @ 7a219348a): L1-L5 ------------

def test_default_mounts_source_reads_the_file_not_the_path_string():
    # L1 (operator compute-node run): MountsSource("/proc/mounts") stored
    # the path STRING as mounts text -> lustre_paths()==[] -> check()
    # returned [] for every script in production. A positional argument
    # is a PATH and is read; the ambiguous positional-text form is gone.
    with pytest.raises(RuntimeError):
        hpc_io_check.check("for i in 1 2; do x; done",
                           env={"HOME": "/pfs/home/u"},
                           mounts_path=hpc_io_check.MountsSource(
                               "/nonexistent/mounts"))


def test_mounts_source_rejects_the_ambiguous_two_argument_form():
    # L1 contract: passing mounts TEXT positionally was the bug shape
    # (a bare string is a PATH under the new contract, so this would
    # silently try to open a file named like the mounts text). Only one
    # form may be given; two forms together are rejected.
    with pytest.raises(TypeError):
        hpc_io_check.MountsSource("/x", text=LUSTRE_MOUNTS)


def test_default_wiring_finds_lustre_from_real_mounts_file(tmp_path,
                                                            monkeypatch):
    # L1 continued: with the constructor fixed, check() called with NO
    # mounts_path must read /proc/mounts (monkeypatched to a real file).
    mounts = tmp_path / "mounts"
    mounts.write_text("mgmt@o2ib:fs /pfs lustre rw 0 0\n")
    script = "for i in 1 2; do du -sh /pfs/data/ws; done\n"
    monkeypatch.setattr(hpc_io_check, "DEFAULT_MOUNTS", str(mounts))
    findings = hpc_io_check.check(script, env={"HOME": "/pfs/home/u"})
    assert any(f.rule == "tree_walk_on_lustre" for f in findings)


def test_absolute_path_interpreter_in_loop_is_flagged():
    # L3: the 14.11. shape — an absolute-path venv python in a loop.
    script = ("for i in 1 2 3; do\n"
              "  /pfs/data6/envs/run/bin/python -c 'import os'\n"
              "done\n")
    findings = _check(script)
    assert any(f.rule == "lustre_interpreter_in_loop" for f in findings)


def test_python_dash_m_pytest_from_home_is_refused():
    # L4: the 16.11. shape — `python -m pytest` from a HOME checkout.
    script = ("cd /pfs/home/ka_user/delfin\n"
              ".venv/bin/python -m pytest -q tests/\n")
    findings = _check(script)
    assert any(f.rule == "pytest_from_lustre_checkout" for f in findings)


def test_stage_in_copy_to_tmpdir_is_not_a_tree_walk():
    # L5 FALSE REFUSAL: the rule ASKS for a stage-in. One job folder
    # copied TO a node-local destination is the fix, not the problem.
    script = "cp -r /pfs/data6/ws/runs/water $TMPDIR/\n"
    assert _check(script) == []


def test_stage_in_rsync_to_tmpdir_is_not_a_tree_walk():
    script = "rsync -a /pfs/data6/ws/runs/water/ $TMPDIR/water/\n"
    assert _check(script) == []


def test_walk_over_home_root_still_refused():
    # L5 guard direction: walks over HOME/workspace/mount roots stay
    # refused even though $TMPDIR destinations are now allowed.
    script = "cp -r /pfs/home/ka_user/projects $TMPDIR/projects\n"
    findings = _check(script)
    assert any(f.rule == "tree_walk_on_lustre" for f in findings)


def test_bad_pytest_home_case_text_is_refused():
    # Live-run round 2, L4: the .gate/io1_cases/bad_pytest_home.sh text
    # verbatim — cd onto a Lustre checkout, then a RELATIVE interpreter
    # running `python -m pytest`. The pytest line carries no Lustre path
    # itself; the Lustre-ness comes from the cd above it.
    script = (
        "#!/bin/bash\n"
        "cd /pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin"
        "/.delfin/worktrees/runde3-s1\n"
        ".venv/bin/python -m pytest -q tests/\n"
    )
    findings = _check(script)
    assert any(f.rule == "pytest_from_lustre_checkout" for f in findings)

