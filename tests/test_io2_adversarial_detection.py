"""IO2 adversarial tests (reviewer nacht-s2): does the check SEE the walks it
refuses, and does it stay silent only where silence is correct?

Package IO2 puts ``hpc_io_check.check`` in front of two more submit routes
(the dashboard's ``backend_slurm._sbatch`` and an ``sbatch`` the agent types
into bash). A check that is wired into a route but cannot see the script it
is judging adds a step without adding protection: the second route then
reports "clean" for exactly the scripts the first mail was about.

Every test here names the behaviour it protects. Each one is red on
``7abca4483`` and states the measured current value in its docstring, so a
reader can tell a fix from a re-tuning of the test.

Reference: ``delfin/hpc_io_check.py`` (the module under test);
``tests/test_io1_hpc_io_check.py`` (IO1's own suite, unchanged by this file).
"""
from pathlib import Path

import pytest

from delfin import hpc_io_check

REPO = Path(__file__).resolve().parents[1]
TEMPLATES = REPO / "delfin" / "submit_templates"

# /pfs is the mount point the rules are about; "user" is a placeholder, not
# an account of anybody.
LUSTRE_MOUNTS = """\
/dev/sda1 / ext4 rw,relatime 0 0
10.1.2.3@o2ib:/pfs /pfs lustre rw,noatime 0 0
tmpfs /tmp tmpfs rw 0 0
"""
LUSTRE_HOME = {"HOME": "/pfs/home/user"}
MOUNTS = hpc_io_check.MountsSource(text=LUSTRE_MOUNTS)


def findings_of(script, env=None):
    """Every finding, as (line, severity, rule) -- the shape a route
    reports to the user."""
    return [(f.line, f.severity, f.rule) for f in hpc_io_check.check(
        script, env=LUSTRE_HOME if env is None else env, mounts_path=MOUNTS)]


def refusals_of(script, env=None):
    return [f for f in findings_of(script, env=env) if f[1] == "refuse"]


# --- 1. A bare `pytest` line, the shape a submit script actually has -------

def test_a_bare_pytest_alone_on_a_line_is_refused():
    """A line that ENDS with `pytest` is the commonest spelling in a job
    script. ``_INTERPRETER`` (hpc_io_check.py:67) needs a character after the
    word, so today such a line is never looked at: measured ``[]``."""
    assert refusals_of("cd /pfs/home/user/work\npytest\n"), (
        "a bare `pytest` as the last line walks the Lustre checkout just as "
        "surely as `pytest tests/` does")


def test_a_bare_pytest_in_the_middle_of_a_script_is_refused():
    """Measured today: ``[]``. The following line must not be what decides
    whether the walk is seen."""
    assert refusals_of("cd /pfs/home/user/work\npytest\necho done\n")


# --- 2. The node-local `cd` exemption is per script, not per line ----------

def test_a_later_cd_back_onto_lustre_is_not_excused_by_an_earlier_tmpdir_cd():
    """``here`` is searched over the WHOLE text (hpc_io_check.py:279) and the
    exemption fires on any ``cd $VAR`` anywhere (:283). The last `cd` is the
    one that decides the cwd, so this script walks Lustre. Measured ``[]``."""
    script = (
        "cd $TMPDIR\n"
        "cd /pfs/home/user/work\n"
        "pytest tests/\n"
    )
    assert refusals_of(script), "the last cd wins; this still walks Lustre"


def test_the_genuine_node_local_fix_stays_allowed():
    """The exemption exists for a reason: a checkout really unpacked to the
    node must not be refused. Measured today: allowed (and it must stay)."""
    assert not refusals_of("cd $TMPDIR\ntar -xf env.tar\npytest tests/\n")


# --- 3. A shell variable must not hide the walk ---------------------------

def test_a_variable_holding_a_lustre_path_is_refused():
    """``_expand_env_paths`` (:108) expands $HOME only, so every other
    variable -- the shape DELFIN's own template uses for its staging copies
    -- hides the path completely. Measured ``[]``; the same line with the
    path written out refuses."""
    script = 'DELFIN_DIR=/pfs/home/user/software/delfin\nfind "$DELFIN_DIR/delfin" -type f\n'
    assert refusals_of(script), (
        "a walk of a Lustre path must be refused whether the path is "
        "written out or carried in a variable assigned in the script")


def test_a_variable_from_the_environment_is_refused():
    """The same hole one step out: the path arrives in the env dict the
    route already passes. Measured ``[]``."""
    env = dict(LUSTRE_HOME, DELFIN_DIR="/pfs/home/user/software/delfin")
    assert refusals_of('find "$DELFIN_DIR/delfin" -type f\n', env=env)


# --- 4. A wrapper in front of the command must not hide it ----------------

@pytest.mark.parametrize("script", [
    "bash -c 'find /pfs/home/user -type f'\n",
    "sh -c 'du -sh /pfs/home/user'\n",
    "command find /pfs/home/user -type f\n",
    "nice find /pfs/home/user -type f\n",
    "exec find /pfs/home/user -type f\n",
    "( find /pfs/home/user -type f )\n",
])
def test_a_wrapper_does_not_hide_the_walk(script):
    """The rule anchors the command word at the start of the line
    (:304), so any wrapper hides the walk. Measured ``[]`` for every case
    here; the un-wrapped form refuses."""
    assert refusals_of(script)


def test_the_control_that_these_tests_are_measured_against():
    """If this one ever stops firing, every test above is meaningless."""
    assert refusals_of("find /pfs/home/user -type f\n")


# --- 5. DELFIN's own submit script -----------------------------------------

@pytest.mark.parametrize("name,substitutions", [
    # submit_delfin.sh:441 walks the checkout with `find ... -type f`.
    ("submit_delfin.sh", {
        "${DELFIN_DIR}": "/pfs/home/user/software/delfin",
        "$DELFIN_DIR": "/pfs/home/user/software/delfin",
    }),
    # submit_turbomole.sh:11 copies the run back with `rsync -a`; the job
    # directory it copies INTO is the submit directory.
    ("submit_turbomole.sh", {
        "$JOB_DIR": "/pfs/home/user/run",
        "${JOB_DIR}": "/pfs/home/user/run",
    }),
])
def test_a_template_whose_paths_are_written_out_is_still_judged(name, substitutions):
    """Phase 2 asks that DELFIN's generated script pass clean. It passes with
    the variable and refuses with the path written out -- so "clean" today is
    a property of the blindness above, not of the template. Measured: the
    rendered ``submit_delfin.sh`` gives 3 refusals, the shipped text 0."""
    text = (TEMPLATES / name).read_text(encoding="utf-8")
    rendered = text
    for needle, path in substitutions.items():
        rendered = rendered.replace(needle, path)
    assert refusals_of(rendered), (
        "a submit script must be judged on the paths it will use, not on "
        "the variable names that carry them")
