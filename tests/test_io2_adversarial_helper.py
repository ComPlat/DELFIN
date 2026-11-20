"""IO2 adversarial tests (reviewer nacht-s2): the shared submit-time helper.

``hpc_io_check.check_submit(script, env)`` is the one place both submit
routes and ``cluster_submit`` are meant to agree. A helper that decides
"path or text?" wrongly fails open, and a fail-open here is worse than no
check at all: the route reports "clean" and the person believes it was
looked at.

Every test states the measured value on the reviewed commit, so a fix can
be told from a re-tuned test.

Reviewed: builder nacht-s1, commit 9bdabb8b9 (IO2 phase 1).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin import hpc_io_check

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


# --- 1. A path that names nothing must not be read as script text ---------

def test_a_str_path_that_does_not_exist_is_an_error_not_a_pass(lustre, tmp_path):
    """The helper's docstring promises "Raises OSError when a script path
    cannot be read -- the check never fails open". Measured: that holds for a
    ``Path`` (the builder's own test) and NOT for a ``str``. A str naming no
    file falls through to "this is script text" (check_submit, the
    ``Path(script).is_file()`` branch) and the path itself is checked as a
    one-line script, which contains no command -- so it returns
    ``([], [])``: a clean verdict for a script nobody read.

    This is the same shape as IO1's L1, which the QS already made a finding.
    """
    with pytest.raises(OSError):
        hpc_io_check.check_submit(str(tmp_path / "missing.sh"),
                                  env={"HOME": "/home/x"})


def test_the_two_spellings_of_one_missing_path_agree(lustre, tmp_path):
    """A Path raises; the same name as a str must not take another route to
    a verdict. Measured: Path raises OSError, str returns ([], [])."""
    missing = tmp_path / "missing.sh"
    outcomes = []
    for spelling in (missing, str(missing)):
        try:
            outcomes.append(hpc_io_check.check_submit(
                spelling, env={"HOME": "/home/x"}))
        except OSError:
            outcomes.append("OSError")
    assert outcomes[0] == outcomes[1], (
        "the same script, named as Path and as str, must get the same verdict")


# --- 2. What the helper does correctly must stay correct ------------------

def test_a_real_path_is_read_and_refused(lustre, tmp_path):
    script = tmp_path / "run.sh"
    script.write_text(WALK, encoding="utf-8")
    for spelling in (script, str(script)):
        refused, warned = hpc_io_check.check_submit(
            spelling, env={"HOME": "/home/x"})
        assert refused and "tree_walk_on_lustre" in refused[0], spelling


def test_script_text_with_newlines_is_never_taken_for_a_path(lustre):
    refused, _ = hpc_io_check.check_submit(WALK, env={"HOME": "/home/x"})
    assert refused and "tree_walk_on_lustre" in refused[0]


def test_a_clean_script_is_clean(lustre):
    refused, warned = hpc_io_check.check_submit(
        "#!/bin/bash\ntar -xf env.tar -C $TMPDIR\n", env={"HOME": "/home/x"})
    assert refused == [] and warned == []


# --- 3. A refusal must be actionable --------------------------------------

def test_a_refusal_says_which_line(lustre):
    """A person fixing the script is told the file and the line today
    (cluster_tools.submit formats ``name:line [rule]``). The helper drops
    both. Measured: the refused entry is ``[tree_walk_on_lustre] ...`` with
    no line number, so the same helper cannot report what IO1 reported."""
    script = "#!/bin/bash\necho one\nfind /pfs/data6 -type f\n"
    refused, _ = hpc_io_check.check_submit(
        script, env={"HOME": "/home/x"})
    assert refused
    assert "3" in refused[0], (
        "the refusal must name the line the walk is on (it is line 3)")


# --- 4. The package's headline claim: the routes cannot drift -------------

def test_the_helper_reports_exactly_what_cluster_tools_reports(lustre, tmp_path):
    """The commit's stated purpose is that the dashboard backend, the bash
    gate and ``cluster_submit`` "cannot drift apart in what they call a
    refusal". Measured: they already differ -- cluster_tools.submit prefixes
    ``name:line`` (cluster_tools.py:403-404) and check_submit does not, so
    the same script yields two different strings and a test comparing them
    would be comparing two formats. Whichever form is chosen, both must
    produce it."""
    script = "#!/bin/bash\necho one\nfind /pfs/data6 -type f\n"
    findings = hpc_io_check.check(script, env={"HOME": "/home/x"})
    as_cluster_tools_reports = [
        f"run.sh:{f.line} [{f.rule}] {f.why} instead: {f.instead}"
        for f in findings if f.severity == "refuse"]
    path = tmp_path / "run.sh"
    path.write_text(script, encoding="utf-8")
    refused, _ = hpc_io_check.check_submit(path, env={"HOME": "/home/x"})
    assert refused == as_cluster_tools_reports


# --- 5. A caller that forgets `env` gets a weaker check -------------------

def test_forgetting_env_expands_home_from_os_environ(lustre, monkeypatch):
    """The mechanism, not this host's HOME: ``env`` defaults to None and
    ``check_submit`` must then expand $HOME from ``os.environ`` so the rule
    still works. Pinning *this host's* HOME on /pfs would assert machine
    state instead (H4 of the first round, corrected: I had written it that
    way and the builder was right to push back)."""
    script = '#!/bin/bash\nfind "$HOME/work" -type f\n'

    # A Lustre HOME in the process env: no env passed, the rule must fire.
    monkeypatch.setenv("HOME", "/pfs/home/u")
    assert hpc_io_check.check_submit(script)[0], (
        "env omitted must expand $HOME from os.environ and refuse a Lustre "
        "home")

    # A node-local HOME in the process env: nothing to refuse, and that is
    # correct -- it proves the refusal above came from the HOME, not luck.
    monkeypatch.setenv("HOME", "/home/u")
    assert not hpc_io_check.check_submit(script)[0]


def test_an_explicit_env_wins_over_os_environ(lustre, monkeypatch):
    """The other half: a caller that knows the cluster HOME passes it, and
    that value is used rather than the process's."""
    script = '#!/bin/bash\nfind "$HOME/work" -type f\n'
    monkeypatch.setenv("HOME", "/home/u")
    assert hpc_io_check.check_submit(script, env={"HOME": "/pfs/home/u"})[0]
    monkeypatch.setenv("HOME", "/pfs/home/u")
    assert not hpc_io_check.check_submit(script, env={"HOME": "/home/u"})[0]
