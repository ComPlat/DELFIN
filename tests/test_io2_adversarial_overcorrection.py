"""IO2 adversarial tests (reviewer nacht-s2), round 4: the fix must not
over-correct.

A check that is wired into every submit route is a gate. Closing a
false-negative (a walk that slipped through) is only half the job: if the
same change starts refusing ordinary submit scripts, it blocks real work
and the next person will weaken the rule to get their job in. These are the
ordinary scripts and the newly-caught scripts, side by side.

Measured on builder nacht-s1's 19042b858 (W1-W3, D1, D2).
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/pfs/home/user"}


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


def refusals(script):
    return [(f.line, f.rule) for f in hpc_io_check.check(
        script, env=ENV) if f.severity == "refuse"]


# --- 9. Ordinary submit scripts must keep working ------------------------

@pytest.mark.parametrize("script", [
    # The canonical fix the rule itself prescribes.
    "#!/bin/bash\ncd $TMPDIR\ntar -xf env.tar\npython run.py\n",
    # Stage the job's own folder to the node, then run there.
    '#!/bin/bash\ncp -a inputs "$TMPDIR"/\ncd $TMPDIR\nmpirun orca in.inp\n',
    # One-directory reads are explicitly allowed by the rule's own advice.
    "#!/bin/bash\nfind . -maxdepth 1 -name '*.inp'\n",
    "#!/bin/bash\ncd /pfs/data6/work\nfind . -maxdepth 1 -name '*.out'\n",
    # A stage-in of ONE job folder (IO1's L5) stays allowed.
    '#!/bin/bash\nrsync -a /pfs/data6/work/run42/ "$TMPDIR"/\n',
    # Wrappers on a harmless command.
    "#!/bin/bash\nsleep 5\nnice -n 10 echo hello\n",
    "#!/bin/bash\nenv | sort\n",
    # The loop rule is about interpreters on Lustre, not about reading a
    # queue file to decide what to launch.
    '#!/bin/bash\nfor f in *.inp; do python run.py "$f"; done\n',
])
def test_an_ordinary_submit_script_is_not_refused(lustre, script):
    """Each of these is what a working job script looks like, and each must
    pass. This is the half of the rule that keeps it usable."""
    assert not refusals(script)


# --- 10. What the round-2 fix closed stays closed ------------------------

@pytest.mark.parametrize("script", [
    "nice -n 10 find /pfs/data6 -type f\n",
    "nice -19 du -sh /pfs/data6\n",
    "env X=1 find /pfs/data6 -type f\n",
    "env A=1 B=2 du -sh /pfs/data6\n",
    "env -i find /pfs/data6 -type f\n",
    "bash -c 'find /pfs/data6 -type f' zero\n",
    "sh -c 'du -sh /pfs/data6' \"$0\"\n",
])
def test_a_wrapper_with_arguments_still_does_not_hide_the_walk(lustre, script):
    """W1/W2 on 19042b858: the wrapper spellings that carry their own
    argument. All were ``[]`` before this commit."""
    assert refusals(script)


# --- 11. A chained cd, in both directions --------------------------------

def test_a_chain_that_ends_on_lustre_is_refused(lustre):
    """W3, false-negative direction: the LAST cd is the cwd."""
    assert refusals("cd /tmp && cd /pfs/data6/w\npytest tests/\n")


def test_a_chain_that_ends_node_local_is_allowed(lustre):
    """W3, false-positive direction: the script really does end on
    node-local /tmp, so refusing it would be wrong. This was refused before
    19042b858 -- the more instructive half of the bug."""
    assert not refusals("cd /pfs/data6/w && cd /tmp\npytest tests/\n")
    assert not refusals("cd /tmp && cd /pfs/data6/w && cd /tmp\npytest tests/\n")


def test_a_semicolon_chain_ends_on_lustre(lustre):
    assert refusals("cd /tmp; cd /pfs/data6/w\npytest tests/\n")
