"""On a shared cluster, "cancel everything" is never the agent's call.

The host-process guard already stops the agent killing the session it is
running in — that was the pkill incident. One step out from that, nothing
looked: `killall -u <user>` ends that person's editor, shell and running
calculations, and `scancel -u <user>` ends every job they have queued,
including ones this session knows nothing about and ones a colleague is
waiting on.

Hours of other people's compute, and no version of "the agent meant well"
gives it back. Cancelling ONE named job stays reachable and goes through
the ordinary gate — the point is the difference between a job and a user.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def perms():
    return A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))


@pytest.mark.parametrize("cmd,label", [
    ("scancel -u ka_ew7404", "every job of a user"),
    ("scancel --user=ka_ew7404", "same, long form"),
    ("scancel -u $USER --state=PENDING", "with extra flags"),
    ("killall -u ka_ew7404", "every process of a user"),
    ("pkill -u ka_ew7404", "same via pkill"),
    ("pkill -U 1001", "numeric uid"),
    ("scontrol suspend 12345", "suspending someone's job"),
])
def test_user_wide_destruction_is_denied(perms, cmd, label):
    assert perms.matches_bash_deny(cmd), label


@pytest.mark.parametrize("cmd,label", [
    ("scancel 123456", "cancelling one named job"),
    ("squeue -u $USER", "looking at the queue"),
    ("sbatch job.sh", "submitting"),
    ("sacct -j 123456", "checking accounting"),
    ("scontrol show job 123456", "inspecting a job"),
])
def test_ordinary_cluster_work_still_runs(perms, cmd, label):
    """A guard that stops the work is not a guard, and this is a chemistry
    framework whose whole purpose is submitting and watching jobs."""
    assert not perms.matches_bash_deny(cmd), label


def test_the_deny_survives_a_compound_command(perms):
    """Hiding it behind a safe first segment is the oldest trick there is."""
    assert perms.matches_bash_deny("squeue -u $USER && scancel -u $USER")


def test_the_host_process_guard_is_still_there():
    """This is one step out from that incident, not a replacement for it.

    The ancestry is supplied rather than read from /proc: in a test process
    there is no Voila ancestor, so asserting against the live chain would
    pass or fail for reasons that have nothing to do with the guard.
    """
    ancestry = [(4242, "voila", "voila --port=8866 dashboard.ipynb")]
    assert A._kill_targets_host_process("pkill -f voila", ancestry)
    assert A._kill_targets_host_process("kill -9 4242", ancestry)
    assert not A._kill_targets_host_process("kill -9 999999", ancestry)
